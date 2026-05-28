"""
resp_workflow
=============

Full RESP charge derivation workflow:

    capped MOL2  --xyz-->  xTB pre-opt  -->  ORCA opt  -->  ORCA SP (GBW)
                                                                   |
                                                       orca_2mkl -molden
                                                                   |
                                                                   v
                                                       Multiwfn RESP fit
                                                                   |
                                              antechamber -c rc -cf RESP.chg
                                                                   |
                                                                   v
                                                          final RESP MOL2

Requires the following executables (resolved via CLI / env var / PATH):
    orca, orca_2mkl, Multiwfn_noGUI, xtb

Environment variables honoured:
    NSAA_ORCA_EXE         : explicit path to orca
    NSAA_ORCA_2MKL_EXE    : explicit path to orca_2mkl
    NSAA_MULTIWFN_EXE     : explicit path to Multiwfn_noGUI
    NSAA_XTB_EXE          : explicit path to xtb
    NSAA_RESP_MULTIPLICITY: integer overriding the guessed spin multiplicity
    NSAA_RESP_NPROCS      : number of ORCA processes (default 4)
"""
from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional


# =========================
# EXECUTABLE RESOLVER
# =========================
def _resolve_executable(*, cli_value: Optional[str], env_var: str, program_name: str) -> str:
    if cli_value:
        exe = Path(cli_value)
        if not exe.exists():
            raise FileNotFoundError(f"{program_name} not found at {exe}")
        return str(exe.resolve())

    env_val = os.environ.get(env_var)
    if env_val:
        exe = Path(env_val)
        if not exe.exists():
            raise FileNotFoundError(f"{program_name} not found at {exe}")
        return str(exe.resolve())

    found = shutil.which(program_name)
    if found:
        return str(Path(found).resolve())

    raise RuntimeError(f"{program_name} not found. Provide via CLI, ENV or PATH.")


# =========================
# DATA STRUCTURES
# =========================
@dataclass(frozen=True)
class Mol2Atom:
    atom_id: int
    name: str
    x: float
    y: float
    z: float
    atom_type: str
    subst_id: int
    subst_name: str
    charge: float


_TWO_LETTER_ELEMENTS = {
    "CL": "Cl", "BR": "Br", "NA": "Na", "MG": "Mg",
    "ZN": "Zn", "FE": "Fe", "CA": "Ca", "CU": "Cu", "MN": "Mn",
}

_ATOMIC_NUMBERS = {
    "H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
    "NA": 11, "MG": 12, "P": 15, "S": 16, "CL": 17,
    "K": 19, "CA": 20, "MN": 25, "FE": 26, "CU": 29,
    "ZN": 30, "BR": 35, "I": 53,
}


# =========================
# BASIC HELPERS
# =========================
def _parse_mol2_atoms(mol2_path: str) -> List[Mol2Atom]:
    atoms: List[Mol2Atom] = []
    in_atoms = False
    with open(mol2_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if s.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if s.startswith("@<TRIPOS>") and in_atoms:
                break
            if not in_atoms or not s:
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            charge = float(parts[8]) if len(parts) >= 9 else 0.0
            atoms.append(
                Mol2Atom(
                    atom_id=int(parts[0]),
                    name=parts[1],
                    x=float(parts[2]),
                    y=float(parts[3]),
                    z=float(parts[4]),
                    atom_type=parts[5],
                    subst_id=int(parts[6]),
                    subst_name=parts[7],
                    charge=charge,
                )
            )
    if not atoms:
        raise ValueError(f"No atoms parsed from MOL2: {mol2_path}")
    return atoms


def _guess_element(atom_name: str, atom_type: str) -> str:
    at_alpha = re.sub(r"[^A-Za-z]", "", str(atom_type or "")).upper()

    if at_alpha:
        for token, proper in _TWO_LETTER_ELEMENTS.items():
            if at_alpha.startswith(token):
                return proper
        symbol = at_alpha[0].upper()
        if symbol in _ATOMIC_NUMBERS:
            return symbol

    nm_alpha = re.sub(r"[^A-Za-z]", "", str(atom_name or "")).upper()
    if nm_alpha:
        symbol = nm_alpha[0].upper()
        if symbol in _ATOMIC_NUMBERS:
            return symbol

    return "C"


def _atomic_number(symbol: str) -> int:
    return _ATOMIC_NUMBERS.get(symbol.upper(), 0)


def _electron_count(mol2_path: str, charge: int) -> int:
    atoms = _parse_mol2_atoms(mol2_path)
    return sum(_atomic_number(_guess_element(a.name, a.atom_type)) for a in atoms) - int(charge)


def _guess_multiplicity(mol2_path: str, charge: int) -> int:
    electrons = _electron_count(mol2_path, charge)
    return 1 if electrons % 2 == 0 else 2


def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return int(raw)
    except Exception as exc:
        raise ValueError(f"Environment variable {name} must be an integer, got: {raw}") from exc


def _run(cmd: list[str], *, cwd: Path, log_path: Path) -> subprocess.CompletedProcess:
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(cwd))
    with log_path.open("a", encoding="utf-8") as log:
        log.write("\n\n=== CMD ===\n")
        log.write(f"CWD: {cwd}\n")
        log.write("CMD: " + " ".join(shlex.quote(x) for x in cmd) + "\n")
        log.write("=== STDOUT ===\n")
        log.write(result.stdout)
        log.write("\n=== STDERR ===\n")
        log.write(result.stderr)
        log.write("\n=== RETURN CODE ===\n")
        log.write(str(result.returncode) + "\n")
    return result


def _write_xyz_from_mol2(mol2_path: Path, xyz_path: Path) -> list[Mol2Atom]:
    atoms = _parse_mol2_atoms(str(mol2_path))

    lines = [str(len(atoms)), f"Generated from {mol2_path.name}"]
    for atom in atoms:
        element = _guess_element(atom.name, atom.atom_type)
        lines.append(f"{element:<2} {atom.x: .10f} {atom.y: .10f} {atom.z: .10f}")

    xyz_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return atoms


def _read_xyz_coordinates(xyz_path: Path) -> list[tuple[float, float, float]]:
    lines = xyz_path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 3:
        raise ValueError(f"XYZ file is too short: {xyz_path}")
    try:
        natoms = int(lines[0].strip())
    except Exception as exc:
        raise ValueError(f"Invalid XYZ header in {xyz_path}") from exc

    coords: list[tuple[float, float, float]] = []
    for line in lines[2:2 + natoms]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Malformed XYZ atom line in {xyz_path}: {line}")
        coords.append((float(parts[1]), float(parts[2]), float(parts[3])))

    if len(coords) != natoms:
        raise ValueError(
            f"XYZ atom count mismatch in {xyz_path}: expected {natoms}, got {len(coords)}"
        )
    return coords


def _run_xtb(xyz: Path, work_dir: Path, xtb_exe: str, charge: int, log: Path) -> Path:
    cmd = [xtb_exe, xyz.name, "--opt", "--chrg", str(charge)]

    env = os.environ.copy()
    if "XTBPATH" not in env:
        xtb_share = Path(xtb_exe).parent.parent / "share" / "xtb"
        if xtb_share.exists():
            env["XTBPATH"] = str(xtb_share)

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(work_dir),
        env=env,
    )

    with log.open("a", encoding="utf-8") as fh:
        fh.write("\n\n=== CMD ===\n")
        fh.write(f"CWD: {work_dir}\n")
        fh.write("CMD: " + " ".join(shlex.quote(x) for x in cmd) + "\n")
        fh.write("=== STDOUT ===\n")
        fh.write(result.stdout)
        fh.write("\n=== STDERR ===\n")
        fh.write(result.stderr)
        fh.write("\n=== RETURN CODE ===\n")
        fh.write(str(result.returncode) + "\n")

    if result.returncode != 0:
        raise RuntimeError(
            "xTB failed.\n"
            f"STDOUT:\n{result.stdout}\n\n"
            f"STDERR:\n{result.stderr}"
        )

    out = work_dir / "xtbopt.xyz"
    if not out.exists():
        raise FileNotFoundError("xtbopt.xyz not found")

    return out


def _write_orca_opt_input(xyz: Path, out: Path, charge: int, mult: int, nprocs: int):
    text = f"""! HF 6-31G* TightSCF Opt TightOpt

%pal
  nprocs {nprocs}
end

%geom
  MaxIter 300
end

* xyzfile {charge} {mult} {xyz.name}
"""
    out.write_text(text, encoding="utf-8")


def _write_orca_sp_input(xyz: Path, out: Path, charge: int, mult: int, nprocs: int):
    text = f"""! HF 6-31G* TightSCF

%pal
  nprocs {nprocs}
end

* xyzfile {charge} {mult} {xyz.name}
"""
    out.write_text(text, encoding="utf-8")


def _parse_multiwfn_resp_output(raw_output: Path) -> list[float]:
    lines = raw_output.read_text(encoding="utf-8", errors="replace").splitlines()
    start_index = None
    for i, line in enumerate(lines):
        if "Successfully converged!" in line:
            start_index = i + 3
            break
    if start_index is None:
        raise RuntimeError(
            "Could not find 'Successfully converged!' in Multiwfn output. RESP fitting may have failed."
        )

    charges: list[float] = []
    for line in lines[start_index:]:
        s = line.strip()
        if not s:
            continue
        if "Sum of charges:" in s:
            break
        if ")" not in s:
            continue
        _, right = s.split(")", 1)
        try:
            charge = float(right.strip())
        except Exception:
            continue
        charges.append(charge)

    if not charges:
        raise RuntimeError("No RESP charges parsed from Multiwfn output.")
    return charges


def _write_antechamber_charge_file(charges: list[float], charge_file: Path) -> None:
    """Write one RESP charge per line in the same atom order as the MOL2 template."""
    charge_file.write_text(
        "\n".join(f"{q:.10f}" for q in charges) + "\n",
        encoding="utf-8",
    )


def _run_antechamber_resp_rc(
    *,
    input_mol2: Path,
    output_mol2: Path,
    charge_file: Path,
    resname: str,
    net_charge: int,
    log_path: Path,
) -> None:
    antechamber_exe = shutil.which("antechamber")
    if not antechamber_exe:
        raise RuntimeError("antechamber not found in PATH")

    cmd = [
        antechamber_exe,
        "-fi", "mol2",
        "-i", str(input_mol2),
        "-fo", "mol2",
        "-o", str(output_mol2),
        "-c", "rc",
        "-cf", str(charge_file),
        "-at", "amber",
        "-bk", resname,
        "-nc", str(net_charge),
        "-s", "2",
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(output_mol2.parent),
    )

    with log_path.open("a", encoding="utf-8") as log:
        log.write("\n\n=== ANTECHAMBER RESP RC ===\n")
        log.write("CMD: " + " ".join(shlex.quote(x) for x in cmd) + "\n")
        log.write("=== STDOUT ===\n")
        log.write(result.stdout)
        log.write("\n=== STDERR ===\n")
        log.write(result.stderr)
        log.write("\n=== RETURN CODE ===\n")
        log.write(str(result.returncode) + "\n")

    if result.returncode != 0:
        raise RuntimeError(
            f"Antechamber RESP RC failed for {resname}. See log file: {log_path}"
        )

    if not output_mol2.exists() or output_mol2.stat().st_size == 0:
        raise RuntimeError(
            f"Antechamber did not create final RESP MOL2 for {resname}: {output_mol2}"
        )


# =========================
# ROBUST MULTIWFN DRIVER
# =========================
def _expect_or_raise(child, patterns, step: str, timeout: int):
    try:
        return child.expect(patterns, timeout=timeout)
    except Exception as exc:
        before = getattr(child, "before", "") or ""
        raise RuntimeError(
            f"Multiwfn failed during step: {step}\n"
            f"{before[-2000:]}"
        ) from exc


def _run_multiwfn_resp(
    *,
    molden_path: Path,
    fitting_dir: Path,
    log_path: Path,
    multiwfn_exe: str,
) -> list[float]:
    import pexpect

    if not molden_path.exists():
        raise FileNotFoundError(f"Molden file not found: {molden_path}")

    fitting_dir.mkdir(parents=True, exist_ok=True)

    conformer_list = fitting_dir / "resp_conformers.txt"
    conformer_list.write_text(f"{molden_path.name}\n", encoding="utf-8")

    raw_output = fitting_dir / "MultiWfn_raw_outputs.txt"

    with raw_output.open("w", encoding="utf-8") as logfile:
        child = pexpect.spawn(
            multiwfn_exe,
            [molden_path.name],
            encoding="utf-8",
            logfile=logfile,
            cwd=str(fitting_dir),
            timeout=120,
        )

        child.expect([r".*Main.*menu.*", r".*Select.*function.*", r".*300.*"], timeout=120)
        child.sendline("7")

        child.expect([r".*ESP.*", r".*electrostatic.*", r".*20.*"], timeout=120)
        child.sendline("18")

        child.expect([r".*RESP.*", r".*fitting.*", r".*1.*"], timeout=120)
        child.sendline("1")

        try:
            idx = child.expect(
                [r".*vdW radius.*", r".*Input.*file.*"],
                timeout=10,
            )
            if idx == 0:
                child.sendline("")
        except Exception:
            pass

        child.sendline(conformer_list.name)

        converged_count = 0

        while True:
            try:
                idx = child.expect(
                    [
                        r".*Successfully converged!.*",
                        r".*\(y/n\).*",
                        pexpect.EOF,
                    ],
                    timeout=300,
                )

                if idx == 0:
                    converged_count += 1
                elif idx == 1:
                    child.sendline("n")
                    break
                elif idx == 2:
                    break
            except Exception:
                break

        child.close(force=True)

    if converged_count < 1:
        raise RuntimeError("Multiwfn did not reach RESP convergence")

    return _parse_multiwfn_resp_output(raw_output)


# =========================
# MAIN WORKFLOW
# =========================
def run_resp_charge_workflow(
    *,
    capped_file: str,
    charged_file: str,
    resname: str,
    net_charge: int,
    residue_dir: str,
    orca_path: Optional[str] = None,
    orca_2mkl_path: Optional[str] = None,
    multiwfn_path: Optional[str] = None,
    xtb_path: Optional[str] = None,
) -> dict:
    residue_path = Path(residue_dir)
    qm = residue_path / "resp_qm"
    fit = residue_path / "resp_fit"
    qm.mkdir(parents=True, exist_ok=True)
    fit.mkdir(parents=True, exist_ok=True)

    log = residue_path / f"{resname}_resp.log"
    log.write_text("", encoding="utf-8")
    final_mol2 = Path(charged_file).resolve()

    orca = _resolve_executable(cli_value=orca_path, env_var="NSAA_ORCA_EXE", program_name="orca")
    orca2mkl = _resolve_executable(cli_value=orca_2mkl_path, env_var="NSAA_ORCA_2MKL_EXE",
                                   program_name="orca_2mkl")
    multiwfn = _resolve_executable(cli_value=multiwfn_path, env_var="NSAA_MULTIWFN_EXE",
                                   program_name="Multiwfn_noGUI")
    xtb = _resolve_executable(cli_value=xtb_path, env_var="NSAA_XTB_EXE", program_name="xtb")

    print(f"[RESP] ORCA: {orca}")
    print(f"[RESP] ORCA_2MKL: {orca2mkl}")
    print(f"[RESP] Multiwfn: {multiwfn}")
    print(f"[RESP] xTB: {xtb}")

    mol2 = Path(capped_file)
    xyz = qm / "resp_input.xyz"
    atoms = _write_xyz_from_mol2(mol2, xyz)

    multiplicity = _env_int(
        "NSAA_RESP_MULTIPLICITY",
        _guess_multiplicity(str(mol2), int(net_charge)),
    )
    nprocs = _env_int("NSAA_RESP_NPROCS", 4)

    xtb_xyz = _run_xtb(xyz=xyz, work_dir=qm, xtb_exe=xtb, charge=int(net_charge), log=log)

    opt_inp = qm / "resp_opt.inp"
    _write_orca_opt_input(xyz=xtb_xyz, out=opt_inp, charge=int(net_charge),
                          mult=multiplicity, nprocs=nprocs)

    opt_result = _run([orca, str(opt_inp.resolve())], cwd=qm, log_path=log)
    if opt_result.returncode != 0:
        raise RuntimeError("ORCA optimization crashed")

    opt_xyz = qm / "resp_opt.xyz"
    if not opt_xyz.exists():
        raise FileNotFoundError("Optimized XYZ not found (resp_opt.xyz)")

    sp_inp = qm / "resp_sp.inp"
    _write_orca_sp_input(xyz=opt_xyz, out=sp_inp, charge=int(net_charge),
                         mult=multiplicity, nprocs=nprocs)

    sp_result = _run([orca, str(sp_inp.resolve())], cwd=qm, log_path=log)
    if sp_result.returncode != 0:
        raise RuntimeError("ORCA SP crashed")
    if "ORCA TERMINATED NORMALLY" not in sp_result.stdout:
        raise RuntimeError("ORCA SP did not terminate normally")

    gbw_path = qm / "resp_sp.gbw"
    if not gbw_path.exists():
        raise FileNotFoundError("SP GBW not found (resp_sp.gbw)")

    mkl_result = _run([orca2mkl, "resp_sp", "-molden"], cwd=qm, log_path=log)
    if mkl_result.returncode != 0:
        raise RuntimeError("orca_2mkl failed")

    molden = qm / "resp_sp.molden.input"
    if not molden.exists():
        raise RuntimeError("Molden file not generated")
    if molden.stat().st_size == 0:
        raise RuntimeError("Molden file is empty")

    copied_molden = fit / molden.name
    copied_molden.write_text(
        molden.read_text(encoding="utf-8", errors="replace"),
        encoding="utf-8",
    )

    if not copied_molden.exists():
        raise RuntimeError("Molden file missing before RESP")

    resp_charges = _run_multiwfn_resp(
        molden_path=copied_molden,
        fitting_dir=fit,
        log_path=log,
        multiwfn_exe=multiwfn,
    )

    if len(resp_charges) != len(atoms):
        raise RuntimeError(
            f"RESP returned {len(resp_charges)} charges, but molecule has {len(atoms)} atoms."
        )

    opt_coords = _read_xyz_coordinates(opt_xyz)

    if len(opt_coords) != len(atoms):
        raise RuntimeError(
            f"Optimized XYZ has {len(opt_coords)} atoms, but capped MOL2 has {len(atoms)} atoms."
        )

    resp_charge_file = residue_path / f"{resname}_resp.chg"
    _write_antechamber_charge_file(resp_charges, resp_charge_file)

    _run_antechamber_resp_rc(
        input_mol2=mol2,
        output_mol2=final_mol2,
        charge_file=resp_charge_file,
        resname=resname,
        net_charge=int(net_charge),
        log_path=log,
    )

    print("RESP pipeline completed successfully")

    return {
        "charge_method": "resp",
        "charge_backend": "xtb+orca+multiwfn+antechamber_rc",
        "net_charge": int(net_charge),
        "multiplicity": int(multiplicity),
        "files": {
            "capped_mol2": str(mol2),
            "input_xyz": str(xyz),
            "xtb_xyz": str(xtb_xyz),
            "opt_input": str(opt_inp),
            "opt_xyz": str(opt_xyz),
            "sp_input": str(sp_inp),
            "sp_gbw": str(gbw_path),
            "sp_molden": str(copied_molden),
            "resp_charge_file": str(resp_charge_file),
            "final_mol2": str(final_mol2),
            "resp_log": str(log),
        },
    }
