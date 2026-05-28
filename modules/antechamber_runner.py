"""
modules.antechamber_runner
==========================

Stage-2 of the NSAA workflow. For every charged single-residue MOL2 produced
by :mod:`modules.residue_processor`, run the canonical AMBER toolchain::

    antechamber     -> .ac     (atom typing under the requested force field)
    prepgen_writer  -> .mc     (main-chain control file for prepgen)
    prepgen         -> .prepin
    parmchk2 (x2)   -> .frcmod (one per force field family)
    tleap           -> .lib + .prmtop + .inpcrd

When ``generate_gmx=True`` the prmtop/inpcrd pair for each residue is
additionally converted to GROMACS ``.top`` and ``.gro`` via ParmEd
(see :mod:`modules.gromacs_converter`).

Returns
-------
dict
    ``{"successful": [resname, ...], "failed": [(resname, reason), ...]}``
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

from modules.prepgen_writer import write_prepgen_mc_file
from modules.mol2_utils import (
    classify_residue_net_charge,
    normalize_resname,
    fix_backbone_atom_types_in_ac,
    fix_backbone_atom_types_in_prepin,
)


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def _run(cmd, cwd: Path, log_file: Path) -> bool:
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(cwd))
    with log_file.open("a", encoding="utf-8") as f:
        f.write("\n\n=== CMD ===\n")
        f.write(" ".join(cmd) if isinstance(cmd, list) else str(cmd))
        f.write("\n=== STDOUT ===\n")
        f.write(result.stdout)
        f.write("\n=== STDERR ===\n")
        f.write(result.stderr)
        f.write("\n=== RETURN CODE ===\n")
        f.write(str(result.returncode) + "\n")
    return result.returncode == 0


def _read_json_if_exists(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _read_residue_meta(residue_dir: Path) -> dict:
    residue_meta = _read_json_if_exists(residue_dir / "residue_meta.json")
    if residue_meta:
        return residue_meta

    capping_meta = _read_json_if_exists(residue_dir / "residue_capping_meta.json")
    if capping_meta:
        return {
            "head_name": capping_meta.get("requested_head_name"),
            "tail_name": capping_meta.get("requested_tail_name"),
            "main_chain": None,
            "pre_head_type": "C",
            "post_tail_type": "N",
            "applied_caps": capping_meta.get("applied_caps", []),
        }
    return {}


def _read_net_charge_for_residue(residue_dir: Path, mol2_path: Path) -> Tuple[int, str]:
    data = _read_residue_meta(residue_dir)
    if "net_charge" in data:
        return int(data["net_charge"]), "meta"

    classified, source = classify_residue_net_charge(
        str(mol2_path),
        resname=normalize_resname(mol2_path.stem),
    )
    if classified is not None:
        return int(classified), source

    return 0, "fallback_0"


# ----------------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------------
def run_antechamber_for_all(
    mol2_files: List[str],
    backbone: str = "ff19SB",
    sidechain: str = "gaff2",
    charge: str = "abcg2",
    generate_gmx: bool = False,
) -> Dict[str, List]:
    """
    Run the AMBER toolchain for every charged residue MOL2.

    If ``generate_gmx`` is True, also produce GROMACS .top / .gro per
    residue via ParmEd.
    """
    successful: List[str] = []
    failed: List[Tuple[str, str]] = []

    amberhome = os.environ.get("AMBERHOME")
    if not amberhome:
        print("AMBERHOME not set. Cannot run AMBER toolchain.")
        for input_mol2 in mol2_files:
            resname = normalize_resname(Path(input_mol2).stem)
            failed.append((resname, "AMBERHOME not set"))
        return {"successful": successful, "failed": failed}

    backbone_map = {
        "ff19SB": "parm19.dat",
        "ff14SB": "parm10.dat",
        "ff99SB": "parm99.dat",
    }
    sidechain_map = {
        "gaff": "gaff.dat",
        "gaff2": "gaff2.dat",
    }

    if backbone not in backbone_map:
        raise ValueError(f"Unsupported backbone: {backbone}")
    if sidechain not in sidechain_map:
        raise ValueError(f"Unsupported sidechain: {sidechain}")

    backbone_parm_path = os.path.join(amberhome, "dat", "leap", "parm", backbone_map[backbone])
    gaff_parm_path = os.path.join(amberhome, "dat", "leap", "parm", sidechain_map[sidechain])

    # Lazy-import so the AMBER stage still runs if parmed is missing and
    # the user did not actually request --gmx.
    if generate_gmx:
        try:
            from modules.gromacs_converter import amber_to_gromacs  # noqa: F401
        except Exception as exc:
            print(f"[WARN] GROMACS conversion requested but ParmEd not "
                  f"available: {exc}. Skipping GROMACS output.")
            generate_gmx = False

    for input_mol2 in mol2_files:
        mol2_path = Path(input_mol2).resolve()
        residue_dir = mol2_path.parent
        resname = normalize_resname(mol2_path.stem)
        meta = _read_residue_meta(residue_dir)

        log_file = residue_dir / f"{resname}.log"
        log_file.write_text("", encoding="utf-8")

        ac_output = residue_dir / f"{resname}.ac"
        mc_file = residue_dir / f"{resname}.mc"
        prepin_file = residue_dir / f"{resname}.prepin"
        lib_file = residue_dir / f"{resname}.lib"
        prmtop_file = residue_dir / f"{resname}.prmtop"
        inpcrd_file = residue_dir / f"{resname}.inpcrd"

        backbone_frcmod_output = residue_dir / f"{resname}_{backbone}.frcmod"
        gaff_frcmod_output = residue_dir / f"{resname}_{sidechain}.frcmod"

        net_charge, _ = _read_net_charge_for_residue(residue_dir, mol2_path)
        print(f"[{resname}] using net charge = {net_charge}")

        # 1) antechamber: MOL2 -> AC
        antechamber_cmd = [
            "antechamber",
            "-i", mol2_path.name,
            "-fi", "mol2",
            "-o", ac_output.name,
            "-fo", "ac",
            "-at", "amber",
            "-nc", str(net_charge),
        ]
        if not _run(antechamber_cmd, residue_dir, log_file):
            print(f"[{resname}] FAILED at AC generation")
            failed.append((resname, "AC generation"))
            continue

        fix_backbone_atom_types_in_ac(str(ac_output))

        # 2) write the prepgen main-chain (.mc) control file
        applied_caps = tuple(str(x).upper() for x in meta.get("applied_caps", []))
        try:
            write_prepgen_mc_file(
                str(mol2_path),
                str(mc_file),
                head_name=meta.get("head_name", "N"),
                tail_name=meta.get("tail_name", "C"),
                main_chain=meta.get("main_chain"),
                charge=float(net_charge),
                central_resname=resname,
                cap_resnames=applied_caps if applied_caps else ("ACE", "NME"),
                pre_head_type=str(meta.get("pre_head_type", "C")),
                post_tail_type=str(meta.get("post_tail_type", "N")),
                infer_mainchain_from_connectivity=True,
            )
        except Exception as exc:
            print(f"[{resname}] FAILED at MC file generation: {exc}")
            failed.append((resname, f"MC file generation: {exc}"))
            continue

        # 3) prepgen: AC + MC -> PREPIN
        prepgen_cmd = [
            "prepgen",
            "-i", ac_output.name,
            "-o", prepin_file.name,
            "-m", mc_file.name,
            "-rn", resname,
        ]
        if not _run(prepgen_cmd, residue_dir, log_file):
            print(f"[{resname}] FAILED at prepgen")
            failed.append((resname, "prepgen"))
            continue

        if fix_backbone_atom_types_in_prepin(str(prepin_file)):
            print(f"[{resname}] corrected backbone atom types in PREPIN file")

        # 4) parmchk2 (backbone + sidechain frcmods)
        parmchk_backbone_cmd = [
            "parmchk2", "-i", ac_output.name, "-f", "ac",
            "-o", backbone_frcmod_output.name, "-a", "Y",
            "-p", backbone_parm_path,
        ]
        parmchk_sidechain_cmd = [
            "parmchk2", "-i", ac_output.name, "-f", "ac",
            "-o", gaff_frcmod_output.name, "-a", "Y",
            "-p", gaff_parm_path,
        ]

        if not _run(parmchk_backbone_cmd, residue_dir, log_file):
            print(f"[{resname}] FAILED at parmchk2 (backbone)")
            failed.append((resname, "parmchk2 backbone"))
            continue
        if not _run(parmchk_sidechain_cmd, residue_dir, log_file):
            print(f"[{resname}] FAILED at parmchk2 (sidechain)")
            failed.append((resname, "parmchk2 sidechain"))
            continue

        # 5) tleap -> .lib + .prmtop + .inpcrd
        # The unit loaded here is the capped+charged MOL2 (ACE-X-NME), so
        # the prmtop/inpcrd describe a capped tripeptide-style system that
        # the user can simulate directly as a single-residue smoke test.
        # Full peptide topologies are built separately by peptide_assembler.
        leap_script = residue_dir / "leap.in"
        leap_script.write_text(
            (
                f"source leaprc.protein.{backbone}\n"
                f"source leaprc.{sidechain}\n"
                f"loadamberparams {backbone_frcmod_output.name}\n"
                f"loadamberparams {gaff_frcmod_output.name}\n"
                f"{resname} = loadmol2 {mol2_path.name}\n"
                f"saveoff {resname} {lib_file.name}\n"
                f"saveamberparm {resname} {prmtop_file.name} {inpcrd_file.name}\n"
                f"quit\n"
            ),
            encoding="utf-8",
        )

        if not _run(["tleap", "-f", leap_script.name], residue_dir, log_file):
            print(f"[{resname}] FAILED at tleap")
            failed.append((resname, "tleap"))
            continue

        # 6) optional GROMACS conversion (per-residue capped system)
        if generate_gmx:
            try:
                from modules.gromacs_converter import amber_to_gromacs
                top_out = residue_dir / f"{resname}.top"
                gro_out = residue_dir / f"{resname}.gro"
                amber_to_gromacs(
                    str(prmtop_file), str(inpcrd_file),
                    str(top_out), str(gro_out),
                )
                print(f"[{resname}] GROMACS .top / .gro written")
            except Exception as exc:
                # GROMACS conversion failure does NOT fail the residue —
                # AMBER outputs are still valid. We just log and warn.
                with log_file.open("a", encoding="utf-8") as f:
                    f.write(f"\n=== GROMACS CONVERSION FAILED ===\n{exc}\n")
                print(f"[{resname}] WARN: GROMACS conversion failed: {exc}")

        print(f"\033[1m[{resname}] parametrization complete\033[0m")
        successful.append(resname)

    return {"successful": successful, "failed": failed}