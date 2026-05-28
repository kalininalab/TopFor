"""
modules.peptide_assembler
=========================

Build full-peptide AMBER (and optionally GROMACS) topology/coordinate
files AFTER the per-residue pipeline has produced .prepin + .frcmod for
every non-standard residue.

Pipeline:
    1. Convert the original peptide MOL2 -> PDB (PyMOL subprocess).
    2. Build a tleap script that:
         * sources leaprc.protein.<backbone> and leaprc.<sidechain>,
         * loadamberparams every NSAA frcmod,
         * loadamberprep every NSAA .prepin    <-- intentional, NOT loadoff
                                                  on the .lib (those carry
                                                  ACE/NME caps embedded in
                                                  the residue unit and
                                                  cannot be polymerised).
         * loadpdb the converted peptide PDB,
         * saveamberparm -> peptide.prmtop / peptide.inpcrd.
    3. If generate_gmx: convert prmtop/inpcrd -> top/gro via ParmEd.

Caveats
-------
* The peptide PDB residue names must match the resnames used during
  single-residue parametrization. Since the splitter reads resnames
  straight from the input MOL2, this normally just works.
* tleap is sensitive to atom-name conventions. If a peptide atom name
  differs from the prep template, inspect the assembly log.
* This module assumes the input peptide does NOT already carry ACE/NME
  caps; tleap will leave both termini as charged amine / carboxylate
  (AMBER's standard treatment). If your peptide has explicit caps,
  rename them to ACE / NME and they will be picked up from the
  force field.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Iterable


def _run(cmd: list[str], cwd: Path, log_file: Path) -> bool:
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(cwd))
    with log_file.open("a", encoding="utf-8") as f:
        f.write("\n\n=== CMD ===\n")
        f.write(" ".join(cmd))
        f.write("\n=== STDOUT ===\n")
        f.write(result.stdout or "")
        f.write("\n=== STDERR ===\n")
        f.write(result.stderr or "")
        f.write(f"\n=== RETURN CODE === {result.returncode}\n")
    return result.returncode == 0


def _mol2_to_pdb(mol2_path: Path, pdb_path: Path, log_file: Path) -> bool:
    script = Path(__file__).parent / "mol2_to_pdb.py"
    if not script.exists():
        raise FileNotFoundError(f"PyMOL converter script missing: {script}")
    return _run(
        ["python", str(script), str(mol2_path), str(pdb_path)],
        cwd=pdb_path.parent,
        log_file=log_file,
    )


def _collect_residue_artifacts(
    out_base: Path,
    resnames: Iterable[str],
    backbone: str,
    sidechain: str,
) -> dict:
    """For each resname, locate its .prepin and both .frcmod files."""
    prepins: list[Path] = []
    frcmods: list[Path] = []
    missing: list[str] = []

    for resname in resnames:
        residue_dir = out_base / resname
        if not residue_dir.is_dir():
            missing.append(resname)
            continue

        prepin = residue_dir / f"{resname}.prepin"
        bb_frcmod = residue_dir / f"{resname}_{backbone}.frcmod"
        sc_frcmod = residue_dir / f"{resname}_{sidechain}.frcmod"

        if not prepin.exists():
            missing.append(resname)
            continue

        prepins.append(prepin)
        if bb_frcmod.exists():
            frcmods.append(bb_frcmod)
        if sc_frcmod.exists():
            frcmods.append(sc_frcmod)

    return {"prepins": prepins, "frcmods": frcmods, "missing": missing}


def assemble_peptide(
    peptide_input: str,
    successful_resnames: Iterable[str],
    out_base: str,
    *,
    backbone: str = "ff19SB",
    sidechain: str = "gaff2",
    generate_gmx: bool = False,
) -> dict:
    """
    Generate prmtop / inpcrd (and optionally .top / .gro) for the
    full peptide.

    Outputs land in <out_base>/peptide/.

    Returns
    -------
    dict with keys: status ("ok"|"failed"), prmtop, inpcrd, top, gro,
    log, and reason (when failed).
    """
    peptide_path = Path(peptide_input).resolve()
    out_base_p = Path(out_base).resolve()
    asm_dir = out_base_p / "peptide"
    asm_dir.mkdir(parents=True, exist_ok=True)

    log_file = asm_dir / "peptide_assembly.log"
    log_file.write_text("", encoding="utf-8")

    failed_payload = {
        "status": "failed",
        "prmtop": None,
        "inpcrd": None,
        "top": None,
        "gro": None,
        "log": str(log_file),
    }

    # ---- 1) PDB feed for tleap ----
    suffix = peptide_path.suffix.lower()
    peptide_pdb = asm_dir / "peptide_input.pdb"

    if suffix == ".pdb":
        shutil.copyfile(peptide_path, peptide_pdb)
    elif suffix == ".mol2":
        if not _mol2_to_pdb(peptide_path, peptide_pdb, log_file):
            return {**failed_payload, "reason": "MOL2 -> PDB conversion failed"}
    else:
        return {**failed_payload,
                "reason": f"Unsupported peptide input suffix: {suffix}"}

    # ---- 2) Locate per-residue artifacts ----
    artifacts = _collect_residue_artifacts(
        out_base_p, successful_resnames, backbone, sidechain,
    )

    # ---- 3) Build & run tleap ----
    prmtop_out = asm_dir / "peptide.prmtop"
    inpcrd_out = asm_dir / "peptide.inpcrd"
    leap_script = asm_dir / "peptide_leap.in"

    lines: list[str] = []
    lines.append(f"source leaprc.protein.{backbone}")
    lines.append(f"source leaprc.{sidechain}")

    for frcmod in artifacts["frcmods"]:
        lines.append(f"loadamberparams {frcmod.resolve()}")
    for prepin in artifacts["prepins"]:
        lines.append(f"loadamberprep {prepin.resolve()}")

    lines.append(f"peptide = loadpdb {peptide_pdb.resolve()}")
    lines.append(f"saveamberparm peptide {prmtop_out.name} {inpcrd_out.name}")
    lines.append("quit")
    leap_script.write_text("\n".join(lines) + "\n", encoding="utf-8")

    if not _run(["tleap", "-f", leap_script.name], asm_dir, log_file):
        return {**failed_payload, "reason": "tleap failed during peptide assembly"}

    if not prmtop_out.exists() or not inpcrd_out.exists():
        return {**failed_payload,
                "reason": "tleap finished but prmtop/inpcrd were not produced"}

    result: dict = {
        "status": "ok",
        "prmtop": str(prmtop_out),
        "inpcrd": str(inpcrd_out),
        "top": None,
        "gro": None,
        "log": str(log_file),
        "missing_residues": artifacts["missing"],
    }

    # ---- 4) Optional GROMACS conversion ----
    if generate_gmx:
        try:
            from modules.gromacs_converter import amber_to_gromacs

            top_out = asm_dir / "peptide.top"
            gro_out = asm_dir / "peptide.gro"
            amber_to_gromacs(
                str(prmtop_out), str(inpcrd_out),
                str(top_out), str(gro_out),
            )
            result["top"] = str(top_out)
            result["gro"] = str(gro_out)
        except Exception as exc:
            with log_file.open("a", encoding="utf-8") as f:
                f.write(f"\n=== GROMACS CONVERSION FAILED ===\n{exc}\n")
            result["gmx_error"] = str(exc)

    return result