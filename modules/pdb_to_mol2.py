"""
pdb_to_mol2
===========

Stand-alone PyMOL script that converts a PDB file to a MOL2 file.

Invoked as a subprocess (so it runs inside the PyMOL Python interpreter):

    python pdb_to_mol2.py <input.pdb> <output.mol2>

Important
---------
This conversion intentionally preserves the PDB protonation state.  Hydrogen
removal/re-addition is deferred to ``capping.py``, where the requested net
charge is already known and the pipeline can correct PyMOL protonation mistakes.

Exit codes
----------
    0  success
    1  bad CLI arguments or missing input
    2  PyMOL not importable in the current interpreter
    3  PyMOL load/save failed or produced an invalid/empty MOL2
"""
from __future__ import annotations

import sys
from pathlib import Path

try:
    import pymol  # type: ignore
except Exception:
    pymol = None


def main() -> int:
    if pymol is None:
        print("ERROR: PyMOL is not available in this Python environment. "
              "Install/enable PyMOL to use pdb_to_mol2.")
        return 2

    pymol.finish_launching(["pymol", "-cq"])

    if len(sys.argv) < 3:
        print("Usage: python pdb_to_mol2.py <input_pdb_path> <output_mol2_path>")
        return 1

    input_pdb = Path(sys.argv[1]).resolve()
    output_mol2 = Path(sys.argv[2]).resolve()

    if not input_pdb.exists():
        print(f"ERROR: {input_pdb} not found!")
        return 1

    output_mol2.parent.mkdir(parents=True, exist_ok=True)

    try:
        pymol.cmd.reinitialize()
        pymol.cmd.load(str(input_pdb), "prot")

        atom_count = int(pymol.cmd.count_atoms("prot"))
        if atom_count <= 0:
            print(f"ERROR: PyMOL loaded no atoms from {input_pdb}")
            return 3

        # Do NOT remove/re-add hydrogens here.  At this point the requested
        # residue net charge has not been resolved yet, and PyMOL h_add can
        # change protonation states (especially carboxylates).  capping.py
        # performs the hydrogen reset later, after the target charge is known.
        pymol.cmd.save(str(output_mol2), "prot")
    except Exception as exc:
        print(f"ERROR: PyMOL PDB->MOL2 conversion failed: {exc}")
        return 3

    if not output_mol2.exists() or output_mol2.stat().st_size == 0:
        print(f"ERROR: PyMOL did not create a non-empty MOL2: {output_mol2}")
        return 3

    try:
        text = output_mol2.read_text(encoding="utf-8", errors="replace")
    except Exception as exc:
        print(f"ERROR: Could not read generated MOL2 {output_mol2}: {exc}")
        return 3


    if "@<TRIPOS>ATOM" not in text or "@<TRIPOS>BOND" not in text:
        print(
            "ERROR: Generated file is not a complete Tripos MOL2 "
            f"(missing ATOM/BOND section): {output_mol2}"
        )
        return 3

    atom_lines = []
    in_atoms = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if stripped.startswith("@<TRIPOS>") and in_atoms:
            break
        if in_atoms and stripped:
            atom_lines.append(line)

    if not atom_lines:
        print(f"ERROR: Generated MOL2 contains no atom records: {output_mol2}")
        return 3

    print(
        f"Conversion successful: {output_mol2} "
        f"({len(atom_lines)} atoms; protonation preserved from PDB)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
