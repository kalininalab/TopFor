"""
pdb_to_mol2
===========

Stand-alone PyMOL script that converts a PDB file to a MOL2 file.

Invoked as a subprocess (so it runs inside the PyMOL Python interpreter):

    python pdb_to_mol2.py <input.pdb> <output.mol2>

Exit codes
----------
    0  success
    1  bad CLI arguments or missing input
    2  PyMOL not importable in the current interpreter
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

    pymol.cmd.reinitialize()
    pymol.cmd.load(str(input_pdb), "prot")
    pymol.cmd.remove("hydro")
    pymol.cmd.h_add("prot")
    pymol.cmd.save(str(output_mol2), "prot")
    print(f"Conversion successful: {output_mol2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
