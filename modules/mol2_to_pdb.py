"""
mol2_to_pdb
===========

Stand-alone PyMOL subprocess script. Mirrors ``pdb_to_mol2.py`` but runs
in the opposite direction. Used by the peptide assembler so we can hand
the original peptide structure to tleap via ``loadpdb`` (more reliable
than ``loadmol2`` for multi-residue peptides).

Usage:
    python mol2_to_pdb.py <input.mol2> <output.pdb>

Exit codes:
    0  success
    1  bad arguments / missing input
    2  PyMOL not importable
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
        print("ERROR: PyMOL is not available in this Python environment.")
        return 2

    pymol.finish_launching(["pymol", "-cq"])

    if len(sys.argv) < 3:
        print("Usage: python mol2_to_pdb.py <input.mol2> <output.pdb>")
        return 1

    input_mol2 = Path(sys.argv[1]).resolve()
    output_pdb = Path(sys.argv[2]).resolve()

    if not input_mol2.exists():
        print(f"ERROR: {input_mol2} not found!")
        return 1

    pymol.cmd.reinitialize()
    pymol.cmd.load(str(input_mol2), "pep")
    pymol.cmd.save(str(output_pdb), "pep")
    print(f"Conversion successful: {output_pdb}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())