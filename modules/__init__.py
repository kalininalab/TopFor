"""
nsaa_paramgen.modules
=====================

Internal modules for the NSAA (Non-Standard Amino Acid) parameter generation tool.

Modules
-------
mol2_utils
    Low-level MOL2 parsing, charge classification, atom-type fixing, and
    molecule validation helpers.

peptide_splitter
    Extract individual non-standard residues out of a full peptide MOL2 file
    and write per-residue MOL2 + split-metadata files.

capping
    Stand-alone PyMOL script that attaches ACE/NME caps to a single residue.
    Invoked as a subprocess so it runs inside the PyMOL Python interpreter.

pdb_to_mol2
    Stand-alone PyMOL script that converts a PDB file to a MOL2 file.

prepgen_writer
    Build a prepgen ``.mc`` file (HEAD_NAME / TAIL_NAME / MAIN_CHAIN /
    OMIT_NAME / PRE_HEAD_TYPE / POST_TAIL_TYPE / CHARGE) from a capped MOL2.

resp_workflow
    Full RESP workflow: xTB pre-optimization -> ORCA optimization + SP ->
    Multiwfn RESP fitting -> antechamber ``-c rc`` for charge assignment.

antechamber_runner
    Run antechamber / prepgen / parmchk2 / tleap to produce ``.ac``,
    ``.prepin``, ``.frcmod``, ``.lib`` files for each charged residue.

residue_processor
    High-level driver that ties everything together: prepares the input MOL2,
    caps it, assigns partial charges, validates the result, and writes a
    single consolidated ``residue_meta.json`` per residue.
"""

__version__ = "2.0.0"
__all__ = [
    "mol2_utils",
    "peptide_splitter",
    "prepgen_writer",
    "resp_workflow",
    "antechamber_runner",
    "residue_processor",
]
