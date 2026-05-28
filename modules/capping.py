"""
capping
=======

Stand-alone PyMOL script that attaches ACE / NME caps to a single residue.

Invoked as a subprocess (so it runs inside the PyMOL Python interpreter)::

    python capping.py <residue_folder> [--head N|NONE] [--tail C|NONE]
                                       [--preserve-oxt]
                                       [--net-charge <int>]

Important behaviour for terminal residues (used by peptide mode and the
``-t/--terminal`` flag in main.py):

    * If ``--preserve-oxt`` is given, the OXT atom is NOT stripped from the
      C terminus and NME capping is skipped automatically (the residue is
      already a chemically complete carboxylate). This is required for
      C-terminal residues of linear peptides where OXT participates in
      hydrogen bonding.

    * If ``--head NONE`` is given, ACE capping is skipped.
    * If ``--tail NONE`` is given, NME capping is skipped.

Charge-aware protonation
------------------------
If ``--net-charge <int>`` is provided, after PyMOL's ``h_add`` step the
capped residue is audited: PyMOL's ``h_add`` does not respect MOL2 formal
charges, so a carboxylate O can silently pick up an H and turn a COO-
side chain into a neutral COOH. When the resulting electron count is odd
for the declared net charge, this script finds carboxylate -OH groups
(an O bonded to an H AND to a C that also has another O neighbour) and
strips an H from one of them, restoring the deprotonated state. This is
exactly the behaviour you want for Asp / Glu (-1) at physiological pH.

The script writes:
    * ``residue_capped.mol2`` - the capped (or uncapped) MOL2.
    * ``residue_capping_meta.json`` - small metadata blob that the residue
      processor later merges into the consolidated ``residue_meta.json``.

Exit codes
----------
    0  success
    1  bad arguments / missing input
    2  PyMOL not importable
    3  capping ran but no output file produced
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

try:
    import pymol  # type: ignore
except Exception:
    pymol = None


MAX_ATOM_NAME_LEN = 4

# Element -> atomic number table. Two-letter elements must be checked first.
_TWO_LETTER_ELEMENTS = ("CL", "BR", "NA", "MG", "FE", "ZN", "CU", "MN", "CA")
_ATOMIC_Z = {
    "H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
    "NA": 11, "MG": 12, "P": 15, "S": 16, "CL": 17,
    "K": 19, "CA": 20, "MN": 25, "FE": 26, "CU": 29,
    "ZN": 30, "BR": 35, "I": 53,
}


# ----------------------------------------------------------------------------
# Atom-name renaming helpers
# ----------------------------------------------------------------------------
def _all_atom_names(selection: str = "all") -> set[str]:
    atoms = pymol.cmd.get_model(selection).atom
    return {str(atom.name).strip() for atom in atoms}


def _make_unique_name(base: str, used: set[str]) -> str:
    base = (base or "X")[:MAX_ATOM_NAME_LEN]
    if base not in used:
        return base

    for i in range(1, 10000):
        suffix = str(i)
        prefix = base[: MAX_ATOM_NAME_LEN - len(suffix)]
        candidate = f"{prefix}{suffix}"
        if candidate not in used:
            return candidate

    raise RuntimeError(f"Could not generate a unique atom name from base '{base}'")


def _safe_unique_rename(selection: str, prefix: str) -> None:
    atoms = pymol.cmd.get_model(selection).atom
    if not atoms:
        return

    used = _all_atom_names("all")

    for i, atom in enumerate(atoms, start=1):
        old_name = str(atom.name).strip()
        used.discard(old_name)
        newname = _make_unique_name(f"{prefix}{i}", used)
        used.add(newname)
        pymol.cmd.alter(f"{selection} and id {atom.id}", f"name='{newname}'")


def _normalize_name_arg(value):
    if value is None:
        return None
    s = str(value).strip()
    if not s or s.upper() in {"NONE", "NULL", "0", ""}:
        return None
    return s


# ----------------------------------------------------------------------------
# Charge-aware protonation audit
# ----------------------------------------------------------------------------
def _element_of(atom) -> str:
    """Best-effort element symbol from a ChemPy atom."""
    s = (getattr(atom, "symbol", "") or "").strip().upper()
    if not s:
        n = (getattr(atom, "name", "") or "").strip().upper()
        s = n[:2] if n else ""
    for two in _TWO_LETTER_ELEMENTS:
        if s.startswith(two):
            return two
    return s[:1] if s else ""


def _electron_count(target_charge: int) -> int:
    """Sum of atomic numbers minus target charge."""
    model = pymol.cmd.get_model("all")
    z_sum = sum(_ATOMIC_Z.get(_element_of(a), 0) for a in model.atom)
    return z_sum - int(target_charge)


def _find_carboxylate_OH_index():
    """
    Find an O atom in a carboxylic-acid -OH context:
    - bonded to >= 1 H
    - bonded to a C that also has another O neighbour (the carbonyl partner)

    Returns the *PyMOL ``index``* of the H atom to remove, or None if no
    candidate is found.
    """
    model = pymol.cmd.get_model("all")

    # Build bond adjacency in 0-based model-list indices.
    adj: dict[int, list[int]] = {}
    for bond in model.bond:
        a, b = bond.index[0], bond.index[1]
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)

    atoms = list(model.atom)

    for o_idx, o_atom in enumerate(atoms):
        if _element_of(o_atom) != "O":
            continue

        h_neighbors = [
            n for n in adj.get(o_idx, [])
            if _element_of(atoms[n]) == "H"
        ]
        if not h_neighbors:
            continue

        c_neighbors = [
            n for n in adj.get(o_idx, [])
            if _element_of(atoms[n]) == "C"
        ]
        for c_idx in c_neighbors:
            other_os = [
                n for n in adj.get(c_idx, [])
                if _element_of(atoms[n]) == "O" and n != o_idx
            ]
            if other_os:
                # Carboxylate -OH found; pick its first H to remove.
                return atoms[h_neighbors[0]].index

    return None


def _fix_protonation_to_match_charge(
    target_charge: int,
    max_iters: int = 6,
) -> tuple[int, bool]:
    """
    Strip Hs from carboxylate -OH groups one at a time until the electron
    count is closed-shell at ``target_charge`` (or no candidates remain).

    Returns (hs_removed, parity_ok).
    """
    removed = 0
    for _ in range(max_iters):
        if _electron_count(target_charge) % 2 == 0:
            return removed, True
        h_index = _find_carboxylate_OH_index()
        if h_index is None:
            return removed, False
        pymol.cmd.remove(f"index {h_index}")
        removed += 1

    return removed, _electron_count(target_charge) % 2 == 0


# ----------------------------------------------------------------------------
# Argument parsing
# ----------------------------------------------------------------------------
def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="capping.py",
        description="Attach ACE / NME caps to a residue MOL2 using PyMOL.",
    )
    p.add_argument("residue_folder", help="Folder containing residue.mol2")
    p.add_argument(
        "--head",
        default="N",
        help="Head atom name (default 'N'). Use 'NONE' to skip ACE capping.",
    )
    p.add_argument(
        "--tail",
        default="C",
        help="Tail atom name (default 'C'). Use 'NONE' to skip NME capping.",
    )
    p.add_argument(
        "--preserve-oxt",
        action="store_true",
        help="Keep the OXT atom intact and skip C-terminal NME capping. "
             "Use this for C-terminal residues of linear peptides.",
    )
    p.add_argument(
        "--net-charge",
        type=int,
        default=None,
        help="Target net charge. When given, the capped residue is audited "
             "after h_add: if the electron count is odd, carboxylate -OH "
             "groups are deprotonated until parity matches. This restores "
             "the deprotonated COO- form for Asp/Glu and other carboxylate "
             "side chains that PyMOL's h_add may silently neutralize.",
    )
    return p


def _parse_cli(argv: list[str]) -> argparse.Namespace:
    """
    Accept both the new (--head/--tail) flag style AND the historical
    positional style ``<folder> <head> <tail>`` so we don't break any
    existing call sites or shell scripts.
    """
    parser = _build_argparser()

    if len(argv) >= 2 and not argv[1].startswith("-"):
        if len(argv) >= 3 and not argv[2].startswith("-"):
            folder = argv[1]
            head = argv[2]
            tail = argv[3] if len(argv) >= 4 and not argv[3].startswith("-") else "C"
            extras = []
            start = 4 if len(argv) >= 4 and not argv[3].startswith("-") else 3
            for token in argv[start:]:
                extras.append(token)
            return parser.parse_args([folder, "--head", head, "--tail", tail, *extras])

    return parser.parse_args(argv[1:])


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main() -> int:
    if pymol is None:
        print("ERROR: PyMOL is not available in this Python environment. "
              "Install/enable PyMOL to use capping.")
        return 2

    pymol.finish_launching(["pymol", "-cq"])

    if len(sys.argv) < 2:
        _build_argparser().print_help()
        return 1

    args = _parse_cli(sys.argv)

    residue_folder = Path(args.residue_folder).resolve()
    head_name = _normalize_name_arg(args.head)
    tail_name = _normalize_name_arg(args.tail)
    preserve_oxt = bool(args.preserve_oxt)
    target_charge = args.net_charge

    input_file = residue_folder / "residue.mol2"
    output_file = residue_folder / "residue_capped.mol2"
    meta_file = residue_folder / "residue_capping_meta.json"

    if not input_file.exists():
        print(f"ERROR: {input_file} not found!")
        return 1

    pymol.cmd.reinitialize()
    pymol.cmd.load(str(input_file), "prot")

    # Hydrogens are always reset; PyMOL re-adds them after capping.
    pymol.cmd.remove("hydro")

    oxt_atom_count = pymol.cmd.count_atoms("name OXT")
    oxt_preserved = preserve_oxt and oxt_atom_count > 0

    if not preserve_oxt:
        pymol.cmd.remove("name OXT")

    if oxt_preserved and tail_name:
        tail_name = None

    applied_caps: list[str] = []
    found_head = False
    found_tail = False

    if head_name:
        head_sel = f"name {head_name}"
        pymol.cmd.select("pk1", head_sel)
        found_head = pymol.cmd.count_atoms("pk1") > 0
        if found_head:
            pymol.cmd.editor.attach_amino_acid("pk1", "ace")
            _safe_unique_rename("resn ACE", "AC")
            applied_caps.append("ACE")

    if tail_name:
        tail_sel = f"name {tail_name} and not resn ACE"
        pymol.cmd.select("pk1", tail_sel)
        found_tail = pymol.cmd.count_atoms("pk1") > 0
        if found_tail:
            pymol.cmd.editor.attach_amino_acid("pk1", "nme")
            _safe_unique_rename("resn NME", "NM")
            applied_caps.append("NME")

    pymol.cmd.h_add("prot")

    # ----- charge-aware carboxylate deprotonation -----
    # PyMOL's h_add does not respect MOL2 formal charges, so a carboxylate O
    # can silently pick up an H. If the resulting electron parity is wrong
    # for the declared net charge, find carboxylate -OH groups and strip Hs
    # one at a time until parity matches. The complementary case
    # (protonating a basic N like Arg/Lys) is handled post-cap in
    # mol2_utils.add_proton_to_basic_n_for_charge.
    protonation_fix = {
        "target_charge": target_charge,
        "hs_removed": 0,
        "parity_ok": None,
    }
    if target_charge is not None:
        pre_electrons = _electron_count(target_charge)
        hs_removed, parity_ok = _fix_protonation_to_match_charge(target_charge)
        post_electrons = _electron_count(target_charge)
        protonation_fix["hs_removed"] = hs_removed
        protonation_fix["parity_ok"] = parity_ok
        protonation_fix["pre_electrons"] = pre_electrons
        protonation_fix["post_electrons"] = post_electrons
        if hs_removed > 0:
            print(
                f"Deprotonated {hs_removed} carboxylate -OH "
                f"(electrons {pre_electrons} -> {post_electrons} at charge {target_charge})"
            )

    pymol.cmd.save(str(output_file), "prot")

    meta = {
        "requested_head_name": head_name,
        "requested_tail_name": tail_name,
        "preserve_oxt_requested": preserve_oxt,
        "oxt_preserved": oxt_preserved,
        "has_head": bool(found_head),
        "has_tail": bool(found_tail),
        "applied_caps": applied_caps,
        "protonation_fix": protonation_fix,
    }
    meta_file.write_text(json.dumps(meta, indent=2), encoding="utf-8")

    if output_file.exists() and output_file.stat().st_size > 0:
        print(f"Capping successful: {output_file}")
        print(json.dumps(meta))
        return 0

    print(f"ERROR: PyMOL failed to create {output_file}!")
    return 3


if __name__ == "__main__":
    raise SystemExit(main())

