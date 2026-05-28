"""
mol2_utils
==========

Low-level helpers shared by every module in the pipeline.

Includes:
    * MOL2 parsing (atoms + bonds).
    * Atom-name and residue-name normalization.
    * Net-charge classification (standard residue table -> sidechain heuristic).
    * Backbone atom-type correction (CA -> CX, etc.) in MOL2 / AC / PREPIN.
    * Cap atom-name normalization (PyMOL "AC*"/"NM*" -> canonical names).
    * Total-charge renormalization to an integer.
    * Main-chain inference (shortest path between head and tail).

This module has NO side effects on import. Every helper is a pure function
that takes file paths or in-memory records as arguments.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Dict, Iterable, List, Optional, Set, Tuple


# ----------------------------------------------------------------------------
# Data classes
# ----------------------------------------------------------------------------
@dataclass(frozen=True)
class AtomRec:
    atom_id: int
    name: str
    atype: str
    charge: float
    subst_name: str = ""


@dataclass(frozen=True)
class BondRec:
    a1: int
    a2: int
    btype: str


# ----------------------------------------------------------------------------
# Residue-name normalization
# ----------------------------------------------------------------------------
def normalize_resname(raw: str) -> str:
    """Strip whitespace, uppercase, drop dotted/underscored suffixes."""
    s = str(raw or "").strip().upper()
    if not s:
        return "UNK"

    s = s.split(".")[0]
    s = s.split("_")[0]
    s = re.sub(r"[^A-Z0-9]", "", s)
    if not s:
        return "UNK"

    if len(s) > 4 and s[:4].isalpha() and s[4:].isdigit():
        return s[:4]
    if len(s) > 3 and s[:3].isalpha() and s[3:].isdigit():
        return s[:3]
    if len(s) > 3 and s[:2].isalpha() and s[2].isalnum() and s[3:].isdigit():
        return s[:3]

    if s.endswith(tuple("0123456789")):
        m = re.match(r"^(.*?)(\d+)$", s)
        if m:
            prefix = m.group(1)
            if len(prefix) in (3, 4):
                return prefix

    return s


# ----------------------------------------------------------------------------
# Element / atomic-number tables
# ----------------------------------------------------------------------------
_ELEMENT_BY_ATYPE_PREFIX = {
    "CL": "Cl", "BR": "Br", "NA": "Na", "MG": "Mg", "ZN": "Zn",
    "FE": "Fe", "CA": "Ca", "CU": "Cu", "MN": "Mn",
    "P": "P", "S": "S", "O": "O", "N": "N", "C": "C",
    "H": "H", "F": "F", "I": "I",
}

_ATOMIC_NUMBERS = {
    "H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
    "NA": 11, "MG": 12, "P": 15, "S": 16, "CL": 17,
    "K": 19, "CA": 20, "MN": 25, "FE": 26, "CU": 29,
    "ZN": 30, "BR": 35, "I": 53,
}

_BOND_ORDER_MAP = {
    "1": 1.0, "2": 2.0, "3": 3.0,
    "AM": 1.0, "AR": 1.5, "DU": 1.0, "UN": 1.0, "NC": 1.0,
}


STANDARD_RESIDUE_CHARGES: Dict[str, int] = {
    "ALA": 0, "ARG": 1, "ASN": 0, "ASP": -1, "CYS": 0,
    "GLN": 0, "GLU": -1, "GLY": 0, "HIS": 0, "HID": 0,
    "HIE": 0, "HIP": 1, "ILE": 0, "LEU": 0, "LYS": 1,
    "MET": 0, "PHE": 0, "PRO": 0, "SER": 0, "THR": 0,
    "TRP": 0, "TYR": 0, "VAL": 0,
}

CAP_RESNAMES = {"ACE", "NME"}


# ----------------------------------------------------------------------------
# Element guessing
# ----------------------------------------------------------------------------
def _guess_element(atom_name: str, atom_type: str) -> str:
    at = str(atom_type or "").strip()
    nm = str(atom_name or "").strip()

    # Tripos MOL2 atom types use ``<element>.<subtype>``: ``N.am``, ``N.pl3``,
    # ``C.3``, ``O.co2``, ``Cl``, ``Na``, etc. The substring *before* the
    # first dot is the element symbol verbatim — we look that up directly
    # rather than doing longest-prefix matching, which would otherwise treat
    # ``N.am`` (amide nitrogen) as ``Na`` (sodium) because ``NAM`` starts
    # with ``NA``.
    if at:
        element_part = at.split(".", 1)[0].strip().upper()
        if element_part in _ELEMENT_BY_ATYPE_PREFIX:
            return _ELEMENT_BY_ATYPE_PREFIX[element_part]
        # Falls through to name-based heuristic if the atom type prefix
        # doesn't directly match a known element symbol.
        if element_part and element_part[:1] in _ELEMENT_BY_ATYPE_PREFIX:
            return _ELEMENT_BY_ATYPE_PREFIX[element_part[:1]]

    name_alpha = re.sub(r"[^A-Za-z]", "", nm)
    if len(name_alpha) >= 2 and name_alpha[:2].upper() in {
        "CL", "BR", "NA", "MG", "ZN", "FE", "CA", "CU", "MN"
    }:
        return _ELEMENT_BY_ATYPE_PREFIX[name_alpha[:2].upper()]

    if name_alpha:
        return name_alpha[0].upper()

    return "C"


def _bond_order_value(bond_type: str) -> float:
    return _BOND_ORDER_MAP.get(str(bond_type or "1").strip().upper(), 1.0)


def get_mass(atom_type: str) -> float:
    masses = {
        "H": 1.008, "HO": 1.008, "HC": 1.008, "H1": 1.008, "H2": 1.008,
        "H3": 1.008, "HN": 1.008,
        "C": 12.011, "CA": 12.011, "CT": 12.011, "CM": 12.011, "C2": 12.011,
        "C3": 12.011, "CX": 12.011,
        "N": 14.007, "NA": 14.007, "N2": 14.007, "N3": 14.007, "N4": 14.007,
        "NB": 14.007,
        "O": 15.999, "O2": 15.999, "OH": 15.999, "OS": 15.999,
        "P": 30.974,
        "S": 32.06, "SH": 32.06,
        "F": 18.998, "CL": 35.45, "BR": 79.904, "I": 126.90,
    }
    return masses.get(str(atom_type).upper(), 12.011)


# ----------------------------------------------------------------------------
# Low-level MOL2 parsers
# ----------------------------------------------------------------------------
def _parse_mol2_atoms(mol2_path: str) -> List[AtomRec]:
    atoms: List[AtomRec] = []
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
            subst_name = parts[7] if len(parts) >= 8 else ""
            atoms.append(AtomRec(int(parts[0]), parts[1], parts[5], charge, subst_name))
    if not atoms:
        raise ValueError(f"No atoms parsed from MOL2: {mol2_path}")
    return atoms


def _parse_mol2_bonds_with_types(mol2_path: str) -> List[BondRec]:
    bonds: List[BondRec] = []
    in_bonds = False
    with open(mol2_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if s.startswith("@<TRIPOS>BOND"):
                in_bonds = True
                continue
            if s.startswith("@<TRIPOS>") and in_bonds:
                break
            if not in_bonds or not s:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            bonds.append(BondRec(int(parts[1]), int(parts[2]), parts[3]))
    return bonds


def _parse_mol2_bonds(mol2_path: str) -> List[Tuple[int, int]]:
    return [(b.a1, b.a2) for b in _parse_mol2_bonds_with_types(mol2_path)]


def extract_charges(mol2_path: str) -> List[float]:
    atoms = sorted(_parse_mol2_atoms(mol2_path), key=lambda a: a.atom_id)
    return [a.charge for a in atoms]


def get_atomtypes(mol2_path: str):
    atoms = sorted(_parse_mol2_atoms(mol2_path), key=lambda a: a.atom_id)
    atom_dicts = [{"id": a.atom_id, "name": a.name, "type": a.atype} for a in atoms]
    unique_types = sorted({a.atype for a in atoms})
    atomtypes_block = [f"; {t} (parameters from forcefield)" for t in unique_types]
    return atom_dicts, atomtypes_block


# ----------------------------------------------------------------------------
# Graph helpers
# ----------------------------------------------------------------------------
def _build_graph(num_atoms: int, bonds: List[Tuple[int, int]]) -> Dict[int, Set[int]]:
    g: Dict[int, Set[int]] = {i: set() for i in range(1, num_atoms + 1)}
    for a, b in bonds:
        if a in g and b in g:
            g[a].add(b)
            g[b].add(a)
    return g


def _build_neighbor_bonds(bonds: Iterable[BondRec]) -> Dict[int, List[Tuple[int, BondRec]]]:
    out: Dict[int, List[Tuple[int, BondRec]]] = {}
    for b in bonds:
        out.setdefault(b.a1, []).append((b.a2, b))
        out.setdefault(b.a2, []).append((b.a1, b))
    return out


def _shortest_path_ids(graph: Dict[int, Set[int]], start: int, goal: int) -> Optional[List[int]]:
    if start == goal:
        return [start]

    queue: List[int] = [start]
    seen = {start}
    prev: Dict[int, int] = {}
    while queue:
        cur = queue.pop(0)
        for nxt in graph.get(cur, set()):
            if nxt in seen:
                continue
            seen.add(nxt)
            prev[nxt] = cur
            if nxt == goal:
                out = [goal]
                walker = goal
                while walker != start:
                    walker = prev[walker]
                    out.append(walker)
                out.reverse()
                return out
            queue.append(nxt)
    return None


def _atomic_number(element: str) -> int:
    return _ATOMIC_NUMBERS.get(element.upper(), 0)


# ----------------------------------------------------------------------------
# Public main-chain inference
# ----------------------------------------------------------------------------
def infer_main_chain_from_mol2(
    mol2_path: str,
    head_name: Optional[str],
    tail_name: Optional[str],
) -> Optional[List[str]]:
    """
    Infer MAIN_CHAIN atom names as the shortest path between ``head_name`` and
    ``tail_name``, excluding the head/tail atoms themselves.

    Examples
    --------
    alpha amino acid  : N - CA - C        -> ["CA"]
    beta amino acid   : N - CA - CB - C   -> ["CA", "CB"]
    longer linkers    : N - ... - C       -> ["CA", "CB", "CG", ...]

    Returns ``None`` if either terminus is missing or no path of length >= 3
    exists between them (e.g. carbocyclic residues where the path can't be
    inferred from connectivity alone).
    """
    if not head_name or not tail_name:
        return None

    try:
        atoms = _parse_mol2_atoms(mol2_path)
        bonds = _parse_mol2_bonds(mol2_path)
    except Exception:
        return None

    by_id: Dict[int, AtomRec] = {a.atom_id: a for a in atoms}
    name_to_id: Dict[str, int] = {}
    for a in atoms:
        name_to_id.setdefault(a.name, a.atom_id)

    if head_name not in name_to_id or tail_name not in name_to_id:
        return None

    g = _build_graph(len(atoms), bonds)
    path = _shortest_path_ids(g, name_to_id[head_name], name_to_id[tail_name])
    if not path or len(path) < 3:
        return None

    names = [by_id[i].name for i in path if i in by_id]
    return [n for n in names if n != head_name and n != tail_name] or None


def has_atom_named(mol2_path: str, atom_name: str) -> bool:
    """Return True if any atom in the MOL2 file has the given name (case-sensitive)."""
    try:
        atoms = _parse_mol2_atoms(mol2_path)
    except Exception:
        return False
    return any(a.name == atom_name for a in atoms)


# ----------------------------------------------------------------------------
# Central-residue inference and net-charge classification
# ----------------------------------------------------------------------------
def _infer_central_residue_name(atoms: List[AtomRec]) -> str:
    non_caps = [normalize_resname(a.subst_name) for a in atoms
                if normalize_resname(a.subst_name) not in CAP_RESNAMES]
    if non_caps:
        counts: Dict[str, int] = {}
        for r in non_caps:
            counts[r] = counts.get(r, 0) + 1
        return max(counts.items(), key=lambda kv: kv[1])[0]
    if atoms:
        return normalize_resname(atoms[0].subst_name)
    return "UNK"


def _atom_type_upper(atom: AtomRec) -> str:
    return str(atom.atype or "").strip().upper()


def _element_upper(atom: AtomRec) -> str:
    return _guess_element(atom.name, atom.atype).upper()


def _oxygen_atom_type(atom: AtomRec) -> str:
    return _atom_type_upper(atom)


def _is_hydrogen(atom: AtomRec) -> bool:
    return _element_upper(atom) == "H"


def _neighbor_atoms(atom_id: int, by_id: Dict[int, AtomRec],
                    nbrs: Dict[int, List[Tuple[int, BondRec]]]) -> List[AtomRec]:
    return [by_id[nid] for nid, _ in nbrs.get(atom_id, []) if nid in by_id]


def _oxygen_bond_info(carbon_id: int, by_id: Dict[int, AtomRec],
                      nbrs: Dict[int, List[Tuple[int, BondRec]]]) -> List[Tuple[AtomRec, BondRec]]:
    out: List[Tuple[AtomRec, BondRec]] = []
    for nid, bond in nbrs.get(carbon_id, []):
        nbr = by_id[nid]
        if _element_upper(nbr) == "O":
            out.append((nbr, bond))
    return out


def _is_amide_n(n_atom_id: int, by_id: Dict[int, AtomRec],
                nbrs: Dict[int, List[Tuple[int, BondRec]]]) -> bool:
    for cid, _ in nbrs.get(n_atom_id, []):
        carbon = by_id[cid]
        if _element_upper(carbon) != "C":
            continue
        oxy = _oxygen_bond_info(carbon.atom_id, by_id, nbrs)
        if len(oxy) >= 1:
            return True
    return False


def _carboxyl_group_charge_from_carbon(
    carbon: AtomRec,
    by_id: Dict[int, AtomRec],
    nbrs: Dict[int, List[Tuple[int, BondRec]]],
) -> Optional[int]:
    if _element_upper(carbon) != "C":
        return None

    oxy = _oxygen_bond_info(carbon.atom_id, by_id, nbrs)
    if len(oxy) != 2:
        return None

    bond_orders = sorted(_bond_order_value(b.btype) for _, b in oxy)
    otypes = [_oxygen_atom_type(o) for o, _ in oxy]

    if all("CO2" in ot for ot in otypes):
        return -1

    if any(ot in {"O.3", "OH", "OS"} for ot in otypes):
        return 0

    if bond_orders == [1.0, 2.0]:
        return 0

    if bond_orders == [1.0, 1.0]:
        if any(ot in {"O.2", "O2"} for ot in otypes) and any(ot in {"O.3", "OH", "OS"} for ot in otypes):
            return 0
        if all(ot in {"O.2", "O2", "O.CO2", "OCO2"} or "CO2" in ot for ot in otypes):
            return -1

    return None


def _nitrogen_formal_charge(
    atom: AtomRec,
    by_id: Dict[int, AtomRec],
    nbrs: Dict[int, List[Tuple[int, BondRec]]],
) -> int:
    if _element_upper(atom) != "N":
        return 0

    at = _atom_type_upper(atom)
    neigh = _neighbor_atoms(atom.atom_id, by_id, nbrs)
    heavy = [n for n in neigh if not _is_hydrogen(n)]
    bond_order_sum = sum(_bond_order_value(b.btype) for _, b in nbrs.get(atom.atom_id, []))

    if _is_amide_n(atom.atom_id, by_id, nbrs):
        return 0

    if "N.4" in at or at == "N4":
        return 1

    if any(_element_upper(n) == "O" for n in heavy) and bond_order_sum < 4.0:
        return 0

    if len(neigh) >= 4 or bond_order_sum >= 4.0:
        return 1

    return 0


def _guanidinium_centers(
    candidate_ids: Iterable[int],
    by_id: Dict[int, AtomRec],
    nbrs: Dict[int, List[Tuple[int, BondRec]]],
) -> Set[int]:
    out: Set[int] = set()
    for atom_id in candidate_ids:
        atom = by_id[atom_id]
        if _element_upper(atom) != "C":
            continue
        n_neighbors = [by_id[nid] for nid, _ in nbrs.get(atom_id, [])
                       if _element_upper(by_id[nid]) == "N"]
        if len(n_neighbors) >= 3:
            out.add(atom_id)
    return out


def _infer_backbone_atom_ids(
    atoms: List[AtomRec],
    bonds: List[BondRec],
    central_resname: str,
) -> Set[int]:
    graph = _build_graph(len(atoms), [(b.a1, b.a2) for b in bonds])
    by_id = {a.atom_id: a for a in atoms}
    central_atoms = [a for a in atoms if normalize_resname(a.subst_name) == central_resname]
    central_ids = {a.atom_id for a in central_atoms}
    if not central_ids:
        return set()

    head_ids: List[int] = []
    tail_ids: List[int] = []
    name_to_ids: Dict[str, List[int]] = {}
    for a in central_atoms:
        name_to_ids.setdefault(a.name.upper(), []).append(a.atom_id)

    for b in bonds:
        a1 = by_id.get(b.a1)
        a2 = by_id.get(b.a2)
        if a1 is None or a2 is None:
            continue
        r1 = normalize_resname(a1.subst_name)
        r2 = normalize_resname(a2.subst_name)
        if r1 == central_resname and r2 == "ACE":
            head_ids.append(a1.atom_id)
        elif r2 == central_resname and r1 == "ACE":
            head_ids.append(a2.atom_id)
        if r1 == central_resname and r2 == "NME":
            tail_ids.append(a1.atom_id)
        elif r2 == central_resname and r1 == "NME":
            tail_ids.append(a2.atom_id)

    if not head_ids and "N" in name_to_ids:
        head_ids.extend(name_to_ids["N"])
    if not tail_ids and "C" in name_to_ids:
        tail_ids.extend(name_to_ids["C"])

    backbone_ids: Set[int] = set()
    if head_ids and tail_ids:
        path = _shortest_path_ids(graph, head_ids[0], tail_ids[0])
        if path:
            backbone_ids.update(i for i in path if i in central_ids)

    if not backbone_ids:
        for nm in ("N", "CA", "C"):
            for atom_id in name_to_ids.get(nm, []):
                if atom_id in central_ids:
                    backbone_ids.add(atom_id)

    for atom_id in list(backbone_ids):
        atom = by_id[atom_id]
        if atom.name.upper() != "C" and _element_upper(atom) != "C":
            continue
        for nbr_id in graph.get(atom_id, set()):
            if nbr_id not in central_ids:
                continue
            nbr = by_id[nbr_id]
            if _element_upper(nbr) == "O":
                backbone_ids.add(nbr_id)

    return backbone_ids


def _sidechain_charge_for_amino_acid_like(
    atoms: List[AtomRec],
    bonds: List[BondRec],
    central_resname: str,
) -> tuple[int, List[str]]:
    by_id = {a.atom_id: a for a in atoms}
    nbrs = _build_neighbor_bonds(bonds)
    central_atoms = [a for a in atoms if normalize_resname(a.subst_name) == central_resname]
    central_ids = {a.atom_id for a in central_atoms}
    backbone_ids = _infer_backbone_atom_ids(atoms, bonds, central_resname)
    sidechain_ids = {i for i in central_ids if i not in backbone_ids}

    charge = 0
    reasons: List[str] = []

    guan_centers = _guanidinium_centers(sidechain_ids, by_id, nbrs)
    if guan_centers:
        charge += len(guan_centers)
        reasons.extend(f"sidechain_guanidinium@C{cid}" for cid in sorted(guan_centers))

    guan_nitrogen_ids: Set[int] = set()
    for cid in guan_centers:
        for nid, _ in nbrs.get(cid, []):
            if nid in sidechain_ids and _element_upper(by_id[nid]) == "N":
                guan_nitrogen_ids.add(nid)

    for atom_id in sorted(sidechain_ids):
        atom = by_id[atom_id]
        if _element_upper(atom) == "C":
            group_charge = _carboxyl_group_charge_from_carbon(atom, by_id, nbrs)
            if group_charge == -1:
                charge -= 1
                reasons.append(f"sidechain_carboxylate@C{atom_id}")

    for atom_id in sorted(sidechain_ids):
        if atom_id in guan_nitrogen_ids:
            continue
        atom = by_id[atom_id]
        if _element_upper(atom) != "N":
            continue
        q = _nitrogen_formal_charge(atom, by_id, nbrs)
        if q > 0:
            charge += q
            reasons.append(f"sidechain_ammonium@N{atom_id}")

    return charge, reasons


def estimate_net_charge(mol2_path: str) -> tuple[Optional[int], str]:
    atoms = _parse_mol2_atoms(mol2_path)
    bonds = _parse_mol2_bonds_with_types(mol2_path)
    by_id = {a.atom_id: a for a in atoms}
    nbrs = _build_neighbor_bonds(bonds)

    charge = 0
    reasons: List[str] = []

    guan_centers = _guanidinium_centers((a.atom_id for a in atoms), by_id, nbrs)
    if guan_centers:
        charge += len(guan_centers)
        reasons.extend(f"guanidinium@C{cid}" for cid in sorted(guan_centers))

    guan_nitrogen_ids: Set[int] = set()
    for cid in guan_centers:
        for nid, _ in nbrs.get(cid, []):
            if _element_upper(by_id[nid]) == "N":
                guan_nitrogen_ids.add(nid)

    for atom in atoms:
        if _element_upper(atom) != "C":
            continue
        gc = _carboxyl_group_charge_from_carbon(atom, by_id, nbrs)
        if gc == -1:
            charge -= 1
            reasons.append(f"carboxylate@C{atom.atom_id}")

    for atom in atoms:
        if atom.atom_id in guan_nitrogen_ids:
            continue
        q = _nitrogen_formal_charge(atom, by_id, nbrs)
        if q > 0:
            charge += q
            reasons.append(f"ammonium@N{atom.atom_id}")

    if not reasons:
        return 0, "heuristic_neutral_no_strong_ionic_motifs"

    return charge, f"heuristic({';'.join(reasons)})"


def detect_formal_charge_from_mol2(mol2_path: str) -> Optional[int]:
    charge, _ = estimate_net_charge(mol2_path)
    return charge


def classify_residue_net_charge(mol2_path: str, resname: Optional[str] = None) -> tuple[Optional[int], str]:
    norm_resname = normalize_resname(resname or "")
    if norm_resname in STANDARD_RESIDUE_CHARGES:
        return STANDARD_RESIDUE_CHARGES[norm_resname], "known_residue_table"

    atoms = _parse_mol2_atoms(mol2_path)
    bonds = _parse_mol2_bonds_with_types(mol2_path)
    central_resname = _infer_central_residue_name(atoms)

    if central_resname in STANDARD_RESIDUE_CHARGES:
        return STANDARD_RESIDUE_CHARGES[central_resname], "known_residue_table_from_mol2"

    backbone_ids = _infer_backbone_atom_ids(atoms, bonds, central_resname)
    central_atoms = [a for a in atoms if normalize_resname(a.subst_name) == central_resname]
    amino_acid_like = bool(backbone_ids) or bool({a.name.upper() for a in central_atoms} & {"N", "CA", "C"})

    if amino_acid_like:
        charge, reasons = _sidechain_charge_for_amino_acid_like(atoms, bonds, central_resname)
        if reasons:
            return charge, f"amino_acid_sidechain({';'.join(reasons)})"
        return 0, "amino_acid_like_neutral_backbone"

    charge, reason = estimate_net_charge(mol2_path)
    return charge, reason


def classify_residue_net_charge_from_full_mol2(
    mol2_path: str,
    *,
    target_subst_id: Optional[int] = None,
    target_resname: Optional[str] = None,
) -> tuple[Optional[int], str]:
    atoms = _parse_mol2_atoms(mol2_path)
    bonds = _parse_mol2_bonds_with_types(mol2_path)

    if target_subst_id is not None:
        selected_atoms = [a for a in atoms if str(a.subst_name).strip()
                          and normalize_resname(a.subst_name) == normalize_resname(target_resname or a.subst_name)]
        exact = []
        with open(mol2_path, "r", encoding="utf-8", errors="replace") as f:
            in_atoms = False
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
                if int(parts[6]) == int(target_subst_id):
                    exact.append(AtomRec(
                        int(parts[0]), parts[1], parts[5],
                        float(parts[8]) if len(parts) >= 9 else 0.0,
                        parts[7],
                    ))
        selected_atoms = exact if exact else selected_atoms
    else:
        selected_atoms = [a for a in atoms if normalize_resname(a.subst_name) == normalize_resname(target_resname or "")]

    if not selected_atoms:
        return None, "target_residue_not_found"

    central_resname = normalize_resname(target_resname or selected_atoms[0].subst_name)
    central_ids = {a.atom_id for a in selected_atoms}
    relevant_atoms = [a for a in atoms if a.atom_id in central_ids]
    relevant_bonds = [b for b in bonds if b.a1 in central_ids and b.a2 in central_ids]

    backbone_ids = _infer_backbone_atom_ids(relevant_atoms, relevant_bonds, central_resname)
    amino_acid_like = bool(backbone_ids) or bool({a.name.upper() for a in selected_atoms} & {"N", "CA", "C"})

    if amino_acid_like:
        charge, reasons = _sidechain_charge_for_amino_acid_like(relevant_atoms, relevant_bonds, central_resname)
        if reasons:
            return charge, f"full_context_amino_acid_sidechain({';'.join(reasons)})"
        return 0, "full_context_amino_acid_neutral_backbone"

    charge, reason = estimate_net_charge_for_subgraph(relevant_atoms, relevant_bonds)
    return charge, f"full_context_{reason}"


def estimate_net_charge_for_subgraph(atoms: List[AtomRec], bonds: List[BondRec]) -> tuple[int, str]:
    by_id = {a.atom_id: a for a in atoms}
    nbrs = _build_neighbor_bonds(bonds)
    charge = 0
    reasons: List[str] = []

    guan_centers = _guanidinium_centers((a.atom_id for a in atoms), by_id, nbrs)
    if guan_centers:
        charge += len(guan_centers)
        reasons.extend(f"guanidinium@C{cid}" for cid in sorted(guan_centers))

    guan_nitrogen_ids: Set[int] = set()
    for cid in guan_centers:
        for nid, _ in nbrs.get(cid, []):
            if _element_upper(by_id[nid]) == "N":
                guan_nitrogen_ids.add(nid)

    for atom in atoms:
        if _element_upper(atom) != "C":
            continue
        gc = _carboxyl_group_charge_from_carbon(atom, by_id, nbrs)
        if gc == -1:
            charge -= 1
            reasons.append(f"carboxylate@C{atom.atom_id}")

    for atom in atoms:
        if atom.atom_id in guan_nitrogen_ids:
            continue
        q = _nitrogen_formal_charge(atom, by_id, nbrs)
        if q > 0:
            charge += q
            reasons.append(f"ammonium@N{atom.atom_id}")

    if not reasons:
        return 0, "heuristic_neutral_no_strong_ionic_motifs"
    return charge, f"heuristic({';'.join(reasons)})"


def electron_count_for_charge(mol2_path: str, charge: int) -> int:
    atoms = _parse_mol2_atoms(mol2_path)
    total_atomic_number = sum(_atomic_number(_guess_element(a.name, a.atype)) for a in atoms)
    return total_atomic_number - int(charge)


def adjust_charge_for_even_electrons(mol2_path: str, proposed_charge: int) -> tuple[int, bool]:
    """
    If ``proposed_charge`` gives an odd electron count, search for the
    *nearest* charge that yields an even (closed-shell) electron count.

    A radical amino-acid residue is chemically improbable, and sqm/AM1-BCC
    will fail outright on odd electrons. When the peptide splitter or the
    full-context classifier is off by one (e.g. because it misjudged the
    protonation state of an OXT-bearing C terminus after PyMOL re-added
    hydrogens), we'd rather correct it here than crash antechamber.

    Returns
    -------
    (corrected_charge, was_adjusted)
        ``was_adjusted`` is True only when the proposed charge was modified.
    """
    proposed = int(proposed_charge)
    if electron_count_for_charge(mol2_path, proposed) % 2 == 0:
        return proposed, False

    # Try ±1 first (most common case: protonation off by one), then ±2.
    # Prefer moving toward 0 when tied: i.e. for proposed = -1, try 0 before -2.
    for delta in (+1, -1, +2, -2):
        candidate = proposed + delta
        if electron_count_for_charge(mol2_path, candidate) % 2 == 0:
            return candidate, True

    # Nothing in ±2 worked (unlikely): hand back the original.
    return proposed, False


# ----------------------------------------------------------------------------
# Basic-N protonation (Arg / Lys / N-terminal amine, etc.)
# ----------------------------------------------------------------------------
def _parse_mol2_atom_coords(mol2_path: str) -> Dict[int, Tuple[float, float, float]]:
    """Return ``{atom_id: (x, y, z)}`` parsed from the ATOM block."""
    out: Dict[int, Tuple[float, float, float]] = {}
    in_atoms = False
    with open(mol2_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if s.startswith("@<TRIPOS>"):
                in_atoms = False
                continue
            if not in_atoms or not s:
                continue
            parts = s.split()
            if len(parts) < 6:
                continue
            try:
                atom_id = int(parts[0])
                x = float(parts[2])
                y = float(parts[3])
                z = float(parts[4])
                out[atom_id] = (x, y, z)
            except (ValueError, IndexError):
                continue
    return out


def _new_h_position(
    target_xyz: Tuple[float, float, float],
    neighbor_xyzs: List[Tuple[float, float, float]],
    bond_length: float = 1.01,
) -> Tuple[float, float, float]:
    """
    Place an H atom roughly opposite the centroid of existing neighbours.

    Geometry doesn't need to be perfect — sqm / RESP will optimize. We just
    need a non-overlapping position so the bond table is consistent.
    """
    tx, ty, tz = target_xyz
    if not neighbor_xyzs:
        return tx + bond_length, ty, tz

    cx = sum(p[0] for p in neighbor_xyzs) / len(neighbor_xyzs)
    cy = sum(p[1] for p in neighbor_xyzs) / len(neighbor_xyzs)
    cz = sum(p[2] for p in neighbor_xyzs) / len(neighbor_xyzs)
    dx, dy, dz = tx - cx, ty - cy, tz - cz
    norm = (dx * dx + dy * dy + dz * dz) ** 0.5
    if norm < 1e-6:
        # Centroid coincides with target — fall back to +X offset.
        return tx + bond_length, ty, tz
    dx, dy, dz = dx / norm, dy / norm, dz / norm
    return tx + bond_length * dx, ty + bond_length * dy, tz + bond_length * dz


def _select_basic_n_to_protonate(
    atoms: List[AtomRec],
    bonds: List[BondRec],
) -> Optional[AtomRec]:
    """
    Pick the most chemically-reasonable basic N to protonate, applying the
    priority order:
        0. Terminal guanidinium N - bonded only to the guanidinium C, no
           other heavy-atom neighbour (Arg NH1 / NH2). This is where the
           +1 form differs from neutral guanidine; protonating here yields
           the symmetric -C(NH2)2+ ion that AMBER parameterises.
        1. Primary amine - exactly one heavy-atom neighbour (a C), 1-2 Hs
           (Lys NZ, N-terminal -NH2 -> -NH3+).
        2. Other basic N with room for an H (secondary amine).
        4. Interior guanidinium N - bonded to two heavy atoms (Arg NE).
           Same H count in neutral and +1 forms, so protonating here
           creates a non-canonical -NH2- that parmchk2 cannot type.

    Amide N's (N bonded to a C that has a C=O) are skipped - they should
    not be protonated to NH2+.
    """
    by_id = {a.atom_id: a for a in atoms}
    nbrs = _build_neighbor_bonds(bonds)

    guanidinium_C_ids = _guanidinium_centers(
        (a.atom_id for a in atoms), by_id, nbrs,
    )

    candidates: List[Tuple[int, int, AtomRec]] = []
    # (priority, current_h_count, atom)

    for atom in atoms:
        if _element_upper(atom) != "N":
            continue
        if _is_amide_n(atom.atom_id, by_id, nbrs):
            continue

        neighbor_pairs = nbrs.get(atom.atom_id, [])
        if len(neighbor_pairs) >= 4:
            continue

        neighbor_ids = [nid for nid, _ in neighbor_pairs]
        h_count = sum(1 for nid in neighbor_ids if _is_hydrogen(by_id[nid]))
        heavy_neighbors = [nid for nid in neighbor_ids if not _is_hydrogen(by_id[nid])]
        c_neighbors = [nid for nid in heavy_neighbors if _element_upper(by_id[nid]) == "C"]
        if not c_neighbors or h_count >= 3:
            continue

        is_guanidinium = any(cid in guanidinium_C_ids for cid in c_neighbors)
        is_terminal = len(heavy_neighbors) == 1

        if is_guanidinium and is_terminal:
            priority = 0           # Arg NH1 / NH2 - canonical protonation site
        elif is_terminal and len(c_neighbors) == 1:
            priority = 1           # Lys NZ, N-terminal -NH2
        elif is_guanidinium:
            priority = 4           # Arg NE - H count doesn't change for +1
        else:
            priority = 2           # secondary amine

        candidates.append((priority, h_count, atom))

    if not candidates:
        return None

    candidates.sort(key=lambda t: (t[0], t[1], t[2].atom_id))
    return candidates[0][2]


def _format_mol2_atom_line(
    atom_id: int,
    name: str,
    x: float,
    y: float,
    z: float,
    atype: str,
    subst_id: int,
    subst_name: str,
    charge: float,
) -> str:
    """Emit a MOL2 ATOM line in a tolerant, space-separated layout."""
    return (
        f"{atom_id:>7d} {name:<8s}"
        f"{x:>10.4f}{y:>10.4f}{z:>10.4f} "
        f"{atype:<6s} {subst_id:>4d} {subst_name:<6s} {charge:>10.4f}\n"
    )


def _format_mol2_bond_line(bond_id: int, a1: int, a2: int, btype: str) -> str:
    return f"{bond_id:>6d}{a1:>6d}{a2:>6d} {btype}\n"


def _splice_atom_and_bond_into_mol2(
    mol2_path: str,
    new_atom_line: str,
    new_bond_line: str,
    new_atom_count: int,
    new_bond_count: int,
) -> None:
    """
    Insert ``new_atom_line`` at the end of the @<TRIPOS>ATOM block and
    ``new_bond_line`` at the end of the @<TRIPOS>BOND block, and update the
    counts in the @<TRIPOS>MOLECULE block. The file is rewritten in place.
    """
    with open(mol2_path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    # Locate section start indices.
    section_at: Dict[str, int] = {}
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith("@<TRIPOS>"):
            name = s[len("@<TRIPOS>"):].strip().upper()
            if name not in section_at:
                section_at[name] = i

    mol_start = section_at.get("MOLECULE")
    atom_start = section_at.get("ATOM")
    bond_start = section_at.get("BOND")
    if mol_start is None or atom_start is None or bond_start is None:
        raise ValueError(f"MOL2 file is missing required sections: {mol2_path}")

    # Update counts: typically the 2nd non-comment line of MOLECULE.
    counts_line_idx = None
    for j in range(mol_start + 1, min(mol_start + 6, len(lines))):
        parts = lines[j].split()
        if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
            counts_line_idx = j
            break
    if counts_line_idx is None:
        raise ValueError(f"Could not find counts line in {mol2_path}")

    counts_parts = lines[counts_line_idx].split()
    counts_parts[0] = str(new_atom_count)
    counts_parts[1] = str(new_bond_count)
    lines[counts_line_idx] = " ".join(counts_parts) + "\n"

    # Find end of ATOM block (start of the next section after ATOM).
    next_after_atom = min(
        (s for s in section_at.values() if s > atom_start),
        default=len(lines),
    )
    # Find end of BOND block (start of the next section after BOND).
    next_after_bond = min(
        (s for s in section_at.values() if s > bond_start),
        default=len(lines),
    )

    # Insert in descending order so the lower index doesn't shift.
    if next_after_bond >= next_after_atom:
        lines.insert(next_after_bond, new_bond_line)
        lines.insert(next_after_atom, new_atom_line)
    else:
        lines.insert(next_after_atom, new_atom_line)
        lines.insert(next_after_bond, new_bond_line)

    with open(mol2_path, "w", encoding="utf-8") as fh:
        fh.writelines(lines)


def add_proton_to_basic_n_for_charge(
    mol2_path: str,
    target_charge: int,
) -> dict:
    """
    If ``mol2_path`` has an odd electron count at ``target_charge``, find a
    basic N atom (guanidinium / primary amine / secondary amine; never an
    amide N) and append an H atom + bond to the file to restore the
    protonated state.

    This is the Arg / Lys / N-terminal complement to the carboxylate
    deprotonation pass done inside capping.py. PyMOL's ``h_add`` does not
    respect MOL2 formal charges, so a -NH3+ (Lys) or -C(NH2)2+ (Arg
    guanidinium) can come out of capping as the neutral -NH2 / =NH-NH2
    form, and the only way to recover the chemically correct charge state
    is to put the missing proton back.

    Returns a dict:
        {
            "hs_added": int,
            "target_n_name": Optional[str],
            "parity_ok_after": bool,
            "pre_electrons": int,
            "post_electrons": int,
        }
    """
    pre = electron_count_for_charge(mol2_path, target_charge)
    info: dict = {
        "hs_added": 0,
        "target_n_name": None,
        "parity_ok_after": (pre % 2 == 0),
        "pre_electrons": pre,
        "post_electrons": pre,
    }
    if pre % 2 == 0:
        return info

    atoms = _parse_mol2_atoms(mol2_path)
    bonds = _parse_mol2_bonds_with_types(mol2_path)
    target = _select_basic_n_to_protonate(atoms, bonds)
    if target is None:
        return info

    coords = _parse_mol2_atom_coords(mol2_path)
    by_id = {a.atom_id: a for a in atoms}
    nbrs = _build_neighbor_bonds(bonds)
    neighbor_xyzs = [
        coords[nid] for nid, _ in nbrs.get(target.atom_id, []) if nid in coords
    ]
    target_xyz = coords.get(target.atom_id, (0.0, 0.0, 0.0))
    hx, hy, hz = _new_h_position(target_xyz, neighbor_xyzs)

    # Generate a unique H name following AMBER conventions. For an N named
    # ``NH2`` whose siblings include ``HH21`` already, pick ``HH22`` so the
    # naming matches the standard library: parmchk2 doesn't care about names
    # (only types), but matching the canonical scheme keeps the MOL2 easy
    # to inspect by hand.
    used_names = {a.name.strip() for a in atoms}
    base_suffix = target.name.strip().lstrip("N")
    h_name: Optional[str] = None
    if base_suffix:
        for digit in range(1, 10):
            candidate = f"H{base_suffix}{digit}"
            if len(candidate) <= 4 and candidate not in used_names:
                h_name = candidate
                break
        if h_name is None:
            # All numbered candidates taken; fall back to the bare "H<suffix>".
            bare = f"H{base_suffix}"[:4]
            if bare not in used_names:
                h_name = bare
    if h_name is None:
        h_name = "H"
        counter = 1
        while h_name in used_names and counter < 10000:
            h_name = f"H{counter}"[:4]
            counter += 1

    new_atom_id = max(a.atom_id for a in atoms) + 1
    new_atom_line = _format_mol2_atom_line(
        atom_id=new_atom_id,
        name=h_name,
        x=hx, y=hy, z=hz,
        atype="H",
        subst_id=1,  # MOL2 lib substructure index; processor renormalises later
        subst_name=(target.subst_name or "RES")[:6],
        charge=0.0,
    )
    new_bond_id = (len(bonds) + 1) if bonds else 1
    new_bond_line = _format_mol2_bond_line(
        bond_id=new_bond_id, a1=target.atom_id, a2=new_atom_id, btype="1"
    )

    _splice_atom_and_bond_into_mol2(
        mol2_path,
        new_atom_line=new_atom_line,
        new_bond_line=new_bond_line,
        new_atom_count=len(atoms) + 1,
        new_bond_count=len(bonds) + 1,
    )

    info["hs_added"] = 1
    info["target_n_name"] = target.name
    info["post_electrons"] = electron_count_for_charge(mol2_path, target_charge)
    info["parity_ok_after"] = (info["post_electrons"] % 2 == 0)
    return info


# ----------------------------------------------------------------------------
# Validation
# ----------------------------------------------------------------------------
def validate_molecule(
    mol2_path: str,
    expected_charge: Optional[int] = None,
    check_parity: bool = True,
) -> List[str]:
    """
    Run basic structural checks on a MOL2 file.

    Parameters
    ----------
    mol2_path:
        Path to the MOL2 file to validate.
    expected_charge:
        If given (and ``check_parity`` is True), warn when the implied
        electron count is odd for this charge.
    check_parity:
        Set to False for pre-cap residue MOL2s in peptide mode, where
        dangling-bond extraction makes the parity check uninformative
        (an internal residue is intrinsically open-shell until capped).
    """
    warnings: List[str] = []
    atoms = _parse_mol2_atoms(mol2_path)
    bonds = _parse_mol2_bonds_with_types(mol2_path)
    if not atoms:
        raise ValueError(f"No atoms found in {mol2_path}")

    atom_names = [a.name.strip() for a in atoms]
    if len(atom_names) != len(set(atom_names)):
        warnings.append("Duplicate atom names detected in MOL2.")

    long_names = [a.name for a in atoms if len(a.name.strip()) > 4]
    if long_names:
        warnings.append(f"Amber atom-name limit exceeded for: {', '.join(sorted(set(long_names)))}")

    graph = _build_graph(len(atoms), [(b.a1, b.a2) for b in bonds])
    seen: Set[int] = set()
    stack = [atoms[0].atom_id]
    while stack:
        cur = stack.pop()
        if cur in seen:
            continue
        seen.add(cur)
        stack.extend(n for n in graph.get(cur, set()) if n not in seen)
    if len(seen) != len(atoms):
        warnings.append("Molecule graph is disconnected.")

    if expected_charge is not None and check_parity:
        total_atomic_number = sum(_atomic_number(_guess_element(a.name, a.atype)) for a in atoms)
        electrons = total_atomic_number - int(expected_charge)
        if electrons % 2 != 0:
            warnings.append(
                f"Expected charge {expected_charge} gives an odd electron count "
                f"({electrons}); sqm/AM1-BCC may fail."
            )

    return warnings


# ----------------------------------------------------------------------------
# Charge renormalization to an integer total
# ----------------------------------------------------------------------------
def renormalize_mol2_partial_charges_to_integer(mol2_path: str, target_integer_charge: int) -> bool:
    charges = extract_charges(mol2_path)
    current_sum = float(sum(charges))
    diff = float(target_integer_charge) - current_sum

    if abs(diff) <= 1e-6:
        return False

    per_atom = diff / float(len(charges))
    lines = Path(mol2_path).read_text(encoding="utf-8", errors="replace").splitlines()
    new_lines: List[str] = []

    in_atoms = False
    atom_seen = 0
    for line in lines:
        s = line.strip()
        if s.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            new_lines.append(line)
            continue
        if s.startswith("@<TRIPOS>") and in_atoms and not s.startswith("@<TRIPOS>ATOM"):
            in_atoms = False
            new_lines.append(line)
            continue

        if in_atoms and s:
            parts = line.split()
            if len(parts) >= 9:
                atom_seen += 1
                q = float(parts[-1]) + per_atom
                parts[-1] = f"{q:.6f}"
                new_lines.append(" ".join(parts))
                continue

        new_lines.append(line)

    if atom_seen != len(charges):
        raise RuntimeError(
            f"Charge renorm mismatch: saw {atom_seen} atom lines but parsed {len(charges)} charges."
        )

    Path(mol2_path).write_text("\n".join(new_lines) + "\n", encoding="utf-8")
    return True


# ----------------------------------------------------------------------------
# Bond / angle / dihedral helpers (for ITP-like exports)
# ----------------------------------------------------------------------------
def get_bonds_angles_dihedrals(mol2_path: str):
    atoms = sorted(_parse_mol2_atoms(mol2_path), key=lambda a: a.atom_id)
    bonds = _parse_mol2_bonds(mol2_path)
    g = _build_graph(len(atoms), bonds)

    bond_lines = [f"{a:>5} {b:>5} 1" for a, b in sorted({tuple(sorted(x)) for x in bonds})]

    angles_set: Set[Tuple[int, int, int]] = set()
    for j in range(1, len(atoms) + 1):
        neigh = sorted(g[j])
        for idx_i in range(len(neigh)):
            for idx_k in range(idx_i + 1, len(neigh)):
                i = neigh[idx_i]
                k = neigh[idx_k]
                angles_set.add((i, j, k) if i < k else (k, j, i))
    angle_lines = [f"{i:>5} {j:>5} {k:>5} 1" for (i, j, k) in sorted(angles_set)]

    dihed_set: Set[Tuple[int, int, int, int]] = set()
    for j, k in sorted({tuple(sorted(b)) for b in bonds}):
        for i in g[j]:
            if i == k:
                continue
            for l in g[k]:
                if l == j or l == i:
                    continue
                tup = (i, j, k, l)
                rev = (l, k, j, i)
                dihed_set.add(tup if tup < rev else rev)

    dihed_lines = [f"{i:>5} {j:>5} {k:>5} {l:>5} 1" for (i, j, k, l) in sorted(dihed_set)]
    return bond_lines, angle_lines, dihed_lines


def get_impropers(mol2_path: str):
    return []


def get_pairs(dihedrals: List[str]) -> List[str]:
    pairs: Set[Tuple[int, int]] = set()
    for d in dihedrals:
        parts = d.split()
        if len(parts) < 4:
            continue
        i = int(parts[0])
        l = int(parts[3])
        a, b = (i, l) if i < l else (l, i)
        pairs.add((a, b))
    return [f"{a:>5} {b:>5} 1" for (a, b) in sorted(pairs)]


# ----------------------------------------------------------------------------
# Backbone atom-type fixes (CA -> CX, etc.)
# ----------------------------------------------------------------------------
def fix_backbone_atom_types(mol2_path: str) -> bool:
    lines = Path(mol2_path).read_text(encoding="utf-8", errors="replace").splitlines()
    new_lines: List[str] = []
    in_atoms = False
    changed = False

    for line in lines:
        s = line.strip()

        if s.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            new_lines.append(line)
            continue

        if s.startswith("@<TRIPOS>") and in_atoms:
            in_atoms = False
            new_lines.append(line)
            continue

        if in_atoms and s:
            parts = line.split()
            if len(parts) >= 6:
                atom_name = parts[1].strip().upper()
                atom_type = parts[5]
                new_type = atom_type

                if atom_name == "CA":
                    new_type = "CX"
                elif atom_name == "N":
                    new_type = "N"
                elif atom_name == "C":
                    new_type = "C"
                elif atom_name == "O":
                    new_type = "O"

                if new_type != atom_type:
                    parts[5] = new_type
                    changed = True

                new_lines.append(" ".join(parts))
                continue

        new_lines.append(line)

    if changed:
        Path(mol2_path).write_text("\n".join(new_lines) + "\n", encoding="utf-8")
    return changed


def fix_backbone_atom_types_in_ac(ac_path: str) -> bool:
    lines = Path(ac_path).read_text(encoding="utf-8", errors="replace").splitlines()
    new_lines: List[str] = []
    changed = False

    atom_line_re = re.compile(
        r"^(ATOM\s+\d+\s+)(\S+)(\s+\S+\s+\d+\s+[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+\s+)(\S+)\s*$"
    )

    for line in lines:
        m = atom_line_re.match(line)
        if not m:
            new_lines.append(line)
            continue

        prefix, atom_name, middle, atom_type = m.groups()
        atom_name_u = atom_name.upper().strip()
        new_type = atom_type

        if atom_name_u == "CA":
            new_type = "CX"
        elif atom_name_u == "N":
            new_type = "N"
        elif atom_name_u == "C":
            new_type = "C"
        elif atom_name_u == "O":
            new_type = "O"

        if new_type != atom_type:
            changed = True

        new_lines.append(f"{prefix}{atom_name}{middle}{new_type}")

    if changed:
        Path(ac_path).write_text("\n".join(new_lines) + "\n", encoding="utf-8")
    return changed


def fix_backbone_atom_types_in_prepin(prepin_path: str) -> bool:
    lines = Path(prepin_path).read_text(encoding="utf-8", errors="replace").splitlines()
    new_lines: List[str] = []
    changed = False

    for line in lines:
        stripped = line.strip()
        if not stripped:
            new_lines.append(line)
            continue

        parts = line.split()
        if len(parts) < 3 or not parts[0].isdigit():
            new_lines.append(line)
            continue

        atom_name = parts[1].upper()
        old_type = parts[2]
        new_type = old_type

        if atom_name == "CA":
            new_type = "CX"
        elif atom_name == "N":
            new_type = "N"
        elif atom_name == "C":
            new_type = "C"
        elif atom_name == "O":
            new_type = "O"

        if new_type != old_type:
            parts[2] = new_type
            changed = True
            new_lines.append(" ".join(parts))
        else:
            new_lines.append(line)

    if changed:
        Path(prepin_path).write_text("\n".join(new_lines) + "\n", encoding="utf-8")
    return changed


def normalize_cap_atom_names_in_mol2(mol2_path: Path) -> None:
    """
    Convert PyMOL cap naming (AC*, NM*) -> canonical Amber-friendly naming.

    Only affects ACE/NME residues. Kept for backward compatibility.
    """
    mol2_path = Path(mol2_path)
    lines = mol2_path.read_text(encoding="utf-8", errors="replace").splitlines()
    new_lines = []

    in_atoms = False

    ace_map = {
        "AC1": "CC1", "AC2": "OC2", "AC3": "CC3",
        "AC4": "HC4", "AC5": "HC5", "AC6": "HC6",
    }
    nme_map = {
        "NM1": "NM1", "NM2": "CM2", "NM3": "HM3",
        "NM4": "HM4", "NM5": "HM5", "NM6": "HM6",
    }

    for line in lines:
        s = line.strip()

        if s.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            new_lines.append(line)
            continue

        if s.startswith("@<TRIPOS>") and in_atoms:
            in_atoms = False
            new_lines.append(line)
            continue

        if in_atoms and s:
            parts = line.split()
            if len(parts) >= 8:
                name = parts[1]
                resn = parts[7].upper()

                if resn == "ACE" and name in ace_map:
                    parts[1] = ace_map[name]
                elif resn == "NME" and name in nme_map:
                    parts[1] = nme_map[name]

                line = " ".join(parts)

        new_lines.append(line)

    mol2_path.write_text("\n".join(new_lines) + "\n", encoding="utf-8")
