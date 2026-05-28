"""
peptide_splitter
================

Extract individual non-standard residues out of a full peptide MOL2 file.

For each non-standard residue encountered, the splitter:
    1. Selects all atoms with the matching ``subst_id`` (which preserves OXT
       atoms on C-terminal residues automatically, because they share that
       residue's substructure id).
    2. Writes a stand-alone ``RESNAME.mol2`` file using a remapped atom-id
       space so the file parses cleanly on its own.
    3. Writes ``RESNAME.split.json`` describing the residue's role in the
       polymer: HEAD / TAIL atom names, MAIN_CHAIN, external peptide bonds,
       polymer topology (internal / N-terminal / C-terminal / cyclic) and
       whether an OXT atom is present.

The ``topology`` and ``has_oxt`` fields are what the rest of the pipeline
uses to decide:
    * which caps to attach (ACE / NME / both / none),
    * whether to preserve OXT during capping.
"""
from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

from .mol2_utils import (
    classify_residue_net_charge_from_full_mol2,
    normalize_resname,
)

STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
}

HIS_VARIANTS = {"HID", "HIE", "HIP"}


@dataclass(frozen=True)
class AtomRec:
    atom_id: int
    name: str
    x: float
    y: float
    z: float
    atom_type: str
    subst_id: int
    subst_name: str
    charge: float


@dataclass(frozen=True)
class BondRec:
    bond_id: int
    a1: int
    a2: int
    bond_type: str


def _parse_sections(mol2_path: str) -> Dict[str, List[str]]:
    sections: Dict[str, List[str]] = {}
    current: Optional[str] = None

    with open(mol2_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("@<TRIPOS>"):
                current = line.strip()
                sections[current] = []
                continue
            if current:
                sections[current].append(line.rstrip("\n"))

    return sections


def _parse_atoms(lines: List[str]) -> List[AtomRec]:
    atoms: List[AtomRec] = []
    for line in lines:
        if not line.strip():
            continue
        p = line.split()
        if len(p) < 8:
            continue

        charge = float(p[8]) if len(p) > 8 else 0.0
        atoms.append(
            AtomRec(
                atom_id=int(p[0]),
                name=p[1],
                x=float(p[2]),
                y=float(p[3]),
                z=float(p[4]),
                atom_type=p[5],
                subst_id=int(p[6]),
                subst_name=p[7],
                charge=charge,
            )
        )
    return atoms


def _parse_bonds(lines: List[str]) -> List[BondRec]:
    bonds: List[BondRec] = []
    for line in lines:
        if not line.strip():
            continue
        p = line.split()
        if len(p) < 4:
            continue

        bonds.append(
            BondRec(
                bond_id=int(p[0]),
                a1=int(p[1]),
                a2=int(p[2]),
                bond_type=p[3],
            )
        )
    return bonds


def _guess_element(atom_name: str, atom_type: str) -> str:
    at = "".join(ch for ch in str(atom_type or "") if ch.isalpha()).upper()
    nm = "".join(ch for ch in str(atom_name or "") if ch.isalpha()).upper()

    for token in ("CL", "BR", "NA", "MG", "ZN", "FE", "CU", "MN", "CA"):
        if at.startswith(token):
            return token
        if nm.startswith(token):
            return token

    if at:
        return at[0]
    if nm:
        return nm[0]
    return "C"


def _build_graph(bonds: Iterable[BondRec]) -> Dict[int, Set[int]]:
    g: Dict[int, Set[int]] = {}
    for b in bonds:
        g.setdefault(b.a1, set()).add(b.a2)
        g.setdefault(b.a2, set()).add(b.a1)
    return g


def _shortest_path(graph: Dict[int, Set[int]], start: int, goal: int) -> Optional[List[int]]:
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


def _dedupe_keep_order(items: Iterable[str]) -> List[str]:
    out: List[str] = []
    seen: Set[str] = set()
    for x in items:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def _is_carbonyl_like(atom: AtomRec, central_ids: Set[int],
                     by_id: Dict[int, AtomRec], graph: Dict[int, Set[int]]) -> bool:
    if _guess_element(atom.name, atom.atom_type) != "C":
        return False

    for nbr_id in graph.get(atom.atom_id, set()):
        if nbr_id not in central_ids:
            continue
        nbr = by_id[nbr_id]
        if _guess_element(nbr.name, nbr.atom_type) == "O":
            return True

    return atom.name.upper() == "C"


def _score_head_candidate(atom: AtomRec, central_ids: Set[int],
                          by_id: Dict[int, AtomRec], graph: Dict[int, Set[int]]) -> int:
    score = 0
    if atom.name.upper() == "N":
        score += 10
    if _guess_element(atom.name, atom.atom_type) == "N":
        score += 5

    for nbr_id in graph.get(atom.atom_id, set()):
        if nbr_id not in central_ids:
            continue
        nbr = by_id[nbr_id]
        if _guess_element(nbr.name, nbr.atom_type) == "C":
            score += 2
            if nbr.name.upper() in {"CA", "CB", "C1", "C2", "C3"}:
                score += 1

    return score


def _score_tail_candidate(atom: AtomRec, central_ids: Set[int],
                          by_id: Dict[int, AtomRec], graph: Dict[int, Set[int]]) -> int:
    score = 0
    if atom.name.upper() == "C":
        score += 10
    if _guess_element(atom.name, atom.atom_type) == "C":
        score += 5
    if _is_carbonyl_like(atom, central_ids, by_id, graph):
        score += 6
    return score


def _infer_polymer_connection_atoms(
    residue_atoms: List[AtomRec],
    all_atoms_by_id: Dict[int, AtomRec],
    all_bonds: List[BondRec],
) -> tuple[Optional[str], Optional[str], List[dict]]:
    central_ids = {a.atom_id for a in residue_atoms}
    all_graph = _build_graph(all_bonds)
    by_id = all_atoms_by_id

    external_bonds: List[dict] = []
    head_candidates: List[tuple[int, str]] = []
    tail_candidates: List[tuple[int, str]] = []

    for bond in all_bonds:
        in_a1 = bond.a1 in central_ids
        in_a2 = bond.a2 in central_ids
        if in_a1 == in_a2:
            continue

        central_id = bond.a1 if in_a1 else bond.a2
        external_id = bond.a2 if in_a1 else bond.a1

        central_atom = by_id[central_id]
        external_atom = by_id[external_id]
        central_elem = _guess_element(central_atom.name, central_atom.atom_type)
        external_elem = _guess_element(external_atom.name, external_atom.atom_type)

        external_bonds.append(
            {
                "central_atom": central_atom.name,
                "central_subst": normalize_resname(central_atom.subst_name),
                "external_atom": external_atom.name,
                "external_subst": normalize_resname(external_atom.subst_name),
                "bond_type": bond.bond_type,
            }
        )

        if central_elem == "N" and external_elem == "C":
            head_candidates.append(
                (_score_head_candidate(central_atom, central_ids, by_id, all_graph), central_atom.name)
            )
        elif central_elem == "C" and external_elem == "N":
            tail_candidates.append(
                (_score_tail_candidate(central_atom, central_ids, by_id, all_graph), central_atom.name)
            )

    head_name = max(head_candidates, default=(None, None))[1]
    tail_name = max(tail_candidates, default=(None, None))[1]
    return head_name, tail_name, external_bonds


def _infer_main_chain(
    residue_atoms: List[AtomRec],
    all_bonds: List[BondRec],
    head_name: Optional[str],
    tail_name: Optional[str],
) -> Optional[List[str]]:
    if not head_name or not tail_name:
        return None

    central_ids = {a.atom_id for a in residue_atoms}
    central_bonds = [b for b in all_bonds if b.a1 in central_ids and b.a2 in central_ids]
    graph = _build_graph(central_bonds)
    by_name: Dict[str, List[int]] = {}
    by_id: Dict[int, AtomRec] = {a.atom_id: a for a in residue_atoms}

    for atom in residue_atoms:
        by_name.setdefault(atom.name, []).append(atom.atom_id)

    if head_name not in by_name or tail_name not in by_name:
        return None

    path = _shortest_path(graph, by_name[head_name][0], by_name[tail_name][0])
    if not path or len(path) < 3:
        return None

    names = [by_id[i].name for i in path]
    names = [x for x in names if x not in {head_name, tail_name}]
    names = _dedupe_keep_order(names)
    return names or None


def _has_oxt(residue_atoms: List[AtomRec]) -> bool:
    """Detect whether the residue carries an OXT atom (terminal carboxylate)."""
    return any(a.name.strip().upper() == "OXT" for a in residue_atoms)


def _has_n_terminal_protons(residue_atoms: List[AtomRec],
                            head_name: Optional[str],
                            all_atoms_by_id: Dict[int, AtomRec],
                            all_bonds: List[BondRec]) -> bool:
    """
    True when the head N atom carries 2 or 3 hydrogens — indicating a free
    N-terminus (NH2 or NH3+) rather than an amide nitrogen in a peptide bond.
    """
    if not head_name:
        return False

    target_n = next((a for a in residue_atoms if a.name == head_name), None)
    if target_n is None:
        return False

    graph = _build_graph(all_bonds)
    h_neighbors = 0
    for nbr_id in graph.get(target_n.atom_id, set()):
        nbr = all_atoms_by_id.get(nbr_id)
        if nbr is None:
            continue
        if _guess_element(nbr.name, nbr.atom_type) == "H":
            h_neighbors += 1

    return h_neighbors >= 2


def _write_residue_mol2(
    output_file: Path,
    resname: str,
    residue_atoms: List[AtomRec],
    residue_bonds: List[BondRec],
) -> None:
    atom_ids = {a.atom_id for a in residue_atoms}
    id_map = {old: new for new, old in enumerate(sorted(atom_ids), start=1)}

    with output_file.open("w", encoding="utf-8") as f:
        f.write("@<TRIPOS>MOLECULE\n")
        f.write(f"{resname}\n")
        f.write(f"{len(residue_atoms)} {len(residue_bonds)} 1\n")
        f.write("SMALL\n")
        f.write("USER_CHARGES\n\n")

        f.write("@<TRIPOS>ATOM\n")
        for atom in sorted(residue_atoms, key=lambda a: a.atom_id):
            f.write(
                f"{id_map[atom.atom_id]:>6} {atom.name:<6} "
                f"{atom.x:>10.4f} {atom.y:>10.4f} {atom.z:>10.4f} "
                f"{atom.atom_type:<6} 1 {resname:<6} {atom.charge:>10.6f}\n"
            )

        f.write("@<TRIPOS>BOND\n")
        for i, bond in enumerate(residue_bonds, start=1):
            f.write(f"{i:>6} {id_map[bond.a1]:>4} {id_map[bond.a2]:>4} {bond.bond_type}\n")

        f.write("@<TRIPOS>SUBSTRUCTURE\n")
        f.write(f"1 {resname} 1\n")


def _candidate_rank(meta: dict, residue_atoms: List[AtomRec],
                    residue_bonds: List[BondRec]) -> tuple[int, int, int, int]:
    has_head = 1 if meta.get("head_name") else 0
    has_tail = 1 if meta.get("tail_name") else 0
    has_both = 1 if has_head and has_tail else 0
    return (has_both, has_head + has_tail, len(residue_atoms), len(residue_bonds))


def extract_nonstandard_residues_from_mol2(input_mol2: str, output_dir: str) -> List[str]:
    sections = _parse_sections(input_mol2)
    atoms = _parse_atoms(sections.get("@<TRIPOS>ATOM", []))
    bonds = _parse_bonds(sections.get("@<TRIPOS>BOND", []))

    if not atoms:
        raise ValueError(f"No atoms parsed from MOL2 file: {input_mol2}")

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    std = STANDARD_RESIDUES | HIS_VARIANTS
    residues: Dict[int, List[AtomRec]] = {}
    atoms_by_id: Dict[int, AtomRec] = {a.atom_id: a for a in atoms}

    for atom in atoms:
        residues.setdefault(atom.subst_id, []).append(atom)

    selected_by_resname: Dict[str, dict] = {}

    for subst_id in sorted(residues):
        residue_atoms = residues[subst_id]
        raw_name = residue_atoms[0].subst_name
        resname = normalize_resname(raw_name)

        if resname in std:
            continue

        central_ids = {a.atom_id for a in residue_atoms}
        residue_bonds = [b for b in bonds if b.a1 in central_ids and b.a2 in central_ids]

        head_name, tail_name, external_bonds = _infer_polymer_connection_atoms(
            residue_atoms,
            atoms_by_id,
            bonds,
        )
        main_chain = _infer_main_chain(residue_atoms, bonds, head_name, tail_name)

        is_polymer_internal = bool(head_name and tail_name)
        is_n_terminal_like = bool(tail_name and not head_name)
        is_c_terminal_like = bool(head_name and not tail_name)
        is_cyclic_like = bool(head_name and tail_name and len(external_bonds) >= 2)

        topology = "isolated"
        if is_cyclic_like:
            topology = "cyclic_like"
        elif is_polymer_internal:
            topology = "internal"
        elif is_n_terminal_like:
            topology = "n_term_like"
        elif is_c_terminal_like:
            topology = "c_term_like"

        has_oxt_atom = _has_oxt(residue_atoms)
        has_n_term_protons = _has_n_terminal_protons(
            residue_atoms, head_name, atoms_by_id, bonds
        )

        full_charge, full_charge_source = classify_residue_net_charge_from_full_mol2(
            input_mol2,
            target_subst_id=subst_id,
            target_resname=resname,
        )

        # Capping plan: skip the corresponding cap on whichever terminus is
        # already chemically satisfied (free amine / OXT carboxylate).
        # For an isolated residue we leave both caps off as well.
        preserve_oxt = topology in {"c_term_like", "isolated"} and has_oxt_atom
        skip_head_cap = topology in {"n_term_like", "isolated"}
        skip_tail_cap = topology in {"c_term_like", "isolated"}

        # If we are preserving OXT (peptide C-terminal), suppress tail_name so
        # the capping step won't try to attach NME and won't strip OXT.
        effective_tail_name = None if skip_tail_cap else tail_name
        effective_head_name = None if skip_head_cap else head_name

        meta = {
            "resname": resname,
            "head_name": effective_head_name,
            "tail_name": effective_tail_name,
            "raw_head_name": head_name,
            "raw_tail_name": tail_name,
            "main_chain": main_chain,
            "pre_head_type": "C",
            "post_tail_type": "N",
            "source_input": str(Path(input_mol2).resolve()),
            "source_subst_id": residue_atoms[0].subst_id,
            "source_subst_name": residue_atoms[0].subst_name,
            "external_bonds": external_bonds,
            "is_polymer_internal": is_polymer_internal,
            "is_n_terminal_like": is_n_terminal_like,
            "is_c_terminal_like": is_c_terminal_like,
            "is_cyclic_like": is_cyclic_like,
            "topology": topology,
            "has_oxt": has_oxt_atom,
            "has_n_terminal_protons": has_n_term_protons,
            "preserve_oxt": preserve_oxt,
            "skip_head_cap": skip_head_cap,
            "skip_tail_cap": skip_tail_cap,
            "full_context_net_charge": int(full_charge) if full_charge is not None else None,
            "full_context_charge_source": full_charge_source,
        }

        rank = _candidate_rank(meta, residue_atoms, residue_bonds)

        candidate = {
            "rank": rank,
            "subst_id": int(subst_id),
            "residue_atoms": residue_atoms,
            "residue_bonds": residue_bonds,
            "meta": meta,
        }

        prev = selected_by_resname.get(resname)
        if prev is None:
            selected_by_resname[resname] = candidate
            continue

        prev_rank = prev["rank"]
        prev_subst_id = int(prev["subst_id"])

        if rank > prev_rank:
            selected_by_resname[resname] = candidate
            continue

        if rank == prev_rank and int(subst_id) < prev_subst_id:
            selected_by_resname[resname] = candidate

    extracted: List[str] = []

    for resname in sorted(selected_by_resname):
        selected = selected_by_resname[resname]
        residue_atoms = selected["residue_atoms"]
        residue_bonds = selected["residue_bonds"]
        meta = selected["meta"]

        output_file = out / f"{resname}.mol2"
        meta_file = out / f"{resname}.split.json"

        _write_residue_mol2(output_file, resname, residue_atoms, residue_bonds)
        meta_file.write_text(json.dumps(meta, indent=2), encoding="utf-8")
        extracted.append(str(output_file))

    return extracted


def extract_nonstandard_residues(input_path: str, output_dir: str) -> List[str]:
    suffix = Path(input_path).suffix.lower()
    if suffix == ".mol2":
        return extract_nonstandard_residues_from_mol2(input_path, output_dir)

    raise ValueError("Peptide splitting currently supports only MOL2 input.")


def main() -> None:
    if len(sys.argv) != 3:
        print("Usage: peptide_splitter.py peptide.mol2 output_dir")
        sys.exit(1)

    files = extract_nonstandard_residues(sys.argv[1], sys.argv[2])
    print(f"Generated {len(files)} unique NSAA residues")


if __name__ == "__main__":
    main()
