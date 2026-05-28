"""
prepgen_writer
==============

Build the ``RESNAME.mc`` file consumed by Amber's ``prepgen`` from a capped
MOL2 file.

The ``.mc`` file tells prepgen which atoms are the polymer head / tail, which
atoms belong to the residue's main chain, which atoms should be omitted from
the final library (i.e. the ACE/NME cap atoms), and what the net charge is.

Rules implemented:
    * HEAD_NAME      - only written if the head atom is present.
    * TAIL_NAME      - only written if the tail atom is present.
    * MAIN_CHAIN     - preferred = shortest path between head and tail,
                       excluding head/tail. Falls back to backbone-element
                       candidates and finally to a single anchor atom.
    * OMIT_NAME      - atoms belonging to applied caps.
    * PRE_HEAD_TYPE  - only written if HEAD_NAME is present.
    * POST_TAIL_TYPE - only written if TAIL_NAME is present.
    * CHARGE         - the residue's net charge.

This module is intentionally low-level: it does not call prepgen itself, it
only writes the ``.mc`` input file.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from collections import Counter, deque
from typing import Dict, List, Tuple, Optional, Set

from .mol2_utils import normalize_resname


@dataclass(frozen=True)
class Mol2Atom:
    atom_id: int
    name: str
    atom_type: str
    subst_name: str


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

            atoms.append(
                Mol2Atom(
                    atom_id=int(parts[0]),
                    name=parts[1],
                    atom_type=parts[5],
                    subst_name=parts[7],
                )
            )

    return atoms


def _parse_mol2_bonds(mol2_path: str) -> List[Tuple[int, int]]:
    bonds: List[Tuple[int, int]] = []
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

            bonds.append((int(parts[1]), int(parts[2])))

    return bonds


def _shortest_path(graph: Dict[int, Set[int]], start: int, goal: int) -> Optional[List[int]]:
    if start == goal:
        return [start]

    q = deque([start])
    prev: Dict[int, int] = {}
    seen = {start}

    while q:
        v = q.popleft()

        for n in graph.get(v, set()):
            if n in seen:
                continue

            seen.add(n)
            prev[n] = v

            if n == goal:
                path = [goal]
                cur = goal

                while cur != start:
                    cur = prev[cur]
                    path.append(cur)

                path.reverse()
                return path

            q.append(n)

    return None


def _looks_like_backbone_candidate(atom: Mol2Atom) -> bool:
    name = atom.name.upper()
    atype = atom.atom_type.upper()

    return (
        name.startswith("C")
        or name.startswith("N")
        or atype.startswith("C")
        or atype.startswith("N")
    )


def _dedupe_preserve_order(items: List[str]) -> List[str]:
    out: List[str] = []
    seen: Set[str] = set()

    for x in items:
        if x not in seen:
            seen.add(x)
            out.append(x)

    return out


def _select_anchor_atom(atoms: List[Mol2Atom]) -> Optional[str]:
    for a in atoms:
        if a.name.upper().startswith("C"):
            return a.name

    for a in atoms:
        if a.atom_type.upper().startswith("C"):
            return a.name

    for a in atoms:
        if not a.atom_type.upper().startswith("H") and not a.name.upper().startswith("O"):
            return a.name

    for a in atoms:
        if not a.atom_type.upper().startswith("H"):
            return a.name

    return None


def _collect_backbone_candidates(
    atoms: List[Mol2Atom],
    *,
    exclude_names: Set[str],
) -> List[str]:
    names: List[str] = []

    for atom in atoms:
        if atom.name in exclude_names:
            continue
        if _looks_like_backbone_candidate(atom):
            names.append(atom.name)

    return _dedupe_preserve_order(names)


def write_prepgen_mc_file(
    file_path: str,
    output_file: str,
    *,
    head_name: str | None = "N",
    tail_name: str | None = "C",
    main_chain: Optional[List[str]] = None,
    charge: float = 0.0,
    central_resname: Optional[str] = None,
    cap_resnames: Tuple[str, ...] = ("ACE", "NME"),
    pre_head_type: str = "C",
    post_tail_type: str = "N",
    infer_mainchain_from_connectivity: bool = True,
) -> str:
    """
    Write a prepgen ``.mc`` file with Amber-compatible semantics.

    Parameters
    ----------
    file_path : str
        Path to the input MOL2 file (typically the capped MOL2).
    output_file : str
        Path where the ``.mc`` file will be written.
    head_name : str | None
        Name of the polymer head atom (typically ``"N"``). Pass ``None`` if
        there is no head (i.e. the residue is a free N-terminus).
    tail_name : str | None
        Name of the polymer tail atom (typically ``"C"``). Pass ``None`` if
        there is no tail (i.e. the residue is a free C-terminus).
    main_chain : list[str] | None
        Explicit MAIN_CHAIN atom names. If ``None``, inferred from connectivity.
    charge : float
        Residue net charge written to the ``CHARGE`` line.
    central_resname : str | None
        Name of the central (non-cap) residue. Auto-detected if ``None``.
    cap_resnames : tuple[str, ...]
        Residue names treated as caps for OMIT_NAME generation.
    pre_head_type : str
        Atom type the head atom connects to outside the residue.
    post_tail_type : str
        Atom type the tail atom connects to outside the residue.
    infer_mainchain_from_connectivity : bool
        Enable shortest-path MAIN_CHAIN inference.

    Returns
    -------
    str
        Path to the written ``.mc`` file.
    """
    atoms = _parse_mol2_atoms(file_path)
    bonds = _parse_mol2_bonds(file_path)

    norm_cap_set = {normalize_resname(c) for c in cap_resnames}

    non_cap_atoms = [a for a in atoms if normalize_resname(a.subst_name) not in norm_cap_set]
    subst_counts = Counter(normalize_resname(a.subst_name) for a in non_cap_atoms)

    if central_resname is None:
        if subst_counts:
            central_resname = subst_counts.most_common(1)[0][0]
        else:
            central_resname = normalize_resname(atoms[0].subst_name) if atoms else "UNK"

    central_resname = normalize_resname(central_resname)

    if subst_counts:
        central_atoms = [a for a in atoms if normalize_resname(a.subst_name) == central_resname]
    else:
        central_atoms = list(atoms)

    if not central_atoms:
        raise ValueError(f"No central residue atoms found for residue '{central_resname}'.")

    omit_atom_names = [
        a.name for a in atoms
        if normalize_resname(a.subst_name) in norm_cap_set and a not in central_atoms
    ]

    id_to_atom = {a.atom_id: a for a in central_atoms}

    name_to_ids: Dict[str, List[int]] = {}
    for a in central_atoms:
        name_to_ids.setdefault(a.name, []).append(a.atom_id)

    graph: Dict[int, Set[int]] = {a.atom_id: set() for a in central_atoms}
    for a1, a2 in bonds:
        if a1 in graph and a2 in graph:
            graph[a1].add(a2)
            graph[a2].add(a1)

    has_head = bool(head_name) and head_name in name_to_ids
    has_tail = bool(tail_name) and tail_name in name_to_ids

    mainchain_names: List[str] = []

    # 1) Explicit user-supplied main chain always wins (if any atoms exist).
    if main_chain:
        mainchain_names = [x for x in main_chain if x in name_to_ids]

    # 2) Preferred Amber polymer rule: shortest path between head and tail.
    elif infer_mainchain_from_connectivity and has_head and has_tail:
        head_id = name_to_ids[head_name][0]
        tail_id = name_to_ids[tail_name][0]
        path = _shortest_path(graph, head_id, tail_id)

        if path and len(path) >= 3:
            path_names = [id_to_atom[i].name for i in path]
            mainchain_names = [n for n in path_names if n != head_name and n != tail_name]

    # 3) If only one terminus is present, keep a conservative fallback so
    #    prepgen still gets useful guidance.
    if not mainchain_names and (has_head or has_tail):
        exclude_names: Set[str] = set()
        if has_head and head_name:
            exclude_names.add(head_name)
        if has_tail and tail_name:
            exclude_names.add(tail_name)
        mainchain_names = _collect_backbone_candidates(central_atoms, exclude_names=exclude_names)

    # 4) Final fallback for carbocyclic / unusual residues.
    if not mainchain_names:
        anchor = _select_anchor_atom(central_atoms)
        if anchor:
            mainchain_names = [anchor]

    mainchain_names = _dedupe_preserve_order(mainchain_names)
    omit_atom_names = _dedupe_preserve_order(omit_atom_names)

    out = Path(output_file)
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w", encoding="utf-8") as fh:
        if has_head and head_name:
            fh.write(f"HEAD_NAME {head_name}\n")

        if has_tail and tail_name:
            fh.write(f"TAIL_NAME {tail_name}\n")

        for mc in mainchain_names:
            fh.write(f"MAIN_CHAIN {mc}\n")

        for name in omit_atom_names:
            fh.write(f"OMIT_NAME {name}\n")

        if has_head:
            fh.write(f"PRE_HEAD_TYPE {pre_head_type}\n")

        if has_tail:
            fh.write(f"POST_TAIL_TYPE {post_tail_type}\n")

        fh.write(f"CHARGE {float(charge)}\n")

    return str(out)


# Backward-compatible alias (older code calls process_mol2_file).
process_mol2_file = write_prepgen_mc_file
