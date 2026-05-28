"""
nsaa-paramgen
=============

Command-line entry point for the NSAA AMBER (+ optional GROMACS)
parameter generation pipeline.

Typical usage
-------------

Single residue::

    python main.py -i MLY.mol2

Peptide mode (extract NSAAs, then process each, then build whole-peptide
prmtop/inpcrd)::

    python main.py -p peptide.mol2

Generate GROMACS outputs as well::

    python main.py -p peptide.mol2 --gmx

Output
------
Per residue: .lib, .frcmod (x2), .prepin, .prmtop, .inpcrd, (+ .top/.gro
if --gmx). Peptide mode adds a ``peptide/`` folder with peptide.prmtop,
peptide.inpcrd, and (if --gmx) peptide.top and peptide.gro.
"""
from __future__ import annotations

import argparse
import glob
import json
import shutil
from pathlib import Path
from typing import Iterable, List

from modules.residue_processor import NonStandardAminoAcidProcessor, TERMINAL_CHOICES
from modules.antechamber_runner import run_antechamber_for_all
from modules.peptide_assembler import assemble_peptide


# ----------------------------------------------------------------------------
# Residue map loader
# ----------------------------------------------------------------------------
def load_residue_map(path: str | None) -> dict:
    if not path:
        return {}

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Mapping file not found: {p}")

    data = json.loads(p.read_text(encoding="utf-8")) or {}
    out: dict[str, dict] = {}
    for k, v in data.items():
        if isinstance(v, dict):
            out[str(k).upper()] = v
    return out


# ----------------------------------------------------------------------------
# Argparser
# ----------------------------------------------------------------------------
def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "nsaa-paramgen: NSAA AMBER (+ optional GROMACS) parameter "
            "generation (ff19SB/ff14SB + gaff/gaff2)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--input", nargs="+",
                       help="Single residue input(s) (.mol2 or .pdb).")
    group.add_argument("-p", "--peptide",
                       help="Peptide input; extracts NSAAs, parametrizes "
                            "each, then assembles whole-peptide prmtop/inpcrd.")
    group.add_argument("-b", "--batch", nargs="+",
                       help="Batch input: file list / directory / glob.")

    parser.add_argument("-t", "--terminal", choices=TERMINAL_CHOICES, default="none")
    parser.add_argument("--backbone", "-bb",
                        choices=["ff14SB", "ff19SB", "ff99SB"], default="ff19SB")
    parser.add_argument("--sidechain", "-sc",
                        choices=["gaff", "gaff2"], default="gaff2")
    parser.add_argument("--charge", "-c",
                        choices=["gas", "bcc", "resp", "abcg2"], default="abcg2")
    parser.add_argument("--gmx", "-gmx", action="store_true",
                        help="Also generate GROMACS .top/.gro for each "
                             "residue (and for the whole peptide in -p mode). "
                             "Requires ParmEd.")
    parser.add_argument("--map", default=None)
    parser.add_argument("--default-net-charge", type=int, default=0)
    parser.add_argument("--net-charge", "-nc", type=int, default=None)
    parser.add_argument("--out", "-o", default=".")

    return parser


# ----------------------------------------------------------------------------
# Batch argument expansion
# ----------------------------------------------------------------------------
def _expand_one_batch_item(batch_arg: str) -> List[str]:
    p = Path(batch_arg)

    if p.exists() and p.is_file():
        if p.suffix.lower() in {".mol2", ".pdb"}:
            return [str(p)]
        items: List[str] = []
        with p.open("r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if s and not s.startswith("#"):
                    items.append(s)
        return items

    if p.exists() and p.is_dir():
        files = sorted(str(x) for x in p.glob("*.mol2"))
        files.extend(sorted(str(x) for x in p.glob("*.pdb")))
        return files

    return sorted(glob.glob(batch_arg))


def expand_batch_arguments(batch_args: Iterable[str]) -> List[str]:
    seen: set[str] = set()
    out: List[str] = []
    for arg in batch_args:
        for item in _expand_one_batch_item(arg):
            if item not in seen:
                seen.add(item)
                out.append(item)
    return out


# ----------------------------------------------------------------------------
# Failure handling
# ----------------------------------------------------------------------------
def _move_failed_folders(out_base: Path, failed_resnames: Iterable[str]) -> None:
    failed_root = out_base / "failed"
    moved_any = False

    for resname in failed_resnames:
        src = out_base / resname
        if not src.exists() or not src.is_dir():
            continue
        if not moved_any:
            failed_root.mkdir(parents=True, exist_ok=True)
            moved_any = True
        dst = failed_root / resname
        if dst.exists():
            shutil.rmtree(dst, ignore_errors=True)
        try:
            shutil.move(str(src), str(dst))
        except Exception as exc:
            print(f"[WARN] could not move failed residue folder {src} -> {dst}: {exc}")


def _write_resname_list(path: Path, resnames: Iterable[str]) -> None:
    seen: set[str] = set()
    ordered: List[str] = []
    for r in resnames:
        if r and r not in seen:
            seen.add(r)
            ordered.append(r)
    path.write_text("\n".join(ordered) + ("\n" if ordered else ""), encoding="utf-8")


def _short_reason(raw: str) -> str:
    text = (raw or "").strip()
    for prefix in ("processor:", "amber:", "peptide:", "top-level:"):
        if text.lower().startswith(prefix):
            text = text[len(prefix):].lstrip()
            break
    first_line = text.splitlines()[0].strip() if text else ""
    return first_line or "unknown error"


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    residue_map = load_residue_map(args.map)

    out_base = Path(args.out).resolve()
    out_base.mkdir(parents=True, exist_ok=True)

    charged_files: List[str] = []
    successful_resnames: List[str] = []
    failed_resnames: List[str] = []
    failure_records: List[tuple[str, str]] = []

    def build_processor(path: str) -> NonStandardAminoAcidProcessor:
        return NonStandardAminoAcidProcessor(
            input_file=path,
            charge_model=args.charge,
            residue_map=residue_map,
            default_net_charge=args.default_net_charge,
            net_charge_override=args.net_charge,
            output_base=str(out_base),
            terminal_mode=args.terminal,
        )

    def process_one(path: str) -> None:
        proc = build_processor(path)
        charged, ok, fails = proc.process_single_residue()
        charged_files.extend(charged)
        successful_resnames.extend(ok)
        for resname, reason in fails:
            failed_resnames.append(resname)
            failure_records.append((resname, reason))

    # ---------- dispatch ----------
    if args.input:
        input_paths = list(args.input)
        if len(input_paths) > 1:
            print(f"\nProcessing {len(input_paths)} residue files\n")
        for fp in input_paths:
            print(f"\n--- Processing {fp} ---")
            try:
                process_one(fp)
            except Exception as exc:
                resname = Path(fp).stem
                print(f"[FAILED] {fp}: {exc}")
                failed_resnames.append(resname)
                failure_records.append((resname, f"top-level: {exc}"))

    elif args.peptide:
        try:
            proc = build_processor(args.peptide)
            charged, ok, fails = proc.process_peptide()
            charged_files.extend(charged)
            successful_resnames.extend(ok)
            for resname, reason in fails:
                failed_resnames.append(resname)
                failure_records.append((resname, reason))
        except Exception as exc:
            print(f"[FAILED] peptide processing: {args.peptide}: {exc}")
            failure_records.append((Path(args.peptide).stem, f"peptide: {exc}"))

    elif args.batch:
        batch_items = expand_batch_arguments(args.batch)
        if not batch_items:
            print(f"Batch mode: nothing matched: {args.batch}")
            return

        print(f"\nBatch mode: {len(batch_items)} residues detected\n")
        for fp in batch_items:
            print(f"\n--- Processing {fp} ---")
            try:
                process_one(fp)
            except Exception as exc:
                resname = Path(fp).stem
                print(f"[FAILED] {fp}: {exc}")
                failed_resnames.append(resname)
                failure_records.append((resname, f"top-level: {exc}"))

        print("\nBatch processing finished")

    # ---------- AMBER toolchain stage ----------
    peptide_result: dict | None = None
    if charged_files:
        amber_result = run_antechamber_for_all(
            charged_files,
            backbone=args.backbone,
            sidechain=args.sidechain,
            charge=args.charge,
            generate_gmx=args.gmx,
        )

        amber_ok = list(amber_result.get("successful", []))
        amber_failed = list(amber_result.get("failed", []))

        amber_ok_set = set(amber_ok)
        amber_failed_set = {r for r, _ in amber_failed}

        full_success = [r for r in successful_resnames if r in amber_ok_set]

        for resname, reason in amber_failed:
            if resname not in failed_resnames:
                failed_resnames.append(resname)
            failure_records.append((resname, f"amber: {reason}"))

        unseen = [r for r in successful_resnames
                  if r not in amber_ok_set and r not in amber_failed_set]
        for resname in unseen:
            failed_resnames.append(resname)
            failure_records.append((resname, "amber: not processed"))

        successful_resnames = full_success

        # ---------- peptide-level assembly (-p only) ----------
        if args.peptide and successful_resnames:
            print("\n--- Assembling whole-peptide topology ---")
            try:
                peptide_result = assemble_peptide(
                    peptide_input=args.peptide,
                    successful_resnames=successful_resnames,
                    out_base=str(out_base),
                    backbone=args.backbone,
                    sidechain=args.sidechain,
                    generate_gmx=args.gmx,
                )
                if peptide_result.get("status") == "ok":
                    print(f"Peptide prmtop : {peptide_result['prmtop']}")
                    print(f"Peptide inpcrd : {peptide_result['inpcrd']}")
                    if peptide_result.get("top"):
                        print(f"Peptide top    : {peptide_result['top']}")
                        print(f"Peptide gro    : {peptide_result['gro']}")
                    missing = peptide_result.get("missing_residues") or []
                    if missing:
                        print(f"[WARN] residues without prepin/frcmod skipped "
                              f"during assembly: {missing}")
                else:
                    print(f"[FAILED] peptide assembly: "
                          f"{peptide_result.get('reason', 'unknown')}")
                    print(f"  see {peptide_result.get('log')}")
            except Exception as exc:
                print(f"[FAILED] peptide assembly: {exc}")
                peptide_result = {"status": "failed", "reason": str(exc)}
    else:
        print("\nNo charged .mol2 files were generated. Skipping AMBER toolchain.")

    # ---------- summary files ----------
    successful_resnames = [r for r in successful_resnames if r not in set(failed_resnames)]

    _write_resname_list(out_base / "successful_residues.txt", successful_resnames)
    _write_resname_list(out_base / "failed_residues.txt", failed_resnames)

    if failed_resnames:
        _move_failed_folders(out_base, failed_resnames)

    short_reasons: dict[str, str] = {}
    for r, reason in failure_records:
        if r in short_reasons:
            continue
        short_reasons[r] = _short_reason(reason)

    print("\n" + "=" * 60)
    print(f"Successful residues : {len(successful_resnames)}")
    print(f"Failed residues     : {len(failed_resnames)}")
    print("=" * 60)
    if successful_resnames:
        print("Successful:")
        for r in successful_resnames:
            print(f"  + {r}")
    if failed_resnames:
        print("Failed:")
        for r in failed_resnames:
            print(f"  - {r}: {short_reasons.get(r, 'unknown error')}")
    if peptide_result and peptide_result.get("status") == "ok":
        print(f"\nPeptide topology written to: {Path(peptide_result['prmtop']).parent}")
    print(f"\nSee {out_base / 'successful_residues.txt'}")
    print(f"See {out_base / 'failed_residues.txt'}")
    if failed_resnames:
        print(f"Failed residue folders moved into {out_base / 'failed'}/")


if __name__ == "__main__":
    main()