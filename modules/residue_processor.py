"""
modules.residue_processor
=========================

Stage-1 of the NSAA workflow. For a single residue (or, in peptide mode, every
non-standard residue extracted from a peptide), this module:

    1. Copies / converts the input to ``residue.mol2``.
    2. Decides the head / tail / OXT-preservation policy by combining
       (a) any user-supplied residue map entry,
       (b) the peptide-splitter metadata,
       (c) the manual ``--terminal`` CLI flag.
    3. Runs ``capping.py`` (PyMOL) with the right flags.
    4. Runs antechamber (or the RESP workflow) to assign partial charges.
    5. Writes a SINGLE consolidated ``residue_meta.json`` per residue, with
       a real ``main_chain`` list inferred from connectivity when none was
       provided.

The class is intentionally stateless beyond its constructor arguments so that
``main.py`` can build one per input file and easily track per-residue
success / failure.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any, Optional

from modules.resp_workflow import run_resp_charge_workflow
from modules.peptide_splitter import extract_nonstandard_residues
from modules.mol2_utils import (
    add_proton_to_basic_n_for_charge,
    adjust_charge_for_even_electrons,
    classify_residue_net_charge,
    electron_count_for_charge,
    fix_backbone_atom_types,
    infer_main_chain_from_mol2,
    has_atom_named,
    normalize_resname,
    renormalize_mol2_partial_charges_to_integer,
    validate_molecule,
)


TERMINAL_CHOICES = ("none", "n", "c", "both")


class NonStandardAminoAcidProcessor:
    """
    Single-residue processor.

    Parameters
    ----------
    input_file:
        Path to a ``.mol2`` or ``.pdb`` file.
    charge_model:
        One of ``"abcg2"``, ``"bcc"``, ``"gas"``, or ``"resp"``.
    residue_map:
        Dict keyed by residue name. Values are dicts with optional keys
        ``head``, ``tail``, ``mainchain``, ``net_charge``, ``pre_head_type``,
        ``post_tail_type``. See ``examples/residue_map.json``.
    default_net_charge:
        Fallback net charge when nothing else specifies one.
    net_charge_override:
        CLI override; takes precedence over everything but the residue map.
    output_base:
        Base directory under which a per-residue folder will be created.
    terminal_mode:
        One of ``"none" | "n" | "c" | "both"``. For single-residue mode this
        lets the user mark the residue as a free amine ("n"), free carboxyl
        / OXT-bearing C terminus ("c"), or both ("both" = isolated /
        free amino acid).
    """

    def __init__(
        self,
        input_file: str,
        charge_model: str = "abcg2",
        residue_map: Optional[dict] = None,
        default_net_charge: int = 0,
        net_charge_override: Optional[int] = None,
        output_base: str = ".",
        terminal_mode: str = "none",
    ):
        self.input_file = input_file
        self.charge_model = str(charge_model).strip().lower()
        self.residue_map = residue_map or {}
        self.default_net_charge = default_net_charge
        self.net_charge_override = net_charge_override
        self.output_base = Path(output_base).resolve()

        tm = str(terminal_mode or "none").strip().lower()
        if tm not in TERMINAL_CHOICES:
            raise ValueError(
                f"terminal_mode must be one of {TERMINAL_CHOICES}, got: {terminal_mode!r}"
            )
        self.terminal_mode = tm

    # ------------------------------------------------------------------
    # Low-level helpers
    # ------------------------------------------------------------------
    def _run(self, cmd: list[str], cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
        return subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=str(cwd) if cwd else None,
        )

    @staticmethod
    def _write_subprocess_log(
        log_path: Path,
        *,
        stage: str,
        cmd: list[str],
        result: subprocess.CompletedProcess,
    ) -> None:
        """Append a full subprocess invocation to a per-residue log file."""
        with log_path.open("a", encoding="utf-8") as f:
            f.write(f"\n\n=== {stage} ===\n")
            f.write("CMD: " + " ".join(cmd) + "\n")
            f.write("=== STDOUT ===\n")
            f.write(result.stdout or "")
            f.write("\n=== STDERR ===\n")
            f.write(result.stderr or "")
            f.write(f"\n=== RETURN CODE === {result.returncode}\n")

    def _get_residue_name(self, file_path: str) -> str:
        p = Path(file_path)

        if p.suffix.lower() == ".mol2":
            lines = p.read_text(encoding="utf-8", errors="ignore").splitlines()

            in_sub = False
            for line in lines:
                if line.startswith("@<TRIPOS>SUBSTRUCTURE"):
                    in_sub = True
                    continue
                if line.startswith("@<TRIPOS>") and in_sub:
                    break
                if in_sub:
                    parts = line.split()
                    if len(parts) >= 2:
                        return normalize_resname(parts[1])

            in_atoms = False
            for line in lines:
                if line.startswith("@<TRIPOS>ATOM"):
                    in_atoms = True
                    continue
                if line.startswith("@<TRIPOS>") and in_atoms:
                    break
                if in_atoms:
                    parts = line.split()
                    if len(parts) >= 8:
                        return normalize_resname(parts[7])

        return normalize_resname(p.stem)

    def _normalize_nullable_name(self, value: object, default: Optional[str]) -> Optional[str]:
        if value is None:
            return default
        s = str(value).strip()
        if not s:
            return default
        if s.upper() in {"NONE", "NULL", "0"}:
            return None
        return s

    def _read_split_meta(self, input_path: Path) -> dict:
        candidates = [
            input_path.with_suffix(".split.json"),
            input_path.parent / f"{input_path.stem}.split.json",
        ]
        for meta_path in candidates:
            if not meta_path.exists():
                continue
            try:
                data = json.loads(meta_path.read_text(encoding="utf-8"))
                if isinstance(data, dict):
                    return data
            except Exception:
                continue
        return {}

    # ------------------------------------------------------------------
    # Configuration assembly
    # ------------------------------------------------------------------
    def _apply_terminal_mode(
        self,
        *,
        head_name: Optional[str],
        tail_name: Optional[str],
        preserve_oxt: bool,
    ) -> tuple[Optional[str], Optional[str], bool]:
        """
        Override head/tail/OXT-preservation according to the manual --terminal flag.

        terminal_mode meanings:
            * "none"  : default, no override (used by peptide mode)
            * "n"     : residue is N-terminal -> skip ACE cap, keep N protons
            * "c"     : residue is C-terminal -> skip NME cap, preserve OXT
            * "both"  : free amino acid       -> skip both caps, preserve OXT
        """
        if self.terminal_mode == "none":
            return head_name, tail_name, preserve_oxt

        if self.terminal_mode in {"n", "both"}:
            head_name = None
        if self.terminal_mode in {"c", "both"}:
            tail_name = None
            preserve_oxt = True

        return head_name, tail_name, preserve_oxt

    def _get_residue_cfg(self, resname: str, input_path: Optional[Path] = None) -> dict:
        raw = dict(self.residue_map.get(resname, {}))
        split_meta = self._read_split_meta(input_path) if input_path else {}

        # Peptide-splitter-determined effective head/tail (after topology rules).
        default_head = split_meta.get("head_name", "N")
        default_tail = split_meta.get("tail_name", "C")

        head_raw = raw.get("head", default_head)
        tail_raw = raw.get("tail", default_tail)

        head_name = self._normalize_nullable_name(head_raw, default_head)
        tail_name = self._normalize_nullable_name(tail_raw, default_tail)

        # OXT-preservation flag from peptide splitter (only true for c_term_like
        # or isolated residues that actually have an OXT atom).
        preserve_oxt = bool(split_meta.get("preserve_oxt", False))

        # Manual --terminal flag (single-residue mode) overrides everything above.
        head_name, tail_name, preserve_oxt = self._apply_terminal_mode(
            head_name=head_name,
            tail_name=tail_name,
            preserve_oxt=preserve_oxt,
        )

        main_chain = raw.get("mainchain")
        if not isinstance(main_chain, list):
            main_chain = split_meta.get("main_chain")
            if not isinstance(main_chain, list):
                main_chain = None

        pre_head_type = str(raw.get("pre_head_type", split_meta.get("pre_head_type", "C")))
        post_tail_type = str(raw.get("post_tail_type", split_meta.get("post_tail_type", "N")))

        return {
            "head_name": head_name,
            "tail_name": tail_name,
            "main_chain": main_chain,
            "pre_head_type": pre_head_type,
            "post_tail_type": post_tail_type,
            "preserve_oxt": preserve_oxt,
            "split_meta": split_meta,
        }

    def _prepare_input_mol2(self, input_path: Path, residue_dir: Path) -> Path:
        residue_file = residue_dir / "residue.mol2"
        suffix = input_path.suffix.lower()

        if suffix == ".mol2":
            residue_file.write_text(
                input_path.read_text(encoding="utf-8", errors="replace"),
                encoding="utf-8",
            )
            return residue_file

        if suffix == ".pdb":
            convert_script = Path(__file__).parent / "pdb_to_mol2.py"
            cmd = ["python", str(convert_script), str(input_path), str(residue_file)]
            result = self._run(cmd, cwd=residue_dir)
            if result.returncode != 0:
                log_path = residue_dir / f"{residue_dir.name}.log"
                self._write_subprocess_log(
                    log_path, stage="PDB_TO_MOL2", cmd=cmd, result=result,
                )
                raise RuntimeError(
                    f"PDB conversion failed (see {log_path})"
                )
            return residue_file

        raise ValueError(f"Unsupported format: {input_path.suffix}")

    def _resolve_input_net_charge(
        self,
        resname: str,
        input_mol2_path: str,
        split_meta: dict,
    ) -> tuple[int, str]:
        if resname in self.residue_map and "net_charge" in self.residue_map[resname]:
            return int(self.residue_map[resname]["net_charge"]), "map"

        if self.net_charge_override is not None:
            return int(self.net_charge_override), "cli"

        full_context = split_meta.get("full_context_net_charge")
        classified, classified_source = classify_residue_net_charge(
            input_mol2_path, resname=resname,
        )

        # Reconciliation:
        #   * Non-zero splitter value wins (it sees external bonds & full
        #     peptide context that the single-residue classifier can't).
        #   * If the splitter returned exactly 0 but the single-residue
        #     classifier found a clear charged group, prefer the classifier.
        #     This catches the case where a peptide MOL2 carries a residue
        #     in (e.g.) the -NH2 / COOH form on disk but the side-chain
        #     pattern is unmistakably Lys / Arg / Asp / Glu - the classifier
        #     is more chemistry-aware than the full-context sum.
        if full_context is not None and int(full_context) != 0:
            return (
                int(full_context),
                str(split_meta.get("full_context_charge_source", "split_full_context")),
            )

        if classified is not None and int(classified) != 0:
            if full_context is not None:
                print(
                    f"[{resname}] charge classifier suggests {classified} "
                    f"(peptide-context returned {full_context}); using classifier"
                )
                return int(classified), f"{classified_source}_over_neutral_split"
            return int(classified), classified_source

        if full_context is not None:
            return (
                int(full_context),
                str(split_meta.get("full_context_charge_source", "split_full_context")),
            )

        if classified is not None:
            return int(classified), classified_source

        return int(self.default_net_charge), "default"

    def _correct_parity_on_capped(
        self,
        *,
        resname: str,
        capped_mol2_path: str,
        net_charge: int,
        charge_source: str,
    ) -> tuple[int, str]:
        """
        After PyMOL capping, verify that the *capped* residue is closed-shell
        at the resolved charge.

        Capping.py has already tried to fix protonation by deprotonating
        carboxylate -OH groups (the Asp / Glu / free-C-terminus case). If
        parity is still off, we try the opposite direction here: append an
        H to a basic N (the Arg guanidinium / Lys NZ / N-terminal amine
        case). Only if both directions fail do we fall back to adjusting
        the declared net charge.

        User-pinned charges (residue map / CLI override) are honoured
        verbatim; we just warn loudly if they leave the capped residue
        open-shell.
        """
        if charge_source in {"map", "cli"}:
            if electron_count_for_charge(capped_mol2_path, net_charge) % 2 != 0:
                # Try protonation even for user-pinned charges, since the
                # whole point of pinning is to assert the chemistry. If a
                # basic N is available, fix it; otherwise warn.
                info = add_proton_to_basic_n_for_charge(capped_mol2_path, net_charge)
                if info["hs_added"] > 0:
                    print(
                        f"[{resname}] protonated basic N ({info['target_n_name']}) "
                        f"to match pinned charge {net_charge} "
                        f"(electrons {info['pre_electrons']} -> {info['post_electrons']})"
                    )
                else:
                    print(
                        f"[{resname}] WARNING: user-pinned charge {net_charge} leaves "
                        f"the capped residue open-shell; sqm/AM1-BCC will likely fail"
                    )
            return net_charge, charge_source

        if electron_count_for_charge(capped_mol2_path, net_charge) % 2 == 0:
            # Capping's deprotonation audit already closed the shell.
            return net_charge, charge_source

        # Try protonating a basic N (Arg guanidinium / Lys NZ / etc.).
        info = add_proton_to_basic_n_for_charge(capped_mol2_path, net_charge)
        if info["hs_added"] > 0 and info["parity_ok_after"]:
            print(
                f"[{resname}] protonated basic N ({info['target_n_name']}) "
                f"to match charge {net_charge} "
                f"(electrons {info['pre_electrons']} -> {info['post_electrons']})"
            )
            return net_charge, charge_source

        # Neither deprotonation nor protonation worked. Fall back to
        # adjusting the declared charge to the nearest closed-shell value.
        corrected, was_adjusted = adjust_charge_for_even_electrons(
            capped_mol2_path, net_charge,
        )
        if not was_adjusted:
            return net_charge, charge_source

        if electron_count_for_charge(capped_mol2_path, corrected) % 2 != 0:
            raise RuntimeError(
                f"capped residue is open-shell at any nearby charge "
                f"(tried {net_charge}, +/-1, +/-2). No carboxylate -OH or "
                f"basic N was available to fix protonation. Pin the "
                f"correct charge with --map or --net-charge."
            )

        print(
            f"[{resname}] post-cap parity check: no deprotonatable -OH or "
            f"protonatable basic N; adjusting charge {net_charge} -> {corrected}"
        )
        return corrected, f"{charge_source}+parity_corrected"

    # ------------------------------------------------------------------
    # Metadata writing (consolidated JSON)
    # ------------------------------------------------------------------
    def _write_residue_meta(
        self,
        residue_dir: Path,
        *,
        resname: str,
        cfg: dict,
        capping_meta: dict,
        net_charge: int,
        charge_source: str,
        validation_warnings: list[str],
        charge_backend_meta: Optional[dict[str, Any]] = None,
        input_mol2_for_inference: Optional[Path] = None,
    ) -> None:
        """
        Write the single consolidated ``residue_meta.json`` for this residue.

        Also infers ``main_chain`` from connectivity when it is missing, so the
        JSON no longer shows ``"main_chain": null`` for typical alpha-amino
        acids. Finally, deletes ``residue_capping_meta.json`` so each residue
        folder ends up with exactly one metadata file.
        """
        split_meta = cfg.get("split_meta", {}) or {}

        # ------ main_chain inference (if missing) ------
        main_chain = cfg.get("main_chain")
        if (
            (not isinstance(main_chain, list) or not main_chain)
            and input_mol2_for_inference is not None
            and cfg.get("head_name")
            and cfg.get("tail_name")
        ):
            inferred = infer_main_chain_from_mol2(
                str(input_mol2_for_inference),
                head_name=cfg.get("head_name"),
                tail_name=cfg.get("tail_name"),
            )
            if inferred:
                main_chain = inferred

        meta = {
            "resname": resname,
            # Resolved capping policy used by this pipeline run:
            "head_name": cfg.get("head_name"),
            "tail_name": cfg.get("tail_name"),
            "main_chain": main_chain,
            "pre_head_type": cfg.get("pre_head_type", "C"),
            "post_tail_type": cfg.get("post_tail_type", "N"),
            "preserve_oxt": bool(cfg.get("preserve_oxt", False)),
            "terminal_mode": self.terminal_mode,
            # Echoed straight from the capping metadata, so we know exactly
            # what PyMOL did to the residue:
            "requested_head_name": capping_meta.get("requested_head_name", cfg.get("head_name")),
            "requested_tail_name": capping_meta.get("requested_tail_name", cfg.get("tail_name")),
            "oxt_preserved": bool(capping_meta.get("oxt_preserved", False)),
            "has_head": bool(capping_meta.get("has_head", False)),
            "has_tail": bool(capping_meta.get("has_tail", False)),
            "applied_caps": list(capping_meta.get("applied_caps", [])),
            # Charge bookkeeping:
            "net_charge": int(net_charge),
            "net_charge_source": charge_source,
            "charge_model": self.charge_model,
            "charge_backend_meta": charge_backend_meta or {},
            "validation_warnings": validation_warnings,
            # Peptide-mode provenance (will be absent / null for single-residue mode):
            "source_input": split_meta.get("source_input"),
            "source_subst_id": split_meta.get("source_subst_id"),
            "source_subst_name": split_meta.get("source_subst_name"),
            "external_bonds": split_meta.get("external_bonds", []),
            "is_polymer_internal": bool(split_meta.get("is_polymer_internal", False)),
            "is_n_terminal_like": bool(split_meta.get("is_n_terminal_like", False)),
            "is_c_terminal_like": bool(split_meta.get("is_c_terminal_like", False)),
            "is_cyclic_like": bool(split_meta.get("is_cyclic_like", False)),
            "topology": split_meta.get("topology"),
            "has_oxt": split_meta.get("has_oxt"),
            "has_n_terminal_protons": split_meta.get("has_n_terminal_protons"),
            "raw_head_name": split_meta.get("raw_head_name"),
            "raw_tail_name": split_meta.get("raw_tail_name"),
            "full_context_net_charge": split_meta.get("full_context_net_charge"),
            "full_context_charge_source": split_meta.get("full_context_charge_source"),
        }

        (residue_dir / "residue_meta.json").write_text(
            json.dumps(meta, indent=2),
            encoding="utf-8",
        )

        # Delete the small capping-meta JSON so each residue folder ends up
        # with a single consolidated metadata file.
        capping_meta_file = residue_dir / "residue_capping_meta.json"
        if capping_meta_file.exists():
            try:
                capping_meta_file.unlink()
            except Exception:
                pass

    # ------------------------------------------------------------------
    # Charge assignment
    # ------------------------------------------------------------------
    def _assign_charges_with_antechamber(
        self,
        *,
        capped_file: Path,
        charged_file: Path,
        resname: str,
        net_charge: int,
    ) -> dict[str, Any]:
        cmd = [
            "antechamber",
            "-fi", "mol2",
            "-i", str(capped_file),
            "-bk", resname,
            "-fo", "mol2",
            "-o", str(charged_file),
            "-c", self.charge_model,
            "-at", "amber",
            "-nc", str(net_charge),
        ]
        charge_result = self._run(cmd, cwd=charged_file.parent)
        if charge_result.returncode != 0:
            log_path = charged_file.parent / f"{resname}_charge.log"
            self._write_subprocess_log(
                log_path, stage="ANTECHAMBER_CHARGE", cmd=cmd, result=charge_result,
            )
            raise RuntimeError(
                f"antechamber charge assignment failed for {resname} (see {log_path})"
            )

        return {
            "charge_method": self.charge_model,
            "charge_backend": "antechamber",
            "net_charge": int(net_charge),
            "files": {
                "input_mol2": str(capped_file),
                "output_mol2": str(charged_file),
            },
        }

    def _assign_charges_with_resp(
        self,
        *,
        capped_file: Path,
        charged_file: Path,
        resname: str,
        net_charge: int,
    ) -> dict[str, Any]:
        return run_resp_charge_workflow(
            capped_file=str(capped_file),
            charged_file=str(charged_file),
            resname=resname,
            net_charge=int(net_charge),
            residue_dir=str(charged_file.parent),
        )

    def _assign_charges(
        self,
        *,
        capped_file: Path,
        charged_file: Path,
        resname: str,
        net_charge: int,
    ) -> dict[str, Any]:
        if self.charge_model in {"bcc", "gas", "abcg2"}:
            return self._assign_charges_with_antechamber(
                capped_file=capped_file,
                charged_file=charged_file,
                resname=resname,
                net_charge=net_charge,
            )
        if self.charge_model == "resp":
            return self._assign_charges_with_resp(
                capped_file=capped_file,
                charged_file=charged_file,
                resname=resname,
                net_charge=net_charge,
            )
        raise ValueError(f"Unsupported charge model: {self.charge_model}")

    # ------------------------------------------------------------------
    # Main per-residue pipeline
    # ------------------------------------------------------------------
    def _process_input_path(self, input_path: Path) -> tuple[str, list[str]]:
        """
        Returns ``(resname, [charged_mol2_path, ...])`` on success.
        Raises on failure (caller catches and records the failed resname).
        """
        resname = self._get_residue_name(str(input_path))
        residue_dir = self.output_base / resname
        residue_dir.mkdir(parents=True, exist_ok=True)

        cfg = self._get_residue_cfg(resname, input_path=input_path)
        input_mol2 = self._prepare_input_mol2(input_path, residue_dir)

        net_charge, charge_source = self._resolve_input_net_charge(
            resname,
            str(input_mol2),
            cfg.get("split_meta", {}) or {},
        )

        # NOTE: parity is NOT checked on the input MOL2. In peptide mode the
        # extraction leaves dangling peptide bonds, so an internal residue is
        # intrinsically open-shell until capped. Checking parity here would
        # produce misleading warnings; we check it after capping instead.
        input_warnings = validate_molecule(
            str(input_mol2), expected_charge=net_charge, check_parity=False,
        )

        # ------ run capping.py (PyMOL) ------
        cap_script = Path(__file__).parent / "capping.py"
        head_arg = cfg["head_name"] if cfg["head_name"] else "NONE"
        tail_arg = cfg["tail_name"] if cfg["tail_name"] else "NONE"

        cap_cmd = [
            "python", str(cap_script), str(residue_dir),
            "--head", str(head_arg),
            "--tail", str(tail_arg),
            "--net-charge", str(int(net_charge)),
        ]

        # Determine whether to ask capping to preserve OXT. We pass the flag
        # whenever:
        #   * the peptide splitter said so (peptide mode), OR
        #   * the user asked for it via --terminal c|both (single-residue mode), OR
        #   * the input actually has an OXT atom AND we already decided to skip
        #     the tail cap (free C terminus).
        wants_preserve_oxt = bool(cfg.get("preserve_oxt", False))
        if not wants_preserve_oxt and cfg["tail_name"] is None:
            if has_atom_named(str(input_mol2), "OXT"):
                wants_preserve_oxt = True

        if wants_preserve_oxt:
            cap_cmd.append("--preserve-oxt")
            # propagate into cfg so the consolidated meta JSON reflects it
            cfg["preserve_oxt"] = True

        cap_result = self._run(cap_cmd, cwd=residue_dir)
        if cap_result.returncode != 0:
            log_path = residue_dir / f"{resname}_capping.log"
            self._write_subprocess_log(
                log_path, stage="CAPPING", cmd=cap_cmd, result=cap_result,
            )
            raise RuntimeError(
                f"Capping failed for {resname} (see {log_path})"
            )

        capped_file = residue_dir / "residue_capped.mol2"
        if not capped_file.exists() or capped_file.stat().st_size == 0:
            raise RuntimeError(f"Capped MOL2 was not created for {resname}: {capped_file}")

        capping_meta_file = residue_dir / "residue_capping_meta.json"
        try:
            capping_meta = json.loads(capping_meta_file.read_text(encoding="utf-8"))
        except Exception:
            capping_meta = {}

        # ------ post-cap parity correction ------
        # PyMOL h_add inside capping.py does not respect MOL2 formal charges,
        # so the capped residue may end up with one more / fewer proton than
        # the resolved charge implies. We verify parity on the file antechamber
        # will actually consume and adjust by +/-1 (then +/-2) if needed.
        net_charge, charge_source = self._correct_parity_on_capped(
            resname=resname,
            capped_mol2_path=str(capped_file),
            net_charge=net_charge,
            charge_source=charge_source,
        )

        validation_warnings = (
            input_warnings
            + validate_molecule(str(capped_file), expected_charge=net_charge)
        )
        print(f"[{resname}] net charge = {net_charge}")
        print(f"[{resname}] charge model = {self.charge_model}")
        if cfg.get("preserve_oxt"):
            print(f"[{resname}] OXT preserved (terminal residue, NME cap skipped)")
        for w in validation_warnings:
            print(f"[{resname}] warning: {w}")

        # ------ charge assignment ------
        charged_file = residue_dir / f"{resname}.mol2"
        charge_backend_meta = self._assign_charges(
            capped_file=capped_file,
            charged_file=charged_file,
            resname=resname,
            net_charge=int(net_charge),
        )

        renormalize_mol2_partial_charges_to_integer(str(charged_file), int(net_charge))
        fix_backbone_atom_types(str(charged_file))

        # ------ consolidated metadata ------
        self._write_residue_meta(
            residue_dir,
            resname=resname,
            cfg=cfg,
            capping_meta=capping_meta,
            net_charge=net_charge,
            charge_source=charge_source,
            validation_warnings=validation_warnings,
            charge_backend_meta=charge_backend_meta,
            input_mol2_for_inference=input_mol2,
        )

        return resname, [str(charged_file)]

    # ------------------------------------------------------------------
    # Public top-level entry points
    # ------------------------------------------------------------------
    def process_single_residue(self) -> tuple[list[str], list[str], list[tuple[str, str]]]:
        """
        Process the single input residue.

        Returns
        -------
        (charged_files, successful_resnames, failed_resnames)
            * ``charged_files`` is forwarded to the AMBER stage.
            * ``successful_resnames`` lists residues that survived the
              charge-assignment / capping stage.
            * ``failed_resnames`` is a list of ``(resname, reason)`` tuples.
        """
        path = Path(self.input_file).resolve()
        try:
            resname, charged = self._process_input_path(path)
            return charged, [resname], []
        except Exception as exc:
            resname = self._get_residue_name(str(path))
            print(f"[FAILED] {resname}: {exc}")
            return [], [], [(resname, f"processor: {exc}")]

    def process_peptide(self) -> tuple[list[str], list[str], list[tuple[str, str]]]:
        """
        Extract every non-standard residue from the peptide and process each
        one independently. Per-residue failures are recorded but do not abort
        the rest of the batch.
        """
        peptide_path = Path(self.input_file).resolve()
        split_dir = self.output_base / "_extracted_nonstandard_residues"
        split_dir.mkdir(parents=True, exist_ok=True)

        extracted_paths = extract_nonstandard_residues(str(peptide_path), str(split_dir))

        charged_files: list[str] = []
        successful: list[str] = []
        failed: list[tuple[str, str]] = []

        for residue in extracted_paths:
            try:
                resname, charged = self._process_input_path(Path(residue).resolve())
                charged_files.extend(charged)
                successful.append(resname)
            except Exception as exc:
                resname = normalize_resname(Path(residue).stem)
                print(f"[FAILED] {resname}: {exc}")
                failed.append((resname, f"processor: {exc}"))
                continue

        return charged_files, successful, failed
