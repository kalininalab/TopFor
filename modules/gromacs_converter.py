"""
modules.gromacs_converter
=========================

Convert AMBER topology/coordinate files to GROMACS format using ParmEd.

ParmEd is the standard library for inter-format conversion of MD topologies
and is bundled with AmberTools (also installable via `pip install parmed`
or `conda install -c conda-forge parmed`). No MD is run here; this is a
parameter-file conversion only.

Public API
----------
    amber_to_gromacs(prmtop, inpcrd, top_out, gro_out) -> dict
    split_top_into_itp(top_path, itp_out)              -> str
    ParmEdNotAvailable                                 -> exception

NOTE: please confirm the ParmEd `save(..., format="gromacs")` call against
your installed parmed version (`python -c "import parmed; print(parmed.__version__)"`);
the public API has been stable for years but you should verify before
relying on it for production runs.
"""
from __future__ import annotations

from pathlib import Path


class ParmEdNotAvailable(RuntimeError):
    """Raised when the `parmed` Python package cannot be imported."""


def _require_parmed():
    try:
        import parmed  # noqa: F401
    except Exception as exc:
        raise ParmEdNotAvailable(
            "ParmEd is required for GROMACS conversion but could not be "
            "imported. Install via `pip install parmed` or activate an "
            "AmberTools environment that exposes its Python bindings."
        ) from exc


def amber_to_gromacs(
    prmtop: str,
    inpcrd: str,
    top_out: str,
    gro_out: str,
) -> dict:
    """
    Convert an AMBER (prmtop, inpcrd) pair to GROMACS (top, gro).

    Returns
    -------
    dict
        ``{"top": <path>, "gro": <path>}`` on success.

    Raises
    ------
    ParmEdNotAvailable, FileNotFoundError, RuntimeError
    """
    _require_parmed()
    import parmed as pmd

    prmtop_p = Path(prmtop)
    inpcrd_p = Path(inpcrd)
    if not prmtop_p.exists():
        raise FileNotFoundError(f"prmtop not found: {prmtop_p}")
    if not inpcrd_p.exists():
        raise FileNotFoundError(f"inpcrd not found: {inpcrd_p}")

    Path(top_out).parent.mkdir(parents=True, exist_ok=True)
    Path(gro_out).parent.mkdir(parents=True, exist_ok=True)

    try:
        structure = pmd.load_file(str(prmtop_p), str(inpcrd_p))
        # Explicit format keyword for the .top write; .gro is extension-detected.
        structure.save(str(top_out), format="gromacs", overwrite=True)
        structure.save(str(gro_out), overwrite=True)
    except Exception as exc:
        raise RuntimeError(
            f"ParmEd failed to convert {prmtop_p.name}: {exc}"
        ) from exc

    return {"top": str(top_out), "gro": str(gro_out)}


def split_top_into_itp(top_path: str, itp_out: str) -> str:
    """
    Best-effort textual extraction of the ``[ moleculetype ]`` ... block
    (up to but excluding ``[ system ]``) from a GROMACS .top into a
    stand-alone .itp.

    The ``[ atomtypes ]`` block is intentionally NOT included — most users
    want to keep atomtypes in the master .top so they aren't duplicated
    across multiple includes.

    Please inspect the resulting .itp; ParmEd's layout is conventional but
    other tools organize sections differently.
    """
    lines = Path(top_path).read_text(encoding="utf-8").splitlines()
    start = end = None
    for i, line in enumerate(lines):
        s = line.strip().lower()
        if start is None and s.startswith("[ moleculetype"):
            start = i
            continue
        if start is not None and s.startswith("[ system"):
            end = i
            break

    if start is None:
        raise RuntimeError(f"No [ moleculetype ] block found in {top_path}")
    if end is None:
        end = len(lines)

    body = "\n".join(lines[start:end]).rstrip() + "\n"
    Path(itp_out).parent.mkdir(parents=True, exist_ok=True)
    Path(itp_out).write_text(body, encoding="utf-8")
    return str(itp_out)