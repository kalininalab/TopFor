#!/usr/bin/env bash
# =============================================================
# setup_env.sh  –  run ONCE after activating the topfor env
# =============================================================
# Usage:
#   conda activate topfor
#   bash setup_env.sh
#
# What it does:
#   1. Writes the required AMBERHOME export into the conda env's
#      activation hook so it is set automatically every time you
#      run  `conda activate topfor`.
#   2. Verifies every tool TopFor needs is reachable.
#   3. Prints a summary and an optional next-step for RESP.
# =============================================================

set -e

# ── Colour helpers ──────────────────────────────────────────
RED='\033[0;31m'
GRN='\033[0;32m'
YLW='\033[1;33m'
BLD='\033[1m'
RST='\033[0m'

ok()   { echo -e "  ${GRN}✓${RST}  $*"; }
warn() { echo -e "  ${YLW}!${RST}  $*"; }
fail() { echo -e "  ${RED}✗${RST}  $*"; ERRORS=$((ERRORS+1)); }

ERRORS=0

# ── Guard: must be inside the topfor conda env ──────────────
if [[ -z "$CONDA_PREFIX" ]]; then
    echo -e "${RED}Error:${RST} No conda environment is active."
    echo "  Run:  conda activate topfor"
    exit 1
fi

ENV_NAME=$(basename "$CONDA_PREFIX")
echo -e "\n${BLD}TopFor environment setup${RST}  (env: $ENV_NAME)"
echo "────────────────────────────────────────────────"

# ── 1. Persist AMBERHOME into the activation hook ───────────
ACTIVATE_DIR="$CONDA_PREFIX/etc/conda/activate.d"
DEACTIVATE_DIR="$CONDA_PREFIX/etc/conda/deactivate.d"
mkdir -p "$ACTIVATE_DIR" "$DEACTIVATE_DIR"

cat > "$ACTIVATE_DIR/topfor_env_vars.sh" << 'EOF'
export AMBERHOME="$CONDA_PREFIX"
EOF

cat > "$DEACTIVATE_DIR/topfor_env_vars.sh" << 'EOF'
unset AMBERHOME
EOF

# Apply immediately for this session
export AMBERHOME="$CONDA_PREFIX"
ok "AMBERHOME set to $AMBERHOME (will auto-set on future activations)"

# ── 2. Verify required tools (always needed) ────────────────
echo -e "\n${BLD}Checking required tools:${RST}"

check_exe() {
    local name="$1"
    local purpose="$2"
    if command -v "$name" &>/dev/null; then
        local ver
        ver=$("$name" --version 2>&1 | head -1 || true)
        ok "$name  –  $purpose  ($ver)"
    else
        fail "$name  –  $purpose  NOT FOUND"
    fi
}

check_exe "antechamber" "AMBER atom typing + charge assignment"
check_exe "prepgen"     "residue template (.prepin) generation"
check_exe "parmchk2"    "missing parameter detection (.frcmod)"
check_exe "tleap"       "AMBER library (.lib) generation"

# pymol has no version flag on PATH; test via python import
if python -c "import pymol" &>/dev/null 2>&1; then
    PYMOL_VER=$(python -c "import pymol; print(pymol.__version__)" 2>/dev/null || echo "?")
    ok "pymol  –  residue capping + PDB→MOL2 conversion  ($PYMOL_VER)"
else
    fail "pymol  –  residue capping + PDB→MOL2 conversion  NOT importable"
fi

# pexpect – needed for RESP; warn rather than fail if absent
if python -c "import pexpect" &>/dev/null 2>&1; then
    PEXPECT_VER=$(python -c "import pexpect; print(pexpect.__version__)" 2>/dev/null || echo "?")
    ok "pexpect  –  Multiwfn automation (RESP pathway)  ($PEXPECT_VER)"
else
    warn "pexpect  –  NOT importable (only needed for --charge resp)"
fi

# ── 3. Check optional RESP tools ────────────────────────────
echo -e "\n${BLD}Checking optional RESP tools:${RST}"

check_resp_exe() {
    local envvar="$1"
    local name="$2"
    local purpose="$3"
    local resolved="${!envvar:-}"
    if [[ -n "$resolved" ]]; then
        if [[ -x "$resolved" ]]; then
            ok "$name  –  $purpose  (via $envvar)"
        else
            fail "$name  –  $purpose  (${envvar}=${resolved} not executable)"
        fi
    elif command -v "$name" &>/dev/null; then
        ok "$name  –  $purpose  (found on PATH)"
    else
        warn "$name  –  $purpose  NOT FOUND  (set $envvar or add to PATH)"
    fi
}

check_resp_exe "NSAA_XTB_EXE"       "xtb"            "semi-empirical pre-optimiser"
check_resp_exe "NSAA_ORCA_EXE"      "orca"           "HF geometry optimisation"
check_resp_exe "NSAA_ORCA_2MKL_EXE" "orca_2mkl"      "GBW→Molden wavefunction converter"
check_resp_exe "NSAA_MULTIWFN_EXE"  "Multiwfn_noGUI" "RESP electrostatic-potential fitting"

# ── 4. Parameter file sanity check ──────────────────────────
echo -e "\n${BLD}Checking AMBER parameter files:${RST}"

for parm in parm19.dat parm10.dat gaff2.dat gaff.dat; do
    PARM_PATH="$AMBERHOME/dat/leap/parm/$parm"
    if [[ -f "$PARM_PATH" ]]; then
        ok "$parm"
    else
        fail "$parm  NOT FOUND at $PARM_PATH"
    fi
done

# ── 5. Summary ──────────────────────────────────────────────
echo ""
echo "────────────────────────────────────────────────"
if [[ $ERRORS -eq 0 ]]; then
    echo -e "${GRN}${BLD}All required tools verified.${RST}"
    echo -e "You can now run:  ${BLD}python main.py -p peptide.mol2${RST}"
else
    echo -e "${RED}${BLD}$ERRORS required tool(s) missing.${RST}"
    echo "  Ensure ambertools and pymol-open-source installed:"
    echo "    conda env update -f environment.yml --prune"
fi

echo ""
echo -e "${YLW}RESP pathway (--charge resp):${RST}"
echo "  Requires ORCA and Multiwfn, which must be installed manually."
echo "  See environment.yml for download links and env-var instructions."
echo ""
