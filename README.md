# TopFor

AMBER parameter generation for **noncanonical amino acids (ncAAs)**.

Given a residue MOL2 / PDB (or a whole peptide), `TopFor` produces a
complete AMBER-ready parameter set:

```
residue.mol2  ->  residue_capped.mol2  ->  RES.mol2 (with charges)
                                      ->  RES.ac
                                      ->  RES.mc
                                      ->  RES.prepin
                                      ->  RES_ff19SB.frcmod
                                      ->  RES_gaff2.frcmod
                                      ->  RES.lib
```

It also writes a single consolidated `residue_meta.json` per residue
documenting every decision the pipeline made (head/tail atoms, OXT
preservation, net charge source, applied caps, etc.).

---

## Installation

The detailed installation and environment setup guide is available here:

[Installation and environment setup guide](Installation/environment_setup.md)

The pipeline shells out to a number of standard tools. Make sure these are
on `PATH` (or set the relevant environment variables for the RESP backend,
see below):

| Tool                | Required for                          |
|---------------------|---------------------------------------|
| `antechamber`       | always                                |
| `prepgen`           | always (AMBER toolchain stage)        |
| `parmchk2`          | always (AMBER toolchain stage)        |
| `tleap`             | always (AMBER toolchain stage)        |
| `pymol`             | always (capping + PDB->MOL2)          |
| `xtb`               | only when `--charge resp`             |
| `orca`, `orca_2mkl` | only when `--charge resp`             |
| `Multiwfn_noGUI`    | only when `--charge resp`             |
| Python 3.10+        | always                                |

`AMBERHOME` must be set so the parameter database can be located.

---
## Quick start

### Single residue

```bash
topfor -i MVA.mol2
```

### Several residues at once (shell glob)

`-i` accepts multiple files, so the usual `*.mol2` shell glob works:

```bash
topfor -i *.mol2
```

### Free amino acid (e.g. for a fragment that already has OXT and a free amine)

```bash
topfor -i MVA.mol2 --terminal both
```

### C-terminal residue of a peptide (keep OXT, skip NME cap)

```bash
topfor -i MVA.mol2 --terminal c
```

### N-terminal residue of a peptide (keep free amine, skip ACE cap)

```bash
topfor -i MVA.mol2 --terminal n
```

### Peptide mode

Extract every non-standard residue from a peptide and parameterize each.
Terminal residues are detected automatically: the C-terminal residue keeps its OXT, the N-terminal residue keeps its free amine, and internal residues are capped with ACE/NME for charge derivation.

```bash
topfor -p peptide.pdb
```

### Batch mode

`-b` accepts any combination of:

- a directory (all `.mol2` and `.pdb` inside it are processed),
- a plain-text list file (one path per line, `#` for comments),
- a glob pattern (quoted, so the shell does not pre-expand it).

```bash
topfor -b residues/
topfor -b mol2_files.txt
topfor -b "data/*.mol2"
topfor -b residues/ "extra/*.mol2"     # multiple args concatenated
```

### Choosing a charge model

```bash
topfor -i MVA.mol2 -c abcg2     # default
topfor -i MVA.mol2 -c bcc
topfor -i MVA.mol2 -c gas
topfor -i MVA.mol2 -c resp      # needs xTB + ORCA + Multiwfn
```

### Choosing force fields

```bash
topfor -i MVA.mol2 -bb ff19SB -sc gaff2   # default
topfor -i MVA.mol2 -bb ff14SB -sc gaff
```

### Residue map

Hand-curated overrides for tricky residues live in a JSON file passed via
`--map`. A working example is at `examples/residue_map.json`:

```json
{
  "MVA": {
    "head": "N",
    "tail": "C",
    "mainchain": ["CA"],
    "net_charge": 0,
    "pre_head_type": "C",
    "post_tail_type": "N"
  },
  "PCA": {
    "head": "NONE",
    "tail": "C",
    "net_charge": 0
  }
}
```

Keys are residue names (case-insensitive on input, upper-cased internally).
Values are dicts with any of the following optional fields:

| Field            | Meaning                                                          |
|------------------|------------------------------------------------------------------|
| `head`           | Head atom name. Use `"NONE"` to skip ACE capping.                |
| `tail`           | Tail atom name. Use `"NONE"` to skip NME capping.                |
| `mainchain`      | Ordered list of atom names between head and tail.                |
| `net_charge`     | Integer net charge for this residue. Overrides auto-detection.   |
| `pre_head_type`  | AMBER atom-type symbol expected upstream of `head` (default "C").|
| `post_tail_type` | AMBER atom-type symbol expected downstream of `tail` (default "N"). |

Anything missing falls back to topology detection + connectivity inference,
so the map only needs to contain real exceptions.

Run the tool with the map:

```bash
topfor -i MVA.mol2 --map examples/residue_map.json
```

---

## Output structure

For each input residue, a folder is created under `--out/-o` (default: current
directory):

```
out/
├── MVA/
│   ├── residue.mol2            # exact copy / PyMOL conversion of the input
│   ├── residue_capped.mol2     # post-PyMOL capping (ACE/NME or OXT preserved)
│   ├── MVA.mol2                # with partial charges
│   ├── MVA.ac
│   ├── MVA.mc
│   ├── MVA.prepin
│   ├── MVA_ff19SB.frcmod
│   ├── MVA_gaff2.frcmod
│   ├── MVA.lib
│   ├── MVA.log                 # AMBER toolchain log
│   └── residue_meta.json       # single consolidated metadata file
├── failed/                     # created only if there were failures
│   └── BAD/                    # full working folder of any failed residue
└── successful_residues.txt
└── failed_residues.txt
```

`successful_residues.txt` and `failed_residues.txt` each contain one residue
name per line, suitable for downstream scripting.

The single consolidated `residue_meta.json` includes (among others):

- `RES`, `head_name`, `tail_name`, `main_chain`, `net_charge`,
  `charge_model`, `topology`, `applied_caps`
- `preserve_oxt` / `oxt_preserved` — whether OXT was kept for a terminal residue
- `terminal_mode` — the `--terminal` flag that was active
- `is_polymer_internal` / `is_n_terminal_like` / `is_c_terminal_like` —
  peptide-mode topology classification
- `charge_backend_meta` — paths and provenance from the charge-assignment stage
- `validation_warnings` — any sanity-check messages from the molecule validator

---

## Why `--terminal` matters

C-terminal residues of a linear peptide carry an **OXT** atom that is needed
for proper hydrogen bonding and to match the chemical reality of the system.
A bare residue extracted from such a peptide must therefore **keep** its
OXT and must **not** be NME-capped, adding NME on top of OXT produces a
geometrically broken topology, and removing OXT to attach NME destroys
chemistry that the parameter set is supposed to represent.

In peptide mode (`-p`), the splitter detects this automatically and tells
the capping step to preserve OXT. In single-residue mode (`-i`), the user
must say so explicitly with `--terminal c` (C-terminal) or `--terminal both`
(free amino acid). For N-terminal residues, `--terminal n` skips ACE capping
so the residue's free amine is preserved.

---

## Module map

```
topfor/
├── examples/
│   └── residue_map.json           map schema + worked example
├── Installation/
│   ├── environment.yml
│   ├── environment_setup.md
│   └── setup_env.sh
├── modules/
│   ├── __init__.py
│   ├── residue_processor.py       Stage 1 (capping + charges)
│   ├── peptide_splitter.py        peptide -> per-residue extractor
│   ├── capping.py                 PyMOL ACE/NME capping (OXT-aware)
│   ├── pdb_to_mol2.py             PyMOL PDB -> MOL2 helper
│   ├── prepgen_writer.py          prepgen .mc control file writer
│   ├── antechamber_runner.py      Stage 2 (AMBER toolchain)
│   ├── resp_workflow.py           xTB + ORCA + Multiwfn RESP backend
│   └── mol2_utils.py              MOL2 parsing, validation, helpers                      
├── test/ 
│   ├── MVA.mol2                      Working example of N-methylated Valine
│   └── MVA/                          Parameters and topology folder  
├── main.py
├── README.md   
└── topfor
```
