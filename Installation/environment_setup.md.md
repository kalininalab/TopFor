# TopFor Environment Setup

This file explains how to set up the computer environment needed to run **TopFor**.

Do not worry if you are new to Conda, AmberTools, PyMOL, ORCA, Multiwfn, or terminal commands. Follow the steps slowly and in the same order.

---

## 1. What this setup does

TopFor needs some other tools to work.

For normal TopFor runs, you need:

- **Python** — runs the TopFor code
- **AmberTools** — makes AMBER files like `.prepin`, `.frcmod`, and `.lib`
- **PyMOL** — helps with capping and PDB to MOL2 conversion

For RESP charge calculation, you also need:

- **xTB** — makes the starting structure better before ORCA
- **ORCA** — runs the quantum chemistry calculation
- **orca_2mkl** — comes with ORCA and converts ORCA files for Multiwfn
- **Multiwfn_noGUI** — fits the RESP charges
- **pexpect** — lets Python control Multiwfn automatically

Most normal tools are installed by:

```bash
environment.yml
```

The setup check is done by:

```bash
setup_env.sh
```

The `topfor` command is provided by the file named:

```bash
topfor
```

That `topfor` file is a small launcher. It runs `main.py` for you, so users do **not** need to type `python main.py`.

---

## 2. Files you should have

Your TopFor folder should look like this:

```text
TopFor/
├── topfor
├── Installation/
│   ├── environment.yml
│   ├── README.md
│   └── setup_env.sh
└── other topfor files
    
```
---

## 3. Install Conda first

Before creating the TopFor environment, install Conda.

You can use one of these:

- Miniconda
- Anaconda
- Miniforge
- Mambaforge

After installing it, open a terminal and check:

```bash
conda --version
```

If it prints a version number, Conda is working.

---

## 4. Create the TopFor environment

Go to the `Installation/` folder:

```bash
cd path/to/TopFor/Installation
```

Create the Conda environment:

```bash
conda env create -f environment.yml
```

This may take some time because AmberTools and PyMOL are large.

Now activate the environment:

```bash
conda activate topfor
```
---

## 5. Run the setup script

While you are still inside the `Installation/` folder, run:

```bash
bash setup_env.sh
```

This script does these checks:

1. It sets `AMBERHOME`.
2. It checks AmberTools.
3. It checks PyMOL.
4. It checks optional RESP tools.
5. It checks AMBER parameter files.

You only need to run this once after creating the environment.

After that, every time you run:

```bash
conda activate topfor
```

`AMBERHOME` should be set automatically.

---

## 6. Check if the normal setup is ready

For normal TopFor usage, this is enough:

```text
✓ antechamber
✓ prepgen
✓ parmchk2
✓ tleap
✓ pymol
```

If the script says:

```text
All required tools verified.
```

then the normal TopFor setup is ready.

Warnings about ORCA or Multiwfn are okay if you are **not** using RESP charges.

---

## 7. Make the `topfor` command work

TopFor has a file named:

```bash
topfor
```

This file is the command launcher.

### Step 1: go back to the main TopFor folder

If you are inside `Installation/`, go one folder top:

```bash
cd ..
```

Now you should be inside the main `TopFor/` folder.

### Step 2: make the launcher executable

Run:

```bash
chmod +x topfor
```

### Step 3: add TopFor to PATH for this terminal

Run this from the main `TopFor/` folder:

```bash
export PATH="$(pwd):$PATH"
```

Now you can run:

```bash
topfor --help
```
---

## 8. Important note about RESP charges

Normal TopFor charge methods work with the Conda environment.

RESP is different.

If you use:

```bash
-c resp
```

or:

```bash
--charge resp
```

then TopFor needs ORCA and Multiwfn.

These are not installed fully by Conda, so you must install them manually.

For RESP, TopFor checks these commands or paths:

| Tool | Needed for |
|---|---|
| `xtb` | first quick structure improvement |
| `orca` | quantum chemistry calculation |
| `orca_2mkl` | converts ORCA output for Multiwfn |
| `Multiwfn_noGUI` | RESP charge fitting |
| `pexpect` | Python automation of Multiwfn |

`xTB` and `pexpect` should already come from the Conda environment.

ORCA and Multiwfn must be installed manually.

---

## 9. Install ORCA for RESP

Go to the ORCA website and download ORCA.

Website:

```text
https://www.faccts.de/customer
```

Download the Linux shared OpenMPI build.

Example file name:

```text
orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg.tar.xz
```

Make a software folder:

```bash
mkdir -p "$HOME/software"
```

Move the ORCA `.tar.xz` file into that folder, or use the real path to where it was downloaded.

Extract it:

```bash
tar -xf orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg.tar.xz -C "$HOME/software"
```

After extraction, you should have a folder like:

```text
$HOME/software/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
```

Inside that folder, these files should exist:

```text
orca
orca_2mkl
```

Make sure they are executable:

```bash
chmod +x "$HOME/software/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca"
chmod +x "$HOME/software/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca_2mkl"
```

ORCA may also need OpenMPI.

On Ubuntu or Debian, you can install it with:

```bash
sudo apt install openmpi-bin libopenmpi-dev
```
---

## 10. Install Multiwfn for RESP

Go to the Multiwfn download page.

Website:

```text
http://sobereva.com/multiwfn/download.html
```

Download the Linux noGUI version.

Example file name:

```text
Multiwfn_2026.3.27_bin_Linux_noGUI.zip
```

Extract it into your software folder:

```bash
unzip Multiwfn_2026.3.27_bin_Linux_noGUI.zip -d "$HOME/software"
```

After extraction, you should have a folder like:

```text
$HOME/software/Multiwfn_2026.3.27_bin_Linux_noGUI
```

Make `Multiwfn_noGUI` executable:

```bash
chmod +x "$HOME/software/Multiwfn_2026.3.27_bin_Linux_noGUI/Multiwfn_noGUI"
```

---

## 11. Check xTB

xTB should already be installed by the Conda environment.

Check it:

```bash
conda activate topfor
which xtb
xtb --version
```

If `xtb` is not found, update the environment:

```bash
cd path/to/TopFor/Installation
conda env update -f environment.yml --prune
```

If you installed xTB manually, set this variable:

```bash
export NCAA_XTB_EXE="/path/to/xtb"
```

Use the real path to your `xtb` file.

---

## 12. Set RESP paths

TopFor needs to know where ORCA and Multiwfn are installed.

Use the real folders on your computer.

Example:

```bash
export ORCA_DIR="$HOME/software/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg"
export MULTIWFN_DIR="$HOME/software/Multiwfn_2026.3.27_bin_Linux_noGUI"

export NCAA_ORCA_EXE="$ORCA_DIR/orca"
export NCAA_ORCA_2MKL_EXE="$ORCA_DIR/orca_2mkl"
export NCAA_MULTIWFN_EXE="$MULTIWFN_DIR/Multiwfn_noGUI"

export Multiwfnpath="$MULTIWFN_DIR"

export PATH="$ORCA_DIR:$MULTIWFN_DIR:$PATH"
export LD_LIBRARY_PATH="$ORCA_DIR/lib:$LD_LIBRARY_PATH"
```

Now test:

```bash
"$NCAA_ORCA_EXE" --version
test -x "$NCAA_ORCA_2MKL_EXE" && echo "orca_2mkl found"
test -x "$NCAA_MULTIWFN_EXE" && echo "Multiwfn_noGUI found"
```

Then run the setup check again:

```bash
cd path/to/TopFor/Installation
bash setup_env.sh
```

If ORCA and Multiwfn are found, RESP is ready.

---