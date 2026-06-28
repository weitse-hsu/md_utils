# md_utils
[![CircleCI](https://dl.circleci.com/status-badge/img/gh/weitse-hsu/md_utils/tree/main.svg?style=shield)](https://dl.circleci.com/status-badge/redirect/gh/weitse-hsu/md_utils/tree/main)
[![codecov](https://codecov.io/gh/weitse-hsu/md_utils/graph/badge.svg?token=PgncVdVkDa)](https://codecov.io/gh/weitse-hsu/md_utils)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

A collection of utility functions and command-line tools for GROMACS-based MD simulations.

## Installation

```bash
pip install git+https://github.com/weitse-hsu/general_utils.git
pip install .
```

## Command-line tools

### Simulation

| Command | Description |
|---|---|
| `prep_simulation` | Prepare a GROMACS production simulation from input files |
| `process_gmx_traj` | Process a GROMACS trajectory |
| `create_index_grps` | Create GROMACS index groups from a selections file |

For generating GROMACS-compatible topologies, use the standalone [acpype](https://github.com/alanwilter/acpype) package:

```bash
conda install -c conda-forge acpype  # or: pip install acpype
```

### Structure

| Command | Description |
|---|---|
| `identify_pocket` | Identify pocket residues around a ligand |
| `identify_ss_residues` | Identify secondary structure domains in a protein structure |
| `map_residues` | Map residue numbers between a reference PDB and a simulation GRO file |

Run any command with `--help` for full usage details.

## Authors

- Wei-Tse Hsu, University of Oxford (wei-tse.hsu@bioch.ox.ac.uk)
