# Covalent Ligand Dataset for Protein-Ligand Docking

This repository contains tools and datasets for generating and curating covalent protein-ligand complexes from the Protein Data Bank (PDB). The dataset is designed for benchmarking covalent docking methods and structure prediction tools by removing redundancy and train-test data leakage.

## Table of Contents

- [Overview](#overview)
- [Dataset Description](#dataset-description)
- [JSON File Formats](#json-file-formats)
- [Installation](#installation)
- [Usage](#usage)

## Overview

This project provides a curated dataset of covalent protein-ligand complexes extracted from the PDB, along with tools to:

1. **Extract covalent ligand structures** from PDB with quality filters
2. **Remove train-test leakage** based on sequence and ligand similarity
3. **Cluster and deduplicate** structures to create non-redundant datasets
4. **Generate cross-docking pairs** for evaluating cross-docking scinario
5. **Prepare input files** for various docking and structure prediction tools (GalaxyCDock, AutoDock4, DOCK6)

## Dataset Description

### Training Datasets

1. **Scarpino Set** (`train_data/scarpino_set.json`)
   - Covalent complexes from the Scarpino et al. benchmark
   - Only Cystein-targetting ligands
   - Pre-2018 structures
   - Used as training/reference set

2. **PDBbind Data** (`train_data/pdbbind_data.json`)
   - Non-covalent protein-ligand complexes from PDBbind
   - Used as training set of GalaxyDock-DL, which is utilized for rescoring in GalaxyCDock web server.

### Test Datasets

1. **Covalent Test Set** (`test_data/cov_set.json`)
   - Newly curated covalent complexes (pre-2025)
   - Filtered for quality and non-redundancy within the test set
   - Remove redundancy from the Scarpino set
   - Extended Target residues including Cys, Ser, Lys, Tyr, and His
   - Main benchmark set

2. **Cross-Docking Set** (`cross_set/cross_covset.json`)
   - Pairs of structures with same protein but different ligands
   - Used for cross-docking validation
   - Includes ligand similarity cutoff of maximum 0.5

### Directory Structure

``
cdock_dataset/
├── README.md                          # This file
├── requirements.txt                   # Python dependencies
├── check_binding_site.py             # Filter structures by binding site quality
├── clustering.py                     # Clustering and leakage removal
├── gen_testset.py                    # Generate test set from PDB
├── gen_input.py                      # Generate input files for docking tools
├── test_data/
│   ├── cov_set.json                  # Main covalent test set
│   ├── covset_seq.fasta              # Protein sequences for test set
│   └── cov_rcsb_19950101_20250919.txt # PDB IDs list (downloaded from RCSB PDB)
├── train_data/
│   ├── scarpino_set.json             
│   ├── pdbbind_data.json             
│   ├── train_seq.fasta               # Training set sequences of the Scarpino and PDBbind sets.
│   └── scarpino.txt                  # Scarpino PDB IDs
└── cross_set/
    ├── cross_covset.json             # Cross-docking pairs of the main test set.
    └── cross_scarpino.json           # Scarpino cross-docking pairs
```

## JSON File Formats

### Covalent Ligand Dataset Format

Files: `cov_set.json`, `scarpino_set.json`

Each entry represents a covalent protein-ligand complex:

```json
{
  "pdb": "1b12",                    // PDB ID (4-character code)
  "date": "1999-11-24",             // Deposition date (YYYY-MM-DD)
  "exp_type": "X-RAY DIFFRACTION",  // Experimental method
  "resolution": 1.95,                // Resolution in Angstroms
  "uniprot_id": "P00803",           // UniProt accession of protein
  
  // Protein residue information (covalently bonded residue)
  "res_name": "SER",                // Residue name (3-letter code)
  "res_atom": "OG",                 // Atom name forming covalent bond
  "res_atom_alt_id": "?",           // Alternate location indicator
  "res_num": "15",                  // Residue number (label)
  "res_num_auth": "90",             // Residue number (author)
  "res_chain": "A",                 // Chain ID (label)
  "res_chain_auth": "A",            // Chain ID (author)
  
  // Ligand information
  "lig_name": "1PN",                // Ligand name (CCD code)
  "lig_atom": "C7",                 // Atom name forming covalent bond
  "lig_atom_alt_id": "?",           // Alternate location indicator
  "lig_num": ".",                   // Ligand number (label)
  "lig_num_auth": "1001",           // Ligand number (author)
  "lig_chain": "E",                 // Chain ID (label)
  "lig_chain_auth": "A",            // Chain ID (author)
  "lig_atom_num": 20,               // Number of heavy atoms in ligand
  "lig_n_tor": 7,                   // Number of rotatable bonds
  "smiles": "C[C@@H](OC(C)=O)..."   // SMILES string of ligand
}
```

### PDBbind Dataset Format

File: `pdbbind_data.json`

Non-covalent complexes with binding affinity data:

```json
{
  "pdb": "2tpi",                    // PDB ID
  "lig_name": "2-mer",              // Ligand name
  "resolution": "2.10",             // Resolution (Angstroms)
  "year": 1982,                     // Deposition year
  "smiles": "CC[C@H](C)[C@H]...",   // SMILES string
  "affinity": 4.309803919971486,    // -log(affinity) value
  "affinity_type": "Kd",            // Affinity type (Kd, Ki, IC50)
  "affinity_sign": "=",             // Equality sign (=, >, <, ~)
  "uniprot": "P00760"               // UniProt ID (legacy field name)
}
```


### Cross-Docking Set Format

Files: `cross_covset.json`, `cross_scarpino.json`

Pairs for cross-docking experiments:

```json
{
  "lig_id": "5ydm_DUW_B",           // Ligand structure ID (pdb_ligname_chain)
  "rec_id": "5ydl_DUV_A",           // Receptor structure ID
  "uniprot_id": "A0A0E3JLZ0",       // UniProt ID (same protein)
  "tanimoto": 0.3611111111,         // Ligand Tanimoto similarity
  "source": "additional"            // Source ("existing" or "additional")
}
```

**Key Fields:**

- **lig_id**: Structure to use as ligand (format: `{pdb}_{lig_name}_{chain}`)
- **rec_id**: Structure to use as receptor
- **tanimoto**: Ligand similarity (0-1, higher = more similar)
- **source**: 
  - "existing": Both structures in original covalent dataset. The receptor structure is in the bound-form with a covalent ligand.
  - "additional": The receptor structure is in the bound-form with non-covalent ligand

**Usage**: Dock ligand from `lig_id` structure into protein from `rec_id` structure

## Installation

### Requirements

- Python 3.8+
- MMseqs2 (for sequence clustering)
- OpenBabel (for molecular format conversion)

### Setup

1. Clone the repository:
```bash
git clone https://github.com/yourusername/cdock_dataset.git
cd cdock_dataset
```

2. Install Python dependencies:
```bash
pip install -r requirements.txt
```

3. Install external tools:

**MMseqs2** (sequence clustering):
```bash
# Ubuntu/Debian
sudo apt-get install mmseqs2

# macOS
brew install mmseqs2

# From source
git clone https://github.com/soedinglab/MMseqs2.git
cd MMseqs2
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..
make && make install
```

**OpenBabel** (molecular format conversion):
```bash
# Ubuntu/Debian
sudo apt-get install openbabel

# macOS
brew install open-babel

# Conda
conda install -c conda-forge openbabel
```

## Usage

### 1. Generate Covalent Test Set from PDB

Extract covalent ligands from PDB files:

```bash
python gen_testset.py \
  --input-list test_data/cov_rcsb_19950101_20250919.txt \
  --output test_data/cov_set_raw.json \
  --cif-dir /path/to/cif/files \
  --min-date 2018-01-01 \
  --max-resolution 2.5
```

**Options:**
- `--input-list`: Text file with PDB IDs (one per line) downloaded from RCSB PDB
- `--output`: Output JSON file path
- `--cif-dir`: Directory containing mmCIF files
- `--min-date`: Minimum deposition date (YYYY-MM-DD)
- `--max-resolution`: Maximum resolution in Angstroms

### 2. Filter by Binding Site Quality

Remove structures with poor binding site geometry:

```bash
python check_binding_site.py \
  --input test_data/cov_set_raw.json \
  --output test_data/cov_set_filtered.json \
  --cif-dir /path/to/cif/files \
  --cutoff 0.65
```

**Options:**
- `--cutoff`: Minimum contact ratio (default: 0.65)
  - Higher values = stricter quality filter
  - Recommended range: 0.6-0.75

### 3. Remove Train-Test Leakage

Filter test set to avoid similar proteins/ligands in training set:

```bash
python clustering.py remove-leakage \
  --test-fasta test_data/covset_seq.fasta \
  --train-fasta train_data/train_seq.fasta \
  --output test_data/cov_set_no_leak.json \
  --seq-threshold 0.3 \
  --lig-threshold 0.5
```

**Options:**
- `--seq-threshold`: Maximum sequence identity (default: 0.3)
- `--lig-threshold`: Maximum Tanimoto similarity (default: 0.5)

### 4. Remove Redundancy

Create non-redundant dataset by clustering:

```bash
python clustering.py remove-redundancy \
  --input-fasta test_data/covset_seq.fasta \
  --input-json test_data/cov_set_no_leak.json \
  --output test_data/cov_set_final.json \
  --seq-threshold 0.3 \
  --lig-threshold 0.5
```

**Options:**
- `--seq-threshold`: Sequence identity for clustering (default: 0.3)
- `--lig-threshold`: Ligand similarity for clustering (default: 0.5)

### 5. Generate Input Files for Docking

Prepare structure files for docking programs:

```bash
python gen_input.py \
  --input test_data/cov_set_final.json \
  --cif-dir /path/to/cif/files \
  --output-dir docking_inputs
```

This generates:
- Protein PDB/CIF files
- Ligand MOL2 files
- Configuration files for AutoDock4, DOCK6, etc.

## Workflow

### Complete Pipeline for Dataset Curation

```bash
# 1. Extract covalent complexes from PDB
python gen_testset.py \
  --input-list test_data/pdb_list.txt \
  --output test_data/cov_set_raw.json \
  --cif-dir cif_files/ \
  --min-date 2018-01-01

# 2. Filter by binding site quality
python check_binding_site.py \
  --input test_data/cov_set_raw.json \
  --output test_data/cov_set_filtered.json \
  --cif-dir cif_files/ \
  --cutoff 0.65

# 3. Remove train-test leakage
python clustering_refactored.py remove-leakage \
  --test-fasta test_data/covset_seq.fasta \
  --train-fasta train_data/train_seq.fasta \
  --output test_data/cov_set_no_leak.json

# 4. Remove redundancy
python clustering_refactored.py remove-redundancy \
  --input-fasta test_data/covset_seq.fasta \
  --input-json test_data/cov_set_no_leak.json \
  --output test_data/cov_set_final.json

# 5. Generate docking inputs
python gen_input.py \
  --input test_data/cov_set_final.json \
  --cif-dir cif_files/ \
  --output-dir docking_inputs/
```

## Dataset Statistics

### Covalent Test Set (`cov_set.json`)

- **Total structures**: 647 complexes
- **Date range**: 1995-2025
- **Resolution**: ≤ 2.5 Å
- **Experiment types**: X-ray crystallography, Cryo-EM
- **Covalent residues**: CYS (55%), SER (24%), LYS (18%), others (3%)
- **Unique proteins**: 346 (based on UniProt ID)


## Quality Filters

### Structural Quality
- X-ray or Cryo-EM, resolution ≤ 2.5 Å
- Single covalent bond per ligand
- Binding site contact ratio ≥ 0.65

### Chemical Properties (Ligands)
- Heavy atom count: 8-50 atoms
- Molecular weight: ≤ 700 Da
- Rotatable bonds: ≤ 20
- cLogP: 0-7.5
- H-bond donors: ≤ 5
- H-bond acceptors: ≤ 10
- Polar surface area: ≤ 200 Ų
- Allowed elements: H, C, N, O, F, P, S, Cl, Br, I

### Redundancy Removal
- Sequence identity: < 30%
- Ligand Tanimoto similarity: < 0.5
- Same protein + same ligand → keep best resolution

## File Naming Conventions

### Structure IDs

Format: `{pdb}_{lig_name}_{lig_chain}`

Examples:
- `5ydm_DUW_B`: PDB 5ydm, ligand DUW, chain B
- `1b12_1PN_A`: PDB 1b12, ligand 1PN, chain A

### Generated Files

From `gen_input.py`:
- `{pdb}_{lig}_{chain}_prot.pdb`: Protein structure
- `{pdb}_{lig}_{chain}_lig.mol2`: Ligand structure
- `{pdb}_{lig}_{chain}_lig_ca.mol2`: Ligand + CA-CB of covalent residue (for GalaxyCDock input)
- `{pdb}_{lig}_{chain}_comp.cif`: Complete complex