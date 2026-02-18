"""Configuration constants for covalent ligand dataset generation.

This module defines various constants and filters used throughout the dataset
generation pipeline, including:
- Experimental quality thresholds
- Chemical property filters
- Biologically relevant/irrelevant residues
- Covalent bond definitions
"""

import datetime

# ============================================================================
# Experimental Quality Filters
# ============================================================================

# Minimum deposition date for structures (ISO format: YYYY-MM-DD)
DEPOSIT_DATE_CUTOFF = datetime.date(2018, 1, 1)

# Allowed experimental methods
EXP_TYPES = ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"]

# Maximum resolution for structures (Angstroms)
RESOLUTION_CUTOFF = 2.5

# Maximum covalent bond length (Angstroms)
BOND_LENGTH_CUTOFF = 2.0

# Binding site distance cutoff (Angstroms)
BINDING_SITE_CUTOFF = 6.0

# Minimum protein chain length (to exclude peptides)
MIN_CHAIN_LENGTH = 20


# ============================================================================
# Ligand Property Filters (Drug-like properties)
# ============================================================================

# Ligand size constraints
MIN_HEAVY_ATOMS = 8
MAX_HEAVY_ATOMS = 50

# Molecular weight limit (Da)
MAX_MOLECULAR_WEIGHT = 700

# Rotatable bonds (flexibility)
MAX_ROTATABLE_BONDS = 20

# Lipophilicity (cLogP range)
MIN_CLOGP = 0.0
MAX_CLOGP = 7.5

# Hydrogen bonding capacity
MAX_H_BOND_DONORS = 5
MAX_H_BOND_ACCEPTORS = 10

# Polar surface area (Ų)
MAX_POLAR_SURFACE_AREA = 200

# Allowed chemical elements in ligands
ALLOWED_ATOMS = [
    "H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"
]


# ============================================================================
# Amino Acid Definitions
# ============================================================================

# Standard canonical amino acids
CANONICAL_RESIDUES = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL"
]

# Non-canonical amino acids (modified residues)
NON_CANONICAL_RESIDUES = [
    "CXM",  # S-carboxymethyl cysteine
    "CME",  # S-methyl cysteine
    "PTR",  # O-phosphotyrosine
    "HYP",  # Hydroxyproline
    "CSO",  # S-hydroxycysteine
    "MSE",  # Selenomethionine
    "MLZ",  # N-methyl lysine
    "TPO",  # Phosphothreonine
    "MEN",  # N-methyl asparagine
    "3CT",  # 3-carboxytyrosine
    "SNN",  # 
    "FTR",  # 
    "9IJ",  # 
    "SEP",  # Phosphoserine
    "FME",  # N-formyl methionine
    "AYA",  # N-acetylalanine
    "KCX",  # Lysine carboxylic acid
    "TYS",  # O-sulfotyrosine
    "PCA",  # Pyroglutamic acid
    "PHI",  # Iodophenylalanine
]

# Amino acids that can form covalent bonds with ligands
COVALENT_RESIDUES = [
    "CYS",  # Cysteine (thiol group)
    "SER",  # Serine (hydroxyl group)
    "LYS",  # Lysine (amine group)
    "TYR",  # Tyrosine (phenol group)
    "HIS",  # Histidine (imidazole)
    "THR",  # Threonine (hydroxyl group)
    "ASP",  # Aspartic acid (carboxyl)
    "GLU",  # Glutamic acid (carboxyl)
]


# ============================================================================
# Cofactors and Biologically Irrelevant Molecules
# ============================================================================

# Common cofactors (should be kept in structures)
COFACTORS = [
    "FAD",  # Flavin adenine dinucleotide
    "NAD",  # Nicotinamide adenine dinucleotide
    "NAP",  # NADP+
    "HEM",  # Heme
    "FMN",  # Flavin mononucleotide
    "COA",  # Coenzyme A
    "GSH",  # Glutathione
    "THF",  # Tetrahydrofolate
    "SAM",  # S-adenosylmethionine
    "PLP",  # Pyridoxal phosphate
    "ADP",  # Adenosine diphosphate
    "ATP",  # Adenosine triphosphate
    "GDP",  # Guanosine diphosphate
    "MG",   # Magnesium ion
    "CA",   # Calcium ion
    "ZN",   # Zinc ion
    "FE",   # Iron ion
    "CU",   # Copper ion
    "MN",   # Manganese ion
    "NI",   # Nickel ion
    "CL",   # Chloride ion
    "K",    # Potassium ion
    "NA",   # Sodium ion
    "NAG",  # N-acetyl-D-glucosamine
    "FE2",  # Fe2+ ion
    "ANP",  # Phosphoaminophosphonic acid-adenylate ester
    "SAH",  # S-adenosyl-L-homocysteine
]

# Biologically irrelevant molecules (crystallization artifacts, solvents)
BIOLOGICALLY_IRRELEVANT = [
    "HOH",  # Water
    "EDO",  # Ethylene glycol
    "WAT",  # Water
    "CL",   # Chloride
    "SO4",  # Sulfate
    "DMS",  # Dimethyl sulfoxide
    "PO4",  # Phosphate
    "NO3",  # Nitrate
    "SCN",  # Thiocyanate
    "UNX",  # Unknown
    "UNL",  # Unknown ligand
    "GOL",  # Glycerol
    "ACT",  # Acetate
    "PEG",  # Polyethylene glycol
    "PGE",  # Triethylene glycol
    "MLI",  # Malonate ion
    "CO",   # Cobalt
    "IMD",  # Imidazole
    "YT3",  # 
    "TSL",  # 
    "ACN",  # Acetone
    "MES",  # 2-(N-morpholino)-ethanesulfonic acid
    "MPD",  # 2-methyl-2,4-pentanediol
    "MYR",  # Myristic acid
    "NTK",  # 
    "1PE",  # Pentaethylene glycol
    "BR",   # Bromide
    "FMT",  # Formate
    "TLA",  # L-tartaric acid
    "P6G",  # Hexaethylene glycol
    "PG6",  # Hexaethylene glycol
    "PG4",  # Tetraethylene glycol
    "CIT",  # Citrate
]

# Ligands to exclude (problematic or non-relevant)
EXCLUDE_LIGANDS = [
    "PLM",  # Phospholipid
    "MYR",  # Myristic acid (lipid)
    "DEP",  # Deposit
    "CXM",  # Terminally modified peptide
    "MK7",  # Free ligand form overlapping in PDB (6qw9)
    "W6Z",  # Free ligand form overlapping in PDB (6zam)
]


# ============================================================================
# Atom Connectivity (for covalent bond handling)
# ============================================================================

# Maps atoms to their neighboring backbone/sidechain atoms
# Used to include appropriate context atoms when extracting covalent bonds
RESIDUE_ATOM_NEIGHBORS = {
    "CB": "CA",   # Beta carbon connects to alpha carbon
    "OG": "CB",   # Serine hydroxyl
    "SG": "CB",   # Cysteine thiol
    "ND1": "CE1", # Histidine
    "NE2": "CE1", # Histidine
    "CE": "CD",   # Lysine
    "NZ": "CE",   # Lysine amine (can be linked)
    "OE1": "CD",  # Glutamate
    "OE2": "CD",  # Glutamate
    "OG1": "CB",  # Threonine
    "OH": "CE",   # Tyrosine
    "OD1": "CG",  # Aspartate
    "OD2": "CG",  # Aspartate
}


# ============================================================================
# Clustering and Similarity Thresholds
# ============================================================================

# Sequence similarity threshold for leakage detection
SEQ_IDENTITY_THRESHOLD = 0.3  # 30%

# Sequence coverage threshold for MMseqs2
SEQ_COVERAGE_THRESHOLD = 0.8  # 80%

# Ligand Tanimoto similarity threshold
LIGAND_SIMILARITY_THRESHOLD = 0.5  # 50%

# Binding site contact ratio threshold
CONTACT_RATIO_THRESHOLD = 0.65


# ============================================================================
# Time-based Splits (for model training cutoffs)
# ============================================================================

# AlphaFold 3 training data cutoff
AF3_DATE_CUTOFF = datetime.date(2021, 9, 30)

# Boltz-2 training data cutoff
BOLTZ2_DATE_CUTOFF = datetime.date(2023, 6, 1)


# ============================================================================
# Manual Corrections
# ============================================================================

# Manual UniProt ID corrections for structures with annotation errors
UNIPROT_MANUAL_CORRECTIONS = {
    "2gko": "Q45681",
    # Add more as needed
}

# Ligands with RDKit torsion calculation errors
# Manual torsion counts for problematic ligands
RDKIT_ERROR_LIGAND_TORSIONS = {
    "VS3": 18,
    "GW9": 4,
    # Add more as needed
}


# ============================================================================
# File Paths (can be overridden)
# ============================================================================

# Default directories (relative to project root)
DEFAULT_CIF_DIR = "../../inputs/cif"
DEFAULT_MOL2_DIR = "../../inputs/mol2"
DEFAULT_PDB_DIR = "../../inputs/pdb"
DEFAULT_FASTA_DIR = "../../inputs/fasta"
