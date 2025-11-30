import os
import json
import datetime
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from pathlib import Path
import pandas as pd
import subprocess as sp
import requests
import numpy as np
from bs4 import BeautifulSoup

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
from rdkit.Chem.rdMolDescriptors import (
    CalcExactMolWt, 
    CalcNumRotatableBonds, 
    CalcCrippenDescriptors, 
    CalcNumHBD, 
    CalcNumHBA, 
    CalcTPSA
    )

#============
# get the Cov set
#============
# check deposit date and experiment type
#   - after 20180101
#   - after 19950101 for all_cov_set
#   - X-ray or cryoEM
#   - resolution <= 2.5
# check linkage
#   - only one covalent bond
#   - not linked to peptide : check the seq_len of linked protein
#   - exclude saccharide (linked residue name)
# save as json file
#  - pdb code
#  - deposit date
#  - experiment type
#  - resolution
#  - linked residue chain
#  - linked residue res number
#  - linked residue name
#  - linked residue atom name
#  - linked ligand chain
#  - linked ligand res number
#  - linked ligand name (CCD code)
#  - linked ligand atom name
#  - uniprot ID
#  - lig_atom_num
#  - lig_n_tor

#============
# input preparation
#============
# filtering and clustering
# - uniprotID, seq_id (from TM-score), ligand CCD code
# - remove duplicates (based on seq_id and ligand CCD code)

# make cdock(ad4,dock6) input files
# - protein pbd, ligand mol2 file (two versions: CB-SG/OG or CA-CB-SG/OG)
# make boltz, af3 input json files
# make rf3 mmcif files
 
#============
# dataset information for analysis
#============
# time split
af3_d_cutoff = 20210930
boltz2_d_cutoff = 20230601

# ligand similarity (tanimoto) to training set
# protein similarity (seq_id) to training set

DEPOSIT_DATE_CUTOFF = datetime.date(2018, 1, 1)
EXP_TYPES = ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"]
RES_CUTOFF = 2.5
MIN_CHAIN_LENGTH = 20 # to exclude peptide
BOND_LENGTH_CUTOFF = 2.0 # maximum covalent bond lenght, including disulfide bond 
BINDING_SITE_CUTOFF = 6.0
ALLOWED_ATOMS = ["H","C","N","O","F","P","S","Cl","Br",
                 "I" ## <- included in Scarpino set, unlike the author's description
                 ]
LINKED_RES = ["CYS","SER"]
STANDARD_AA = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL'
}

EXCLUDE_LIGANDS = [
    "PLM", # phospholipid
    "MYR",
    "DEP",
]

## TODO: need to calc torsion angle independent way (no RDKit)
## there might be code somewhere to calculate rotatble bonds from .mol2 file... i need to find it..
RDKIT_ERROR_LIG_TORSION = {
    "VS3":18,
    "GW9":4,
}

CIF_DIR = "../../inputs/cif"

def get_uniprot_id_from_api(pdb: str, entity_id: int):
    url = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb}"
    r = requests.get(url.format(pdb=pdb))
    data = r.json()

    uniprot = data[pdb]['UniProt']
    for k, v in uniprot.items():
        for mapping in v['mappings']:
            if mapping['entity_id'] == entity_id:
                return k


def get_field(
        mmcif_dict, field) -> list[str] | None:
    if field not in mmcif_dict:
        return None
    o = mmcif_dict[field]
    if not isinstance(o, list):
        o = [o]
    if o == '?':
        return ''
    return o


def get_seq_from_mmcif_dict(mmcif_dict, chain_id):
    # chain_id is auth_chain_id
    if '_entity_poly.pdbx_strand_id' not in mmcif_dict:
        return None, None

    entity_ids = get_field(mmcif_dict, "_entity_poly.entity_id")
    strand_ids = get_field(mmcif_dict, "_entity_poly.pdbx_strand_id")
    seq_s = get_field(mmcif_dict, "_entity_poly.pdbx_seq_one_letter_code_can")
                      
    for i, strands in enumerate(strand_ids):
        if chain_id in strands.split(','):
            return seq_s[i], int(entity_ids[i])
        
    return None, None


def get_uniprot_id_from_mmcif_dict(mmcif_dict, entity_id):
    if "_struct_ref.pdbx_db_accession" not in mmcif_dict:
        return None
    
    entity_ids = get_field(mmcif_dict, "_struct_ref.entity_id")
    uniprot_ids = get_field(mmcif_dict, "_struct_ref.pdbx_db_accession")

    for i, eid in enumerate(entity_ids):
        if int(eid) == int(entity_id):
            return uniprot_ids[i]

    return None


def check_ligand_is_saccharide(mmcif_dict, ligand_name):
    try:
        if '_chem_comp.type' in mmcif_dict and '_chem_comp.id' in mmcif_dict:
            comp_ids = get_field(mmcif_dict, '_chem_comp.id')
            comp_types = get_field(mmcif_dict, '_chem_comp.type')

            for i, comp_id in enumerate(comp_ids):
                if comp_id == ligand_name:
                    if 'saccharide' in comp_types[i].lower(): # TODO: compound names for sacharride. need to check
                        return True
    except KeyError:
        return False

    return False


def check_ligand_in_exclude_list(ligand_name):
    if ligand_name in EXCLUDE_LIGANDS:
        return True
    return False


def _get_mol_from_mmcif_data(mmcif_dict, ligand_name):
    mol = Chem.RWMol()
    atom_map = {}

    try:
        atom_ids = get_field(mmcif_dict, '_chem_comp_atom.atom_id')
        comp_ids = get_field(mmcif_dict, '_chem_comp_atom.comp_id')
        type_symbols = get_field(mmcif_dict, '_chem_comp_atom.type_symbol')
        type_symbols = [s.capitalize() for s in type_symbols]

        for i, comp in enumerate(comp_ids):
            if comp == ligand_name:
                atom_name = atom_ids[i]
                atom_symbol = type_symbols[i]
                if atom_symbol == "H":
                    continue
                new_atom = Chem.Atom(atom_symbol)
                atom_idx = mol.AddAtom(new_atom)
                atom_map[atom_name] = atom_idx

        bond_orders = {
            'sing': Chem.BondType.SINGLE,
            'doub': Chem.BondType.DOUBLE,
            'trip': Chem.BondType.TRIPLE,
            'arom': Chem.BondType.AROMATIC
        }

        bond_atom1 = get_field(mmcif_dict, '_chem_comp_bond.atom_id_1')
        bond_atom2 = get_field(mmcif_dict, '_chem_comp_bond.atom_id_2')
        bond_types = get_field(mmcif_dict, '_chem_comp_bond.value_order')
        bond_comp_ids = get_field(mmcif_dict, '_chem_comp_bond.comp_id')

        for i, comp in enumerate(bond_comp_ids):
            if comp == ligand_name:
                atom_name1 = bond_atom1[i]
                atom_name2 = bond_atom2[i]
                order = bond_types[i]
                if atom_name1 in atom_map and atom_name2 in atom_map:
                    idx1 = atom_map[atom_name1]
                    idx2 = atom_map[atom_name2]
                    bond_type = bond_orders.get(order, Chem.BondType.UNSPECIFIED)
                    # exception case : nitro group N+O2 (neutral charge, 1.5 bond)
                    # it raises valence excess error due to N with 4 bonds (C:1, O-:1, O:2)
                    # TODO: handle delocalized bonds properly... 
                    mol.AddBond(idx1, idx2, bond_type)

        Chem.SanitizeMol(mol)

    except KeyError as e:
        return None, f"Missing required mmCIF data: {e}"
    except Exception as e:
        return None, f"Failed to build RDKit Mol: {e}"

    return mol, "Success"


def check_ligand_properties(mmcif_dict, ligand_name):
    mol, mol_err = _get_mol_from_mmcif_data(mmcif_dict, ligand_name)
    if mol is None:
        print(f"Ligand {ligand_name} skipped: {mol_err}")
        return False
    
    num_atoms = mol.GetNumAtoms()
    if not (8 <= num_atoms <= 50):
        print(f"Ligand {ligand_name} skipped: number of atoms {num_atoms} out of range")
        return False
    
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    if any(atom not in ALLOWED_ATOMS for atom in atoms):
        print(f"Ligand {ligand_name} skipped: contains disallowed atoms")
        return False
    
    n_tor = CalcNumRotatableBonds(mol)
    if n_tor > 20:
        print(f"Ligand {ligand_name} skipped: number of rotatable bonds {n_tor} exceeds limit")
        return False
    
    mw = CalcExactMolWt(mol)
    if mw > 700:
        print(f"Ligand {ligand_name} skipped: molecular weight {mw} exceeds limit")
        return False
    
    clogp, _ = CalcCrippenDescriptors(mol)
    if not (0 <= clogp <= 7.5):
        print(f"Ligand {ligand_name} skipped: ClogP {clogp} out of range")
        return False
    
    h_donor = CalcNumHBD(mol)
    if h_donor > 5:
        print(f"Ligand {ligand_name} skipped: number of H-bond donors {h_donor} exceeds limit")
        return False
    
    h_acc = CalcNumHBA(mol)
    if h_acc > 10:
        print(f"Ligand {ligand_name} skipped: number of H-bond acceptors {h_acc} exceeds limit")
        return False
    
    psa = CalcTPSA(mol)
    if psa > 200:
        print(f"Ligand {ligand_name} skipped: polar surface area {psa} exceeds limit")
        return False
    
    return True


def check_covalent_bond_length(structure, bond_info):
    lig_chain = bond_info['lig_chain_auth']
    lig_res_num = int(bond_info['lig_num_auth'])
    lig_atom_name = bond_info['lig_atom']
    lig_atom_alt_id = bond_info['lig_atom_alt_id']

    prot_chain = bond_info['prot_chain_auth']
    prot_res_num = int(bond_info['prot_res_num_auth'])
    prot_atom_name = bond_info['prot_atom']
    prot_atom_alt_id = bond_info['prot_atom_alt_id']

    lig_atom = None
    prot_atom = None

    for model in structure:
        for chain in model:
            if chain.get_id() == lig_chain:
                for residue in chain:
                    if (residue.get_id()[1] == lig_res_num and 
                        residue.get_resname() == bond_info['lig_name']):
                        for atom in residue:
                            if (atom.get_id() == lig_atom_name and 
                                ((lig_atom_alt_id in (' ', '?')) or (atom.get_altloc() == lig_atom_alt_id))):
                                lig_atom = atom
            
            if chain.get_id() == prot_chain:
                for residue in chain:
                    if (residue.get_id()[1] == prot_res_num and 
                        residue.get_resname() == bond_info['prot_res']):
                        for atom in residue:
                            if (atom.get_id() == prot_atom_name and 
                                ((prot_atom_alt_id in (' ', '?')) or (atom.get_altloc() == prot_atom_alt_id))):
                                prot_atom = atom

    if lig_atom is None or prot_atom is None:
        print(f"Could not find ligand or protein atom: {bond_info}")
        return False

    distance = lig_atom - prot_atom
    if np.linalg.norm(distance) > BOND_LENGTH_CUTOFF:
        print(f"Covalent bond length {distance:.2f} exceeds cutoff {BOND_LENGTH_CUTOFF}")
        return False

    return True


def check_single_chain_contact(structure, info):
    ligand_atoms = []
    protein_atoms = []
    for atom in structure.get_atoms():
        resname = atom.get_parent().get_resname()
        auth_chain = atom.get_parent().get_parent().get_id()
        res_n = atom.get_parent().get_id()[1]
        if (resname == info['lig_name'] 
            and auth_chain == info['lig_chain_auth']
            and res_n == int(info['lig_num_auth'])):
            ligand_atoms.append(atom)
        elif atom.get_id() != 'H' and resname in STANDARD_AA:
            # only consider atoms of standard residues to check binding site. No "H"etero atoms.
            protein_atoms.append(atom)

    if not ligand_atoms:
        print("Ligand residues not found in structure")
        return False

    ns = NeighborSearch(protein_atoms)
    binding_site_chains = set()

    for lig_atom in ligand_atoms:
        neighbors = ns.search(lig_atom.get_coord(), BINDING_SITE_CUTOFF)
        for neighbor in neighbors:
            chain = neighbor.get_parent().get_parent().get_id()
            binding_site_chains.add(chain)
    
    if len(binding_site_chains) > 1:
        return False

    return True


def get_num_ligand_atoms(mmcif_dict, 
                         ligand_name,
                         ligand_chain,
                         ligand_num):
    atom_res_ids = get_field(mmcif_dict, '_atom_site.label_comp_id')
    atom_chain_ids = get_field(mmcif_dict, '_atom_site.label_asym_id')
    atom_res_nums = get_field(mmcif_dict, '_atom_site.label_seq_id')
    atom_types = get_field(mmcif_dict, '_atom_site.type_symbol')

    count = 0
    for i in range(len(atom_res_ids)):
        # skip Hydrogen atoms
        if atom_types[i] == 'H':
            continue

        check_res = atom_res_ids[i] == ligand_name
        check_chain = atom_chain_ids[i] == ligand_chain
        check_num = atom_res_nums[i] == str(ligand_num)
        if check_res and check_chain and check_num:
            count += 1
    return count


def filter_nonspecific_binding(info_list):
    # consider a ligand as a nonspecific binder 
    # if it appears in more than 1 residues in the same chain
    lig_chain_map = {}

    for info in info_list:
        key = (info['lig_name'], info['res_chain'])
        lig_chain_map.setdefault(key, [])
        lig_chain_map[key].append(info)

    filtered_list = []
    for infos in lig_chain_map.values():
        if len(infos) == 1:
            filtered_list.append(infos[0])
        else:
            print(f"Nonspecific binding detected for ligand {infos[0]['lig_name']} in chain {infos[0]['res_chain']}")

    return filtered_list


def process_mmcif_file(fn, 
                       check_entity_property=True,
                       check_date=True,
                       other_residue=False):
    print(f"Processing {fn}")
    
    mmcif_dict = MMCIF2Dict(fn)
    pdb_id = Path(fn).stem.lower()

    try:
        deposit_date = mmcif_dict[
            '_pdbx_database_status.recvd_initial_deposition_date'][0]
        deposit_date = datetime.datetime.strptime(deposit_date, "%Y-%m-%d").date()
        exp_type = mmcif_dict['_exptl.method'][0]
        
        if exp_type not in EXP_TYPES:
                return None
        
        if exp_type == 'X-RAY DIFFRACTION':
            resol = float(mmcif_dict['_refine.ls_d_res_high'][0])
        else:
            resol = float(mmcif_dict['_em_3d_reconstruction.resolution'][0])
        
        if check_entity_property:
            if deposit_date < DEPOSIT_DATE_CUTOFF and check_date:
                return None
            
            if resol > RES_CUTOFF:
                return None
    except (KeyError, IndexError, ValueError):
        print(f"experiment type or resolution error: {fn}")
        return None
    
    cov_bonds = []
    if '_struct_conn.id' in mmcif_dict:
        bond_ids = get_field(mmcif_dict, '_struct_conn.id')
        conn_types = get_field(mmcif_dict, '_struct_conn.conn_type_id')
        ptnr1_comp = get_field(mmcif_dict, '_struct_conn.ptnr1_label_comp_id')
        ptnr1_seq = get_field(mmcif_dict, '_struct_conn.ptnr1_label_seq_id')
        ptnr1_atom = get_field(mmcif_dict, '_struct_conn.ptnr1_label_atom_id')
        ptnr1_alt_id = get_field(mmcif_dict, '_struct_conn.pdbx_ptnr1_label_alt_id')
        ptnr1_asym = get_field(mmcif_dict, '_struct_conn.ptnr1_label_asym_id')
        ptnr1_chain_auth = get_field(mmcif_dict, '_struct_conn.ptnr1_auth_asym_id')
        ptnr1_seq_auth = get_field(mmcif_dict, '_struct_conn.ptnr1_auth_seq_id')
        ptnr2_comp = get_field(mmcif_dict, '_struct_conn.ptnr2_label_comp_id')
        ptnr2_seq = get_field(mmcif_dict, '_struct_conn.ptnr2_label_seq_id')
        ptnr2_atom = get_field(mmcif_dict, '_struct_conn.ptnr2_label_atom_id')
        ptnr2_alt_id = get_field(mmcif_dict, '_struct_conn.pdbx_ptnr2_label_alt_id')
        ptnr2_asym = get_field(mmcif_dict, '_struct_conn.ptnr2_label_asym_id')
        ptnr2_chain_auth = get_field(mmcif_dict, '_struct_conn.ptnr2_auth_asym_id')
        ptnr2_seq_auth = get_field(mmcif_dict, '_struct_conn.ptnr2_auth_seq_id')

        for i, conn_type in enumerate(conn_types):
            if conn_type == 'covale':
                bond = {
                    'id': bond_ids[i],
                    'ptnr1_comp': ptnr1_comp[i],
                    'ptnr1_seq': ptnr1_seq[i],
                    'ptnr1_atom': ptnr1_atom[i],
                    'ptnr1_alt_id': ptnr1_alt_id[i],
                    'ptnr1_asym': ptnr1_asym[i],
                    'ptnr2_comp': ptnr2_comp[i],
                    'ptnr2_seq': ptnr2_seq[i],
                    'ptnr2_atom': ptnr2_atom[i],
                    'ptnr2_alt_id': ptnr2_alt_id[i],
                    'ptnr2_asym': ptnr2_asym[i],
                    'ptnr1_chain_auth': ptnr1_chain_auth[i],
                    'ptnr1_seq_auth': ptnr1_seq_auth[i],
                    'ptnr2_chain_auth': ptnr2_chain_auth[i],
                    'ptnr2_seq_auth': ptnr2_seq_auth[i]
                }
                cov_bonds.append(bond)
    
    if not cov_bonds:
        print(f"No covalent bonds found: {fn}")
        return None

    linked_ligands = {}
    for bond in cov_bonds:
        is_p1_protein = bond['ptnr1_comp'] in STANDARD_AA
        is_p2_protein = bond['ptnr2_comp'] in STANDARD_AA

        if is_p1_protein and not is_p2_protein:
            lig_key = (bond['ptnr2_asym'], bond['ptnr2_comp'], bond['ptnr2_seq'])
            linked_ligands.setdefault(lig_key, [])
            bond_info = {
                'chain': bond['ptnr1_asym'],
                'res_name': bond['ptnr1_comp'],
                'res_num': bond['ptnr1_seq'],
                'res_atom': bond['ptnr1_atom'],
                'res_atom_alt_id': bond['ptnr1_alt_id'],
                'chain_auth': bond['ptnr1_chain_auth'],
                'res_num_auth': bond['ptnr1_seq_auth'],
                'lig_chain_auth': bond['ptnr2_chain_auth'],
                'lig_num_auth': bond['ptnr2_seq_auth'],
                'lig_atom': bond['ptnr2_atom'],
                'lig_atom_alt_id': bond['ptnr2_alt_id'],
            }
            linked_ligands[lig_key].append(bond_info)

        elif is_p2_protein and not is_p1_protein:
            lig_key = (bond['ptnr1_asym'], bond['ptnr1_comp'], bond['ptnr1_seq'])
            linked_ligands.setdefault(lig_key, [])
            bond_info = {
                'chain': bond['ptnr2_asym'],
                'res_name': bond['ptnr2_comp'],
                'res_num': bond['ptnr2_seq'],
                'res_atom': bond['ptnr2_atom'],
                'res_atom_alt_id': bond['ptnr2_alt_id'],
                'chain_auth': bond['ptnr2_chain_auth'],
                'res_num_auth': bond['ptnr2_seq_auth'],
                'lig_chain_auth': bond['ptnr1_chain_auth'],
                'lig_num_auth': bond['ptnr1_seq_auth'],
                'lig_atom': bond['ptnr1_atom'],
                'lig_atom_alt_id': bond['ptnr1_alt_id'],
            }
            linked_ligands[lig_key].append(bond_info)

    if not linked_ligands:
        print(f"No linked ligands found: {fn}")
        return None
    
    single_cov_ligands = []
    for (asym, comp, seq), bond_info in linked_ligands.items():
        if len(bond_info) == 1:
            bond_info = bond_info[0]
            single_cov_ligands.append({
                'lig_chain': asym,
                'lig_res_num': seq,
                'lig_name': comp,
                'lig_atom': bond_info['lig_atom'],
                'lig_atom_alt_id': bond_info['lig_atom_alt_id'],
                'lig_chain_auth': bond_info['lig_chain_auth'],
                'lig_num_auth': bond_info['lig_num_auth'],
                'prot_chain': bond_info['chain'],
                'prot_res': bond_info['res_name'],
                'prot_res_num': bond_info['res_num'],
                'prot_chain_auth': bond_info['chain_auth'],
                'prot_res_num_auth': bond_info['res_num_auth'],
                'prot_atom': bond_info['res_atom'],
                'prot_atom_alt_id': bond_info['res_atom_alt_id'],
            })

    if not single_cov_ligands:
        print(f"No single covalent ligands found: {fn}")
        return None

    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_id, fn)
    except Exception as e:
        print(f"Skipped {fn}: failed to parse structure - {e}")
        return None

    valid_ligands = []
    for info in single_cov_ligands:
        linked_res = info['prot_res']
        protein_chain = info['prot_chain']
        protein_seq, entity_id = get_seq_from_mmcif_dict(
            mmcif_dict, protein_chain) 
        
        if  not other_residue and (linked_res not in LINKED_RES):
            print(f"Linked residue: {linked_res} not in allowed list: {fn}")
            continue

        if check_entity_property:
            if check_covalent_bond_length(structure, info) is False:
                print(f"Covalent bond length exceeds cutoff: {fn}")
                continue
            
            if check_ligand_is_saccharide(mmcif_dict, info['lig_name']):
                print(f"Ligand {info['lig_name']} identified as saccharide: {fn}")
                continue

            if check_ligand_in_exclude_list(info['lig_name']):
                print(f"Ligand {info['lig_name']} is in exclude list: {fn}")
                continue

            if check_ligand_properties(mmcif_dict, info['lig_name']) is False:
                print(f"Ligand {info['lig_name']} failed property checks: {fn}")
                continue

            if not protein_seq or len(protein_seq) <= MIN_CHAIN_LENGTH:
                print(f"Linked protein chain {protein_chain} too short or not found: {fn}")
                continue

        try:
            uniprot_id = get_uniprot_id_from_mmcif_dict(mmcif_dict, entity_id)
        except Exception as e:
            uniprot_id = None
        lig_atom_num = get_num_ligand_atoms(mmcif_dict, 
                                            info['lig_name'],
                                            info['lig_chain'],
                                            info['lig_res_num'])
        
        if info['lig_name'] in RDKIT_ERROR_LIG_TORSION.keys():
            n_tor = RDKIT_ERROR_LIG_TORSION[info['lig_name']]
        else:
            n_tor = CalcNumRotatableBonds(
                _get_mol_from_mmcif_data(mmcif_dict, info['lig_name'])[0])

        result = {  # noqa: F811
            'pdb': pdb_id,
            'date': str(deposit_date),
            'exp_type': exp_type,
            'resolution': resol,
            'uniprot_id': uniprot_id,
            'res_name': linked_res,
            'res_atom': info['prot_atom'],
            'res_atom_alt_id': info['prot_atom_alt_id'],
            'res_num': info['prot_res_num'],
            'res_num_auth': info['prot_res_num_auth'],
            'res_chain': protein_chain,
            'res_chain_auth': info['prot_chain_auth'],
            'lig_name': info['lig_name'], # CCD code
            'lig_atom': info['lig_atom'],
            'lig_atom_alt_id': info['lig_atom_alt_id'],
            'lig_num': info['lig_res_num'],
            'lig_num_auth': info['lig_num_auth'],
            'lig_chain': info['lig_chain'],
            'lig_chain_auth': info['lig_chain_auth'],
            'lig_atom_num': lig_atom_num,
            'lig_n_tor': n_tor
        }
        valid_ligands.append(result)

    valid_ligands = filter_nonspecific_binding(valid_ligands)

    final_ligands = []
    for info in valid_ligands:
        if check_entity_property:
            if not check_single_chain_contact(structure, info):
                # multiple chain contact should be considered at last.
                # Otherwise, multiple contact non-specific ligands are removed, 
                # and a single remaining one passes the filter
                print(f"Ligand interacts with multiple chains: {fn}")
                continue
        final_ligands.append(info)

    if not final_ligands:
        return None

    return final_ligands


def write_json(pdb_s, out_fn, 
               check_entity_property=True,
               check_date=True,
               other_residue=False):

    cif_s = []
    for pdb in pdb_s:
        if os.path.exists(f"{CIF_DIR}/{pdb}.cif"):
            cif_s.append(f"{CIF_DIR}/{pdb}.cif")
            continue
        sp.run(['cif_get',pdb])
        sp.run(['mv',f"{pdb}.cif",str(CIF_DIR)])
        cif_s.append(f"{CIF_DIR}/{pdb}.cif")

    # #DEBUG
    # cif_s = [f"{CIF_DIR}/7duq.cif"]

    data = []
    for cif in cif_s:
        cov_ligs = process_mmcif_file(cif, 
                    check_entity_property=check_entity_property,
                    check_date=check_date,
                    other_residue=other_residue)
        if cov_ligs is None:
            continue
        if len(cov_ligs) > 1:
            redundant_ligs = {}
            for info in cov_ligs:
                # same protein and ligand, but different chain
                key = (info['lig_name'], info['res_num'], info['res_name'])
                if key not in redundant_ligs:
                    redundant_ligs[key] = [info]
                else:
                    redundant_ligs[key].append(info)
            
            cov_ligs = [v[0] for k, v in redundant_ligs.items()]

        data.extend(cov_ligs)

    df = pd.DataFrame(data)

    redundant_pdb_list = []
    for g, gg in df.groupby(['lig_name', 'res_name', 'res_num_auth', 'uniprot_id']):
        if len(gg) > 1:
            best_entry = gg.loc[gg['resolution'].idxmin()]
            redundant_entries = gg[gg.index != best_entry.name]
            redundant_pdb_list.extend(redundant_entries['pdb'].tolist())
    
    df = df[~df['pdb'].isin(redundant_pdb_list)].reset_index(drop=True)
    df['smiles'] = df['lig_name'].apply(get_smiles)

    # DEBUG
    df.to_json(out_fn, orient="records", indent=2)


def compute_tanimoto_similarity(smi1, smi2, radius=2, nBits=2048):
    try:
        mol1 = Chem.MolFromSmiles(smi1)
        mol2 = Chem.MolFromSmiles(smi2)
        if mol1 is None or mol2 is None:
            print(f"Failed to read SMILES: {smi1}, {smi2}")
            return 0.0
        
        generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nBits)
        
        fp1 = generator.GetFingerprint(mol1)
        fp2 = generator.GetFingerprint(mol2)
        
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
        return similarity

    except Exception as e:
        print(f"Error computing Tanimoto similarity: {e}")
        return 0.0


def check_file_exists(file_path):
    return (os.path.exists(file_path) and os.path.getsize(file_path) > 0)

def get_smiles(ligand_name):
    # some ligand CIF files failed to download from RCSB (eg. J5G,J5J)
    cif_path = f"{CIF_DIR}/{ligand_name}.cif"

    if check_file_exists(cif_path) is False:
        url = f"https://files.rcsb.org/ligands/download/{ligand_name}.cif"
        sp.run(["wget", "-O", cif_path, url])
        if not check_file_exists(cif_path):
            print(f"Failed to download CIF for ligand {ligand_name}")
            return get_smiles_from_rcsb(ligand_name)

    with open(cif_path) as f:
        for l in f:
            if "SMILES_CANONICAL" in l:
                smiles = l.split()[-1].strip("\"")
                return smiles
    return None


def get_smiles_from_rcsb(ligand_name):
    url = f"https://www.rcsb.org/ligand/{ligand_name}"
    r = requests.get(url)
    if r.status_code != 200:
        print(f"Failed to get SMILES for ligand {ligand_name} from RCSB")
        return None
    
    soup = BeautifulSoup(r.text, 'html.parser')
    can_smiles = []
    for tr in soup.find_all("tr"):
        tds = tr.find_all("td")
        if len(tds) >= 4 and "Canonical SMILES" in tr.text:
            smiles = tds[3].get_text(strip=True)
            can_smiles.append(smiles)
    
    return can_smiles[0] if can_smiles else None



def get_pdbs_for_uniprot(uniprot_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{uniprot_id}"
    r = requests.get(url)
    if r.status_code != 200:
        print(f"Failed to get PDBs for UniProt {uniprot_id}")
        return []
    data = r.json()
    if uniprot_id not in data:
        return []
    pdbs = list(data[uniprot_id]['PDB'].keys())
    return pdbs


MOL2_PATH = Path(__file__).parents[2]/"inputs/mol2"
def make_cross_set():
    # with open("cov_set.json") as f:
    #     covset_data = json.load(f)
    with open("all_cov_set.json") as f:
        covset_data = json.load(f)
    with open("scarpino_set.json") as f:
        scaprino_data = json.load(f)

    data = covset_data + scaprino_data

    df = pd.DataFrame(data)
    # TODO: fix missing uniprot cases.. (only-digit cases or pdb_id)
    df = df.dropna(subset=['uniprot_id'])
    df = df[df['uniprot_id']!= ""]
    invalid_uniprot = np.array([len(i) for i in df['uniprot_id'].values]) <= 4
    df = df[~invalid_uniprot].reset_index(drop=True)

    cross_data = []
    for uni, g in df.groupby('uniprot_id'):
        if len(g) == 1:
            # Need apo/noncov-ligand bound structures.
            print(f"Single entry for uniprot {uni}, {g['pdb'].values[0]}, skipped.")
            pdbs = get_pdbs_for_uniprot(uni)
            if len(pdbs) == 0:
                continue
            cross_data.append((g['pdb'].values[0], '', uni, ''))
            continue
        for pdb in g['pdb']:
            q_smiles = g[g['pdb'] == pdb]['smiles'].values[0]
            closest = None
            tanimoto = -1
            for _, row in g.iterrows():
                if row['pdb'] == pdb:
                    continue
                tani = compute_tanimoto_similarity(
                    q_smiles,
                    row['smiles'])
                if tani > tanimoto and tani < 1.0:
                    # TODO DEBUG: the same molecule of different stereoisomers are listed under different names within the same UniProt group. (duplicated in the cov_set)
                    tanimoto = tani
                    closest = row['pdb']
            if closest is None:
                continue
            cross_data.append((pdb, closest, uni, tanimoto))
    cross_df = pd.DataFrame(cross_data, columns=['lig_pdb', 'rec_pdb', 'uniprot_id', 'tanimoto'])
    
    covset_list = [d['pdb'] for d in covset_data]
    scarpino_list = [d['pdb'] for d in scaprino_data]

    cross_covset = cross_df[cross_df['lig_pdb'].isin(covset_list)]
    cross_scarpino = cross_df[cross_df['lig_pdb'].isin(scarpino_list)]

    cross_covset.to_json("cross_allcovset.json", orient="records", indent=2)
    # cross_scarpino.to_json("cross_scarpino.json", orient="records", indent=2)


def cov_set():
    with open("cov_rcsb_20180101_20250919.txt") as f:
        pdb_s = [line.strip() for line in f if line.strip()]

    write_json(pdb_s, "cov_set.json")


def scarpino_set():
    with open("scarpino.txt") as f:
        pdb_s = [line.strip() for line in f if line.strip()]

    write_json(pdb_s, "scarpino_set.json",
               check_entity_property=False,
               check_date=False)

def all_cov_set():
    with open("cov_rcsb_19950101_20250919.txt") as f:
        pdb_s = [line.strip() for line in f if line.strip()]

    write_json(pdb_s, "all_cov_set.json",
               check_date=False,
               other_residue=True)

if __name__ == "__main__":
    # cov_set()
    # scarpino_set()
    # all_cov_set()
    make_cross_set()
