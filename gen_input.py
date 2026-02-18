import json
import pandas as pd
import subprocess as sp
import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB import StructureBuilder, Atom, Residue, Chain
from pathlib import Path
import requests

BINDING_SITE_CUTOFF = 6.0
COFACTORS = [
    "FAD", "NAD", "NAP", "HEM", "FMN", "COA", "GSH", "THF",
    "SAM", "PLP", "ADP", "ATP", "GDP", "MG", "CA", "ZN", "FE",
    "CU", "MN", "NI", "CL", "K", "NA", "NAG", "FE2", "ANP", "SAH"
]
BIOIRR_RES = [
    "HOH", "EDO", "WAT", "CL", "SO4", "DMS",
    "PO4", "NO3", "SCN", "UNX", "UNL", "GOL",
    "ACT", "PEG", "PGE", "MLI", "CO", "IMD",
    "YT3", "TSL", "ACN", "MES", "MPD", "MYR",
    "NTK", "1PE", "BR", "FMT"," TLA", "P6G", "PG6",
    "PG4", "CIT"
] + [
    "MK7",  # free ligand form of cov_lig JJT (6qw9), overlaped in pdb
    "W6Z", # free ligand form of cov_lig W6Z (6zam), overlaped in pdb
]
CAN_RES = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL"
]
NONCAN_RES = [
    'CXM', 'CME', 'PTR', 'HYP','CSO', 'MSE', 'MLZ', 
    'TPO', 'MEN', '3CT', 'SNN', 'FTR','9IJ', 'SEP', 
    'FME', 'AYA', 'KCX', 'TYS', 'PCA', 'PHI'
]
RES_NEIGHBOR = {
    "CB": "CA", # PDB deposit bug: CB of SER is assigned to ligand (eg.3ncl)
    "OG": "CB",
    "SG": "CB",
    "ND1": "CE1",
    "NE2": "CE1",
    "CE": "CD", # PDB deposit bug: NZ is assigned to ligand (eg. 8h7f)
    "NZ": "CE",
    "OE1": "CD",
    "OE2": "CD",
    "OG1": "CB",
    "OH": "CE",
    "OD1":"CG",
    "OD2":"CG",
}

with open("cov_set.json", "r") as f:
    data = json.load(f)
df = pd.DataFrame(data)

INPUT_DIR = Path(__file__).parents[2]/"inputs"
CIF_DIR = INPUT_DIR/"cif"
MOL2_DIR = INPUT_DIR/"mol2"
PDB_DIR = INPUT_DIR/"pdb"
FASTA_DIR = INPUT_DIR/"fasta"


def get_cofactors(cif, info):
    # get residue names except canonical residue, water and the linked ligand
    parser = MMCIFParser()
    try:
        structure = parser.get_structure(info["pdb"], cif)
    except Exception as e:
        print(f"Error parsing {cif}: {e}")
        return []
    
    lig_atoms = []
    het_atoms = []
    model = structure[0]
    for atom in model.get_atoms():
        resname = atom.get_parent().get_resname()
        chain_auth = atom.get_parent().get_parent().get_id()
        res_n = atom.get_parent().get_id()[1]
        if (chain_auth == info['lig_chain_auth'] 
            and resname == info['lig_name']
            and res_n == int(info['lig_num_auth'])):
            lig_atoms.append(atom)
        elif resname in BIOIRR_RES + CAN_RES + [info['lig_name']]:
            continue
        else:
            het_atoms.append(atom)
    
    if het_atoms == []:
        return []
    ns = NeighborSearch(het_atoms)
    cof_info = set()
    overlap_ligs = []
    for lig_atom in lig_atoms:
        neighbors = ns.search(lig_atom.get_coord(), BINDING_SITE_CUTOFF)
        for neighbor in neighbors:
            res_id = neighbor.get_parent().get_full_id()
            if res_id in overlap_ligs:
                continue
            
            d = np.linalg.norm(lig_atom.get_coord() - neighbor.get_coord())
            if d < 1.0: # clash: overlaped free ligand
                overlap_lig_id = (neighbor.get_parent().get_full_id())
                overlap_ligs.append(overlap_lig_id)
                continue
            resname = neighbor.get_parent().get_resname()
            auth_res_id = neighbor.get_parent().get_id()[1]
            auth_chain_id = neighbor.get_parent().get_parent().get_id()
            cof_info.add((resname, auth_chain_id, auth_res_id))

    return list(cof_info)


def init_builder(name):
    builder = StructureBuilder.StructureBuilder()
    builder.init_structure(name)
    builder.init_model(0)
    return builder.get_structure()


def add_residue(structure,
                chain_id,
                resname,
                het_flag,
                seq_id,
                ins_code,
                ):
    
    model = structure[0]
    if chain_id not in [c.id for c in model]:
        chain = Chain.Chain(chain_id)
        model.add(chain)
    else:
        chain = model[chain_id]

    res_id = (het_flag, seq_id, ins_code)
    residue = Residue.Residue(res_id, resname, '')
    chain.add(residue)
    return residue


def copy_atom(atom):
    return Atom.Atom(
        atom.get_name(),
        atom.coord,
        atom.bfactor,
        atom.occupancy,
        atom.altloc,
        atom.fullname,
        atom.serial_number,
        element=atom.element,
    )


# DEBUG:TODO
def format_cif_atomline(data):
    # line = (
    #     f"{data[0]:<6} {data[1]:<5} {data[2]:<2} {data[3]:<4} "
    #     f"{data[4]} {data[5]:<3} {data[6]} {data[7]} "
    #     f"{data[8]:<4} {data[9]} {data[10]:<7} {data[11]:<7} {data[12]:<7} "
    #     f"{data[13]} {data[14]}  {data[15]} {data[16]:3} {data[17]:3} " 
    #     f"{data[18]} {data[19]:<4} {data[20]}\n"
    # )
    line = " ".join(data) + "\n"
    return line


def cif_label_seq_renumbering(ref_cif, cif, out_cif):
    mmcif_dict = MMCIF2Dict(ref_cif)
    
    auth_seq_ids = get_field(mmcif_dict, "_atom_site.auth_seq_id")
    label_seq_ids = get_field(mmcif_dict, "_atom_site.label_seq_id")
    auth_asym_ids = get_field(mmcif_dict, "_atom_site.auth_asym_id")
    label_asym_ids = get_field(mmcif_dict, "_atom_site.label_asym_id")

    auth_to_label = {}
    for (a_i, a_c, l_i, l_c) in zip(
            auth_seq_ids, auth_asym_ids, label_seq_ids, label_asym_ids):
        auth_to_label[(a_c, a_i)] = (l_c, l_i)

    wrt = []
    with open(cif, "r") as f:
        for l in f:
            if l.startswith("ATOM") or l.startswith("HETATM"):
                s = l.split()
                auth_seq_id = s[15]
                auth_asym_id = s[16]
                key = (auth_asym_id, auth_seq_id)
                label_asym_id, label_seq_id = auth_to_label[key]
                s[8] = label_seq_id
                s[6] = label_asym_id
                new_l = format_cif_atomline(s)
                wrt.append(new_l)
            else:
                wrt.append(l)

    with open(out_cif, "w") as f:
        f.writelines(wrt)

# lig mol2 and lig_ca mol2
def gen_cif(prefix, cif, info, cof_info):
    parser = MMCIFParser()
    try:
        structure = parser.get_structure(prefix, cif)
    except Exception as e:
        print(f"Error parsing {cif}: {e}")
        return

    lig_only = init_builder('lig_only')
    lig = init_builder('lig')
    lig_ca = init_builder('lig_ca')
    lig_ad4 = init_builder('lig_ad4')
    prot = init_builder('prot')
    prot_ala = init_builder('prot_ala')
    full = init_builder('full')

    res_chain_auth = info["res_chain_auth"]
    res_num_auth = int(info["res_num_auth"])
    res_atom = info["res_atom"]

    lig_chain_auth = info["lig_chain_auth"]
    lig_num_auth = int(info["lig_num_auth"])

    model = structure[0]
    for chain in model:
        for residue in chain:
            if residue.get_resname() in BIOIRR_RES:
                continue
            res_auth_id = residue.get_full_id()[3][1]
            res_auth_chain_id = residue.get_full_id()[2]
            res_name = residue.get_resname()

            is_lig = (res_auth_chain_id == lig_chain_auth 
                      and res_auth_id == lig_num_auth)
            is_res = (res_auth_chain_id == res_chain_auth 
                      and res_auth_id == res_num_auth)
            is_receptor = (res_auth_chain_id == res_chain_auth)
            is_receptor_canonical = (res_auth_chain_id == res_chain_auth
                        and res_name in CAN_RES)
            is_cofactor = (
                (res_name, res_auth_chain_id, res_auth_id) in cof_info
            )

            if is_lig:
                r_lig = add_residue(lig, chain.id, res_name, *residue.id)
                r_lig_ca = add_residue(lig_ca, chain.id, res_name, *residue.id)
                r_lig_ad4 = add_residue(lig_ad4, chain.id, res_name, *residue.id)
                r_lig_only = add_residue(lig_only, chain.id, res_name, *residue.id)
                r_full = add_residue(full, chain.id, res_name, *residue.id)

                linked_atom = residue[info["lig_atom"]]
                for atom in residue:
                    if atom.element == "H" and atom - linked_atom < 1.0:                        
                        # In case the H atoms is bonded to the linked atom,
                        # resulting in clash in docking.
                        continue
                    if atom.get_name() in ["CB", "CA", res_atom]:
                        atom.name = atom.get_name() + "9"  # to avoid name clash
                    r_lig.add(copy_atom(atom))
                    r_lig_ca.add(copy_atom(atom))
                    r_lig_ad4.add(copy_atom(atom))
                    r_full.add(copy_atom(atom))
                    r_lig_only.add(copy_atom(atom))
                continue

            elif is_res:
                r_lig = add_residue(lig, chain.id, res_name, *residue.id)
                r_lig_ca = add_residue(lig_ca, chain.id, res_name, *residue.id)
                r_lig_ad4 = add_residue(lig_ad4, chain.id, res_name, *residue.id)
                r_prot = add_residue(prot, chain.id, res_name, *residue.id)
                r_prot_ala = add_residue(prot_ala, chain.id, "ALA", *residue.id)
                r_full = add_residue(full, chain.id, res_name, *residue.id)

                linked_atom = residue[res_atom]
                for atom in residue:
                    if atom.element == "H" and atom - linked_atom < 1.25:
                        continue
                    r_lig_ad4.add(copy_atom(atom))
                    r_prot.add(copy_atom(atom))
                    r_full.add(copy_atom(atom))
                    if atom.get_name() in ["CA", "N", "C", "O", "CB"]:
                        r_prot_ala.add(copy_atom(atom))
                    if atom.get_name() not in ["N", "C", "O"]:
                        r_lig_ca.add(copy_atom(atom))
                    if atom.get_name() in [res_atom, RES_NEIGHBOR[res_atom]]:
                        r_lig.add(copy_atom(atom))
                continue

            if is_receptor:
                r_full = add_residue(full, chain.id, res_name, *residue.id)
                for atom in residue:
                    r_full.add(copy_atom(atom))
            
            if is_receptor_canonical:
                r_prot = add_residue(prot, chain.id, res_name, *residue.id)
                r_prot_ala = add_residue(prot_ala, chain.id, res_name, *residue.id)
                for atom in residue:
                    r_prot.add(copy_atom(atom))
                    r_prot_ala.add(copy_atom(atom))

            if is_cofactor:
                r_prot = add_residue(prot, chain.id, res_name, *residue.id)
                r_prot_ala = add_residue(prot_ala, chain.id, res_name, *residue.id)
                for atom in residue:
                    r_prot.add(copy_atom(atom))
                    r_prot_ala.add(copy_atom(atom))


    lig_tmp_cif = CIF_DIR/f"{prefix}_lig_tmp.cif"
    lig_ca_tmp_cif = CIF_DIR/f"{prefix}_lig_ca_tmp.cif"
    lig_ad4_tmp_cif = CIF_DIR/f"{prefix}_lig_ad4_tmp.cif"
    prot_tmp_cif = CIF_DIR/f"{prefix}_prot_tmp.cif"
    prot_ala_tmp_cif = CIF_DIR/f"{prefix}_prot_ala_tmp.cif"
    full_tmp_cif = CIF_DIR/f"{prefix}_comp_tmp.cif"
    lig_only_tmp_cif = CIF_DIR/f"{prefix}_lig_only_tmp.cif"

    lig_cif = CIF_DIR/f"{prefix}_lig.cif"
    lig_ca_cif = CIF_DIR/f"{prefix}_lig_ca.cif"
    lig_ad4_cif = CIF_DIR/f"{prefix}_lig_ad4.cif"
    prot_cif = CIF_DIR/f"{prefix}_prot.cif"
    prot_ala_cif = CIF_DIR/f"{prefix}_prot_ala.cif"
    full_cif = CIF_DIR/f"{prefix}_comp.cif"
    lig_only_cif = CIF_DIR/f"{prefix}_lig_only.cif"

    io = MMCIFIO()

    io.set_structure(lig)
    io.save(str(lig_tmp_cif))
    cif_label_seq_renumbering(cif, lig_tmp_cif, lig_cif)

    io.set_structure(lig_ca)
    io.save(str(lig_ca_tmp_cif))
    cif_label_seq_renumbering(cif, lig_ca_tmp_cif, lig_ca_cif)

    io.set_structure(lig_ad4)
    io.save(str(lig_ad4_tmp_cif))
    cif_label_seq_renumbering(cif, lig_ad4_tmp_cif, lig_ad4_cif)

    io.set_structure(lig_only)
    io.save(str(lig_only_tmp_cif))
    cif_label_seq_renumbering(cif, lig_only_tmp_cif, lig_only_cif)

    io.set_structure(prot)
    io.save(str(prot_tmp_cif))
    cif_label_seq_renumbering(cif, prot_tmp_cif, prot_cif)

    io.set_structure(prot_ala)
    io.save(str(prot_ala_tmp_cif))
    cif_label_seq_renumbering(cif, prot_ala_tmp_cif, prot_ala_cif)

    io.set_structure(full)
    io.save(str(full_tmp_cif))
    cif_label_seq_renumbering(cif, full_tmp_cif, full_cif)

    cov_bond_info = f"""
# 
_struct_conn.id                            covale1 
_struct_conn.conn_type_id                  covale 
_struct_conn.pdbx_leaving_atom_flag        none 
_struct_conn.pdbx_PDB_id                   ? 
_struct_conn.ptnr1_label_asym_id           {info['res_chain']} 
_struct_conn.ptnr1_label_comp_id           {info['res_name']} 
_struct_conn.ptnr1_label_seq_id            {info['res_num']} 
_struct_conn.ptnr1_label_atom_id           {info['res_atom']} 
_struct_conn.pdbx_ptnr1_label_alt_id       ? 
_struct_conn.pdbx_ptnr1_PDB_ins_code       ? 
_struct_conn.pdbx_ptnr1_standard_comp_id   ? 
_struct_conn.ptnr1_symmetry                1_555 
_struct_conn.ptnr2_label_asym_id           {info['lig_chain']} 
_struct_conn.ptnr2_label_comp_id           {info['lig_name']} 
_struct_conn.ptnr2_label_seq_id            . 
_struct_conn.ptnr2_label_atom_id           {info['lig_atom']} 
_struct_conn.pdbx_ptnr2_label_alt_id       ? 
_struct_conn.pdbx_ptnr2_PDB_ins_code       ? 
_struct_conn.ptnr1_auth_asym_id            {info['res_chain_auth']} 
_struct_conn.ptnr1_auth_comp_id            {info['res_name']} 
_struct_conn.ptnr1_auth_seq_id             {info['res_num_auth']} 
_struct_conn.ptnr2_auth_asym_id            {info['lig_chain_auth']} 
_struct_conn.ptnr2_auth_comp_id            {info['lig_name']} 
_struct_conn.ptnr2_auth_seq_id             {info['lig_num_auth']} 
_struct_conn.ptnr2_symmetry                1_555 
_struct_conn.pdbx_ptnr3_label_atom_id      ? 
_struct_conn.pdbx_ptnr3_label_seq_id       ? 
_struct_conn.pdbx_ptnr3_label_comp_id      ? 
_struct_conn.pdbx_ptnr3_label_asym_id      ? 
_struct_conn.pdbx_ptnr3_label_alt_id       ? 
_struct_conn.pdbx_ptnr3_PDB_ins_code       ? 
_struct_conn.details                       ? 
_struct_conn.pdbx_dist_value               ? 
_struct_conn.pdbx_value_order              ? 
_struct_conn.pdbx_role                     ? 
# 
_struct_conn_type.id          covale 
_struct_conn_type.criteria    ? 
_struct_conn_type.reference   ? 
# 
"""
    
    with open(lig_cif, "a") as f:
        f.write(cov_bond_info)
    with open(lig_ca_cif, "a") as f:
        f.write(cov_bond_info)
    with open(lig_ad4_cif, "a") as f:
        f.write(cov_bond_info)
    with open(full_cif, "a") as f:
        f.write(cov_bond_info)

    chem_comp_info = get_chem_comp_info(cif, info['lig_name'])

    if chem_comp_info:
        with open(lig_cif, "a") as f:
            f.write(chem_comp_info)
        with open(lig_ca_cif, "a") as f:
            f.write(chem_comp_info)
        with open(lig_ad4_cif, "a") as f:
            f.write(chem_comp_info)
        with open(full_cif, "a") as f:
            f.write(chem_comp_info)
        with open(lig_only_cif, "a") as f:
            f.write(chem_comp_info)

    return lig_cif, lig_ca_cif, lig_ad4_cif, prot_cif, prot_ala_cif, full_cif, lig_only_cif


def remove_digit(s):
    return ''.join([i for i in s if not i.isdigit()])


def check_pdb_atom_type(pdb_path):
    # check if atom type is * in pdb file
    # some Cl, Br atoms are marked as *
    wrt = []
    with open(pdb_path, "r") as f:
        for l in f:
            if l.startswith("HETATM"):
                if l[76:78].strip() == "*":
                    atom_type = remove_digit(l[12:16].strip())
                    if len(atom_type) > 2:
                        # in the case of alphabet naming like CLA, CLD, etc.
                        atom_type = atom_type[:2]
                    new_l = l[:76] + f"{atom_type:>2}" + "\n"
                    wrt.append(new_l)
                else:
                    wrt.append(l)
            else:
                wrt.append(l)
    
    with open(pdb_path, "w") as f:
        f.writelines(wrt)


def cif_to_mol2(cif, mol2):
    tmp_pdb = mol2.with_suffix(".pdb")
    cmd = f"obabel {cif} -o pdb -O {tmp_pdb}"
    try:
        sp.run(cmd, shell=True, check=True, stdout=sp.PIPE, stderr=sp.PIPE)
        # to fix obabel cif-to-pdb conversion bug
        wrt = []
        with open(tmp_pdb) as f:
            for l in f:
                if l.startswith("ATOM"):
                    new_l = "HETATM" + l[6:]
                    wrt.append(new_l)
                else:
                    wrt.append(l)
        with open(tmp_pdb, "w") as f:
            f.writelines(wrt)
        check_pdb_atom_type(tmp_pdb)
        cmd = f"obabel {tmp_pdb} -o mol2 -O {mol2}"
        sp.run(cmd, shell=True, check=True, stdout=sp.PIPE, stderr=sp.PIPE)
        return True
    except sp.CalledProcessError as e:
        print(f"Error converting {cif} to {mol2}: {e.stderr.decode().strip()}")
        return False


def cif_to_pdb(cif, pdb):
    cmd = f"obabel {cif} -o pdb -O {pdb}"
    try:
        sp.run(cmd, shell=True, check=True, stdout=sp.PIPE, stderr=sp.PIPE)
    except sp.CalledProcessError as e:
        print(f"Error converting {cif} to {pdb}: {e.stderr.decode().strip()}")
        return False
    
    # correct bugs in PDB generated by obabel
    cif_chain_s = []
    with open(cif) as f:
        for l in f:
            if l.startswith("ATOM") or l.startswith("HETATM"):
                label_chain = l.split()[6]
                if label_chain not in cif_chain_s:
                    cif_chain_s.append(label_chain)

    wrt = []
    chain_occur_order = []
    with open(pdb) as f:
        for l in f:
            if l.startswith("ATOM") or l.startswith("HETATM"):
                chain_id = l[21]
                if chain_id in chain_occur_order:
                    o = chain_occur_order.index(chain_id)
                if chain_id not in chain_occur_order:
                    chain_occur_order.append(chain_id)
                    o = len(chain_id) - 1

                label_chain_id = cif_chain_s[o]
                new_l = l[:21] + label_chain_id + l[22:]
                
                if l.startswith("ATOM") and l[17:20] not in CAN_RES:
                    new_l = "HETATM" + new_l[6:]
                    wrt.append(new_l)
                else:
                    wrt.append(new_l)
    
    with open(pdb,"w") as f:
        f.writelines(wrt)

    check_pdb_atom_type(pdb)
    
    return True


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
    if '_entity_poly.pdbx_strand_id' not in mmcif_dict:
        return None, None

    entity_ids = get_field(mmcif_dict, "_entity_poly.entity_id")
    strand_ids = get_field(mmcif_dict, "_entity_poly.pdbx_strand_id")
    seq_s = get_field(mmcif_dict, "_entity_poly.pdbx_seq_one_letter_code_can")
                      
    for i, strands in enumerate(strand_ids):
        if chain_id in strands.split(','):
            return seq_s[i], int(entity_ids[i])
        
    return None, None


def zero_crd(cif, out_cif):
    # read cif and put all coordinates to zero
    parser = MMCIFParser()
    try:
        structure = parser.get_structure("temp", cif)
    except Exception as e:
        print(f"Error parsing {cif}: {e}")
        return False
    model = structure[0]
    for chain in model:
        for residue in chain:
            for atom in residue:
                atom.set_coord((0.0, 0.0, 0.0))
    
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(out_cif)
    cif_label_seq_renumbering(cif, out_cif, out_cif)
    
    return True


def prep_ad4_lig(pdb, res, chain, res_num):
    # (1) all start with "HETATM"
    # (2) all have same res, chain, and res_num
    # (3) two-charactered label_asym_id (eg. DA in 7m8r) exception
    wrt = []
    with open(pdb) as f:
        for l in f:
            if l.startswith("ATOM") or l.startswith("HETATM"):
                new_l = "HETATM" + l[6:]
                new_l = new_l[:17] + res + new_l[20:]
                chain_id = remove_digit(new_l[21:23].strip())
                if len(chain_id) == 1:
                    new_l = new_l[:21] + chain + new_l[22:]
                elif len(chain_id) == 2:
                    new_l = new_l[:21] + chain + new_l[23:]
                new_l = new_l[:22] + f"{str(res_num):>4}" + new_l[26:]
                wrt.append(new_l)
            else:
                wrt.append(l)
    
    with open(pdb, "w") as f:
        f.writelines(wrt)


def generate_inputs(cif, info):
    prefix = f"{info['pdb']}_{info['lig_name']}_{info['lig_chain']}"
    cof_info = get_cofactors(cif, info)

    cif_s = gen_cif(prefix, cif, info, cof_info)
    
    lig_cif = cif_s[0]
    lig_ca_cif = cif_s[1]
    lig_ad4_cif = cif_s[2]
    prot_cif = cif_s[3]
    prot_ala_cif = cif_s[4]
    full_cif = cif_s[5]
    lig_only_cif = cif_s[6]

    if not lig_cif:
        return
    
    lig_only_mol2 = MOL2_DIR/f"{prefix}_lig_only.mol2"
    cif_to_mol2(lig_only_cif, lig_only_mol2)
    
    lig_mol2 = MOL2_DIR/f"{prefix}_lig.mol2"
    cif_to_mol2(lig_cif, lig_mol2)

    lig_ca_mol2 = MOL2_DIR/f"{prefix}_lig_ca.mol2"
    cif_to_mol2(lig_ca_cif, lig_ca_mol2)

    lig_ad4_pdb = PDB_DIR/f"{prefix}_ad4lig.pdb"
    cif_to_pdb(lig_ad4_cif,lig_ad4_pdb)
    prep_ad4_lig(lig_ad4_pdb, 
                info['res_name'],
                info['res_chain'],
                info['res_num'])

    prot_pdb = PDB_DIR/f"{prefix}_prot.pdb"
    cif_to_pdb(prot_cif, prot_pdb)

    prot_ala_pdb = PDB_DIR/f"{prefix}_prot_ala.pdb"
    cif_to_pdb(prot_ala_cif, prot_ala_pdb)

    full_zero_cif = str(CIF_DIR/f"{prefix}_comp_zero.cif") # rf3 input
    zero_crd(full_cif, full_zero_cif)

    fasta_fn = FASTA_DIR/f"{prefix}.fasta"
    mmcif_dict = MMCIF2Dict(cif)
    prot_seq, entity_id = get_seq_from_mmcif_dict(mmcif_dict, info['res_chain_auth'])
    if prot_seq:
        with open(fasta_fn, "w") as f:
            f.write(f">{prefix}|{info['res_chain_auth']}|{entity_id}\n")
            f.write(prot_seq + "\n")



def read_covset():
    with open("cov_set.json") as f:
        data = json.load(f)
    with open("scarpino_set.json") as f:
        data += json.load(f)
    return data


def get_id(info):
    return f"{info['pdb']}_{info['lig_name']}_{info['lig_chain']}"

def check_cofactors():
    data = read_covset()

    cof_s = []
    for info in data:
        pdb_code = info["pdb"]
        cif = f"{CIF_DIR}/{pdb_code}.cif"
        cof_info = get_cofactors(cif, info)
        cof = list(set(c[0] for c in cof_info))
        id_ = get_id(info)
        cof_s.append(
            {"id": id_, 
             "cofactors": cof
            }
        )
    exit()     
    with open("cofactors.json", "w") as f:
        json.dump(cof_s, f, indent=2)

    cof_set = set()
    for i in cof_s:
        cof_set.update(i["cofactors"])

    print(f"Unique cofactors: {cof_set}")


def generate_cofactors_mol2_from_ccd():
    # ERROR: obabel fails to convert CCD-style cif into mol2
    with open("cofactors.json") as f:
        data = json.load(f)

    cof_set = set()
    for info in data:
        cof_set.update(info["cofactors"])
    
    for cof in cof_set:
        url = f"https://files.rcsb.org/ligands/download/{cof}.cif"
        response = requests.get(url)
        if response.status_code != 200:
            print(f"Failed to download {cof} cif from RCSB.")
            continue
        cif_path = MOL2_DIR/f"{cof}.cif"
        with open(cif_path, "wb") as f:
            f.write(response.content)
        mol2_path = MOL2_DIR/f"{cof}.mol2"
        sp.run(f"obabel -icif {cif_path} --gen3d -omol2 -O {mol2_path}", shell=True)
        sp.run(f"process_ligand.py {mol2_path}", shell=True)
        print(cif_path, mol2_path)
        

def get_chem_comp_info(cif, lig_name):
    chem_comp_info = []
    mmcif_dict = MMCIF2Dict(cif)

    if '_chem_comp.id' in mmcif_dict:
        comp_ids = get_field(mmcif_dict, '_chem_comp.id')
        comp_types = get_field(mmcif_dict, '_chem_comp.type')
        nstd_flag = get_field(mmcif_dict, '_chem_comp.mon_nstd_flag')
        names = get_field(mmcif_dict, '_chem_comp.name')
        synonyms = get_field(mmcif_dict, '_chem_comp.pdbs_synonyms')
        formulae = get_field(mmcif_dict, '_chem_comp.formula')
        formula_weight = get_field(mmcif_dict, '_chem_comp.formula_weight')

        chem_comp_info.append(
            "#\n"
            "loop_\n"\
            "_chem_comp.id\n"\
            "_chem_comp.type\n"\
            "_chem_comp.mon_nstd_flag\n"\
            "_chem_comp.name\n"\
            "_chem_comp.pdbx_synonyms\n"\
            "_chem_comp.formula\n"
            "_chem_comp.formula_weight\n"
        )
        for i in range(len(comp_ids)):
            cid = comp_ids[i]
            if cid != lig_name:
                continue
            name = names[i] if names else "?"
            ctype = comp_types[i] if comp_types else "?"
            nstd = nstd_flag[i] if nstd_flag else "?"
            syn = synonyms[i] if synonyms else "?"
            formula = formulae[i] if formulae else "?"
            weight = formula_weight[i] if formula_weight else "?"
            line = f"{cid:<3} {ctype} {nstd} \"{name}\" {syn} \"{formula}\" {weight}\n"
            chem_comp_info.append(line)
        chem_comp_info.append("#\n")
    
    if '_chem_comp_atom.comp_id' in mmcif_dict:
        comp_ids = get_field(mmcif_dict, '_chem_comp_atom.comp_id')
        atom_ids = get_field(mmcif_dict, '_chem_comp_atom.atom_id')
        element_s = get_field(mmcif_dict, '_chem_comp_atom.type_symbol')
        ar_s = get_field(mmcif_dict, '_chem_comp_atom.pdbx_aromatic_flag')
        stereo_s = get_field(mmcif_dict, '_chem_comp_atom.pdbx_stereo_config')
        ordinal_s = get_field(mmcif_dict, '_chem_comp_atom.pdbx_ordinal')

        chem_comp_info.append(
            "#\n"
            "loop_\n"\
            "_chem_comp_atom.comp_id\n"\
            "_chem_comp_atom.atom_id\n"\
            "_chem_comp_atom.type_symbol\n"\
            "_chem_comp_atom.pdbx_aromatic_flag\n"\
            "_chem_comp_atom.pdbx_stereo_config\n"\
            "_chem_comp_atom.pdbx_ordinal\n"
        )

        for i in range(len(comp_ids)):
            cid = comp_ids[i]
            if cid != lig_name:
                continue
            atom_id = atom_ids[i]
            element = element_s[i] if element_s else "?"
            ar = ar_s[i] if ar_s else "?"
            st = stereo_s[i] if stereo_s else "?"
            ordinal = ordinal_s[i] if ordinal_s else "?"
            line = f"{cid:<3} {atom_id:<4}   {element}  {ar} {st} {ordinal}\n"
            chem_comp_info.append(line)
        chem_comp_info.append("#\n")

    if '_chem_comp_bond.comp_id' in mmcif_dict:
        comp_ids = get_field(mmcif_dict, '_chem_comp_bond.comp_id')
        atom_id_1s = get_field(mmcif_dict, '_chem_comp_bond.atom_id_1')
        atom_id_2s = get_field(mmcif_dict, '_chem_comp_bond.atom_id_2')
        ord_s = get_field(mmcif_dict, '_chem_comp_bond.value_order')
        arom_s = get_field(mmcif_dict, '_chem_comp_bond.pdbx_aromatic_flag')
        stero_s = get_field(mmcif_dict, '_chem_comp_bond.pdbx_stereo_config')

        chem_comp_info.append(
            "#\n"
            "loop_\n"\
            "_chem_comp_bond.comp_id\n"\
            "_chem_comp_bond.atom_id_1\n"\
            "_chem_comp_bond.atom_id_2\n"\
            "_chem_comp_bond.value_order\n"\
            "_chem_comp_bond.pdbx_aromatic_flag\n"\
            "_chem_comp_bond.pdbx_stereo_config\n"\
            "_chem_comp_bond.pdbx_ordinal\n"
        )
        n = 0
        for i in range(len(comp_ids)):
            cid = comp_ids[i]
            if cid != lig_name:
                continue
            atom1 = atom_id_1s[i]
            atom2 = atom_id_2s[i]
            o = ord_s[i]
            ar = arom_s[i]
            st = stero_s[i]
            line = f"{cid} {atom1:<4s}{atom2:<4s} {o} {ar} {st} {n+1}\n"
            n += 1
            chem_comp_info.append(line)
        chem_comp_info.append("#\n")

    return "".join(chem_comp_info)


def extract_lig_from_cif(cif, lig_name, lig_cif):
    parser = MMCIFParser()
    try:
        structure = parser.get_structure("temp", cif)
    except Exception as e:
        print(f"Error parsing {cif}: {e}")
        return None
    
    model = structure[0]
    lig_structure = init_builder('lig')
    chain_res = None
    for chain in model:
        for residue in chain:
            if residue.get_resname() == lig_name:
                if chain_res is None:
                    chain_res = (chain.id, residue.id)
                if (chain.id, residue.id) != chain_res:
                    continue
                r_lig = add_residue(lig_structure, chain.id, residue.get_resname(), *residue.id)
                for atom in residue:
                    r_lig.add(copy_atom(atom))
                break

    io = MMCIFIO()
    io.set_structure(lig_structure)
    io.save(str(lig_cif))

    chem_comp_info = get_chem_comp_info(cif, lig_name)
    if chem_comp_info:
        with open(lig_cif, "a") as f:
            f.write(chem_comp_info)


def check_empty_file(fn):
    with open(fn) as f:
        wrt = f.read()
    if wrt == "":
        return False
    return True


def generate_cofactors_mol2():
    with open("cofactors.json") as f:
        data = json.load(f)

    cof_set = {}
    for d in data:
        id = d["id"]
        cof = d["cofactors"]
        for c in cof:
            if c not in cof_set:
                cof_set[c] = []
            cof_set[c].append(id.split('_')[0])
    
    for cof, pdbs in cof_set.items():
        pdb = pdbs[0]
        cif = CIF_DIR/f"{pdb}.cif"
        cof_cif = CIF_DIR/f"{cof}.cif"
        extract_lig_from_cif(cif, cof, cof_cif)
        
        mol2 = MOL2_DIR/f"{cof}.mol2"
        sp.run(f"obabel {cof_cif} -o mol2 -O {mol2}", shell=True)
        sp.run(f"process_ligand.py {mol2}", shell=True)
        
        if not check_empty_file(mol2):
            # obabel fails to convert single heavy atom cofactors
            # like FE2, NI, H2S, MG, NA
            print(f"{mol2} mol2 generation failed. {cof_cif}")


def generate_single_atom_cofactors():
    cof_s = {
        "FE2":"FE", 
        "NI":"NI", 
        "H2S":"S", 
        "MG":"MG", 
        "NA":"NA",
    }

    for cof in cof_s:
        pdb = MOL2_DIR/f"{cof}.pdb"
        symbol = cof_s[cof]
        with open(pdb, "w") as f:
            line = f"HETATM    1  {symbol:<3} {cof:<3} A   0     108.666   2.339  49.920  1.00  0.00          {symbol:>2}"
            f.write(line + "\n")
        
        mol2 = MOL2_DIR/f"{cof}.mol2"
        sp.run(f"obabel {pdb} -o mol2 -O {mol2}", shell=True)
        sp.run(f"process_ligand.py {mol2}", shell=True)

        if not check_empty_file(mol2):
            print(f"{mol2} mol2 generation failed.")


def main():
    with open("cov_set_nlnr.json") as f:
        data = json.load(f)

    # DEBUG
    # data = [d for d in data if d["pdb"] == "1ppv"] # 7duq: cof, Cl in lig

    generated = []
    for info in data:
        pdb_code = info["pdb"]
        res_type = info["res_name"]
        if res_type not in ["LYS"]: #DEBUG
            continue
        cif = f"{CIF_DIR}/{pdb_code}.cif"
        if not Path(cif).exists():
            print(f"CIF file {cif} does not exist. Skipping.")
            continue
        
        print(f"Processing {pdb_code}...")
        
        # if check_file_generated(info):
        #     print(f"{pdb_code} inputs already generated, skipping.")
        #     generated.append(get_id(info))
        #     continue
        
        generate_inputs(cif, info)

        prefix = get_id(info)

        lig_mol2 = MOL2_DIR/f"{prefix}_lig.mol2"
        lig_ca_mol2 = MOL2_DIR/f"{prefix}_lig_ca.mol2"
        sp.run(["process_ligand.py", str(lig_mol2)])
        sp.run(["process_ligand.py", str(lig_ca_mol2)])

        lig_init = lig_mol2.with_name(f"{prefix}_lig_init.mol2")
        lig_ca_init = lig_ca_mol2.with_name(f"{prefix}_lig_ca_init.mol2")
        if not lig_init.exists() or not lig_ca_init.exists():
            print(f"{pdb_code} ligand processing failed, skipping. {lig_init}, {lig_ca_init}")
            continue
        print(f"Finished processing {prefix}.")
        generated.append(prefix)

    with open("generated_inputs.txt", "w") as f:
        for g in generated:
            f.write(g + "\n")


def check_file_generated(d):
    prefix = f"{d['pdb']}_{d['lig_name']}_{d['lig_chain']}"
    lig_mol2 = MOL2_DIR/f"{prefix}_lig.mol2"
    lig_ca_mol2 = MOL2_DIR/f"{prefix}_lig_ca.mol2"
    prot_pdb = PDB_DIR/f"{prefix}_prot.pdb"
    prot_ala_pdb = PDB_DIR/f"{prefix}_prot_ala.pdb"
    full_zero_cif = CIF_DIR/f"{prefix}_comp_zero.cif"

    files = [lig_mol2, lig_ca_mol2, prot_pdb, prot_ala_pdb, full_zero_cif]
    for fn in files:
        if not fn.exists():
            return False
        else:
            if not check_empty_file(fn):
                return False
    return True


def filter_generated():
    with open("cov_set_nlnr.json") as f:
        data = json.load(f)
    
    with open("scarpino_set.json") as f:
        sca = json.load(f)

    sca_ids = [f"{d['pdb']}_{d['lig_name']}" for d in sca]

    if Path("generated_inputs.txt").exists():
        with open("generated_inputs.txt") as f:
            generated = [l.strip() for l in f]
        
        filtered = [d for d in data 
                    if (f"{d['pdb']}_{d['lig_name']}" in generated) and 
                    (f"{d['pdb']}_{d['lig_name']}" not in sca_ids)]
        
    else:
        filtered = [d for d in data 
                    if check_file_generated(d)]

    with open("cov_set_nlnr.json", "w") as f:
        json.dump(filtered, f, indent=2)


def merge_fasta_help(data, out_fn):    
    ids = [f"{d['pdb']}_{d['lig_name']}_{d['lig_chain']}" for d in data]
    FASTA_DIR = Path(__file__).parents[2]/"inputs"/"fasta"

    wrt = []
    for id_ in ids:
        fn = FASTA_DIR/f"{id_}.fasta"
        if not fn.exists():
            print(f"Missing fasta: {fn}")
            continue
        with open(fn) as f:
            for l in f:
                if l.startswith('>'):
                    l = f">{id_}" + "\n"
                wrt.append(l)

    with open(out_fn, "w") as f:
        f.writelines(wrt)


def merge_fasta_for_clustering():
    with open("cov_set.json") as f:
        covset = json.load(f)
    with open("scarpino_set.json") as f:
        sca = json.load(f)

    merge_fasta_help(sca, "scarpino_seq.fasta")
    merge_fasta_help(covset, "covset_seq.fasta")


def merge_fasta_for_msa():
    with open("cov_set_nlnr.json") as f:
        covset = json.load(f)
    with open("scarpino_set.json") as f:
        sca = json.load(f)
    all_data = covset + sca

    merge_fasta_help(all_data, "msa_gen.fasta")


def get_entity_id(chain_id, strand_ids, entity_ids):
    for i, strands in enumerate(strand_ids):
        if chain_id in strands.split(','):
            entity_id = int(entity_ids[i])
            return entity_id
    return None


def get_uniprot(entity_id, ref_entity_ids, db_codes):
    for i, eid in enumerate(ref_entity_ids):
        if int(eid) == entity_id:
            uniprot = db_codes[i]
            return uniprot
    return None


def gen_receptor_cif(prefix, cif, resn, target_uniprot_id):
    parser = MMCIFParser()
    mmcif_dict = MMCIF2Dict(cif)
    try:
        structure = parser.get_structure("temp", cif)
        model = structure[0]
    except Exception as e:
        print(f"Error parsing {cif}: {e}")
        return None
    
    entity_ids = mmcif_dict.get('_entity_poly.entity_id', [])
    strand_ids = get_field(mmcif_dict, "_entity_poly.pdbx_strand_id")

    ref_entity_ids = mmcif_dict.get('_struct_ref.entity_id', [])
    db_codes = mmcif_dict.get('_struct_ref.pdbx_db_accession', []) # This is usually the UniProt ID
    
    tmp_prot_ala_cif = CIF_DIR/f"{prefix}_prot_ala_tmp.cif"
    tmp_prot_cif = CIF_DIR/f"{prefix}_prot_tmp.cif"

    prot_ala_cif = CIF_DIR/f"{prefix}_prot_ala.cif"
    prot_cif = CIF_DIR/f"{prefix}_prot.cif"

    prot = init_builder('prot')
    prot_ala = init_builder('prot_ala')

    print("number of chains:", len(model))
    for chain in model:
        chain_id = chain.id

        entity_id = get_entity_id(chain_id, strand_ids, entity_ids)
        if entity_id is None:
            continue
        
        uniprot = get_uniprot(entity_id, ref_entity_ids, db_codes)
        
        if (uniprot is None) or (uniprot != target_uniprot_id):
            continue
        
        label_seq_ids = get_field(mmcif_dict, '_atom_site.label_seq_id')
        label_asym_ids = get_field(mmcif_dict, '_atom_site.label_asym_id')
        auth_asym_ids = get_field(mmcif_dict, '_atom_site.auth_asym_id')
        atom_serial_numbers = get_field(mmcif_dict, '_atom_site.id')
        
        # Create mapping from serial number to label_seq_id
        serial_to_label_seq = {}
        if atom_serial_numbers and label_seq_ids:
            for serial, label_seq, label_asym, auth_asym in zip(atom_serial_numbers, label_seq_ids, label_asym_ids, auth_asym_ids):
                serial_to_label_seq[int(serial)] = (label_seq, label_asym, auth_asym)
        
        for residue in chain:
            if residue.get_resname() not in CAN_RES:
                continue
            # Add residue once before processing atoms
            r_prot = add_residue(prot, chain.id, residue.get_resname(), *residue.id)
            r_prot_ala = add_residue(prot_ala, chain.id, "ALA", *residue.id)
            
            for atom in residue:
                serial_num = atom.get_serial_number()
                res_label_n = serial_to_label_seq.get(serial_num, (None, None, None))[0]
                
                atom_name = atom.get_name()
                if res_label_n == str(resn) and atom_name not in ['CA', 'C', 'N', 'O', 'CB']:
                    r_prot.add(copy_atom(atom))
                else:
                    r_prot_ala.add(copy_atom(atom))
                    r_prot.add(copy_atom(atom))
            
        break  # Assuming only one chain matches

    if len(list(prot.get_atoms())) == 0:
        return False
    
    io = MMCIFIO()
    io.set_structure(prot)
    io.save(str(tmp_prot_cif))
    cif_label_seq_renumbering(cif, tmp_prot_cif, prot_cif)

    io.set_structure(prot_ala)
    io.save(str(tmp_prot_ala_cif))
    cif_label_seq_renumbering(cif, tmp_prot_ala_cif, prot_ala_cif)

    return True


def receptor_prep_for_cross():
    with open("cross_covset.json") as f:
        cross_covset = json.load(f)
    with open("cross_scarpino.json") as f:
        cross_sca = json.load(f)

    with open("cov_set_nlnr.json") as f:
        cov = json.load(f)
    with open("scarpino_set.json") as f:
        sca = json.load(f)

    covset = cov + sca
    covset = {get_id(d): d for d in covset}

    cross_set = cross_covset + cross_sca
    data = [
            {
             'pdb': d['rec_id'].split('_')[0],
             'resn':covset[d['lig_id']]['res_num'],
             'uniprot_id': covset[d['lig_id']]['uniprot_id'],
             'rec_id': d['rec_id'],
             }
             for d in cross_set if d['source']=='additional'
    ]
    
    failed = [] # due to not matching uniprot id...
    for info in data:
        pdb = info["pdb"]
        resn = info["resn"]
        uniprot_id = info["uniprot_id"]

        cif = CIF_DIR/f"{pdb}.cif"
        if not cif.exists():
            sp.run(['cif_get',pdb])
            sp.run(['mv', f"{pdb}.cif", str(CIF_DIR/f"{pdb}.cif")])
        
        gen_success = gen_receptor_cif(info['rec_id'], cif, resn, uniprot_id)

        if not gen_success:
            failed.append(info['rec_id'])

        prot_ala_cif = CIF_DIR/f"{info['rec_id']}_prot_ala.cif"
        prot_cif = CIF_DIR/f"{info['rec_id']}_prot.cif"
        if prot_ala_cif.exists() and prot_cif.exists():
            prot_ala_pdb = PDB_DIR/f"{info['rec_id']}_prot_ala.pdb"
            prot_pdb = PDB_DIR/f"{info['rec_id']}_prot.pdb"
            cif_to_pdb(prot_ala_cif, prot_ala_pdb)
            cif_to_pdb(prot_cif, prot_pdb)

    if failed:
        print(f"Failed to generate receptor cif: ({len(failed)})")
        for f in failed:
            print(f)

    filtered = [d for d in cross_set if (d['rec_id'] not in failed) or (d['source'] != 'additional')]
    with open("cross_set.json", "w") as f:
        json.dump(filtered, f, indent=2)


def filter_generated_cross():
    # check whether the prot.cif and prot_ala.cif are empty
    with open("cross_covset.json") as f:
        cov = json.load(f)
    with open("cross_scarpino.json") as f:
        sca = json.load(f)

    filtered_cov = []
    for d in cov:
        rec_id = d['rec_id']
        prot_cif = CIF_DIR/f"{rec_id}_prot.cif"
        prot_ala_cif = CIF_DIR/f"{rec_id}_prot_ala.cif"
        if prot_cif.exists() and prot_ala_cif.exists():
            if prot_cif.stat().st_size < 1024 or prot_ala_cif.stat().st_size < 1024:
                print(f"Filtering out {rec_id} due to empty receptor cif.")
                continue
            else:
                filtered_cov.append(d)

    with open("cross_covset.json", "w") as f:
        json.dump(filtered_cov, f, indent=2)

    filtered_sca = []
    for d in sca:
        rec_id = d['rec_id']
        prot_cif = CIF_DIR/f"{rec_id}_prot.cif"
        prot_ala_cif = CIF_DIR/f"{rec_id}_prot_ala.cif"
        if prot_cif.exists() and prot_ala_cif.exists():
            if prot_cif.stat().st_size < 1024 or prot_ala_cif.stat().st_size < 1024:
                print(f"Filtering out {rec_id} due to empty receptor cif.")
                continue
            else:
                filtered_sca.append(d)

    with open("cross_scarpino.json", "w") as f:
        json.dump(filtered_sca, f, indent=2)


def make_user_provided_ccd():
    af3_fail = "../run_modeling/af3_fail.txt"
    boltz_fail = "../run_modeling/boltz_fail.txt"

    fail_ligs = []
    with open(af3_fail) as f:
        for l in f:
            lig = l.split('_')[1]
            fail_ligs.append(lig)
    with open(boltz_fail) as f:
        for l in f:
            lig = l.split('_')[1]
            fail_ligs.append(lig)
    fail_ligs = set(fail_ligs)

    print(f"Total failed ligands: {len(fail_ligs)}")
    for lig in fail_ligs:
        out_fn = CIF_DIR/f"{lig}.cif"
        url = f"https://files.rcsb.org/ligands/download/{lig}.cif"
        response = requests.get(url)
        if response.status_code != 200:
            print(f"Failed to download {lig} cif from RCSB.")
            continue
        with open(out_fn, "wb") as f:
            f.write(response.content)



if __name__ == "__main__":
    # main()
    # filter_generated()
    # check_cofactors()
    # generate_cofactors_mol2()
    # generate_single_atom_cofactors()
    # merge_fasta_for_msa()
    # receptor_prep_for_cross()
    # filter_generated_cross()
    make_user_provided_ccd()