"""Binding Site Quality Checker for Covalent Ligand Dataset.

This module calculates contact ratios between ligands and protein binding sites
to filter out structures with poor binding site quality. It computes the overlap
volume between ligand atoms and nearby protein atoms within a specified distance.
"""

import json
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import MMCIFParser
import pandas as pd


def calc_sphere_overlap_volume(r1, crd1, r2, crd2):
    """Calculate the overlapping volume between two spheres.
    
    Args:
        r1 (float): Radius of the first sphere
        crd1 (np.ndarray): Coordinates of the first sphere center
        r2 (float): Radius of the second sphere
        crd2 (np.ndarray): Coordinates of the second sphere center
    
    Returns:
        float: Overlapping volume between the two spheres
    """
    d = np.linalg.norm(crd1 - crd2)
    if d >= r1 + r2:
        return 0
    if d <= abs(r1 - r2):
        return (4/3) * np.pi * min(r1, r2)**3
    
    part1 = (r1**2 - r2**2 + d**2) / (2*d)
    part2 = (r2**2 - r1**2 + d**2) / (2*d)
    
    h1 = r1 - part1
    h2 = r2 - part2
    
    vol1 = (1/3) * np.pi * h1**2 * (3*r1 - h1)
    vol2 = (1/3) * np.pi * h2**2 * (3*r2 - h2)
    
    return vol1 + vol2

def calc_contact_ratio(cif, info, 
                       contact_distance=4.0,
                       mode='sphere'):
    """Calculate the contact ratio between ligand and protein binding site.
    
    Args:
        cif (str|Path): Path to the mmCIF file
        info (dict): Dictionary containing ligand and protein chain/residue information
        contact_distance (float): Distance cutoff for defining contact (default: 4.0 Å)
        mode (str): Method for calculating contact ratio ('sphere' or 'box')
    
    Returns:
        float: Contact ratio value (higher = better binding site quality)
    """
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure(info['pdb'], cif)
    except Exception as e:
        print(f"Error parsing {cif}: {e}")
        return 0
    
    model = structure[0]
    
    lig_chain_auth = info['lig_chain_auth']
    lig_num_auth = info['lig_num_auth']
    lig_name = info['lig_name']
    
    res_chain_auth = info['res_chain_auth']
    
    try:
        ligand_residue = None
        for chain in model:
            if chain.id == lig_chain_auth:
                for residue in chain:
                    res_id = residue.id[1]
                    if str(res_id) == str(lig_num_auth) and residue.resname.strip() == lig_name.strip():
                        ligand_residue = residue
                        break
                if ligand_residue:
                    break
        
        if ligand_residue is None:
            print(f"Could not find ligand {lig_name} in chain {lig_chain_auth} at position {lig_num_auth}")
            return 0
        
        ligand_atoms = []
        for atom in ligand_residue:
            ligand_atoms.append(atom)
        
        if len(ligand_atoms) == 0:
            print("No ligand atoms found")
            return 0
        
        protein_atoms = []
        for chain in model:
            if chain.id == res_chain_auth:
                for residue in chain:
                    protein_atoms.extend([atom for atom in residue])

        if len(protein_atoms) == 0:
            print("No protein atoms found")
            return 0
        

        lig_crd = np.array([atom.coord for atom in ligand_atoms])
        max_crd = np.max(lig_crd, axis=0)
        min_crd = np.min(lig_crd, axis=0)
        box_min = min_crd - contact_distance
        box_max = max_crd + contact_distance

        box_size = box_max - box_min
        box_volume = box_size[0] * box_size[1] * box_size[2] # 1A grid space

        atom_in_box = []
        for atom in protein_atoms:
            crd = atom.coord
            if (box_min[0] <= crd[0] <= box_max[0] and
                box_min[1] <= crd[1] <= box_max[1] and
                box_min[2] <= crd[2] <= box_max[2]):
                atom_in_box.append(atom)
        
        if mode == 'box':
            contact_ratio = len(atom_in_box) / box_volume
            return contact_ratio
        
        if mode == 'sphere':
            overlap_v = 0
            for atom in atom_in_box:
                crd = atom.coord
                
                for atom in ligand_atoms:
                    l_crd = atom.coord
                    overlap_v += calc_sphere_overlap_volume(
                        contact_distance, l_crd, 1.5, crd)  
        
            lig_volume = 0
            for atom in ligand_atoms:
                lig_volume += (4/3) * np.pi * (contact_distance**3)

            contact_ratio = overlap_v / lig_volume
            return contact_ratio

    except Exception as e:
        print(f"Error calculating contact ratio for {info['pdb']}: {e}")
        return 0


def test(mode='box'):
    good_pdb_s = [
        "5qh8",
        "5rfj",
        "5rfo",
        "5rfp",
        "6ges",
        "7m05",
    ]
    bad_pdb_s = [
        "7duq",
        "5qip",
        "5rek",
        "8p7v",
    ]
    
    good_contact_ratio = []
    for pdb in good_pdb_s:
        cif = f"../../inputs/cif/{pdb}.cif"
        cr = calc_contact_ratio(cif, data[pdb], mode=mode)
        good_contact_ratio.append(cr)

    bad_contact_ratio = []
    for pdb in bad_pdb_s:
        cif = f"../../inputs/cif/{pdb}.cif"
        cr = calc_contact_ratio(cif, data[pdb], mode=mode)
        bad_contact_ratio.append(cr)
    
    print("Good contact ratios:", np.array(good_contact_ratio).mean())
    print("Bad contact ratios:", np.array(bad_contact_ratio).mean())
    for i in zip(good_pdb_s, good_contact_ratio):
        print(f"Good: {i[0]} : {i[1]:.4f}")
    for i in zip(bad_pdb_s, bad_contact_ratio):
        print(f"Bad: {i[0]} : {i[1]:.4f}")

    plt.figure()
    plt.hist(good_contact_ratio, bins=10, alpha=0.5, label="Good")
    plt.hist(bad_contact_ratio, bins=10, alpha=0.5, label="Bad")
    plt.legend()
    plt.xlabel("Contact Ratio")
    plt.ylabel("Frequency")
    plt.title("Contact Ratio Distribution")
    plt.savefig("contact_ratio.png")


def check_binding_site(input_json, output_json, cif_dir, cutoff=0.65):
    """Filter covalent ligand dataset by binding site contact ratio.
    
    Args:
        input_json (str|Path): Path to input JSON file with ligand data
        output_json (str|Path): Path to output JSON file for filtered data
        cif_dir (str|Path): Directory containing mmCIF files
        cutoff (float): Minimum contact ratio threshold (default: 0.65)
    
    Returns:
        list: Filtered dataset with good binding site quality
    """
    with open(input_json) as f:
        data = json.load(f)
    
    cif_dir = Path(cif_dir)
    filtered = []
    
    for i, item in enumerate(data):
        pdb = item['pdb']
        cif = cif_dir / f"{pdb}.cif"
        
        if not cif.exists():
            print(f"Warning: CIF file not found for {pdb}")
            continue
            
        cr = calc_contact_ratio(cif, item, mode='sphere')
        if cr >= cutoff:
            filtered.append(item)
        
        if (i + 1) % 100 == 0:
            print(f"Processed {i + 1}/{len(data)} structures")

    # Statistics
    residue_types = {}
    for item in filtered:
        res = item["res_name"]
        residue_types[res] = residue_types.get(res, 0) + 1

    print(f"\nFiltered {len(filtered)}/{len(data)} structures with contact ratio >= {cutoff}")
    print(f"Residue type distribution: {residue_types}")

    with open(output_json, "w") as f:
        json.dump(filtered, f, indent=2)
    
    print(f"\nFiltered dataset saved to {output_json}")
    return filtered


def main():
    """Command-line interface for binding site quality checking."""
    parser = argparse.ArgumentParser(
        description="Filter covalent ligand dataset by binding site contact ratio"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input JSON file with covalent ligand data"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output JSON file for filtered data"
    )
    parser.add_argument(
        "--cif-dir", "-c",
        required=True,
        help="Directory containing mmCIF files"
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=0.65,
        help="Minimum contact ratio threshold (default: 0.65)"
    )
    
    args = parser.parse_args()
    check_binding_site(args.input, args.output, args.cif_dir, args.cutoff)


if __name__ == "__main__":
    main()
    