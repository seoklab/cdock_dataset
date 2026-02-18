#!/usr/bin/env python3
"""
JSON Dataset Inspector

A utility tool to inspect and analyze covalent ligand dataset JSON files.
Provides statistics, validation, and export functionality.
"""

import json
import argparse
from pathlib import Path
from collections import Counter
import sys


def load_json(filepath):
    """Load JSON dataset file."""
    try:
        with open(filepath) as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"Error: File not found: {filepath}")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in {filepath}: {e}")
        sys.exit(1)


def print_statistics(data, dataset_type="covalent"):
    """Print comprehensive statistics about the dataset."""
    print(f"\n{'='*70}")
    print(f"Dataset Statistics")
    print(f"{'='*70}\n")
    
    print(f"Total entries: {len(data)}")
    
    if dataset_type == "covalent":
        print_covalent_stats(data)
    elif dataset_type == "pdbbind":
        print_pdbbind_stats(data)
    elif dataset_type == "cross":
        print_cross_docking_stats(data)


def print_covalent_stats(data):
    """Print statistics for covalent dataset."""
    # Basic counts
    pdb_ids = [d['pdb'] for d in data]
    uniprot_ids = [d.get('uniprot_id') for d in data if d.get('uniprot_id')]
    lig_names = [d['lig_name'] for d in data]
    res_names = [d['res_name'] for d in data]
    
    print(f"Unique PDB entries: {len(set(pdb_ids))}")
    print(f"Unique proteins (UniProt): {len(set(uniprot_ids))}")
    print(f"Unique ligands: {len(set(lig_names))}")
    
    # Date range
    dates = [d['date'] for d in data if d.get('date')]
    if dates:
        print(f"\nDate range: {min(dates)} to {max(dates)}")
    
    # Resolution statistics
    resolutions = [d['resolution'] for d in data if d.get('resolution')]
    if resolutions:
        print(f"\nResolution statistics:")
        print(f"  Min: {min(resolutions):.2f} Å")
        print(f"  Max: {max(resolutions):.2f} Å")
        print(f"  Mean: {sum(resolutions)/len(resolutions):.2f} Å")
    
    # Experiment types
    exp_types = Counter(d['exp_type'] for d in data if d.get('exp_type'))
    print(f"\nExperiment types:")
    for exp_type, count in exp_types.most_common():
        print(f"  {exp_type}: {count}")
    
    # Covalent residue distribution
    res_counter = Counter(res_names)
    print(f"\nCovalent residue distribution:")
    for res, count in res_counter.most_common():
        pct = 100 * count / len(data)
        print(f"  {res}: {count} ({pct:.1f}%)")
    
    # Ligand size statistics
    atom_counts = [d['lig_atom_num'] for d in data if d.get('lig_atom_num')]
    if atom_counts:
        print(f"\nLigand size (heavy atoms):")
        print(f"  Min: {min(atom_counts)}")
        print(f"  Max: {max(atom_counts)}")
        print(f"  Mean: {sum(atom_counts)/len(atom_counts):.1f}")
    
    # Rotatable bonds
    n_tors = [d['lig_n_tor'] for d in data if d.get('lig_n_tor') is not None]
    if n_tors:
        print(f"\nRotatable bonds:")
        print(f"  Min: {min(n_tors)}")
        print(f"  Max: {max(n_tors)}")
        print(f"  Mean: {sum(n_tors)/len(n_tors):.1f}")


def print_pdbbind_stats(data):
    """Print statistics for PDBbind dataset."""
    print(f"Unique PDB entries: {len(set(d['pdb'] for d in data))}")
    print(f"Unique ligands: {len(set(d['lig_name'] for d in data))}")
    
    # Year distribution
    years = [d['year'] for d in data if d.get('year')]
    if years:
        print(f"\nYear range: {min(years)} to {max(years)}")
    
    # Resolution
    resolutions = [float(d['resolution']) for d in data 
                  if d.get('resolution') and d['resolution'] != 'NMR']
    if resolutions:
        print(f"\nResolution statistics:")
        print(f"  Min: {min(resolutions):.2f} Å")
        print(f"  Max: {max(resolutions):.2f} Å")
        print(f"  Mean: {sum(resolutions)/len(resolutions):.2f} Å")
    
    # Affinity types
    affinity_types = Counter(d['affinity_type'] for d in data 
                            if d.get('affinity_type'))
    print(f"\nAffinity types:")
    for aff_type, count in affinity_types.most_common():
        print(f"  {aff_type}: {count}")
    
    # Affinity statistics
    affinities = [d['affinity'] for d in data if d.get('affinity')]
    if affinities:
        print(f"\nAffinity statistics (-log scale):")
        print(f"  Min: {min(affinities):.2f}")
        print(f"  Max: {max(affinities):.2f}")
        print(f"  Mean: {sum(affinities)/len(affinities):.2f}")


def print_cross_docking_stats(data):
    """Print statistics for cross-docking dataset."""
    print(f"Total cross-docking pairs: {len(data)}")
    
    # Unique proteins
    uniprot_ids = set(d['uniprot_id'] for d in data if d.get('uniprot_id'))
    print(f"Unique proteins: {len(uniprot_ids)}")
    
    # Source distribution
    sources = Counter(d['source'] for d in data if d.get('source'))
    print(f"\nSource distribution:")
    for source, count in sources.most_common():
        print(f"  {source}: {count}")
    
    # Tanimoto similarity
    tanimotos = [d['tanimoto'] for d in data if d.get('tanimoto')]
    if tanimotos:
        print(f"\nTanimoto similarity:")
        print(f"  Min: {min(tanimotos):.3f}")
        print(f"  Max: {max(tanimotos):.3f}")
        print(f"  Mean: {sum(tanimotos)/len(tanimotos):.3f}")


def validate_dataset(data, dataset_type="covalent"):
    """Validate dataset entries for completeness."""
    print(f"\n{'='*70}")
    print(f"Dataset Validation")
    print(f"{'='*70}\n")
    
    if dataset_type == "covalent":
        required_fields = [
            'pdb', 'date', 'exp_type', 'resolution', 'res_name', 
            'res_atom', 'lig_name', 'lig_atom', 'lig_atom_num', 'lig_n_tor'
        ]
    elif dataset_type == "pdbbind":
        required_fields = ['pdb', 'lig_name', 'resolution', 'year', 'affinity']
    elif dataset_type == "cross":
        required_fields = ['lig_id', 'rec_id', 'uniprot_id', 'tanimoto']
    else:
        print("Unknown dataset type")
        return
    
    missing_counts = {field: 0 for field in required_fields}
    invalid_entries = []
    
    for i, entry in enumerate(data):
        for field in required_fields:
            if field not in entry or entry[field] is None:
                missing_counts[field] += 1
                if i not in invalid_entries:
                    invalid_entries.append(i)
    
    if any(missing_counts.values()):
        print("Found missing or null fields:")
        for field, count in missing_counts.items():
            if count > 0:
                pct = 100 * count / len(data)
                print(f"  {field}: {count} entries ({pct:.1f}%)")
        print(f"\nTotal entries with issues: {len(invalid_entries)}")
    else:
        print("✓ All entries have required fields")
    
    return len(invalid_entries) == 0


def export_summary(data, output_file, dataset_type="covalent"):
    """Export summary to a text file."""
    with open(output_file, 'w') as f:
        # Redirect print to file
        import sys
        old_stdout = sys.stdout
        sys.stdout = f
        
        print_statistics(data, dataset_type)
        validate_dataset(data, dataset_type)
        
        sys.stdout = old_stdout
    
    print(f"\nSummary exported to: {output_file}")


def list_entries(data, limit=10):
    """List first N entries in detail."""
    print(f"\n{'='*70}")
    print(f"Sample Entries (showing first {limit})")
    print(f"{'='*70}\n")
    
    for i, entry in enumerate(data[:limit], 1):
        print(f"Entry {i}:")
        for key, value in entry.items():
            print(f"  {key}: {value}")
        print()


def main():
    parser = argparse.ArgumentParser(
        description="Inspect and analyze covalent ligand dataset JSON files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Show statistics
  python inspect_json.py test_data/cov_set.json --stats
  
  # Validate dataset
  python inspect_json.py test_data/cov_set.json --validate
  
  # List first 5 entries
  python inspect_json.py test_data/cov_set.json --list 5
  
  # Export summary to file
  python inspect_json.py test_data/cov_set.json --stats --export summary.txt
        """
    )
    
    parser.add_argument(
        'input',
        help='Input JSON file to inspect'
    )
    parser.add_argument(
        '--type',
        choices=['covalent', 'pdbbind', 'cross'],
        default='covalent',
        help='Dataset type (default: covalent)'
    )
    parser.add_argument(
        '--stats',
        action='store_true',
        help='Show dataset statistics'
    )
    parser.add_argument(
        '--validate',
        action='store_true',
        help='Validate dataset completeness'
    )
    parser.add_argument(
        '--list',
        type=int,
        metavar='N',
        help='List first N entries in detail'
    )
    parser.add_argument(
        '--export',
        metavar='FILE',
        help='Export summary to text file'
    )
    
    args = parser.parse_args()
    
    # Load data
    data = load_json(args.input)
    
    print(f"\nLoaded: {args.input}")
    print(f"Dataset type: {args.type}")
    
    # Execute requested actions
    if args.stats or args.export:
        print_statistics(data, args.type)
    
    if args.validate or args.export:
        validate_dataset(data, args.type)
    
    if args.list:
        list_entries(data, args.list)
    
    if args.export:
        export_summary(data, args.export, args.type)
    
    # Default: show basic stats if no action specified
    if not any([args.stats, args.validate, args.list, args.export]):
        print_statistics(data, args.type)


if __name__ == '__main__':
    main()
