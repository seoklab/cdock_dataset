"""Clustering and Filtering Module for Covalent Ligand Dataset.

This module provides tools for:
1. Removing sequence and ligand similarity leakage between train and test sets
2. Clustering proteins by sequence similarity and ligands by chemical similarity
3. Removing redundant structures
4. Analyzing dataset overlap and time-based splits
"""

import json
import argparse
from pathlib import Path
import subprocess
import pandas as pd
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import rdFingerprintGenerator


class CovalentDatasetClustering:
    """Handle clustering and filtering operations for covalent ligand datasets."""
    
    def __init__(self, 
                 cov_set_path="test_data/cov_set.json",
                 scarpino_set_path="train_data/scarpino_set.json",
                 pdbbind_data_path="train_data/pdbbind_data.json",
                 fasta_dir="fasta"):
        """Initialize clustering module with dataset paths.
        
        Args:
            cov_set_path (str): Path to covalent test set JSON
            scarpino_set_path (str): Path to Scarpino training set JSON
            pdbbind_data_path (str): Path to PDBbind data JSON
            fasta_dir (str): Directory containing FASTA files
        """
        self.cov_set_path = Path(cov_set_path)
        self.scarpino_set_path = Path(scarpino_set_path)
        self.pdbbind_data_path = Path(pdbbind_data_path)
        self.fasta_dir = Path(fasta_dir)
        
        self.cov_data = None
        self.scarpino_data = None
        self.pdbbind_data = None
        
    def load_datasets(self):
        """Load all dataset JSON files."""
        print("Loading datasets...")
        
        with open(self.cov_set_path) as f:
            self.cov_data = json.load(f)
        
        with open(self.scarpino_set_path) as f:
            self.scarpino_data = json.load(f)
        
        with open(self.pdbbind_data_path) as f:
            self.pdbbind_data = json.load(f)
        
        print(f"Loaded {len(self.cov_data)} covalent structures")
        print(f"Loaded {len(self.scarpino_data)} Scarpino structures")
        print(f"Loaded {len(self.pdbbind_data)} PDBbind structures")
    
    @staticmethod
    def get_structure_id(data_dict):
        """Generate unique ID for a structure.
        
        Args:
            data_dict (dict): Structure information dictionary
        
        Returns:
            str: Unique structure ID (pdb_ligname_chain)
        """
        pdb_code = data_dict['pdb']
        lig_name = data_dict['lig_name']
        lig_chain = data_dict.get('lig_chain', data_dict.get('lig_chain_auth', 'A'))
        return f"{pdb_code}_{lig_name}_{lig_chain}"
    
    @staticmethod
    def run_mmseqs_search(train_fasta, test_fasta, output_prefix='mmseqs_search', 
                         min_seq_id=0.3, coverage=0.8):
        """Run MMseqs2 sequence similarity search.
        
        Args:
            train_fasta (str|Path): Training set FASTA file
            test_fasta (str|Path): Test set FASTA file
            output_prefix (str): Prefix for output files
            min_seq_id (float): Minimum sequence identity threshold
            coverage (float): Minimum coverage threshold
        
        Returns:
            Path: Path to alignment TSV file
        """
        print(f"Running MMseqs2 search with min_seq_id={min_seq_id}, coverage={coverage}")
        
        train_db = f"{output_prefix}_train_DB"
        test_db = f"{output_prefix}_test_DB"
        result_db = f"{output_prefix}_search_res"
        alignment_file = f"{output_prefix}_alignment.tsv"
        tmp_dir = Path("tmp")
        tmp_dir.mkdir(exist_ok=True)

        # Create databases
        subprocess.run(
            f"mmseqs createdb {train_fasta} {train_db}",
            shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        subprocess.run(
            f"mmseqs createdb {test_fasta} {test_db}",
            shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        # Run search
        subprocess.run(
            f"mmseqs search {test_db} {train_db} {result_db} {tmp_dir} "
            f"--min-seq-id {min_seq_id} -c {coverage}",
            shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        # Convert to readable format
        subprocess.run(
            f"mmseqs convertalis {test_db} {train_db} {result_db} {alignment_file}",
            shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        # Move files to tmp directory
        output_dir = tmp_dir / output_prefix
        output_dir.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            f"mv {output_prefix}* {output_dir}",
            shell=True, check=True
        )

        return output_dir / Path(alignment_file).name

    @staticmethod
    def calc_max_ligand_similarity(query_struct, target_structs):
        """Calculate maximum Tanimoto similarity between query ligand and target ligands.
        
        Args:
            query_struct (dict): Query structure with 'smiles' key
            target_structs (list): List of target structure dictionaries
        
        Returns:
            float: Maximum Tanimoto similarity (0-1)
        """
        query_smiles = query_struct.get("smiles")
        if query_smiles is None:
            return 0.0
        
        try:
            query_mol = Chem.MolFromSmiles(query_smiles)
            if query_mol is None:
                return 0.0
        except Exception:
            return 0.0
        
        mgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        query_fp = mgen.GetFingerprint(query_mol)

        max_sim = 0.0
        for target in target_structs:
            target_smiles = target.get("smiles")
            if target_smiles is None:
                continue
            
            # Exact match
            if query_smiles == target_smiles:
                return 1.0
            
            try:
                target_mol = Chem.MolFromSmiles(target_smiles)
                if target_mol is None:
                    continue
                target_fp = mgen.GetFingerprint(target_mol)
                sim = DataStructs.TanimotoSimilarity(query_fp, target_fp)
                max_sim = max(max_sim, sim)
            except Exception:
                continue
        
        return max_sim

    @staticmethod
    def run_mmseqs_clustering(input_fasta, output_prefix="clustering", 
                             min_seq_id=0.3, coverage=0.8):
        """Run MMseqs2 sequence clustering.
        
        Args:
            input_fasta (str|Path): Input FASTA file
            output_prefix (str): Prefix for output files
            min_seq_id (float): Minimum sequence identity for clustering
            coverage (float): Minimum coverage threshold
        
        Returns:
            Path: Path to cluster TSV file
        """
        print(f"Running MMseqs2 clustering with min_seq_id={min_seq_id}")
        
        db_name = f"{output_prefix}_DB"
        cluster_db = f"{output_prefix}_clu"
        tmp_dir = Path("tmp")
        tmp_dir.mkdir(exist_ok=True)
        tsv_file = f"{output_prefix}_cluster.tsv"

        subprocess.run(
            f"mmseqs createdb {input_fasta} {db_name}",
            shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        subprocess.run(
            f"mmseqs cluster {db_name} {cluster_db} {tmp_dir} "
            f"--min-seq-id {min_seq_id} -c {coverage}",
            shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        subprocess.run(
            f"mmseqs createtsv {db_name} {db_name} {cluster_db} {tsv_file}",
            shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        output_dir = tmp_dir / output_prefix
        output_dir.mkdir(parents=True, exist_ok=True)  
        subprocess.run(
            f"mv {output_prefix}* {output_dir}",
            shell=True, check=True
        )
        
        return output_dir / Path(tsv_file).name

    @staticmethod
    def cluster_ligands(structures, tanimoto_threshold=0.5):
        """Cluster ligands by Tanimoto similarity.
        
        Args:
            structures (list): List of structure dictionaries with 'smiles'
            tanimoto_threshold (float): Similarity threshold for clustering
        
        Returns:
            list: List of clusters (each cluster is a list of structures)
        """
        print(f"Clustering {len(structures)} ligands with threshold {tanimoto_threshold}")
        
        mgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        fps = []
        valid_indices = []
        
        for i, struct in enumerate(structures):
            smiles = struct.get("smiles")
            if smiles:
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        fp = mgen.GetFingerprint(mol)
                        fps.append(fp)
                        valid_indices.append(i)
                except Exception:
                    continue

        n = len(fps)
        visited = [False] * n
        clusters = []

        for i in range(n):
            if visited[i]:
                continue
            
            cluster = [structures[valid_indices[i]]]
            visited[i] = True
            
            for j in range(i + 1, n):
                if visited[j]:
                    continue
                
                sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                if sim >= tanimoto_threshold:
                    cluster.append(structures[valid_indices[j]])
                    visited[j] = True
            
            clusters.append(cluster)

        print(f"Created {len(clusters)} ligand clusters")
        return clusters

    def remove_leakage(self, test_fasta, train_fasta, output_json,
                      seq_id_threshold=0.3, tanimoto_threshold=0.5):
        """Remove test-train leakage based on sequence and ligand similarity.
        
        Args:
            test_fasta (str|Path): Test set FASTA file
            train_fasta (str|Path): Training set FASTA file
            output_json (str|Path): Output path for filtered test set
            seq_id_threshold (float): Sequence identity threshold
            tanimoto_threshold (float): Tanimoto similarity threshold
        
        Returns:
            list: Filtered test set without leakage
        """
        print("Removing train-test leakage...")
        
        alignment_file = self.run_mmseqs_search(
            train_fasta, test_fasta,
            output_prefix='leak_detection',
            min_seq_id=seq_id_threshold
        )
        
        aligned_df = pd.read_csv(
            alignment_file,
            sep='\t',
            names=['query', 'target', 'pident', 'alnlen', 'mismatches',
                   'gapopens', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bitscore']
        )

        # Create data mappings
        test_data = {self.get_structure_id(d): d for d in self.cov_data}
        train_data = {self.get_structure_id(d): d for d in self.scarpino_data}
        
        # Add PDBbind data
        for d in self.pdbbind_data:
            train_data[f"{d['pdb']}_{d['lig_name']}_pdbbind"] = d

        leaked_queries = []
        
        for query_id, group in aligned_df.groupby('query'):
            if query_id not in test_data:
                continue
                
            query_struct = test_data[query_id]
            target_structs = [
                train_data[tid] for tid in group['target'] if tid in train_data
            ]
            
            max_tanimoto = self.calc_max_ligand_similarity(query_struct, target_structs)
            
            if max_tanimoto >= tanimoto_threshold:
                leaked_queries.append(query_id)
        
        filtered_test = [
            d for d in self.cov_data
            if self.get_structure_id(d) not in leaked_queries
        ]

        print(f"Original test set size: {len(self.cov_data)}")
        print(f"Filtered test set size (no leakage): {len(filtered_test)}")
        print(f"Removed {len(leaked_queries)} structures with leakage")

        with open(output_json, "w") as f:
            json.dump(filtered_test, f, indent=2)
        
        print(f"Saved filtered test set to {output_json}")
        return filtered_test

    def remove_redundancy(self, input_fasta, input_json, output_json,
                         seq_id_threshold=0.3, tanimoto_threshold=0.5):
        """Remove redundant structures by clustering.
        
        Args:
            input_fasta (str|Path): Input FASTA file
            input_json (str|Path): Input JSON file
            output_json (str|Path): Output JSON file
            seq_id_threshold (float): Sequence clustering threshold
            tanimoto_threshold (float): Ligand clustering threshold
        
        Returns:
            list: Non-redundant dataset
        """
        print("Removing redundancy...")
        
        cluster_file = self.run_mmseqs_clustering(
            input_fasta,
            output_prefix='redundancy',
            min_seq_id=seq_id_threshold
        )
        
        cluster_df = pd.read_csv(cluster_file, sep='\t', names=['rep', 'member'])

        with open(input_json) as f:
            data = json.load(f)
        
        data_dict = {self.get_structure_id(d): d for d in data}
        
        filtered_data = []
        
        for rep_id, group in cluster_df.groupby('rep'):
            members = group['member'].tolist()
            
            if len(members) <= 1:
                if members[0] in data_dict:
                    filtered_data.append(data_dict[members[0]])
                continue
            
            # Cluster by ligand similarity
            member_structs = [data_dict[m] for m in members if m in data_dict]
            lig_clusters = self.cluster_ligands(member_structs, tanimoto_threshold)
            
            # Take representative from each ligand cluster
            for lig_cluster in lig_clusters:
                filtered_data.append(lig_cluster[0])
        
        print(f"Original dataset size: {len(data)}")
        print(f"Non-redundant dataset size: {len(filtered_data)}")

        with open(output_json, "w") as f:
            json.dump(filtered_data, f, indent=2)
        
        print(f"Saved non-redundant dataset to {output_json}")
        return filtered_data


def main():
    """Command-line interface for dataset clustering and filtering."""
    parser = argparse.ArgumentParser(
        description="Cluster and filter covalent ligand datasets"
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Leakage removal command
    leak_parser = subparsers.add_parser('remove-leakage', 
                                        help='Remove train-test leakage')
    leak_parser.add_argument('--test-fasta', required=True, 
                            help='Test set FASTA file')
    leak_parser.add_argument('--train-fasta', required=True,
                            help='Train set FASTA file')
    leak_parser.add_argument('--output', required=True,
                            help='Output JSON file')
    leak_parser.add_argument('--seq-threshold', type=float, default=0.3,
                            help='Sequence identity threshold (default: 0.3)')
    leak_parser.add_argument('--lig-threshold', type=float, default=0.5,
                            help='Ligand Tanimoto threshold (default: 0.5)')
    
    # Redundancy removal command
    red_parser = subparsers.add_parser('remove-redundancy',
                                       help='Remove redundant structures')
    red_parser.add_argument('--input-fasta', required=True,
                           help='Input FASTA file')
    red_parser.add_argument('--input-json', required=True,
                           help='Input JSON file')
    red_parser.add_argument('--output', required=True,
                           help='Output JSON file')
    red_parser.add_argument('--seq-threshold', type=float, default=0.3,
                           help='Sequence clustering threshold (default: 0.3)')
    red_parser.add_argument('--lig-threshold', type=float, default=0.5,
                           help='Ligand clustering threshold (default: 0.5)')
    
    # Common arguments
    for p in [leak_parser, red_parser]:
        p.add_argument('--cov-set', default='test_data/cov_set.json',
                      help='Covalent test set JSON')
        p.add_argument('--scarpino-set', default='train_data/scarpino_set.json',
                      help='Scarpino training set JSON')
        p.add_argument('--pdbbind-data', default='train_data/pdbbind_data.json',
                      help='PDBbind data JSON')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return
    
    clustering = CovalentDatasetClustering(
        cov_set_path=args.cov_set,
        scarpino_set_path=args.scarpino_set,
        pdbbind_data_path=args.pdbbind_data
    )
    clustering.load_datasets()
    
    if args.command == 'remove-leakage':
        clustering.remove_leakage(
            args.test_fasta,
            args.train_fasta,
            args.output,
            args.seq_threshold,
            args.lig_threshold
        )
    
    elif args.command == 'remove-redundancy':
        clustering.remove_redundancy(
            args.input_fasta,
            args.input_json,
            args.output,
            args.seq_threshold,
            args.lig_threshold
        )


if __name__ == "__main__":
    main()
