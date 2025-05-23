#!/usr/bin/env python
# coding: utf-8

# In[1]:


from gget import blast
import json


# In[2]:


def collect_fasta_files(root="fasta_output"):
    """
    Recursively collect all .fasta files under a directory.
    """
    fasta_files = []
    for root_dir, _, files in os.walk(root):
        for file in files:
            if file.endswith(".fasta") and not file.startswith("."):
                fasta_files.append(os.path.join(root_dir, file))
    return fasta_files


def read_fasta_sequence(file_path):
    """
    Read raw DNA/protein sequence from FASTA file (skipping header).
    """
    with open(file_path, "r") as f:
        lines = f.readlines()
        return ''.join(line.strip() for line in lines if not line.startswith(">"))

def run_gget_blast_all(program="blastp", database="nr", limit=10, fasta_dir="fasta_output", output_dir="blast_results_gget"):
    """
    Run gget BLAST search on all FASTA files inside `fasta_dir`.

    Parameters:
        program (str): One of 'blastp', 'blastn', 'blastx', 'tblastn', 'tblastx'
        database (str): Target BLAST DB ('nr', 'nt', 'swissprot', etc.)
        limit (int): Max number of top hits per sequence
        fasta_dir (str): Root directory of .fasta files
        output_dir (str): Where to save resulting .json files

    Returns:
        None
    """
    os.makedirs(output_dir, exist_ok=True)
    fasta_files = collect_fasta_files(fasta_dir)

    if not fasta_files:
        print("❌ No FASTA files found to BLAST.")
        return

    print(f"\n📁 Found {len(fasta_files)} FASTA file(s) to BLAST using gget.")

    for fasta_path in fasta_files:
        gene_name = os.path.basename(fasta_path).replace(".fasta", "")
        print(f"\n🔬 Running BLAST for: {gene_name}")

        try:
            sequence = read_fasta_sequence(fasta_path)
            if not sequence:
                print(f"⚠️ Empty sequence in: {fasta_path}")
                continue

            # gget BLAST
            result = blast(
                sequence,
                program=program,
                database=database,
                limit=limit,
                save=False,
                json=True,
                verbose=False
            )

            if not result or not isinstance(result, list):
                print(f"⚠️ No BLAST hits returned for {gene_name}")
                continue

            # Save result
            out_path = os.path.join(output_dir, f"{gene_name}_{program}_blast.json")
            with open(out_path, "w") as f:
                json.dump(result, f, indent=2)
            print(f"✅ Saved to: {out_path}")

        except Exception as e:
            print(f"❌ Failed to BLAST {gene_name}: {e}")

