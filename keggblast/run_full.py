import os
from keggblast.utils import (
        fetch_kegg_orthology,
        parse_gene_table,
        load_species_data,
        map_species_from_single_input,
        fetch_gene_entry,
        extract_sequence,
        write_fasta_file
    )
from keggblast.blast_ncbi import run_ncbi_blast
from keggblast.json_tools import parse_json_blast_to_table

def run_full_pipeline_single(
    ko_id,
    sequence_type="both",          # "amino", "gene", or "both"
    blast_program="blastp",
    blast_database="nr",
    tax_query=None                 # e.g. "txid4751[ORGN]"
):
    """
    Run the full pipeline for a single species:
    - KO parsing
    - Species matching
    - Gene extraction
    - FASTA download
    - NCBI BLAST
    - JSON-to-DataFrame parsing

    Args:
        ko_id (str): KEGG Orthology ID (e.g. "K09252")
        sequence_type (str): "amino", "gene", or "both"
        blast_program (str): NCBI BLAST program (e.g. "blastp")
        blast_database (str): NCBI database (e.g. "nr")
        tax_query (str): Optional taxonomy filter (e.g. "txid9606[ORGN]")
    """
    print("🚀 Running single-species full pipeline...\n")

    # 1. KO → Gene Table
    ko_entry = fetch_kegg_orthology(ko_id)
    gene_df = parse_gene_table(ko_entry)

    # 2. Match species + extract genes
    species_df = load_species_data()
    matched_name, species_id, gene_list = map_species_from_single_input(species_df, gene_df)

    if not gene_list:
        print("❌ No genes found, aborting.")
        return

    # 3. Fetch sequences and save FASTA
    output_root = "fasta_output"
    sp_dir = os.path.join(output_root, species_id)
    os.makedirs(sp_dir, exist_ok=True)

    for gene_id in gene_list:
        entry = fetch_gene_entry(f"{species_id}:{gene_id}")

        if sequence_type in ("amino", "both"):
            aa_seq = extract_sequence(entry, "AASEQ")
            if aa_seq:
                write_fasta_file(os.path.join(sp_dir, f"{gene_id}_amino.fasta"), gene_id, aa_seq)

        if sequence_type in ("gene", "both"):
            nt_seq = extract_sequence(entry, "NTSEQ")
            if nt_seq:
                write_fasta_file(os.path.join(sp_dir, f"{gene_id}_gene.fasta"), gene_id, nt_seq)

    # 4. Run NCBI BLAST
    run_ncbi_blast_all(
        program=blast_program,
        database=blast_database,
        tax_query=tax_query,
        fasta_dir=sp_dir,
        output_dir="blast_results_ncbi"
    )

    # 5. Parse JSON → DataFrame
    print("\n📊 Parsing BLAST results...")
    df = parse_json_blast_to_table("blast_results_ncbi")
    df.to_csv("single_blast_results.csv", index=False)
    print("✅ Saved BLAST result table to single_blast_results.csv")

def run_full_pipeline_csv(
    csv_path,
    ko_id,
    sequence_type="amino",  # or "gene" or "both"
    blast_program="blastp",
    blast_database="nr",
    tax_query=None
):
    """
    Run the full KO ➡️ FASTA ➡️ BLAST ➡️ table pipeline for a CSV of species.

    Parameters:
        csv_path (str): Path to input CSV with column 'species'
        ko_id (str): KEGG Orthology ID
        sequence_type (str): "amino", "gene", or "both"
        blast_program (str): blastp, blastn, etc.
        blast_database (str): nr, nt, etc.
        tax_query (str): e.g. 'txid4751[ORGN]', or None for global search

    Output:
        Saves FASTA files, JSON BLAST results, and final CSV summary
    """

    # 1. Load KO entry + parse genes
    print(f"📥 Fetching KO: {ko_id}")
    ko_text = fetch_kegg_orthology(ko_id)
    gene_df = parse_gene_table(ko_text)
    species_df = load_species_data()

    # 2. Match species from CSV
    results = map_species_from_csv(csv_path, species_df, gene_df, output_file="species_id_and_genes.csv")

    if results.empty:
        print("❌ No valid species/gene matches.")
        return
        
    print(f"\n📦 Preparing sequences for {len(results)} matched species...")

    for _, row in results.iterrows():
        sp_name = row['species']
        sp_id = row['KEGG Species ID']
        genes = row['Genes'].split(' ') if row['Genes'] != 'none found' else []
    
        for gene_id in genes:
            entry_text = fetch_gene_entry(f"{sp_id}:{gene_id}")
            if sequence_type in ["amino", "both"]:
                aa_seq = extract_sequence(entry_text, "AASEQ")
                if aa_seq:
                    aa_path = f"fasta_output/{sp_id}/{gene_id}_amino.fasta"
                    os.makedirs(os.path.dirname(aa_path), exist_ok=True)
                    write_fasta_file(aa_path, gene_id, aa_seq)
    
            if sequence_type in ["gene", "both"]:
                nt_seq = extract_sequence(entry_text, "NTSEQ")
                if nt_seq:
                    nt_path = f"fasta_output/{sp_id}/{gene_id}_gene.fasta"
                    os.makedirs(os.path.dirname(nt_path), exist_ok=True)
                    write_fasta_file(nt_path, gene_id, nt_seq)

    # 3. Run BLAST
    print("\n🚀 Launching NCBI BLAST...")
    run_ncbi_blast_all(
        program=blast_program,
        database=blast_database,
        tax_query=tax_query,
        fasta_dir="fasta_output",
        output_dir="blast_results_ncbi"
    )

    # 4. Parse JSONs
    print("\n📊 Parsing all BLAST JSONs to table...")
    df = parse_json_blast_to_table("blast_results_ncbi")
    df.to_csv("blast_summary_csv_mode.csv", index=False)
    print("✅ Final BLAST table saved: blast_summary_csv_mode.csv")
