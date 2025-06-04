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
    from keggblast.json_tools import parse_json_blast_to_table
    print("\n📊 Parsing BLAST results...")
    df = parse_json_blast_to_table("blast_results_ncbi")
    df.to_csv("single_blast_results.csv", index=False)
    print("✅ Saved BLAST result table to single_blast_results.csv")
