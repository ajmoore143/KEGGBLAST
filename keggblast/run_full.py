def run_full_pipeline_single(
    ko_id: str,
    species_input: str,
    sequence_type: str = "both",  # "amino", "gene", or "both"
    blast_program: str = "blastp",
    blast_database: str = "nr",
    taxonomy_filter: str = None
):
    """
    Full pipeline for one species: 
    - Fetch KO
    - Match species
    - Get gene list
    - Save FASTAs
    - Run NCBI BLAST

    Args:
        ko_id (str): KEGG Orthology ID (e.g. "K09252")
        species_input (str): Species name or KEGG ID (e.g. "homo sapiens" or "hsa")
        sequence_type (str): "amino", "gene", or "both"
        blast_program (str): NCBI BLAST program (e.g. "blastp")
        blast_database (str): NCBI database (e.g. "nr")
        taxonomy_filter (str): Optional taxonomy filter (e.g. "txid9606[ORGN]")
    """
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

    print(f"\nKO: {ko_id} | Species: {species_input} | Seq: {sequence_type} | BLAST: {blast_program}\n")

    # 1. Fetch KO + gene table
    ko_text = fetch_kegg_orthology(ko_id)
    gene_df = parse_gene_table(ko_text)

    # 2. Load species + match input
    species_df = load_species_data()
    sp_name, sp_id, gene_list = map_species_from_single_input(species_df, gene_df)
    if not gene_list:
        print("❌ No genes found.")
        return

    # 3. Create output dir
    folder = f"fasta_output/{sp_id.lower()}"
    os.makedirs(folder, exist_ok=True)

    for gene in gene_list:
        entry = fetch_gene_entry(f"{sp_id}:{gene}")
        if sequence_type in ["amino", "both"]:
            aa_seq = extract_sequence(entry, "AASEQ")
            if aa_seq:
                write_fasta_file(f"{folder}/{gene}_amino.fasta", gene, aa_seq)
        if sequence_type in ["gene", "both"]:
            nt_seq = extract_sequence(entry, "NTSEQ")
            if nt_seq:
                write_fasta_file(f"{folder}/{gene}_gene.fasta", gene, nt_seq)

    # 4. Collect generated FASTA and run BLAST
    from keggblast.blast_ncbi import collect_fasta_files, read_fasta_sequence
    fasta_files = collect_fasta_files(folder)
    for fasta in fasta_files:
        name = os.path.basename(fasta).replace(".fasta", "")
        seq = read_fasta_sequence(fasta)
        run_ncbi_blast(seq, blast_program, blast_database, name, tax_query=taxonomy_filter)

    
    #5. Parsing JSON into a csv file
    from keggblast.blast_utils import blast_json_to_df

    print("\nParsing BLAST results into tables...")

    blast_dir = "blast_results_ncbi"
    for file in os.listdir(blast_dir):
        if file.endswith(".json"):
            path = os.path.join(blast_dir, file)
            df = blast_json_to_df(path)
            out_csv = path.replace(".json", ".csv")
            df.to_csv(out_csv, index=False)
            print(f"✅ Parsed table saved to: {out_csv}")

    print("All steps complete.")


def run_full_pipeline_csv(
    ko_id: str,
    csv_path: str,
    sequence_type: str = "both",
    blast_program: str = "blastp",
    blast_database: str = "nr",
    taxonomy_filter: str = None
):
    """
    Full pipeline for a CSV of species:
    - Fetch KO + parse table
    - Match species via fuzzy lookup
    - Save FASTAs for all matches
    - Run NCBI BLAST

    Args:
        ko_id (str): KEGG Orthology ID
        csv_path (str): Path to CSV (must contain 'species' column)
        sequence_type (str): "amino", "gene", or "both"
        blast_program (str): NCBI BLAST program
        blast_database (str): NCBI database
        taxonomy_filter (str): Optional taxonomy query (e.g. "txid1239[ORGN]")
    """
    import pandas as pd
    from keggblast.utils import (
        fetch_kegg_orthology, parse_gene_table, load_species_data,
        map_species_from_csv, fetch_gene_entry, extract_sequence, write_fasta_file
    )
    from keggblast.blast_ncbi import run_ncbi_blast, collect_fasta_files, read_fasta_sequence

    # 1. Get gene table
    ko_text = fetch_kegg_orthology(ko_id)
    gene_df = parse_gene_table(ko_text)
    species_df = load_species_data()

    # 2. Match CSV species → KEGG IDs
    match_table = map_species_from_csv(csv_path, species_df, gene_df)

    for _, row in match_table.iterrows():
        sp_name = row['Matched Name']
        sp_id = row['KEGG Species ID']
        genes = row['Genes'].split(';') if row['Genes'] != 'none found' else []
        if sp_id == "no match" or not genes:
            continue

        out_dir = f"fasta_output/{sp_id.lower()}"
        os.makedirs(out_dir, exist_ok=True)

        for gene in genes:
            entry = fetch_gene_entry(f"{sp_id}:{gene}")
            if sequence_type in ["amino", "both"]:
                aa = extract_sequence(entry, "AASEQ")
                if aa:
                    write_fasta_file(f"{out_dir}/{gene}_amino.fasta", gene, aa)
            if sequence_type in ["gene", "both"]:
                nt = extract_sequence(entry, "NTSEQ")
                if nt:
                    write_fasta_file(f"{out_dir}/{gene}_gene.fasta", gene, nt)

    # 3. Run BLAST on saved FASTAs
    all_fastas = collect_fasta_files("fasta_output")
    for fasta in all_fastas:
        gene = os.path.basename(fasta).replace(".fasta", "")
        seq = read_fasta_sequence(fasta)
        run_ncbi_blast(seq, blast_program, blast_database, gene, tax_query=taxonomy_filter)

    #4. Parsing JSON into a csv file
    from keggblast.blast_utils import blast_json_to_df

    print("\nParsing BLAST results into tables...")

    blast_dir = "blast_results_ncbi"
    for file in os.listdir(blast_dir):
        if file.endswith(".json"):
            path = os.path.join(blast_dir, file)
            df = blast_json_to_df(path)
            out_csv = path.replace(".json", ".csv")
            df.to_csv(out_csv, index=False)
            print(f"✅ Parsed table saved to: {out_csv}")

    print("All steps complete.")

