import os
import requests

def extract_genes_for_species(gene_df, species_id, verbose=False):
    """
    Extract clean list of gene symbols for a given species from a KEGG KO gene table.

    Args:
        gene_df (DataFrame): Parsed gene table from a KO entry.
        species_id (str): KEGG species ID to filter by (e.g. 'hsa').
        verbose (bool): If True, prints info about gene extraction.

    Returns:
        list: List of gene symbols (strings).
    """
    if not isinstance(species_id, str) or not species_id.strip():
        raise ValueError("❌ species_id must be a non-empty string.")

    if "Species ID" not in gene_df.columns or "Genes" not in gene_df.columns:
        raise KeyError("❌ gene_df must have columns: 'Species ID' and 'Genes'.")

    subset = gene_df[gene_df['Species ID'].str.lower() == species_id.lower()]

    if subset.empty:
        if verbose:
            print(f"⚠️ No gene entry found for species ID: {species_id}")
        return []

    genes_raw = subset.iloc[0]["Genes"]
    if not isinstance(genes_raw, str):
        return []

    # Split by whitespace (not comma!)
    gene_list = [g.strip().split()[-1] for g in genes_raw.strip().split()]
    
    if verbose:
        print(f"🧬 Extracted {len(gene_list)} gene(s) for {species_id}:")
        print(gene_list)

    return gene_list

def extract_sequence(entry_text, key="AASEQ"):
    """
    Extract amino acid or nucleotide sequence from KEGG entry.

    Args:
        entry_text (str): Full KEGG gene entry text.
        key (str): Either "AASEQ" or "NTSEQ".

    Returns:
        str or None: The extracted sequence in uppercase, or None if not found.
    """
    if key not in ["AASEQ", "NTSEQ"]:
        raise ValueError("❌ Key must be 'AASEQ' or 'NTSEQ'")

    lines = entry_text.splitlines()
    seq_lines = []
    capture = False

    for line in lines:
        if line.startswith(key):
            capture = True
            continue
        if capture:
            if not line.startswith(" " * 12):  # block ends if indent is gone
                break
            seq_lines.append(line.strip())

    if seq_lines:
        sequence = ''.join(seq_lines).upper()
        return sequence
    return None

def write_fasta_file(path, header, sequence):
    """
    Write a sequence to a FASTA file with proper wrapping.

    Args:
        path (str): Full path to output .fasta file.
        header (str): FASTA header line (no '>' included).
        sequence (str): DNA or protein sequence (will be wrapped at 70 chars).

    Raises:
        ValueError: If header or sequence are missing/invalid.
        OSError: If file cannot be written.
    """
    if not header or not isinstance(header, str):
        raise ValueError("❌ FASTA header must be a non-empty string.")

    if not sequence or not isinstance(sequence, str):
        raise ValueError("❌ Sequence must be a non-empty string.")

    # Ensure file extension
    if not path.lower().endswith(".fasta"):
        path += ".fasta"

    try:
        with open(path, 'w') as f:
            f.write(f">{header}\n")
            for i in range(0, len(sequence), 70):
                f.write(sequence[i:i+70] + "\n")
        print(f"✅ FASTA saved: {path}")
    except Exception as e:
        raise OSError(f"❌ Could not write FASTA file: {e}")
