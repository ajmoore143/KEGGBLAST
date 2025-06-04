# KEGGBLAST: KEGG-Powered Gene Discovery + BLAST Automation
KEGGBLAST is a Python toolkit that automates the process of:

- Extracting gene IDs from KEGG Orthology (KO) entries
- Downloading amino acid + nucleotide sequences from KEGG
-  Matching species via fuzzy name search or KEGG ID (e.g. "hamo sapiens" → hsa)
- Saving gene tables and FASTA files
- Running BLAST using gget or NCBI API
- Supports optional taxonomic filters
- Auto-caches KEGG species dictionary for reuse

## Installation
Clone the repo or unzip the archive.

Then in your terminal or Anaconda Prompt, from the outer keggblast/ folder:

```bash
pip install .
```

## FULL EXAMPLE — End-to-End Pipeline
```python

# 1. === Import tools ===
from keggblast.utils import (
    fetch_kegg_orthology, parse_gene_table, load_species_data,
    map_species_from_single_input, extract_genes_for_species,
    fetch_gene_entry, extract_sequence, write_fasta_file
)
from keggblast.blast_gget import run_gget_blast_all
from keggblast.blast_ncbi import run_ncbi_blast_all

# 2. === Get KO entry and parse ===
ko_entry = fetch_kegg_orthology("K09252")
gene_df = parse_gene_table(ko_entry)

# 3. === Match species ===
species_df = load_species_data()
matched_name, species_id, gene_list = map_species_from_single_input(species_df, gene_df)

# 4. === Download FASTA sequences ===
for gene in gene_list:
    entry = fetch_gene_entry(f"{species_id}:{gene}")
    aa_seq = extract_sequence(entry, "AASEQ")
    nt_seq = extract_sequence(entry, "NTSEQ")

    if aa_seq:
        write_fasta_file(f"fasta_output/{species_id}/{gene}_amino.fasta", gene, aa_seq)
    if nt_seq:
        write_fasta_file(f"fasta_output/{species_id}/{gene}_gene.fasta", gene, nt_seq)

# 5. === Run BLAST (choose one backend) ===

## --- Option A: gget BLAST ---
run_gget_blast_all(program="blastp", database="nr")

## --- Option B: NCBI BLAST (with optional taxonomic filter) ---
run_ncbi_blast_all(
    program="blastp",
    database="nr",
    tax_query="txid4751[ORGN]"  # Fungi taxid
)
```

## 📂 Output Structure
```bash
📁 fasta_output/
   └── hsa/
       ├── GENE1_amino.fasta
       ├── GENE1_gene.fasta
       └── GENE2_amino.fasta

📁 blast_results_gget/
   ├── GENE1_blastp_blast.json

📁 blast_results_ncbi/
   ├── GENE1_ncbi_blast.xml
```

## 🧬 FASTA Extraction Logic
For each KO gene found in a matched species:

Downloads both ``AASEQ`` and ``NTSEQ`` blocks

Writes ``.fasta`` files into folders named after species KEGG ID

## BLAST Support

| Tool     | Interface  | Format     | Taxonomic Filtering      |
| -------- | ---------- | ---------- | ------------------------ |
| gget     | Python API | JSON       | ❌ No                     |
| NCBI API | HTTP POST  | XML / Text | ✅ Yes (`txidXXXX[ORGN]`) |


## Auto Species Matching
Handles typos or fuzzy names. For example:

```text
Input: "hamo sapiens"
→ Match: "homo sapiens"
→ KEGG ID: hsa
```
## Author
Developed by Alex Moore
Mentored by Brady Hislop


