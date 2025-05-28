# KEGGBLAST 

KEGGBLAST is a Python toolkit for:

- Extracting genes from KEGG orthology (KO) entries  
- Matching species by name or KEGG ID (fuzzy search + auto ID conversion)  
- Parsing and exporting gene tables  
- Saving amino acid and gene sequences as FASTA  
- Running BLAST via `gget` or NCBI API (supports taxonomic filters)  
- Auto-updating species dictionary from KEGG

---
## Installation

Clone or unzip the project, then from the outer `keggblast/` folder:

'''bash
pip install .

## Quick Example

'''python
from keggblast.utils import fetch_kegg_orthology, parse_gene_table<br>
from keggblast.fasta_tools import fetch_gene_entry, extract_sequence

# Fetch KO entry
ko_entry = fetch_kegg_orthology("K09252")

# Parse gene table from entry
gene_df = parse_gene_table(ko_entry)

# Get gene sequence
entry = fetch_gene_entry("ncr:NCU09026")
aa = extract_sequence(entry, "AASEQ")
nt = extract_sequence(entry, "NTSEQ")

# Run BLAST (via gget or NCBI)
from keggblast.blast_gget import run_gget_blast_all
run_gget_blast_all(program="blastp", database="nr")

from keggblast.blast_ncbi import run_ncbi_blast_all
run_ncbi_blast_all(program="blastn", database="nt", tax_query="txid9606[ORGN]")

## Output Structure
FASTA saved under:

fasta_output/{species_id}/GENE_amino.fasta
fasta_output/{species_id}/GENE_gene.fasta

BLAST results saved to:

blast_results_gget/
blast_results_ncbi/

## Species Matching
Auto-fuzzy matches species names and IDs (e.g. "hamo sapiens" → "homo sapiens" → hsa). You can run on single species or upload CSV files of multiple species.

## Author
Developed by Alex Moore
Mentored by Brady Hislop 
