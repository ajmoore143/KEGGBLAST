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

```python
pip install .
```
