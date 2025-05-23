#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import time
import requests


# In[2]:


def parse_ncbi_blast_text(text):
    """
    Convert plain-text BLAST output to structured list of hits (limited detail).
    Returns list of dicts.
    """
    hits = []
    current = {}

    for line in text.splitlines():
        if line.startswith(">"):
            if current:
                hits.append(current)
            current = {"subject_title": line[1:].strip()}
        elif "Score =" in line:
            parts = line.split(",")
            score = parts[0].split("=")[-1].strip()
            evalue = parts[1].split("=")[-1].strip() if len(parts) > 1 else None
            current["bit_score"] = score
            current["evalue"] = evalue
    if current:
        hits.append(current)

    return hits


# In[3]:


def run_ncbi_blast_all(program="blastp", database="nr", tax_query=None, fasta_dir="fasta_output", output_dir="blast_results_ncbi"):
    """
    Run NCBI BLAST with optional taxonomy filtering and save results as JSON.

    Args:
        program (str): blastp, blastn, blastx, tblastn, tblastx
        database (str): nr, nt, refseq_rna, etc.
        tax_query (str): e.g. 'txid9606[ORGN]'
        fasta_dir (str): Folder with input FASTA files
        output_dir (str): Folder to save results

    Returns:
        None
    """
    os.makedirs(output_dir, exist_ok=True)
    fasta_files = collect_fasta_files(fasta_dir)

    print(f"\n📂 Found {len(fasta_files)} FASTA file(s) for NCBI BLAST.")

    for fasta in fasta_files:
        gene_name = os.path.basename(fasta).replace(".fasta", "")
        print(f"\n🚀 Submitting {gene_name} to NCBI BLAST...")

        try:
            sequence = read_fasta_sequence(fasta)

            params = {
                "CMD": "Put",
                "PROGRAM": program,
                "DATABASE": database,
                "QUERY": sequence
            }
            if tax_query:
                params["ENTREZ_QUERY"] = tax_query

            r = requests.post("https://blast.ncbi.nlm.nih.gov/Blast.cgi", data=params)
            if r.status_code != 200:
                raise Exception("❌ Submission failed")

            rid, rtoe = None, 15
            for line in r.text.splitlines():
                if "RID =" in line:
                    rid = line.split("=")[1].strip()
                if "RTOE =" in line:
                    rtoe = int(line.split("=")[1].strip())
            if not rid:
                raise Exception("❌ RID not found in response.")

            print(f"🧬 RID: {rid} | Estimated wait: {rtoe}s")
            time.sleep(rtoe + 5)

            # Poll for completion
            while True:
                status_check = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi", params={
                    "CMD": "Get",
                    "FORMAT_OBJECT": "SearchInfo",
                    "RID": rid
                })
                if "Status=READY" in status_check.text:
                    print("✅ Search complete.")
                    break
                elif "Status=FAILED" in status_check.text:
                    raise Exception("❌ BLAST job failed.")
                elif "Status=UNKNOWN" in status_check.text:
                    raise Exception("❌ Unknown RID.")
                else:
                    print("⏳ Waiting for result...")
                    time.sleep(30)

            # Get results (TEXT format for easier parsing)
            result = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi", params={
                "CMD": "Get",
                "RID": rid,
                "FORMAT_TYPE": "Text",
                "FORMAT_OBJECT": "Alignment"
            })

            hits = parse_ncbi_blast_text(result.text)

            # Save JSON
            result_path = os.path.join(output_dir, f"{gene_name}_{program}_blast.json")
            with open(result_path, "w") as f:
                json.dump(hits, f, indent=2)
            print(f"💾 Saved: {result_path}")

        except Exception as e:
            print(f"❌ Failed BLAST for {gene_name}: {e}")

