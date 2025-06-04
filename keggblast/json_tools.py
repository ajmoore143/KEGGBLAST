def blast_json_to_table(json_file):
    '''
    Convert BLAST JSON output (from gget or NCBI if json formatted) to a pandas DataFrame.

    Parameters:
        json_file (str): Path to BLAST .json result file

    Returns:
        pd.DataFrame: Table with all BLAST hits and metadata
    '''
    with open(json_file, 'r') as f:
        data = json.load(f)

    # In gget, the data is a list of hits
    if isinstance(data, list):
        df = pd.json_normalize(data)
        return df
    
    # NCBI might wrap results differently (not typical JSON output though)
    elif isinstance(data, dict) and 'BlastOutput2' in data:
        # advanced parsing would go here
        raise NotImplementedError("NCBI JSON parser not implemented yet.")
    
    else:
        raise ValueError("❌ Unknown BLAST JSON structure")

# Example usage:
#df = blast_json_to_table("blast_results_gget/gene1_blastn_blast.json")
#df

def parse_json_blast_to_table(folder_path):
    '''
    Batch-parse all BLAST result JSON files in a folder into one DataFrame.

    Parameters:
        folder_path (str): Directory containing JSON result files

    Returns:
        pd.DataFrame: Combined results
    '''
    import os
    all_rows = []
    
    for file in os.listdir(folder_path):
        if file.endswith(".json"):
            path = os.path.join(folder_path, file)
            try:
                df = blast_json_to_table(path)
                df["__source_file__"] = file
                all_rows.append(df)
            except Exception as e:
                print(f"⚠️ Skipped {file}: {e}")
    
    return pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
