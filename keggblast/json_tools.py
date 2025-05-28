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
