import requests
import pandas as pd
from collections import defaultdict
import os
from datetime import datetime, timedelta
from rapidfuzz import process

def fetch_kegg_orthology(ko_id):
    '''
    Download a KEGG Orthology entry using the KO ID.

    Parameters:
        ko_id (str): KO ID (e.g., 'K09252')

    Returns:
        str: Raw entry text from KEGG

    Raises:
        ValueError: If the KO ID is invalid or the entry is missing.
        requests.RequestException: For network/server issues.
    '''
    if not ko_id or not ko_id.startswith("K"):
        raise ValueError("❌ KO ID must start with 'K' (e.g. 'K09252')")

    url = f"http://rest.kegg.jp/get/{ko_id}"

    try:
        response = requests.get(url, timeout=10)
    except requests.RequestException as e:
        raise requests.RequestException(f"❌ Network error while fetching KO: {e}")

    # Case 1: Server returned a known error (like 400 or 404)
    if response.status_code != 200:
        raise ValueError(f"❌ KO not found: {ko_id} (HTTP {response.status_code})")

    content = response.text.strip()

    # Case 2: Empty or error page
    if not content or "not found" in content.lower() or "<error>" in content.lower():
        raise ValueError(f"❌ KO entry invalid or not found in response body: {ko_id}")

    # Optional sanity check
    if "ENTRY" not in content:
        raise ValueError(f"❌ KO entry returned but missing expected content: {ko_id}")

    return content

def fetch_gene_entry(kegg_gene_id):
    '''
    Fetch a full KEGG gene entry using species:gene format.

    Parameters:
        kegg_gene_id (str): KEGG gene ID like 'hsa:BRCA1'

    Returns:
        str: Full KEGG gene entry text

    Raises:
        ValueError: If the gene ID is invalid or not found.
        requests.RequestException: On network errors.
    '''
    if not kegg_gene_id or ":" not in kegg_gene_id:
        raise ValueError("❌ KEGG gene ID must be in format 'species_id:gene' (e.g., 'hsa:BRCA1')")

    url = f"http://rest.kegg.jp/get/{kegg_gene_id}"

    try:
        response = requests.get(url, timeout=10)
    except requests.RequestException as e:
        raise requests.RequestException(f"❌ Network error while fetching gene: {e}")

    if response.status_code != 200:
        raise ValueError(f"❌ Gene not found: {kegg_gene_id} (HTTP {response.status_code})")

    content = response.text.strip()

    if not content or "not found" in content.lower() or "<error>" in content.lower():
        raise ValueError(f"❌ Gene entry invalid or missing for: {kegg_gene_id}")

    if "ENTRY" not in content:
        raise ValueError(f"❌ Entry for {kegg_gene_id} is malformed or incomplete.")

    return content

def parse_gene_table(entry_text):
    '''
    Extracts the "GENES" section from a KEGG KO entry and returns
    a DataFrame with columns: Species ID, Genes, Number of Genes.

    Parameters:
        entry_text (str): Full KO entry text from KEGG.

    Returns:
        DataFrame: species_id → genes (str) → count
    '''
    in_genes = False
    data = defaultdict(list)

    for line in entry_text.splitlines():
        if line.startswith("GENES"):
            in_genes = True
            line = line[12:]
        elif in_genes and line.startswith(" " * 12):
            line = line[12:]
        elif in_genes:
            break  # exited the GENES block
        else:
            continue  # skip until GENES section starts

        if ":" in line:
            try:
                species_id, genes = line.strip().split(":", 1)
                gene_list = genes.strip().split()
                data["Species ID"].append(species_id.strip())
                data["Genes"].append(" ".join(gene_list))
                data["Number of Genes"].append(len(gene_list))
            except Exception as e:
                print(f"⚠️ Skipped malformed line: {line} -> {e}")

    if not data:
        raise ValueError("❌ No 'GENES' section found or it was empty in KO entry.")

    return pd.DataFrame(data)

def save_gene_table_as_csv(df, filename="kegg_gene_table.csv"):
    '''
    Save parsed KEGG gene table to a CSV file.

    Parameters:
        df (DataFrame): The gene table to save.
        filename (str): Name of output file (default = 'kegg_gene_table.csv').

    Raises:
        ValueError: If input is not a DataFrame or is empty.
        OSError: If file cannot be written to path.
    '''
    if not isinstance(df, pd.DataFrame):
        raise ValueError("❌ Provided object is not a pandas DataFrame.")

    if df.empty:
        raise ValueError("⚠️ Gene table DataFrame is empty. Nothing to save.")

    # Auto-append .csv if missing
    if not filename.lower().endswith(".csv"):
        filename += ".csv"

    try:
        df.to_csv(filename, index=False)
        print(f"✅ Gene table saved to: {filename}")
    except OSError as e:
        raise OSError(f"❌ Failed to save CSV: {e}")
        
# 📁 Config
SPECIES_CSV = "species_cache.csv"
META_FILE = "species_cache_meta.txt"
UPDATE_INTERVAL_DAYS = 7

def update_species_list(species_csv=SPECIES_CSV, meta_file=META_FILE):
    """
    Fetch the latest KEGG organism list and cache it locally as CSV.

    Parameters:
        species_csv (str): Path to save the species list CSV.
        meta_file (str): Path to store metadata (last updated date).

    Returns:
        pd.DataFrame: Parsed species info.
    
    Raises:
        requests.RequestException: If network/server fails.
        ValueError: If the response format is invalid.
    """
    print("🌐 Fetching latest KEGG species list...")
    url = "http://rest.kegg.jp/list/organism"

    try:
        res = requests.get(url, timeout=10)
        res.raise_for_status()
    except requests.RequestException as e:
        raise requests.RequestException(f"❌ Failed to fetch species list: {e}")

    data = []
    for line in res.text.strip().split('\n'):
        parts = line.split('\t')
        if len(parts) >= 3:
            tax_id = parts[0]
            species_id = parts[1]
            full_name = parts[2].strip()

            # Split full name into Latin + Common name
            if '(' in full_name and ')' in full_name:
                latin_name = full_name.split('(')[0].strip().lower()
                common_name = full_name.split('(')[1].split(')')[0].strip().lower()
            else:
                latin_name = full_name.lower()
                common_name = ""

            data.append({
                'Taxonomy ID': tax_id,
                'Species ID': species_id,
                'Species Name': latin_name,
                'Common Name': common_name
            })

    if not data:
        raise ValueError("❌ Parsed species list is empty. Format may have changed.")

    df_species = pd.DataFrame(data)

    try:
        df_species.to_csv(species_csv, index=False)
        with open(meta_file, "w") as f:
            f.write(datetime.now().strftime("%Y-%m-%d"))
        print(f"✅ Saved KEGG species to: {species_csv}")
        print(f"🕒 Metadata timestamp written to: {meta_file}")
    except Exception as e:
        raise OSError(f"❌ Failed to write cache files: {e}")

    return df_species

def is_cache_stale(species_csv="species_cache.csv", meta_file="species_cache_meta.txt", interval_days=7):
    """
    Determine whether the species cache is stale or missing.

    Parameters:
        species_csv (str): Path to species CSV.
        meta_file (str): Path to metadata file.
        interval_days (int): Age threshold for considering cache stale.

    Returns:
        bool: True if the cache needs refreshing, False otherwise.
    """
    if not os.path.exists(species_csv) or not os.path.exists(meta_file):
        print("❌ Cache files not found. Will update KEGG species list.")
        return True

    try:
        with open(meta_file, "r") as f:
            last_updated = datetime.strptime(f.read().strip(), "%Y-%m-%d")
        age = (datetime.now() - last_updated).days
        if age > interval_days:
            print(f"🔁 Cache is {age} days old. Needs refresh.")
            return True
        else:
            print(f"✅ Cache is fresh (updated {age} days ago).")
            return False
    except Exception as e:
        print(f"⚠️ Error reading metadata: {e}")
        return True

def load_species_data(species_csv="species_cache.csv", meta_file="species_cache_meta.txt", interval_days=7):
    """
    Load the species DataFrame from cache, or update it if stale.

    Parameters:
        species_csv (str): Path to species CSV file.
        meta_file (str): Path to metadata file.
        interval_days (int): Age threshold for refresh (days).

    Returns:
        pd.DataFrame: Cached or freshly pulled KEGG species data.
    
    Raises:
        OSError: If file cannot be read or written.
    """
    try:
        if is_cache_stale(species_csv, meta_file, interval_days):
            return update_species_list(species_csv, meta_file)
        else:
            print("📁 Using cached KEGG species list.")
            df = pd.read_csv(species_csv)
            if df.empty:
                raise ValueError("❌ Species CSV is empty or corrupted.")
            return df
    except Exception as e:
        raise OSError(f"❌ Failed to load species data: {e}")

def map_species_from_single_input(species_df, gene_df):
    """
    Prompt user for a single species name or ID, map it to KEGG ID,
    extract gene list from KO table, and display results.

    Args:
        species_df (DataFrame): Loaded KEGG species list.
        gene_df (DataFrame): Parsed KO gene table.

    Returns:
        tuple: (species_name, species_id, gene_list) or None if failed.
    """
    latin_common_names = []
    species_ids = []

    for _, row in species_df.iterrows():
        if isinstance(row['Species Name'], str):
            latin_common_names.append(row['Species Name'].lower())
        if isinstance(row['Common Name'], str):
            latin_common_names.append(row['Common Name'].lower())
        if isinstance(row['Species ID'], str):
            species_ids.append(row['Species ID'].lower())

    while True:
        species_input = input("\n🧠 Enter your species (e.g. 'homo sapiens', 'mouse', or 'hsa'): ").strip().lower()

        # Direct KEGG ID match
        if len(species_input) == 3 and species_input in species_ids:
            for _, row in species_df.iterrows():
                if species_input == row['Species ID'].lower():
                    species_id = row['Species ID']
                    matched_name = row['Species Name']
                    print(f"\n✅ Matched KEGG ID directly: {species_id} ({matched_name})")
                    break
            break

        # Fuzzy match
        top_matches = process.extract(species_input, latin_common_names, limit=3, score_cutoff=80)
        if not top_matches:
            print("❌ No close matches found. Please try again.")
            continue

        print("\n🤖 Possible matches:")
        for idx, (match, score, _) in enumerate(top_matches, 1):
            print(f"[{idx}] {match} — {score:.1f}%")

        choice = input("❓ Choose correct species number or type 'retry': ").strip().lower()
        if choice == 'retry':
            continue
        elif choice.isdigit() and 1 <= int(choice) <= len(top_matches):
            confirmed = top_matches[int(choice) - 1][0]
        else:
            print("❌ Invalid selection.")
            continue

        for _, row in species_df.iterrows():
            if confirmed == str(row['Species Name']).lower() or confirmed == str(row['Common Name']).lower():
                species_id = row['Species ID']
                matched_name = row['Species Name']
                break

        if species_id:
            print(f"\n✅ Final match: {species_id} ({matched_name})")
            break

    # Extract genes
    filtered_df = gene_df[gene_df['Species ID'].str.lower() == species_id.lower()]
    if filtered_df.empty:
        print("⚠️ No genes found for this species in KO entry.")
        return None

    gene_str = filtered_df.iloc[0]['Genes']
    raw_genes = gene_str.split(' ')
    gene_list = [g.strip().split()[-1] for g in raw_genes]

    print(f"\n🧬 Genes for {matched_name}:")
    for g in gene_list:
        print(f"• {g}")
    print(f"\n✅ Total genes: {len(gene_list)}")

    return matched_name, species_id, gene_list

def map_species_from_csv(csv_path, species_df, gene_df, output_file="species_id_and_genes.csv"):
    """
    Match species in a CSV to KEGG IDs and extract genes.

    Args:
        csv_path (str): Path to CSV with a 'species' column.
        species_df (DataFrame): Loaded KEGG species list.
        gene_df (DataFrame): Parsed KO gene table.
        output_file (str): Filename to save the results.

    Returns:
        DataFrame: Updated with KEGG ID, matched name, and genes.
    """
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"❌ File not found: {csv_path}")

    df_input = pd.read_csv(csv_path)
    if 'species' not in df_input.columns:
        raise KeyError("❌ CSV must have a column named 'species'")

    latin_common_names = []
    for _, row in species_df.iterrows():
        if isinstance(row['Species Name'], str):
            latin_common_names.append(row['Species Name'].lower())
        if isinstance(row['Common Name'], str):
            latin_common_names.append(row['Common Name'].lower())

    matched_ids, matched_names, gene_lists = [], [], []

    for name in df_input['species']:
        result = process.extractOne(name.lower(), latin_common_names, score_cutoff=85)
        if result:
            match, score, _ = result
            matched_id, matched_name = None, None
            for _, row in species_df.iterrows():
                if str(row['Species Name']).lower() == match or str(row['Common Name']).lower() == match:
                    matched_id = row['Species ID']
                    matched_name = row['Species Name']
                    break
            if matched_id:
                filtered_df = gene_df[gene_df['Species ID'].str.lower() == matched_id.lower()]
                if not filtered_df.empty:
                    gene_str = filtered_df.iloc[0]['Genes']
                    raw_genes = gene_str.split(' ')
                    genes = [g.strip().split()[-1] for g in raw_genes]
                else:
                    genes = []
                matched_ids.append(matched_id)
                matched_names.append(matched_name)
                gene_lists.append(genes)
            else:
                matched_ids.append("no match")
                matched_names.append("no match")
                gene_lists.append([])
        else:
            matched_ids.append("no match")
            matched_names.append("no match")
            gene_lists.append([])

    df_input['KEGG Species ID'] = matched_ids
    df_input['Matched Name'] = matched_names
    df_input['Genes'] = [';'.join(genes) if genes else 'none found' for genes in gene_lists]
    df_input.to_csv(output_file, index=False)

    print(f"\n✅ CSV match and gene extraction complete! Output saved to: {output_file}")
    return df_input
