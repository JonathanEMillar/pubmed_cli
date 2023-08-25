################################################################################

# PubMed search
# -- Retrieves recent PubMed entries matching a MeSH term
#    -- runs as CLI
#    -- arguments
#       -d = days          [range == 1-90] (default == 1) 
#       -q = query         PubMed search term (defaults to config.yaml)
#       -o = output        prints .csv to local directory (path set in config.yaml)
#       -c = clear cache   clears cach

################################################################################
# Author: JE Millar
# Date 2023-08-25
# Version: 1.1
# Python Version: 3.10.9
# en: GB UTF-8
################################################################################

import argparse
import concurrent.futures
import pandas as pd
import os
import time
import urllib
import configparser
import webbrowser
import logging
import shutil
import yaml
from Bio import Entrez, Medline
from tqdm import tqdm
from joblib import Memory

# --- Configuration -----------------------------------------------------------

# Get the directory of the current script
script_dir = os.path.dirname(os.path.realpath(__file__))

# Build the path to the configuration file
config_path = os.path.join(script_dir, 'config.yaml')

# Attempt to load the configuration file
try:
    with open('config.yaml', 'r') as f:
        config = yaml.safe_load(f)
except FileNotFoundError:
    print("Configuration file not found. Please ensure a valid 'config.yaml' file is present.")
    exit(1)
except yaml.YAMLError as err:
    print("Error occurred while loading the configuration file: ", err)
    exit(1)

# Now you can access your parameters
Entrez.email = config['Email']
Entrez.api_key = config['APIKey']
output_directory = config['OutputDirectory']
query_term = config['QueryTerm']

# --- Logging and Cache -------------------------------------------------------

# Configure logging
logging.basicConfig(filename='pubmed.log', level=logging.INFO)

# Create a cache directory for storing fetched article details
cache_dir = ".article_cache"
os.makedirs(cache_dir, exist_ok=True)
memory = Memory(cache_dir, verbose=0)

# --- Functions ---------------------------------------------------------------

@memory.cache
def fetch_article_details(id):
    retries = 3
    wait_time = 10  # Wait time per request in seconds (10 requests per second)
    for i in range(retries):
        try:
            with Entrez.efetch(db='pubmed', id=id, rettype='medline', retmode='text') as fetch_handle:
                record = Medline.read(fetch_handle)
                return record
        except urllib.error.HTTPError as err:
            if err.code == 429:
                if i < retries - 1:  # It's not the final try yet
                    logging.info(f"Fetching. Output available in {wait_time} seconds...")
                    time.sleep(wait_time)  # Wait before retrying
                    continue
            # If it's not a 429 error or we've hit the max number of retries, raise the error
            raise err
        except Exception as e:
            logging.error(f"Error occurred while fetching details for ID {id}: {e}")
            return None

def open_in_default_browser(url):
    if url.startswith("https://doi.org/") or url.startswith("https://pubmed.ncbi.nlm.nih.gov/"):
        webbrowser.open_new(url)
    else:
        logging.error("Unsupported URL format.")

# --- Arguments ---------------------------------------------------------------

# Argument parser
parser = argparse.ArgumentParser(description='Search PubMed for a specific term and a specific number of days back.')
parser.add_argument('-d', '--days', type=int, default=1, help='Number of days back to search.')
parser.add_argument('-o', '--output', action='store_true', help='Output the results to a CSV file.')
parser.add_argument('-q', '--query', type=str, help='Query term for PubMed search.')
parser.add_argument('-c', '--clearcache', action='store_true', help='Clear the cache before running.')  
args = parser.parse_args()

# If clearcache argument is specified, clear the cache
if args.clearcache:
    shutil.rmtree(cache_dir)  # remove cache directory and its contents
    os.makedirs(cache_dir, exist_ok=True)  # create an empty cache directory
    print("Cache cleared.")

# --- Execution ---------------------------------------------------------------

try:
    # Query term
    query_term = args.query if args.query else query_term

    # Use Entrez esearch to find ids of articles matching the query term
    handle = Entrez.esearch(db='pubmed', term=query_term, reldate=args.days, datetype='pdat')
    record = Entrez.read(handle)
    idlist = record['IdList']

    # Print a space
    print()

    # Print header
    print(f"Searching for recent publications for {query_term}...")

    # Print a space
    print()

    # Print the total number of publications found
    print(f"Total number of publications found: {len(idlist)}\n")

    # Fetch details for each id through efetch using multithreading
    with concurrent.futures.ThreadPoolExecutor() as executor:
        records = list(tqdm(executor.map(fetch_article_details, idlist), total=len(idlist), desc="Fetching Details"))

    # Print a space after the progress bar
    print()

    # Prepare a list for storing the results
    data = []

    # Populate the data list and print the results
    for i, record in enumerate(records, 1):
        if record is None:
            continue

        title = record.get("TI", "?")
        first_author = record.get("AU", ["?"])[0]
        journal = record.get("TA", "?")
        pub_date = record.get("DP", "?")
        pubmed_id = record.get("PMID", "?")
        doi_url = f"https://doi.org/{record.get('LID', '?')}"
        pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}"

        data.append({"Title": title, "First Author": first_author, "Journal": journal, "Publication Date": pub_date, "PubMed ID": pubmed_id, "DOI URL": doi_url})

        print(f"{i}. Title: {title}")
        print("   First Author:", first_author)
        print("   Journal:", journal)
        print("   Publication Date:", pub_date)
        print("   PubMed ID:", pubmed_id)
        print("   DOI URL:", doi_url)
        print("   PubMed URL:", pubmed_url)
        print()

    # If output argument is provided, write to a CSV file
    if args.output:
        df = pd.DataFrame(data)
        output_path = os.path.join(output_directory, 'pubmed_results.csv')
        df.to_csv(output_path, index=False)
        print(f"Results have been written to: {output_path}")
    
    try:
        # Article selection and abstract viewing
        while True:
            selection = input("Enter the number of the article to view its abstract (or 'q' to quit): ")
            if selection == 'q':
                break
            index = int(selection) - 1
            if 0 <= index < len(records):
                selected_record = records[index]
                pubmed_id = selected_record.get("PMID", "?")
                pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}"  # update pubmed_url
                abstract = selected_record.get("AB")
                if abstract:
                    print(f"\nAbstract of Article {index + 1}:\n{abstract}\n")
                else:
                    print(f"No abstract available for Article {index + 1}\n")

                open_prompt = input("Open article URL in default browser? (y/n): ")
                if open_prompt.lower() == "y":
                    open_in_default_browser(pubmed_url)

                print()  # Add an empty line for better formatting

            else:
                print("Invalid selection. Please enter a valid article number.\n")
    except ValueError:
        print("Invalid input. Please enter a valid article number or 'q' to quit.\n")

except Exception as e:
    logging.error(f"An error occurred: {e}")