################################################################################
# PubMed ARDS search
# -- Retrieves recent PubMed entries matching the ARDS MeSH term
#    -- runs as CLI
#    -- arguments
#       -d = days          [range == 1-90] (default == 1) 
#       -o = output        prints .csv to downloads folder
#       -c = clear cache   clears cache
#       -q = query         allows you to specify an altrnative query 
################################################################################
# Author: JE Millar
# Date 2023-07-06
# Version: 1.0
# Python Version: 3.10.9
# en: GB UTF-8
################################################################################

import argparse
import concurrent.futures
import pandas as pd
import os
import time
import urllib
from Bio import Entrez, Medline
from tqdm import tqdm
from joblib import Memory
import webbrowser
import logging

# Configure logging
logging.basicConfig(filename='pubmed.log', level=logging.INFO)

# Create a cache directory for storing fetched article details
cache_dir = ".article_cache"
os.makedirs(cache_dir, exist_ok=True)
memory = Memory(cache_dir, verbose=0)

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

def open_in_safari(url):
    if url.startswith("https://doi.org/"):
        webbrowser.get('safari').open_new(url)
    elif url.startswith("https://pubmed.ncbi.nlm.nih.gov/"):
        webbrowser.get('safari').open_new(url)
    else:
        logging.error("Unsupported URL format.")

# Argument parser
parser = argparse.ArgumentParser(description='Search PubMed for a specific term and a specific number of days back.')
parser.add_argument('-d', '--days', type=int, default=1, help='Number of days back to search.')
parser.add_argument('-o', '--output', action='store_true', help='Output the results to a CSV file.')
parser.add_argument('-q', '--query', type=str, default='Respiratory Distress Syndrome, Adult', help='Query term for PubMed search.')
args = parser.parse_args()

try:
    # Replace 'youremail@domain.com' with your actual email
    Entrez.email = 'youremail@domain.com'

    # Insert your NCBI API Key
    Entrez.api_key = 'NCBI_API_Key'

    # Query term
    query_term = args.query

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
        output_path = os.path.expanduser("~/Downloads/pubmed_results.csv")
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

                open_prompt = input("Open article URL in Safari? (y/n): ")
                if open_prompt.lower() == "y":
                    open_in_safari(pubmed_url)

                print()  # Add an empty line for better formatting

            else:
                print("Invalid selection. Please enter a valid article number.\n")
    except ValueError:
        print("Invalid input. Please enter a valid article number or 'q' to quit.\n")

except Exception as e:
    logging.error(f"An error occurred: {e}")