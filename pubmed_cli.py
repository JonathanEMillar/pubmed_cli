import argparse
import concurrent.futures
import sys
import os
import time
import urllib
import webbrowser
import logging
import shutil
import pandas as pd
import yaml
from Bio import Entrez, Medline
from tqdm import tqdm
from joblib import Memory


# --- Custom exceptions -------------------------------------------------------

class ConfigurationError(Exception):
    pass

class FetchError(Exception):
    pass

# --- Configuration -----------------------------------------------------------

def setup_config():
    """
    Sets up the configuration for the script.

    This function reads the configuration from a 'config.yaml' file in the same directory as the script.
    It sets the email, API key, output directory, and query term for the Entrez module.

    Raises:
        ConfigurationError: If the configuration file is not found or an error occurs while loading it.

    Returns:
        dict: The loaded configuration.
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))
    config_path = os.path.join(script_dir, 'config.yaml')

    try:
        with open('config.yaml', 'r') as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        raise ConfigurationError("Configuration file not found. Please ensure a valid 'config.yaml' file is present.")
        exit(1)
    except yaml.YAMLError as err:
        raise ConfigurationError(f"Error occurred while loading the configuration file: {err}")
        exit(1)

    Entrez.email = config['Email']
    Entrez.api_key = config['APIKey']
    output_directory = config['OutputDirectory']
    query_term = config['QueryTerm']

    return config

config = setup_config()

# --- Logging and Cache -------------------------------------------------------

def setup_logger():
    """
    Sets up the logger for the script.

    This function creates a logger named 'pubmedcli_logger' that writes log messages to a file named 'pubmedcli.log'.
    The logger's level is set to DEBUG, so it logs all messages of level DEBUG and higher.
    The log messages are formatted to include the timestamp, logger name, message level, and message.

    Returns:
        logging.Logger: The created logger.
    """
    logger = logging.getLogger('pubmedcli_logger')
    logger.setLevel(logging.DEBUG)  

    handler = logging.FileHandler('pubmedcli.log')
    handler.setLevel(logging.DEBUG)  

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)

    logger.addHandler(handler)

    logger.info('Logger created')

    return logger

logger = setup_logger()
cache_dir = ".article_cache"
os.makedirs(cache_dir, exist_ok=True)
memory = Memory(cache_dir, verbose=0)

# --- Component functions -----------------------------------------------------

def parse_arguments():
    """
    Parses command-line arguments for the PubMed search script.

    The following command-line arguments are defined:
    -d, --days: The number of days back to search in PubMed. It's an integer and the default value is 1.
    -o, --output: A flag indicating whether to output the results to a CSV file. It's a boolean and the default value is False.
    -q, --query: The query term for the PubMed search. It's a string and the default value is taken as the 'QueryTerm' from the configuration file.
    -c, --clearcache: A flag indicating whether to clear the cache before running the script. It's a boolean and the default value is False.

    Returns:
        argparse.Namespace: An object that holds the values of all defined command-line arguments. 
    """
    parser = argparse.ArgumentParser(description='Search PubMed for a specific term and a specific number of days back.')
    parser.add_argument('-d', '--days', type=int, default=1, help='Number of days back to search.')
    parser.add_argument('-o', '--output', action='store_true', help='Output the results to a CSV file.')
    parser.add_argument('-q', '--query', type=str, default=config['QueryTerm'], help='Query term for PubMed search.')
    parser.add_argument('-c', '--clearcache', action='store_true', help='Clear the cache before running.')  
    args = parser.parse_args()
    return args

def validate_args(args):
    """
    Validates the command-line arguments.

    Args:
        args (argparse.Namespace)

    Raises:
        SystemExit: If the number of days is not a positive integer or the query term is empty.
    """
    if args.days <= 0:
        logger.error("The number of days must be a positive integer.")
        print("Error: The number of days must be a positive integer.")
        sys.exit(1)

    if args.query is None or args.query.strip() == "":
        logger.error("The query term must not be empty.")
        print("Error: The query term must not be empty.")
        sys.exit(1)

def clear_cache_if_needed(args):
    """
    Clears the cache directory if the clear cache argument is specified.

    Args:
        args (argparse.Namespace)
    """
    if args.clearcache:
        shutil.rmtree(cache_dir)  
        os.makedirs(cache_dir, exist_ok=True)  
        print("Cache cleared.")
        
@memory.cache
def fetch_article_details(id):
    """
    Fetches the details of a PubMed article given its ID.

    Args:
        id (str): The PubMed ID of the article.

    Returns:
        dict

    Raises:
        FetchError: If an error occurs while fetching the article details.
    """
    retries = 3
    wait_time = 10  
    for i in range(retries):
        try:
            with Entrez.efetch(db='pubmed', id=id, rettype='medline', retmode='text') as fetch_handle:
                record = Medline.read(fetch_handle)
                logger.info(f"Successfully fetched details for ID {id}")
                return record
        except urllib.error.HTTPError as err:
            if err.code == 429:
                if i < retries - 1:  
                    logger.info(f"Fetching. Output available in {wait_time} seconds...")
                    time.sleep(wait_time)  
                    continue
            logger.error(f"HTTP error occurred while fetching details for ID {id}: {err}")
            raise FetchError(f"HTTP error occurred while fetching details for ID {id}: {err}")
        except urllib.error.URLError as e:
            logger.error(f"Network error occurred while fetching details for ID {id}: {e}")
            raise FetchError(f"Network error occurred while fetching details for ID {id}: {e}")
        except Exception as e:
            logger.error(f"Unexpected error occurred while fetching details for ID {id}: {e}")
            raise FetchError(f"Unexpected error occurred while fetching details for ID {id}: {e}")

def fetch_article_ids(args):
    """
    Fetches the IDs of PubMed articles that match the provided query term and date range.

    Args:
        args (argparse.Namespace)

    Returns:
        list
    """
    handle = Entrez.esearch(db='pubmed', term=args.query, reldate=args.days, datetype='pdat', retmax=1000)
    record = Entrez.read(handle)
    idlist = record['IdList']
    print(f"Found {len(idlist)} articles matching the query term")
    logger.info(f"Found {len(idlist)} articles matching the query term")
    return idlist

def fetch_all_article_details(idlist):
    """
    Fetches the details for a list of PubMed articles using multithreading.

    It uses a ThreadPoolExecutor to fetch the details for multiple articles concurrently, which can significantly 
    speed up the process when dealing with a large number of articles. The progress of the fetch operation is 
    displayed using a tqdm progress bar.

    Args:
        idlist (list): A list of PubMed IDs for the articles to fetch.

    Returns:
        list
    """
    with concurrent.futures.ThreadPoolExecutor() as executor:
        records = list(tqdm(executor.map(fetch_article_details, idlist), total=len(idlist), desc="Fetching Details"))
    return records

def print_and_store_results(records, args):
    """
    Prints the details of each publication and optionally stores the results in a CSV file.

    The details include the title, first author, journal, publication date, PubMed ID, DOI URL, and PubMed URL of each publication.
    If the 'output' command-line argument is specified, the function also stores these details in a CSV file.

    Args:
        records (list): A list of dictionaries, where each dictionary contains the details of a PubMed article.
        args (argparse.Namespace)

    Prints:
        The details of each publication.

    Writes:
        If the 'output' command-line argument is specified, writes the details of each publication to a CSV file.
    """
    data = []

    for i, record in enumerate(records, 1):
        if record is None:
            logger.error(f"Failed to fetch details for article {i}")
            continue

        title = record.get("TI", "?")
        first_author = record.get("AU", ["?"])[0]
        journal = record.get("TA", "?")
        pub_date = record.get("DP", "?")
        pubmed_id = record.get("PMID", "?")
        doi = record.get('LID', '?').replace('[doi]', '')  
        doi_url = f"https://doi.org/{doi}"
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

    if args.output:
        df = pd.DataFrame(data)
        output_path = os.path.join(output_directory, 'pubmed_results.csv')
        df.to_csv(output_path, index=False)
        logger.info(f"Results have been written to: {output_path}")
        print(f"Results have been written to: {output_path}")
        
def interact_with_user(records):
    """
    Interacts with the user to display the abstract of a selected article and optionally open its URL.

    Args:
        records (list): A list of dictionaries, where each dictionary contains the details of a PubMed article.

    Raises:
        ValueError: If the user's input cannot be converted to an integer.
    """
    try:
        while True:
            selection = input("Enter the number of the article to view its abstract (or 'q' to quit): ")
            if selection == 'q':
                break
            index = int(selection) - 1
            if 0 <= index < len(records):
                selected_record = records[index]
                pubmed_id = selected_record.get("PMID", "?")
                pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}"  
                abstract = selected_record.get("AB")
                if abstract:
                    print()
                    print(f"\nAbstract of Article {index + 1}:\n{abstract}\n")
                else:
                    print(f"No abstract available for Article {index + 1}\n")
                open_prompt = input("Open article URL in default browser? (y/n): ")
                if open_prompt.lower() == "y":
                    open_in_default_browser(pubmed_url)
                print()  
            else:
                print("Invalid selection. Please enter a valid article number.\n")
    except ValueError:
        print("Invalid input. Please enter a number or 'q' to quit.")
        
def prompt_for_arguments():
    """
    Prompts the user for new arguments and returns them as an argparse.Namespace object.

    Returns:
        argparse.Namespace
    """
    query = input("Enter the query term: ")
    days = int(input("Enter the number of days: "))
    clearcache = input("Clear cache? (y/n): ").lower() == 'y'
    output = input("Enter output file name (or leave blank for no output file): ")
    args = argparse.Namespace(query=query, days=days, clearcache=clearcache, output=output if output else None)
    return args

def open_in_default_browser(url):
    """
    Opens a URL in the default web browser.

    Args:
        url (str): The URL to open.

    Raises:
        logging.error: If the URL format is unsupported.
    """
    if url.startswith("https://doi.org/") or url.startswith("https://pubmed.ncbi.nlm.nih.gov/"):
        webbrowser.open_new(url)
        logger.info(f"Opened URL {url} in default browser")
    else:
        logger.error("Unsupported URL format.")
        
def process_query(args):
    """
    Processes a PubMed query based on the provided arguments.

    This function performs the following steps:
    1. Validates the provided arguments.
    2. Clears the cache if the 'clearcache' argument is specified.
    3. Fetches the IDs of PubMed articles that match the provided query term and date range.
    4. Fetches the details for each of the fetched articles using multithreading.
    5. Prints the details of each publication and optionally stores the results in a CSV file.
    6. Interacts with the user to display the abstract of a selected article and optionally open its URL.

    Args:
        args (argparse.Namespace): The command-line arguments parsed by argparse.

    Raises:
        ValueError: If the arguments are not valid.
        Exception: If there is an error fetching the article IDs or details.
    """
    validate_args(args)
    clear_cache_if_needed(args)
    print()
    idlist = fetch_article_ids(args)
    print()
    records = fetch_all_article_details(idlist)
    print()
    print_and_store_results(records, args)
    interact_with_user(records)
        
# --- Main function -----------------------------------------------------------

def main():
    """
    Executes the main workflow of the PubMed search script.

    This function performs the following steps:
    1. Parses the command-line arguments for the first run.
    2. Processes the PubMed query based on the provided arguments.
    3. Asks the user if they wish to perform another query.
    4. If the user wishes to perform another query, prompts the user for new arguments and repeats the process.
    5. If the user does not wish to perform another query, breaks the loop and ends the program.

    Raises:
        ValueError: If the arguments are not valid.
        Exception: If there is an error fetching the article IDs or details.
    """
    args = parse_arguments()  # Use command-line arguments for the first run
    while True:
        process_query(args)
        another_query = input("Do you wish to perform another query? (y/n): ")
        if another_query.lower() != 'y':
            break
        args = prompt_for_arguments() 

if __name__ == "__main__":
    main()