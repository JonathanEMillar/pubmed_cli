# PubMed Search Script

This repository contains a Python script that allows you to search PubMed for publications based on a specific query term and a specified number of days back. The script utilizes the Biopython library to interact with the PubMed database.

## Features

- Search PubMed for recent publications using a specific query term
- Specify the number of days back to search
- Fetch article details including title, authors, journal, publication date, PubMed ID, and DOI URL
- View abstracts of selected articles
- Open article URLs in default browser
- Option to output the results to a CSV file

## Prerequisites

Before running the script, make sure you have the following prerequisites installed:

- Python 3.x
- Biopython
- tqdm
- pandas
- argparse
- concurrent.futures
- urllib
- configparser
- webbrowser
- logging
- shutil
- yaml
- joblib 

## Usage

1. Clone this repository to your local machine or download the script file.
2. Install the required dependencies by running the following command:

   ```shell
   pip install biopython tqdm pandas

3. Create a config.yml and save it in the same directory. This should contain the following:

   ```yaml
   Email: [your email]
   APIKey: [your Pubmed API key]
   ```
   You can obtain an API key from the [NCBI API Key website](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/).
6. Open a terminal or command prompt and navigate to the directory where the script is located.
7. Run the script using the following command:

   ```shell
   python pubmed_cli.py -q "Respiratory Distress Syndrome, Adult" -d 7 -o
   ```
   This command will search PubMed for publications containing the term "Respiratory Distress Syndrome, Adult" within the last 7 days and output the results to a CSV file.
   Use the -q option to specify the query term, -d to specify the number of days back to search, and -o to enable CSV output. You can omit the options to use the default values.
8. The script will display the search results in the terminal or command prompt. If the -o option is used, the results will also be written to a CSV file in the ~/Downloads directory.
9. The program caches search results, saving additional API calls when the query is repeated. To clear the cache use the -c flag.
