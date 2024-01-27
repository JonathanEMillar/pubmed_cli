# PubMed command line interface

This repository contains a command line interface designed to search PubMed for a specific query term and for a specified number of prior days. I use this as a tool to regularly view publications in my field. Pubmed cli is The program is written in python and uses the Biopython library to interact with the PubMed database. 

## Features

- Search PubMed for recent publications using a specific query term
- Specify the number of days back to search
- Fetch article details including title, authors, journal, publication date, PubMed ID, and DOI URL
- Option to output the results to a .csv file
- View abstracts of selected articles
- Open article URLs in default browser

## Prerequisites

Before running the script, install the dependencies:

```
pip install -r requirements.txt
```

## Installation

1. Clone this repository to your local machine or download the `pubmed_cli.py` script.
2. Create a config.yaml and save it in the same directory. This should contain the following:

   ```yaml
   Email: [your email]
   APIKey: [your Pubmed API key]
   OutputDirectory: ["/Path/to/output/directory"]
   QueryTerm: ["example query"]
   ```
   You can obtain an API key from the [NCBI API Key website](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/).

## Use

1. Open a terminal and navigate to the directory where the script is located. 
2. Run the script using the following command:

   ```shell
   python pubmed_cli.py -q "Respiratory Distress Syndrome, Adult" -d 7 -o
   ```
   This command will search PubMed for publications containing the term "Respiratory Distress Syndrome, Adult" within the last 7 days and output the results to a CSV file.

   Flags:

    -d, --days: The number of days back to search in PubMed. It's an integer and the default value is 1.

    -o, --output: A flag indicating whether to output the results to a CSV file. It's a boolean and the default value is False.

    -q, --query: The query term for the PubMed search. It's a string and the default value is taken as the 'QueryTerm' from the configuration file.

     -c, --clearcache: A flag indicating whether to clear the cache before running the script. It's a boolean and the default value is False.

   A default search term can be set in the config.yaml.
4. The script will display the search results in the terminal. 
5. The program caches search results, saving additional API calls when the query is repeated. To clear the cache use the -c flag.
6. Individual record can be selected to view the abstract.
7. An option is available to open the slected records PubMed entry in the deafult browser.

## To-Do

- [x] Add ability to try addiitonal queries without restarting the program
- [x] Create unit tests
