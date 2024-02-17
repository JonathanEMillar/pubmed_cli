# first line: 150
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
