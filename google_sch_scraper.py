# We will make smart scrapper using selenium in this notebook
# Along with the associated functions required.

# The purpose of this python file is to serve as library of functions for our main function

# Importing neccesary libraries
from selenium import webdriver 
from bs4 import BeautifulSoup
import time
import os
import subprocess
import numpy as np
import pandas as pd
from Bio import Entrez, SeqIO
from urllib.request import urlretrieve
import requests
from colorama import Fore
import random



# Function to check if there is captcha present on the called url
def check_captcha(html_soup):
    soup_as_str = str(html_soup)
    
    if 'captcha' in soup_as_str:
        return True
    else:
        return False


# Making the function to generate the url given a pubmed id

def get_paper_url(pubmed_id):
    Entrez.email = "your_email@example.com"  # Provide your email address

    # Fetch the summary information for the specified PubMed ID
    handle = Entrez.esummary(db="pubmed", id=pubmed_id)
    record = Entrez.read(handle)
    handle.close()

    # Extract relevant information such as title and authors
    paper_info = {
        "title": record[0]["Title"],
        "authors": record[0]["AuthorList"],
    }
    
    # Now generating the url for google scholar 
    query = f"{paper_info['title']} {', '.join(paper_info['authors'])}"

    url = f"https://scholar.google.com/scholar?q={'+'.join(query.split())}"
    
    return [url, paper_info['title'],paper_info['authors']]

    
# Making function to hold our program execution for upto 10 seconds so that we 
# get blocked by google scholar less frequently

def pause_for_one_minute():
    random_time = random.randint(5,10)
    print(f"Pausing for {random_time} seconds...")
    time.sleep(random_time)
    print("Resuming normal execution.")    

# Function to call Selenium driver and scrape the citations of the paper

def scrape(url):
    # Making driver using Firefox
    driver = webdriver.Firefox()
    
    # Opening the window
    driver.get(url)
    
    # Inspecting the page to see if their is a captcha
    page_source = driver.page_source
    
    # Using beautiful soup to parse the page
    soup = BeautifulSoup(page_source, 'html.parser')
    
    captcha_present = check_captcha(soup)
    
    if captcha_present == True:
        # Now we need to tell the user to solve the captcha
        user_input = input("Please solve the captcha manually and press enter after you do so:")
        
        # After the user presses enter we will take a timeout 
        # after that we will again inspect the page and get the citations
        pause_for_one_minute()
        
        # Now we inspect the page
        page_source = driver.page_source
        
        soup = BeautifulSoup(page_source, 'html.parser')
        soup_str = str(soup)
        # closing the driver and returning the number of citations
        driver.quit()
#         print(soup)
        return soup_str[soup_str.find("Cited by"):].split()[2][:-4]
    
    else:
        # Incase we don't get a captcha to solve we directly extract the citations and return them
        # first we close the driver
        driver.quit()
#         print(soup)
        soup_str = str(soup)
        return soup_str[soup_str.find("Cited by"):].split()[2][:-4]
        
        
        
    
    
    
