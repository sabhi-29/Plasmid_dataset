# Making a separate python file for fetching plasmid list

# Importing neccesary libraries
import os
import subprocess
import numpy as np
import pandas as pd
from Bio import Entrez, SeqIO
from urllib.request import urlretrieve
import requests
from colorama import Fore
from bs4 import BeautifulSoup
from scholarly import ProxyGenerator
import time
import random



def plasmid_list():
    # Reading the file into dataframe
    plsdb_excel =  'https://docs.google.com/spreadsheets/d/e/2PACX-1vQl2cdQK5RSpUZb1n7PxAXl8a0SNw0THRBC7dJBoU0LhehmPF2Q72QnM5TNjlLlQQ/pub?output=xlsx'
    df_plsdb = pd.read_excel(plsdb_excel, header=0, index_col=False, keep_default_na=True)

    # Extracting the description list
    descriptions = [x.split() for x in df_plsdb['Description_NUCCORE']]

    plasmids = []

    for splits in descriptions:
        try:
            plasmids.append(splits[splits.index('plasmid')+1])
        except:
            continue

    # Keeping only the unique values in the list
    # print(len(plasmids))
    plasmids = list(set(plasmids))
    plasmid_final = [x[:-1] for x in plasmids if len(x) > 1]
    # print(len(plasmids))
    return plasmid_final

def extract_plasmid_names(genbank_file):
    plasmid_names = set()
    with open(genbank_file, "r") as file:
        for record in SeqIO.parse(file, "genbank"):
            for feature in record.features:
                if feature.type == 'source' and 'plasmid' in feature.qualifiers:
                    plasmid_names.add(feature.qualifiers['plasmid'][0])
    return plasmid_names



def plasmid_list_orlek():
    # Getting the genbank file from the system
    orlek_gb_path =  'ncbi_enterobac_download_orlek.gb'
    orlek_plasmid_names = list(extract_plasmid_names(orlek_gb_path))
    
    return orlek_plasmid_names
    
    

def get_data():
    plsdb_plasmid = plasmid_list()
    orlek_plasmid = plasmid_list_orlek()

    # Merging the two lists 
    plasmids = plsdb_plasmid + orlek_plasmid
    plasmids = set(plasmids)
    plasmids = [x for x in plasmids if len(x) > 1]
    return plasmids

