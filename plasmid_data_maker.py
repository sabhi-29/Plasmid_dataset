# Keeping the plasmid data extractor function separately

# Making function to create folders to store plasmid details for various hits
# separately, if the plasmid folder has already been created the function returns us the path to the folder

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
import scrapper
import abricate_plasmid



def create_new_folder(folder_name):
    # Check if the folder already exists
    if not os.path.exists(folder_name):
        # The path does not exist so we will be making a new folder
        # Create the folder
        try:
            os.makedirs(folder_name)
            return os.path.abspath(folder_name)
        except:
            # If for some reason we are unable to create a folder due to the name of the plasmid 
            # we will create a folder with a different name and store the hit result there
            if not os.path.exists(f"Miscellaneous_{folder_name}_end"):
                folder_name = f"Miscellaneous_{folder_name}_end"
                os.makedirs(folder_name)
                return os.path.abspath(folder_name)
            else:
#                 print(f"Folder weird_name already exists.")
                return os.path.abspath(f"Miscellaneous_{folder_name}_end")
    else:
#         print(f"Folder '{folder_name}' already exists.")
        return os.path.abspath(folder_name)



############ BETA VERSION OF THE ABOVE FUNCTION (WORK IN PROGRESS)################################

# Writing the function that will pull all the plasmid hits with updated criteria
# Set up a way to find the hits that contain the journalllll of publication and the author, webscrap to find citations
# add the hit to list only if citations > 10

# We also keep a track of which plasmids and accession number pass both of our filters

def plasmid_data_improved(plasmid_name_1, plasmid_email_id, filter_pass_list):

############## REMEMBER TO UNCOMMENT THIS PORTION WHEN RUNNING FINALLY ###################################
    # Making the folder for the plamsmid hits and storing it's path
    folder_name = plasmid_name_1
    if '.' in folder_name:
        folder_name = folder_name.split('.')
        temp_str = ''
        for sub_str in folder_name:
            temp_str += sub_str
            temp_str += '_'
#         folder_name = temp_str
    else:
        temp_str = folder_name
    if '/' in temp_str:
        temp_str_1 = temp_str.split('/')
        temp_str = ''
        for sub_str in temp_str_1:
            temp_str += sub_str
            temp_str += '_'
            
    path = create_new_folder(temp_str)
    print(f"Folder made for {folder_name} as {temp_str}")

##########################################################################################################

    # Setting up our query and accessing Entrez
    query = f"Plasmid+{plasmid_name_1}"
    email = plasmid_email_id
    Entrez.email = email

    # Searching GenBank now
    handle = Entrez.esearch(db = 'nucleotide', term = query, retmax = 100)
    record_id = Entrez.read(handle)['IdList']
    handle.close()
    # print(record_id)
    # Now storing the data we got in form of a dataframe
    # First making the dictionary to store the data

    plasmid_data = {"Accession No.": [], 'Organism': [], 'Topology': [], 'Plasmid': [], 'Nucleotide Sequence': [],
                  'Locus Tag': [], 'Gene Name': [], 'Gene Location': [], 'Gene Product': [], 'Gene Sequence': []}

    # Going through each record we got from the entrez esearch
    for records in record_id:
        # Extracting the genbank record for the id
        
        scrapper.pause_for_one_minute()
        genbank_record = Entrez.efetch(db = 'nucleotide', id = records, rettype = 'gb', retmode = 'text')
        # Reading the record
        gb_record = SeqIO.read(genbank_record, 'genbank')
        # For each record that we retrieve we will perform some filtering to get the nucleotide sequence of our 
        # interest - seeing if the decription contains the plasmid name of interest and if it is complete or not
        
        # First order filtering (FILTER 1) - the record is non-empty and contains the keyword 'complete sequence'
        # along with the plasmid name
        if (len(gb_record.id) != 0 and len(gb_record.description) != 0) and len(gb_record.seq) != 0:
            
            # We will check if the sequence is of our plasmid of interest and complete or not, if yes then we proceed
            # Extracting Accession Number and description
            
            
            accession = [str(gb_record.annotations['accessions'][0])]
            description = gb_record.description

            ##### DEBUGGING PRINT STATEMENT ##############################
            print(f"Current hit: {plasmid_name_1} :- {description}")
            ################################################################


            # We are going to add the record only if it has the plasmid of interest in the description + is a complete sequence
            if plasmid_name_1.lower() in description.lower() and ('complete sequence' in description.lower()):
                
                ##### DEBUGGING PRINT STATEMENT ##############################
#                 print(gb_record.annotations.get("references")[0].pubmed_id)
                ################################################################
                
                # SECOND ORDER FILTERING (FILTER 2.1) - the record must have a PubMed ID
                if not gb_record.annotations.get("references")[0].pubmed_id:
                    print(f"The hit {accession[0]} does not carry a pubmed id....... SKIPPING!")
                    continue
                # Now we are sure that the record contains the plasmid name, is complete and has a pubmed id
                
                
                # Now we know that the sequence has passed FILTER 1 and FILTER 2 we proceed making a dataframe. 
                organism = [str(gb_record.annotations["source"])]
                plasmid_name = [str(description.split("plasmid ")[-1].split(",")[0])]
                locus_tag = []
                topology = [str(gb_record.annotations['topology'])]
                gene_prdt = []
                gene_seq = []
                location = []
                gene_name = []

                ###########################################
                # Getting nucleotide sequence
                handle = Entrez.efetch(db="nucleotide", id=accession[0], rettype="fasta", retmode="text")
                fasta_sequence = handle.read()
                handle.close()
                fasta_seq_split = [x for x in fasta_sequence.split()]
                # To get the nucleotide sequence we get everything after 'complete sequence'
                # Finding complete sequence
                i = 1
                while i < len(fasta_seq_split) and fasta_seq_split[i-1] != 'sequence':        
                    i += 1
                
                if i >= len(fasta_seq_split):
                    print("Invalid fasta file, skipping this hit.")
                    continue
                
                
                
                # Now we are at the start of the nucleotide sequence we store all that follows
                
                nuc_seq = ''
                while i < len(fasta_seq_split):
                    nuc_seq += fasta_seq_split[i]
                    i += 1
                # print(nuc_seq)
                nucleotide_seq = [str(nuc_seq)]
                ###########################################################
                
                # Before proceeding further there is one last filter (FILTER 3)
                # That the hit has to pass in order to be considereddd
                # THAT IS - (THE RESULT FROM PLASMID FINDER MUST BE NON-EMPTY) 
################################################################################################################################                
                # [VERIFY THIS FROM PROF. bCOZ PLASMIDFINDER DATASET IS SMALL]
                
                # First saving the nucleotide sequence of the record as fasta file
                fasta_path_pla_name = plasmid_name[0]
                fasta_path_name = fasta_path_pla_name
                if '/' in fasta_path_pla_name:
                    fasta_path_name = ''
                    for sub_str in fasta_path_pla_name.split('/'):
                        fasta_path_name += sub_str
                        fasta_path_name += '_'
                
           
                fasta_path = f"{path}/{accession[0]}_{fasta_path_name}.fa"
                with open(fasta_path, 'w') as file:
                    file.write(fasta_sequence)
                # Now using abricate.plasmid_finder on this saved file
                pf_df = abricate_plasmid.plasmid_finder(fasta_path)
                if type(pf_df) == int:
                    print("The response was invalid")
                    continue
                
                # Checking if we got a non-empty result
                if pf_df.shape[0] == 0:
                    print(f"Plasmid Finder returns empty for the current hit! Ac. No: {accession[0]}_{plasmid_name[0]}")
                    print("Skipping the current hit!")
                    os.remove(f"{path}/{accession[0]}_{fasta_path_name}.fa_pf.csv")
                    continue
                 
                # Getting Abricate result similiarly
                abr_df = abricate_plasmid.abricate(fasta_path)
                if type(abr_df) == int:
                    print("The response was invalid")
                    continue
                
                # Checking if we got a non-empty result
                if abr_df.shape[0] == 0:
                    print("Empty dataframe from ABRICATE, hence deleting the resulting ABRICATE .csv file")
                    os.remove(f"{path}/{accession[0]}_{fasta_path_name}.fa_resistance.csv")
  
                
###############################################################################################################################                
                # SECOND ORDER FILTERING (FILTER 2.2) - the record must have more than 10 citations
                
                # Storing the pubmed id
                paper_pubmed_id = gb_record.annotations.get("references")[0].pubmed_id
                # Getting the paper title and the authors using previously defined function
                
                url, paper_title, paper_authors = scrapper.get_paper_url(paper_pubmed_id)
            ####################################################################################################
#                 paper_info = get_paper_info(paper_pubmed_id)
#                 paper_title = paper_info['title']
#                 paper_authors = paper_info['authors']
#                 print(f"{paper_title} {paper_authors}")
            ####################################################################################################

                # Now getting the number of citations of the paper
                paper_citations = int(scrapper.scrape(url))
#                 paper_citations = int(get_citations_count_from_google_scholar(paper_title, paper_authors))
                # Checking if the citations > 10
#                 print(f"{paper_title} {paper_authors}")
                print(f"The number of citations of the paper having this hit is {paper_citations}.")
                print(f"The paper information is as follows - {paper_title} by {', '.join(paper_authors)}")
    
    
                
                if paper_citations < 10:
                    print('The hit contains less than 10 citations..... SKIPPING!')
                    continue
                
#                 print(f"Current hit is published in {paper_title} by {', '.join(paper_authors)} and has {paper_citations} citations.")
                
    

    
                # Now we continue, since we successfully know that the nucleotide sequence associated with the pubmed id
                # has > 10 citations as well, we will keep a record of the plasmid name and accession number that
                # surpass both the set filters. We will use it when using Abricate and Plasmid Finder
                
                if str(description.split("plasmid ")[-1].split(",")[0]) in filter_pass_list.keys():
                    filter_pass_list[str(description.split("plasmid ")[-1].split(",")[0])].append(accession[0])
                else:
                    filter_pass_list[str(description.split("plasmid ")[-1].split(",")[0])] = []
                
                
        
                # Now we will access the record and store the bacteria name, accession number,
                # topology, plasmid, nucleotide sequence, locus tag for each gene with gene name present
                # and store their respective names, location, gene product thier sequence

                
                # First we access the features and get all the data of interest
                for features in gb_record.features:
                    if features.type == 'CDS':
                        # We first see if the qualifiers contains gene, if not then we skip that gene
                        if 'gene' in features.qualifiers.keys():
                            # In case translation, locus tag, gene product is missing then we simply leave that cell empty
                            # print(features.qualifiers.keys())
                            # Storing the gene name first
                            gene_name.append(str(features.qualifiers.get('gene')[0]))
                            # Now storing other data - locus tag, gene product, translation, sequence
                            if 'locus_tag' in features.qualifiers.keys():
                                locus_tag.append(str(features.qualifiers.get('locus_tag')[0]))
                            else:
                                locus_tag.append('NaN')
                            if 'product' in features.qualifiers.keys():
                                gene_prdt.append(str(features.qualifiers['product'][0]))
                            else:
                                gene_prdt.append('NaN')
                            if 'translation' in features.qualifiers.keys():
                                gene_seq.append(str(features.qualifiers['translation'][0]))
                            else:
                                gene_seq.append('NaN')
#                             print(gene_name)
                        # Storing the location of the gene on the nucleotide sequenceeeeeeee
                    if features.type == 'gene':
                        if 'gene' in features.qualifiers.keys():
                            location.append(str(features.location))
#                             print(location)
                
                # We will make the dataframe and store it only if all the columns are of the same length
#                 print(len(location), len(locus_tag))
                if (len(locus_tag) == len(location) and len(locus_tag)!= 0):
                    plasmid_data["Accession No."] = accession*len(locus_tag)
                    plasmid_data["Organism"] = organism*len(locus_tag)
                    plasmid_data["Topology"] = topology*len(locus_tag)
                    plasmid_data["Plasmid"] = plasmid_name*len(locus_tag)
                    plasmid_data["Nucleotide Sequence"] = nucleotide_seq*len(locus_tag)
                    plasmid_data["Locus Tag"] = locus_tag
                    plasmid_data["Gene Length"] = [len(x) for x in gene_seq]
                    plasmid_data["Gene Name"] = gene_name
                    plasmid_data["Gene Sequence"] = gene_seq
                    plasmid_data["Gene Product"] = gene_prdt
                    plasmid_data['Gene Location'] = location
                
                
#                 print(plasmid_data)
                # Making the dataframe
                dataframe = pd.DataFrame(plasmid_data)
                
                # Saving the dataframe
                dataframe.to_excel(f'{path}/{accession[0]}_{plasmid_name[0]}_data.xlsx', index = False)
                print(f"DataFrame made for {plasmid_name[0]} at /{accession[0]}")
                
                # And we return the path of the directory where we have saved all the hits for a particular plasmid
                # In our main program we will delete the directory if it is empty so that we don't create many folders

####################################### DONT MODIFY THE CODE BELOW [WORK IN PROGRESS(NOV-29'23)#################################

                return path

                

# Now we have to integrate Abricate and plasmid finder with our resulting dataframes we made
                
# #                 dataframe.head()
# #             else:
# #                 print("The hit does not pass the filter")