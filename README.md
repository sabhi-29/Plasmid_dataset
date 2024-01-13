# Plasmid Annotation Project

## Introduction

Welcome to the Plasmid Annotation project! This project aims to annotate plasmid data, providing valuable information about resistance genes and plasmid finder details. To run this project on your system, follow the instructions below.

## Prerequisites

Before running the main program (`Plasmid_data_version_3.ipynb`), make sure you have the following libraries installed on your system:

- [abricate](https://github.com/tseemann/abricate)
- [prokka](https://github.com/tseemann/prokka)
- [biopython](https://biopython.org/)
- [selenium](https://www.selenium.dev/)

You can install these libraries using the following commands:

```bash
pip install abricate prokka biopython selenium
```

Just so you know, an environment YAML file will be provided in the future to make it easier for the installation process.

## Project Structure

- **Plasmid_data_version_3.ipynb**: The main Jupyter Notebook file to run the project.
- **abricate_plasmid.py**: Python script containing all the functions to run abricate on plasmid sequences.
- **credentials.py**: File containing email and name credentials for accessing the entrez API.
- **google_sch_scraper.py**: Python script for scraping Google Scholar to get the citation information.
- **plasmid_data_maker.py**: Script containing all the functions creating the plasmid data, integrating all the other functions files.
- **plasmid_list.py**: File to download the list of plasmids  integrating PLSDB and orlek datasets.
- **scrapper.py**: General-purpose scraping functions.

## Project Workflow

1. **Initial Data Retrieval**: Download plasmid list from PLSDB and orlek datasets (10k plasmids).
2. **First order Filtering**: Eliminate plasmids without PubMed IDs using the entrez API.
3. **Data Creation**: Create credentials for entrez API and input them into `plasmid_data_improved` function. Make the dataset for all the plasmids the passed the first order filter, after they pass through plasmid finder filter, and two other filters.
4. **Output Files**: Generate dataframes, plasmid finder output CSV, and resistance genes CSV for each plasmid.
5. **Future Enhancement**: The project is a work in progress, aiming to use a large language model for annotating unseen plasmid nucleotide sequences.

## Contributions

This project is open to contributions. If you have any suggestions or improvements, please feel free to let me know. Your feedback is valuable in enhancing the functionality and efficiency of the Plasmid Annotation project.
