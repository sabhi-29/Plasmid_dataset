import subprocess
import pandas as pd




def abricate(nuc_seq_fasta_path):
    try:

        # abricate_cmd = f'abricate {input_file_fas} --db ncbi --csv > {input_file_fas}.csv'
        abricate_cmd = f'abricate {nuc_seq_fasta_path} --db ncbi --csv > {nuc_seq_fasta_path}_resistance.csv'
        result = subprocess.check_output(abricate_cmd, shell = True, universal_newlines = True)
        abr_df = pd.read_csv(f'{nuc_seq_fasta_path}_resistance.csv')
        return abr_df
    except:
        print("There is error handling this file. Please see if the path is correct.")
        return 0


def plasmid_finder(nuc_seq_fasta_path):
    try:

        # abricate_cmd = f'abricate {input_file_fas} --db ncbi --csv > {input_file_fas}.csv'
        plasmid_finder_cmd = f'abricate {nuc_seq_fasta_path} --db plasmidfinder --csv > {nuc_seq_fasta_path}_pf.csv'
        result = subprocess.check_output(plasmid_finder_cmd, shell = True, universal_newlines = True)
        pf_df = pd.read_csv(f'{nuc_seq_fasta_path}_pf.csv')
        return pf_df
    
    except:
        print("There is error handling this file. Please see if the path is correct.")
        return 0
        