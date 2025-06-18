import os
import sys
import numpy as np
import pandas as pd
import OS_Tools
import pdb
import argparse
import csv
import re
from pybiomart import Server
from pyensembl import EnsemblRelease
import primer3

class CommandLine:
    """ Handle the command line, usage and help requests """
    def __init__(self, in_opts=None):
        """ CommandLine constructor: Implements a parser to interpret the command line argv string using argparse. """
        self.parser = argparse.ArgumentParser(description="%(prog)s use Random Forest to classify cells as tranduced or not", 
        add_help=True, 
        prefix_chars='-', 
        usage='python3 %(prog)s [-h] [-v] ...')
     
        ## inputs & outputs
        self.parser.add_argument('-fa', metavar='--fastaFile', action='store', help='ORF fasta file', required=False, type=str, default=None)
        self.parser.add_argument('-bed', metavar='--bedFile', action='store', help='ORF bed file', required=False, type=str, default=None)
        self.parser.add_argument('-title', action='store_true', help='ORF bed file', required=False, default=None)
       
        if in_opts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(in_opts)

def symbol_to_ensembl(species, database, var):
    data = EnsemblRelease(release = 108, species = species)
    ## check if the EnsemblRelease data has been downloaded 
    dbList = subprocess.run(
    ["pyensembl", "list"],
    capture_output=True,
    text=True,
    check=True
)
    if re.search('108', dbList) and re.search(species, dbList):
        x=10
    else: 
        data.download()
        data.index()
    ## initialize ensembl ID column
    database['ensembl_id'] = None
    for i in np.arange(database.shape[0]):
        gene_symbol = database[var][i]
        if isinstance(gene_symbol, str):
            try: 
                database.loc[database[var]==gene_symbol, 'ensembl_id'] = data.gene_ids_of_gene_name(gene_symbol)[0]                    
            except ValueError:
                database.loc[database[var]==gene_symbol, 'ensembl_id'] = None
        else:
            database.loc[database[var]==gene_symbol, 'ensembl_id'] = None
    return database

def main(command=None):
    """ Run specified program(s). """
    if command is None:
        my_command = CommandLine()  # read options from the co(mmand line
    else:
        my_command = CommandLine(command)  # interpret the list passed from the caller of main
    
    fastaFile, bedFile = my_command.args.fa, my_command.args.bed
    # fastaFile = my_command.args.fa
    fastaDat = pd.read_csv(fastaFile, sep=',', header = None)
    bedDat = pd.read_csv(bedFile, sep='\t', header = None)
    bedDat.columns = ['chrom','start','end','name','score','strand']

    headers = []
    sequences = []
    sequence = ""

    # Parse the FASTA file
    for row in fastaDat[0]:
        if row.startswith('>'):  # Header line
            if sequence:  # Save previous sequence before moving to the next
                sequences.append(sequence)
                sequence = ""
            headers.append(row)  # Append the header
        else:
            sequence += row.strip()  # Append the sequence line

    # Append the last sequence (after the loop ends)
    if sequence:
        sequences.append(sequence)

    # Create a new DataFrame with headers and sequences
    fastaDat = pd.DataFrame({
        'Header': headers,
        'Sequence': sequences
    })

    finalDat = pd.DataFrame()
    for row in fastaDat.values:
        gene_symbol = re.sub(">","",row[0])
        seq = row[1]
        if len(seq) >= 18:
            # strand = bedDat[bedDat['name'] == gene_symbol]['strand'].iloc[0]
            if "five_prime_utr" in gene_symbol:
                PRIMER_PICK_LEFT_PRIMER = 1
                PRIMER_PICK_RIGHT_PRIMER = 0
            elif "three_prime_utr" in gene_symbol:
                PRIMER_PICK_LEFT_PRIMER = 0
                PRIMER_PICK_RIGHT_PRIMER = 1
            
            primerParams = {
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_PICK_LEFT_PRIMER': PRIMER_PICK_LEFT_PRIMER,
                'PRIMER_PICK_RIGHT_PRIMER': PRIMER_PICK_RIGHT_PRIMER,
                'PRIMER_PICK_INTERNAL_OLIGO': 0,
                'PRIMER_INTERNAL_MAX_SELF_END': 8,
                'PRIMER_MIN_SIZE': 18,
                'PRIMER_MAX_SIZE': 25,
                'PRIMER_OPT_TM': 60.0,
                'PRIMER_MIN_TM': 57.0,
                'PRIMER_MAX_TM': 63.0,
                'PRIMER_MIN_GC': 20.0,
                'PRIMER_MAX_GC': 80.0,
                'PRIMER_MAX_POLY_X': 100,
                'PRIMER_INTERNAL_MAX_POLY_X': 100,
                'PRIMER_SALT_MONOVALENT': 50.0,
                'PRIMER_DNA_CONC': 50.0,
                'PRIMER_MAX_NS_ACCEPTED': 0,
                'PRIMER_MAX_SELF_ANY': 12,
                'PRIMER_MAX_SELF_END': 8, # the maximum allowable 3'-anchored global alignment score when testing a single primer for self-complementarity

            }
            dat = pd.DataFrame(); i = 0
            while dat.shape[0] == 0:
                if i > 0: 
                    functions = {
                        'PRIMER_MAX_SIZE': lambda x: x + 1 if x < 30 else x,
                        'PRIMER_MIN_TM': lambda x: x - 1 if x > 55 else x,
                        'PRIMER_MAX_TM': lambda x: x + 1 if x < 70 else x
                    }

                    # Apply the functions
                    for key, func in functions.items():
                        if key in primerParams:  # Check if key exists in the data dictionary
                            primerParams[key] = func(primerParams[key])

                t1 = primer3.bindings.design_primers(
                seq_args={'SEQUENCE_ID': gene_symbol, 'SEQUENCE_TEMPLATE': seq}, global_args=primerParams)
                leftDat = pd.DataFrame(t1['PRIMER_LEFT'])
                leftDat['side'] = "left"
                rightDat = pd.DataFrame(t1['PRIMER_RIGHT'])
                rightDat['side'] = 'right'
                dat = pd.concat([dat, leftDat, rightDat])
                dat['gene_symbol'] = gene_symbol
                i+=1

                # Break if 20 iterations are reached
                if i >= 20:
                    break

            finalDat = pd.concat([finalDat, dat])

    finalDat = finalDat[['gene_symbol','side','PENALTY', "SEQUENCE", "COORDS", "TM", "GC_PERCENT", "SELF_ANY_TH", "SELF_END_TH", "HAIRPIN_TH", "END_STABILITY"]]
    spamwriter = csv.writer(sys.stdout, delimiter=',')
    if my_command.args.title:
        spamwriter.writerows(pd.DataFrame(['gene_symbol', 'side', 'PENALTY', "SEQUENCE", "COORDS", "TM", "GC_PERCENT", "SELF_ANY_TH", "SELF_END_TH", "HAIRPIN_TH", "END_STABILITY"]).T.values)
        spamwriter.writerows(finalDat.values)
    else: 
        spamwriter.writerows(finalDat.values)


if __name__ == "__main__":
    main()
