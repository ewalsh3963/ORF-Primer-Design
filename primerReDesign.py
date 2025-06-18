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


class CommandLine:
    """ Handle the command line, usage and help requests """
    def __init__(self, in_opts=None):
        """ CommandLine constructor: Implements a parser to interpret the command line argv string using argparse. """
        self.parser = argparse.ArgumentParser(description="%(prog)s use Random Forest to classify cells as tranduced or not", 
        add_help=True, 
        prefix_chars='-', 
        usage='python3 %(prog)s [-h] [-v] ...')
     
        ## inputs & outputs
        self.parser.add_argument('-gtf', metavar='--gtfFile', action='store', help='GTF file', required=False, type=str, default=None)
        self.parser.add_argument('-fa', metavar='--fastaFile', action='store', help='ORF fasta file', required=False, type=str, default=None)
        self.parser.add_argument('-pd', metavar='--primerDat', action='store', help='ORF bed file', required=False, type=str, default=None)
        self.parser.add_argument('-bed', metavar='--bedFile', action='store', help='ORF bed file', required=False, type=str, default=None)

        if in_opts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(in_opts)

def symbol_to_ensembl(species, database, var):
    data = EnsemblRelease(release = 108, species = species)
    ## check if the EnsemblRelease data has been downloaded 
    dbList = subprocess.run(["pyensembl", "list"], shell = False)
    if re.search('108', dbList) and re.search(species, dbList):
        x=10
    else: 
        # print("Ensembl database NOT found...downloading and indexing now")
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
        my_command = CommandLine()  # read options from the command line
    else:
        my_command = CommandLine(command)  # interpret the list passed from the caller of main
    
    gtfFile, fastaFile, primerDat, bedFile = my_command.args.gtf, my_command.args.fa, my_command.args.pd, my_command.args.bed
    
    ## read in GTF file 
    gtf = pd.read_csv(gtfFile, sep='\t', header = None, comment = "#")
    gtf.columns = ['seqname','source','feature','start','end','score','strand','frame','attribute']
    
    fastaDat = pd.read_csv(fastaFile, sep=',', header = None)
    primerDat = pd.read_csv(primerDat, sep=',', header = 0)
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

    fastaDat['Header'] = fastaDat['Header'].apply(lambda x: re.sub('>','',x))
    # fastaDat['Header'] = fastaDat['Header'].apply(lambda x: re.sub('\(.*','',x))

    primerHits = primerDat['gene_symbol'].to_list()
    primerScans = fastaDat['Header'].to_list()
    primerMisses = [x for x in primerScans if x not in primerHits]

    bedFinal = pd.DataFrame()
    for miss in primerMisses:
        miss =  re.sub(r'\(.*', '', miss)
        tmpDat = bedDat[bedDat['name'] == miss]

        if 'five_prime_utr' in miss:
            if tmpDat['strand'].iloc[0] == '-':
                tmpDat['start'] = tmpDat["start"] - 25
            else:
                tmpDat['end'] = tmpDat["end"] + 25
        elif 'three_prime_utr' in miss:
            if tmpDat['strand'].iloc[0] == '-':
                tmpDat['end'] = tmpDat["end"] + 25
            else:
                tmpDat['start'] = tmpDat["start"] - 25            
        bedFinal = pd.concat([bedFinal, tmpDat])
    bedFinal = bedFinal.sort_values('name') 

    spamwriter = csv.writer(sys.stdout, delimiter='\t')
    spamwriter.writerows(bedFinal.values)        

if __name__ == "__main__":
    main()




# python3 /home/ewalsh/PrimerDesign/primerReDesign.py \
# -gtf /ds-workspace/EW-TempDataStore/ORF/Mus_musculus/Mus_musculus.GRCm39.113.gtf.gz \
# -fa /ds-workspace/EW-TempDataStore/ORF/Mus_musculus/ORF_PrimaryDesign.fa \
# -pd /ds-workspace/EW-TempDataStore/ORF/Mus_musculus/ORF_primers_PrimaryDesign.csv \
# -bed /ds-workspace/EW-TempDataStore/ORF/Mus_musculus/features.bed
