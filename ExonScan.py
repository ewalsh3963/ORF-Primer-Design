import os
import sys
import numpy as np
import pandas as pd
import OS_Tools
import pdb
import argparse
import concurrent.futures
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
        self.parser.add_argument('-fa', metavar='--fastaFile', action='store', help='fasta file', required=False, type=str, default=None)
        self.parser.add_argument('-gtf', metavar='--gtfFile', action='store', help='GTF file', required=False, type=str, default=None)
        self.parser.add_argument('-o', metavar='--ORF_primers', action='store', help='ORF primers file', required=False, type=str)
        self.parser.add_argument('-s', metavar='--species', action='store', help='', required=False, type=str)

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

def CreateBedData(gtf, row):
    feature, feature_type, strand = row['gene'], row['missing_feature'], row['strand']
    feature_gtf = gtf[gtf['attribute'].str.contains('"' + feature + '"', case=True, na=False)]
    if feature_gtf.shape[0] > 0:
        feature_gtf['transcript_biotype'] = feature_gtf['attribute'].apply(lambda x: (match.group(1) if (match := re.search(r'transcript_biotype "([^"]+)"', x)) else None))
        feature_gtf['gene_name'] = feature_gtf['attribute'].apply(lambda x: (match.group(1) if (match := re.search(r'gene_name "([^"]+)"', x)) else None))
        feature_gtf['gene_id'] = feature_gtf['attribute'].apply(lambda x: (match.group(1) if (match := re.search(r'gene_id "([^"]+)"', x)) else None))
        feature_gtf['transcript_id'] = feature_gtf['attribute'].apply(lambda x: (match.group(1) if (match := re.search(r'transcript_id "([^"]+)"', x)) else None))

        if strand == "+" and feature_type == 'three_prime_utr':
            exon_dat = feature_gtf[(feature_gtf['end'] == feature_gtf['end'].max()) & (feature_gtf['feature'] == 'exon')]
            exon_dat['length'] = exon_dat['end'] - exon_dat['start']
            if any(exon_dat['length'] > 25):
                exon_dat['start'] = exon_dat['end'] - 25
                feature_bed = exon_dat[['seqname','start','end','gene_name','score','strand']].drop_duplicates()
                feature_bed['gene_name'] = feature_bed['gene_name'] + "_" + feature_type
            else:
                feature_bed = exon_dat[['seqname','start','end','gene_name','score','strand']].drop_duplicates()
                feature_bed['gene_name'] = feature_bed['gene_name'] + "_" + feature_type
        elif strand == "+" and feature_type == 'five_prime_utr': 
            exon_dat = feature_gtf[(feature_gtf['start'] == feature_gtf['start'].min()) & (feature_gtf['feature'] == 'exon')]
            exon_dat['length'] = exon_dat['end'] - exon_dat['start']
            if any(exon_dat['length'] > 25):
                exon_dat['end'] = exon_dat['start'] + 25
                feature_bed = exon_dat[['seqname','start','end','gene_name','score','strand']].drop_duplicates()
                feature_bed['gene_name'] = feature_bed['gene_name'] + "_" + feature_type
            else:
                feature_bed = exon_dat[['seqname','start','end','gene_name','score','strand']].drop_duplicates()
                feature_bed['gene_name'] = feature_bed['gene_name'] + "_" + feature_type
        elif strand == "-" and feature_type == 'three_prime_utr':
            exon_dat = feature_gtf[(feature_gtf['start'] == feature_gtf['start'].min()) & (feature_gtf['feature'] == 'exon')]
            exon_dat['length'] = exon_dat['end'] - exon_dat['start']
            if any(exon_dat['length'] > 25):
                exon_dat['end'] = exon_dat['start'] + 25
                feature_bed = exon_dat[['seqname','start','end','gene_name','score','strand']].drop_duplicates()
                feature_bed['gene_name'] = feature_bed['gene_name'] + "_" + feature_type
            else:
                feature_bed = exon_dat[['seqname','start','end','gene_name','score','strand']].drop_duplicates()
                feature_bed['gene_name'] = feature_bed['gene_name'] + "_" + feature_type
        elif strand == "-" and feature_type == 'five_prime_utr':
            exon_dat = feature_gtf[(feature_gtf['end'] == feature_gtf['end'].max()) & (feature_gtf['feature'] == 'exon')]
            exon_dat['length'] = exon_dat['end'] - exon_dat['start']
            if any(exon_dat['length'] > 25):
                exon_dat['start'] = exon_dat['end'] - 25
                feature_bed = exon_dat[['seqname','start','end','gene_name','score','strand']].drop_duplicates()
                feature_bed['gene_name'] = feature_bed['gene_name'] + "_" + feature_type
            else:
                feature_bed = exon_dat[['seqname','start','end','gene_name','score','strand']].drop_duplicates()
                feature_bed['gene_name'] = feature_bed['gene_name'] + "_" + feature_type
    else:
        feature_bed = None

    return feature_bed
    
def main(command=None):
    """ Run specified program(s). """
    if command is None:
        my_command = CommandLine()  # read options from the command line
    else:
        my_command = CommandLine(command)  # interpret the list passed from the caller of main
    
    species, gtfFile, orf_primer_file = my_command.args.s, my_command.args.gtf, my_command.args.o
    
    ## read in GTF file 
    gtf = pd.read_csv(gtfFile, sep='\t', header = None, comment = "#")
    gtf.columns = ['seqname','source','feature','start','end','score','strand','frame','attribute']
    
   
    ## read in feature file  
    orfPrimers = pd.read_csv(orf_primer_file, header = 0)
    orfPrimers['gene'] = orfPrimers['gene_symbol'].apply(lambda x: re.sub('_.*','',x))
    orfPrimers['feature_type'] = orfPrimers['gene_symbol'].apply(lambda x: re.sub(r'^.*?_','',x))
    orfPrimers['feature_type'] = orfPrimers['feature_type'].apply(lambda x: re.sub(r'\(.*','',x))
    orfPrimers['strand'] = orfPrimers['gene_symbol'].apply(lambda x: re.sub(r'.*\(','',x))
    orfPrimers['strand'] = orfPrimers['strand'].apply(lambda x: re.sub(r'\)','',x))

    ## what genes do we not have a five primer and three primer primer for 
    dat = orfPrimers.groupby('gene')['feature_type'].nunique().reset_index()
    dat.columns = ['gene','count']
    screen_genes = dat[dat['count'] < 2]['gene'].to_list()
    
    screenDat = orfPrimers[orfPrimers['gene'].isin(screen_genes)][['gene','feature_type','strand']]
    screenDat['missing_feature'] = np.where(screenDat['feature_type'] == "five_prime_utr", 'three_prime_utr', 'five_prime_utr')
    screenDat = screenDat.drop_duplicates()


    # [(CreateBedData, gtf, row) for _, row in screenDat.iterrows()]
    # first_row_index, first_row = next(screenDat.iterrows())
    # CreateBedData(gtf, first_row)

    bedDat = pd.DataFrame()
    with concurrent.futures.ThreadPoolExecutor(max_workers=25) as executor:
        futures = [executor.submit(CreateBedData, gtf, row) for _, row in screenDat.iterrows()]
        for future in concurrent.futures.as_completed(futures):
            bedDat = pd.concat([bedDat, future.result()])

    bedDat['new_strand'] = np.where((bedDat['gene_name'].str.contains('three_prime_utr')) & (bedDat['strand'] == "+"), "-", bedDat['strand'])
    bedDat['new_strand_2'] = np.where((bedDat['gene_name'].str.contains('three_prime_utr')) & (bedDat['strand'] == "-"), "+", bedDat['new_strand'])
    bedDat = bedDat[["seqname", "start", "end", "gene_name", "score", "new_strand_2"]]
    bedDat['start'] = bedDat['start'] - 1
    bedfile = os.path.join("/ds-workspace/EW-TempDataStore/ORF", species, "ExonScan.bed")
    bedDat.to_csv(bedfile, sep='\t', header = None, index = False)

    ## get fasta
    exon_scan_fastafile = os.path.join("/ds-workspace/EW-TempDataStore/ORF", species, "ExonScan.fa")
    cmd = ["/usr/bin/bedtools", "getfasta", "-s", "-name", "-fi", my_command.args.fa, "-bed", bedfile, "-fo", exon_scan_fastafile]
    subprocess.run(cmd)

    #####################################
    ## Design new primers
    fastaDat = pd.read_csv(exon_scan_fastafile, sep=',', header = None)

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

    fastaDat['Header'] = fastaDat['Header'].apply(lambda x: re.sub('^>','',x)) 
    fastaDat['gene'] = fastaDat['Header'].apply(lambda x: re.sub('_.*','',x))
    fastaDat['feature_type'] = fastaDat['Header'].apply(lambda x: re.sub(r'^.*?_','',x))
    fastaDat['feature_type'] = fastaDat['feature_type'].apply(lambda x: re.sub(r'\(.*','',x))
    fastaDat['strand'] = fastaDat['Header'].apply(lambda x: re.sub(r'.*\(','',x))
    fastaDat['strand'] = fastaDat['strand'].apply(lambda x: re.sub(r'\)','',x))
    fastaDat = fastaDat.rename(columns = {'Sequence':'SEQUENCE'})
    primer_dat = pd.concat([fastaDat[['gene','feature_type','strand','SEQUENCE']], orfPrimers[['gene','feature_type','strand','SEQUENCE']]])

    finalDat = primer_dat.drop('strand', axis=1).groupby(['gene', 'feature_type']).sample(n=1).pivot_table(values='SEQUENCE', index=['gene'], columns=['feature_type'], aggfunc='first').reset_index()
    primer_file = os.path.join("/ds-workspace/EW-TempDataStore/ORF", species, "FinalPrimerSet.csv")
    finalDat.to_csv(primer_file, sep=',', header = False, index = False)

    """
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
    outfile = os.path.join("/ds-workspace/EW-TempDataStore/ORF", species, "ORF_primers_exonScan.csv")
    finalDat.to_csv(outfile, sep=',', header = True, index=False)
    """

if __name__ == "__main__":
    main()





