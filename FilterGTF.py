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
        self.parser.add_argument('-f', metavar='--featureFile', action='store', help='feature file', required=False, type=str)
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

def CreateBedData(gtf, feature):
    feature_gtf = gtf[gtf['attribute'].str.contains('"' + feature + '"', case=True, na=False)]
    if feature_gtf.shape[0] > 0:
        feature_gtf['transcript_biotype'] = feature_gtf['attribute'].apply(lambda x: (match.group(1) if (match := re.search(r'transcript_biotype "([^"]+)"', x)) else None))
        feature_gtf['gene_name'] = feature_gtf['attribute'].apply(lambda x: (match.group(1) if (match := re.search(r'gene_name "([^"]+)"', x)) else None))
        feature_gtf['gene_id'] = feature_gtf['attribute'].apply(lambda x: (match.group(1) if (match := re.search(r'gene_id "([^"]+)"', x)) else None))
        feature_gtf['transcript_id'] = feature_gtf['attribute'].apply(lambda x: (match.group(1) if (match := re.search(r'transcript_id "([^"]+)"', x)) else None))
        ## Find the longest coding sequence
        # cds_gtf = feature_gtf[feature_gtf['feature'] == 'CDS']
        cds_gtf = feature_gtf[feature_gtf['feature'] == 'transcript']
        cds_gtf = cds_gtf[cds_gtf['transcript_biotype'] == 'protein_coding']
        if cds_gtf.shape[0] > 0:
            cds_gtf['length'] = abs(cds_gtf['end'] - cds_gtf['start'])
            ## Get the 5' and 3' utr features
            transcript_to_screen = cds_gtf[cds_gtf['length'] == max(cds_gtf['length'])]['transcript_id'].to_list()
            transcript_gtf = feature_gtf[feature_gtf['transcript_id'].isin(transcript_to_screen)]
            transcript_gtf = transcript_gtf[transcript_gtf['feature'].isin(["five_prime_utr","three_prime_utr"])]
            transcript_gtf['length'] = transcript_gtf['end'] - transcript_gtf['start']
            idx = transcript_gtf.groupby('feature')['length'].transform(max) == transcript_gtf['length']    
            transcript_gtf = transcript_gtf[idx]
            transcript_gtf['gene_name'] = transcript_gtf['gene_name'] + "_" + transcript_gtf['feature']
            for index, row in transcript_gtf.iterrows():
                # Modify row values based on conditions or calculations
                if row['strand'] == "+" and row['length'] > 100 and row['feature'] == 'five_prime_utr':
                    diff = row['length'] - 100
                    transcript_gtf.at[index, 'start'] = row['start'] + diff 
                elif row['strand'] == "+" and row['length'] > 100 and row['feature'] == 'three_prime_utr':
                    diff = row['length'] - 100
                    transcript_gtf.at[index, 'end'] = row['end'] - diff 
                elif row['strand'] == "-" and row['length'] > 100 and row['feature'] == 'five_prime_utr':
                    diff = row['length'] - 100
                    transcript_gtf.at[index, 'end'] = row['end'] - diff 
                elif row['strand'] == "-" and row['length'] > 100 and row['feature'] == 'three_prime_utr':
                    diff = row['length'] - 100
                    transcript_gtf.at[index, 'start'] = row['start'] + diff 

            feature_bed = transcript_gtf[['seqname','start','end','gene_name','score','strand']].drop_duplicates()
            # bedDat = pd.concat([bedDat, feature_bed])
        else:
            feature_bed = None
    else:
        feature_bed = None
    return feature_bed
    
def main(command=None):
    """ Run specified program(s). """
    if command is None:
        my_command = CommandLine()  # read options from the command line
    else:
        my_command = CommandLine(command)  # interpret the list passed from the caller of main
    
    species, gtfFile, featureFile = my_command.args.s, my_command.args.gtf, my_command.args.f
    
    ## read in GTF file 
    gtf = pd.read_csv(gtfFile, sep='\t', header = None, comment = "#")
    gtf.columns = ['seqname','source','feature','start','end','score','strand','frame','attribute']
    
    ## read in feature file  
    features = pd.read_csv(featureFile, header = None)
    features.columns = ['feature']
    features = features.drop_duplicates().reset_index().drop('index', axis = 1)
    features = symbol_to_ensembl('Homo_sapiens', features, 'feature')
    ## define the mart/serve to use for querying
    server = Server(host='http://www.ensembl.org', use_cache=False) #server.list_marts()
    mart = server['ENSEMBL_MART_ENSEMBL'] # mart.list_datasets()
    dataset = mart['hsapiens_gene_ensembl'] 

    if species != "Homo_sapiens":
        queryList = features['ensembl_id'].tolist()
        queryList = [i for i in queryList if i is not None] ## remove values that are None

        ## split the query list into chunks of size querySize
        querySize = 100; queryDF = pd.DataFrame()
        queryList = [queryList[i:i + querySize] for i in range(0, len(queryList), querySize)]
        queryFinal = []
        for i in range(0, len(queryList), querySize):
            queryFinal.append(queryList[i:i+querySize])
        queryFinal = queryFinal[0]
        
        for item in queryFinal:
            ## run biomart query 
            ortho_species_biomart = (species[0] + re.sub('.*_','',species)).lower()
            att1 = ortho_species_biomart + '_homolog_ensembl_gene'   
            tmp = dataset.query(attributes=['ensembl_gene_id', att1], filters={'link_ensembl_gene_id': item})
            try: 
                queryDF = pd.concat([queryDF, tmp])
            except ValueError:
                pass
        
        queryDF.columns = ['ensembl_id','ortho_symbol']
        ## merge with symbols 
        features = pd.merge(features, queryDF, how = 'left')
        
        ## get gene symbols for species of interest
        queryList = features['ortho_symbol'].tolist()
        queryList = [i for i in queryList if i is not None and i is not np.nan]  ## remove values that are None
        dataset = mart[ortho_species_biomart + '_gene_ensembl'] 

        ## split the query list into chunks of size querySize
        querySize = 100; queryDF = pd.DataFrame()
        queryList = [queryList[i:i + querySize] for i in range(0, len(queryList), querySize)]
        queryFinal = []
        for i in range(0, len(queryList), querySize):
            queryFinal.append(queryList[i:i+querySize])
        queryFinal = queryFinal[0]
        for item in queryFinal:
            ## run biomart query 
            att1 = "external_gene_name"
            tmp = dataset.query(attributes=['ensembl_gene_id', att1], filters={'link_ensembl_gene_id': item})
            try: 
                queryDF = pd.concat([queryDF, tmp])
            except ValueError:
                pass
        
        queryDF.columns = ['ortho_symbol','ortho_feature']
        features = pd.merge(features, queryDF, how = 'left')

        if species == "Mus_musculus":
            features['ortho_feature'] = np.where(features['ortho_feature'].isna(), features['feature'].apply(lambda x: x[0] + x[1::].lower()), features['ortho_feature'])
        else:
            features['ortho_feature'] = np.where(features['ortho_feature'].isna(), features['feature'], features['ortho_feature'])

        # if species == "Mus_musculus":
        #     features['feature'].apply(lambda x: x[0] + x[1::].lower())
        features_list = features['ortho_feature'].to_list()
    else: 
        features_list = features['ensembl_id'].to_list()
    
    pdb.set_trace()

    # gtf = gtf[gtf['feature'].isin(["start_codon","stop_codon"])]
    bedDat = pd.DataFrame()
    features_list = list(set(features_list))
    features_list = [x for x in features_list if x != None]
    with concurrent.futures.ThreadPoolExecutor(max_workers=25) as executor:
        futures = [executor.submit(CreateBedData, gtf, feature) for feature in features_list]
        for future in concurrent.futures.as_completed(futures):
            bedDat = pd.concat([bedDat, future.result()])

    bedDat = bedDat.sort_values('gene_name') 
    spamwriter = csv.writer(sys.stdout, delimiter='\t')
    spamwriter.writerows(bedDat.values)

if __name__ == "__main__":
    main()



# python3 /home/ewalsh/PrimerDesign/FilterGTF.py \
# -gtf /ds-workspace/EW-TempDataStore/ORF/Mus_musculus/Mus_musculus.GRCm39.113.gtf.gz \
# -f /ds-workspace/EW-TempDataStore/Retrogenix_data/protein_screens/CrossSpeciesORF.csv \
# -s Mus_musculus


# python3 /home/ewalsh/PrimerDesign/FilterGTF.py \
# -gtf /ds-workspace/EW-TempDataStore/ORF/Homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz \
# -f /ds-workspace/EW-TempDataStore/Retrogenix_data/protein_screens/CrossSpeciesORF.csv \
# -s Homo_sapiens


