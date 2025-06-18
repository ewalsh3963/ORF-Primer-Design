import os
import sys
import numpy as np
import pandas as pd
import pdb
import argparse
import csv
import re
from pybiomart import Server
from datetime import datetime
from pyensembl import EnsemblRelease
import requests
from pyensembl import EnsemblRelease
from datetime import datetime
from requests.adapters import HTTPAdapter, Retry


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

class UniProt:
    def __init__(self):
        self.re_next_link = re.compile(r'<(.+)>; rel="next"')
        self.retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
        self.session = requests.Session()
        self.session.mount("https://", HTTPAdapter(max_retries=self.retries))

    def get_next_link(self, headers):
        if "Link" in headers:
            match = self.re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    def get_batch(self, batch_url):
        while batch_url:
            response = self.session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = self.get_next_link(response.headers)

    def query(self, queryList):
        ## split the query list into chunks of size querySize
        querySize = 100; uniprotDat = pd.DataFrame()
        queryList = [queryList[i:i + querySize] for i in range(0, len(queryList), querySize)]
        queryFinal = []
        for i in range(0, len(queryList), querySize):
            queryFinal.append(queryList[i:i+querySize])

        queryFinal = queryFinal[0]

        for item in queryFinal:
            url = 'https://rest.uniprot.org/uniprotkb/search?query={}%&fields=accession,gene_names,xref_ensembl,ft_intramem,ft_transmem,cc_subcellular_location,go_c&size=500&format=tsv'.format(' OR '.join(item))
            for batch, total in self.get_batch(url):
                lines = batch.text.splitlines()
                tmp = pd.DataFrame([x.split('\t') for x  in lines[1:]])
                if tmp.shape[0] > 0:
                    tmp.columns = ['UniProtKD_accension', 'Symbol', 'ensembl_id', 'Intramembrane', 'Transmembrane', "subcell_location","cellular_component"]
                    uniprotDat = pd.concat([uniprotDat, tmp])

        # transmembrane_genes = uniprotDat[(uniprotDat['Transmembrane'] != '')]['Symbol'].unique().tolist()
        # transmembrane_genes = [x.split(' ') for x in transmembrane_genes]
        # transmembrane_genes = [x for sublist in transmembrane_genes for x in sublist]
        return uniprotDat

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
    gtf = gtf[gtf['feature'] == 'gene']
    gtf = gtf[gtf['attribute'].str.contains('protein_coding')]
    gtf[['gene_id','gene_version','gene_name','gene_source','biotype', 'NA']] = gtf['attribute'].str.split(';', expand=True)
    gtf = gtf.drop(["gene_version","gene_source","biotype","NA"], axis = 1)
    
    ## filter the data 
    gtf['gene_id'] = gtf['gene_id'].apply(lambda x: re.sub("gene_id ", "", x))
    gtf['gene_id'] = gtf['gene_id'].apply(lambda x: re.sub("\"", "", x))
    gtf['gene_name'] = gtf['gene_name'].apply(lambda x: re.sub("gene_name ", "", x))
    gtf['gene_name'] = gtf['gene_name'].apply(lambda x: re.sub("\"", "", x))

    ## define the mart/serve to use for querying
    server = Server(host='http://www.ensembl.org', use_cache=False) #server.list_marts()
    mart = server['ENSEMBL_MART_ENSEMBL'] # mart.list_datasets()
    dataset = mart['hsapiens_gene_ensembl']
    queryList = gtf['gene_id'].tolist()
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
        tmp = dataset.query(attributes=['ensembl_gene_id', "uniprot_gn_id", "uniprot_gn_symbol"], filters={'link_ensembl_gene_id': item})
        try: 
            queryDF = pd.concat([queryDF, tmp])
        except ValueError:
            pass
    
    queryDF = queryDF.rename(columns = {'Gene stable ID':'gene_id'})
    gtf = pd.merge(gtf, queryDF, how = 'left')

    myUniProt = UniProt()
    uniprotDat = myUniProt.query(queryList = gtf.dropna()['UniProtKB Gene Name ID'].to_list()) #, species=species)

    filtered_dat = uniprotDat[(uniprotDat['Intramembrane'] != "") & (uniprotDat['Transmembrane'] != '')]
    ## save gene names
    t1 = pd.merge(gtf, filtered_dat, how = "inner", left_on = 'UniProtKB Gene Name ID', right_on = 'UniProtKD_accension')['gene_name']
    t1 = t1.apply(lambda x : re.sub(r"\s+", "", x))
    t1 = t1.drop_duplicates()
    filename = os.path.join("/ds-workspace/EW-TempDataStore/Retrogenix_data/protein_screens", datetime.now().strftime('%Y%m%d') + "_data.csv")
    t1.to_csv(filename, sep=',', header = False, index = False)


if __name__ == "__main__":
    main()
