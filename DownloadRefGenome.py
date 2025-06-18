import os
import sys
import numpy as np
import pandas as pd
import OS_Tools
import argparse


class CommandLine:
    """ Handle the command line, usage and help requests """
    def __init__(self, in_opts=None):
        """ CommandLine constructor: Implements a parser to interpret the command line argv string using argparse. """
        self.parser = argparse.ArgumentParser(description="%(prog)s use Random Forest to classify cells as tranduced or not", 
        add_help=True, 
        prefix_chars='-', 
        usage='python3 %(prog)s [-h] [-v] ...')
     
        ## inputs & outputs
        self.parser.add_argument('-s', metavar='--species', action='store', help='study ID', required=False, type=str, default=None)
        self.parser.add_argument('-o', metavar='--outroot', action='store', help='output root directory', required=False, type=str)

        if in_opts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(in_opts)


def NCBI_get_ref_files(outdir, species, ftp_root="https://ftp.ncbi.nlm.nih.gov"):
    ## get the path to the species reference files from NCBI ftp site
    targetDir = os.path.join(outdir, species)
    OS_Tools.clear_directory(targetDir);
    OS_Tools.ensure_directory(targetDir, critical=False);

    ftp_path = os.path.join(ftp_root, 'genomes', 'refseq', 'vertebrate_mammalian', species.lower(), 'latest_assembly_versions/')
    # command = ["wget -r -nd -q -np -A *_genomic.fna.gz -R index.html*,*_from_genomic.fna.gz -e robots=off -P", os.path.join(targetDir, 'fasta'), ftp_path]
    command = ["wget -r -nd -q -np -A *_genomic.fna.gz -R index.html*,*_from_genomic.fna.gz -e robots=off -P", targetDir, ftp_path]
    subprocess.run(command)

    ftp_path = os.path.join(ftp_root, 'genomes', 'refseq', 'vertebrate_mammalian', species.lower(), 'latest_assembly_versions/')
    # command = ["wget -r -nd -q -np -A *_genomic.gtf.gz -R index.html* -e robots=off -P", os.path.join(targetDir, 'gtf'), ftp_path]
    command = ["wget -r -nd -q -np -A *_genomic.gtf.gz -R index.html* -e robots=off -P", targetDir, ftp_path]
    subprocess.run(command)

def ensembl_get_ref_files(outdir, species, ftp_root="https://ftp.ncbi.nlm.nih.gov"):
    ## get the path to the species reference files from NCBI ftp site
    targetDir = os.path.join(outdir, species)
    OS_Tools.clear_directory(targetDir)

    ftp_path = os.path.join(ftp_root, 'current_fasta', species.lower() + 'dna/')            
    command = ["wget -r -nd -q -np -A *.dna.primary_assembly.fa.gz -R index.html*,*_from_genomic.fna.gz -e robots=off -P", os.path.join(targetDir, 'fasta'), ftp_path]
    subprocess.run(command)
    
    ftp_path = os.path.join(ftp_root, 'current_gtf', species.lower() + '/')
    command = ["wget -r -nd -q -np -A *.gtf.gz -R index.html*,*chr.gtf.gz,*abinitio.gtf.gz,*scaff.gtf.gz  -e robots=off -P", os.path.join(targetDir, 'gtf'), ftp_path]
    subprocess.run(command)

def main(command=None):
    """ Run specified program(s). """
    if command is None:
        my_command = CommandLine()  # read options from the command line
    else:
        my_command = CommandLine(command)  # interpret the list passed from the caller of main
    
    ########################################################################
    ## Load the adata files and integrate them
    species, outdir = my_command.args.s, my_command.args.o
    NCBI_get_ref_files(outdir, species)


if __name__ == "__main__":
    main()