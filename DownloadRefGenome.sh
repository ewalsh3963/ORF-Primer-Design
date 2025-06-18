#!/bin/bash

## Inputs from command line
outdir=$1
species=$2
species_lower=${species,,}

## Fasta file download
ftp_root="https://ftp.ensembl.org/pub"
ftp_path="${ftp_root}/current_fasta/${species_lower}/dna/"

if [ "$species" == "Macaca_fascicularis" ]; then
    # Commands to execute if the variable equals the string
    cmd="wget -r -nd -np -A *.dna.primary_assembly.*.fa.gz -R index.html*,*_from_genomic.fna.gz -e robots=off $ftp_path"
    eval $cmd
else
    # Commands to execute if the variable does NOT equal the string
    cmd="wget -r -nd -np -A *.dna.primary_assembly.fa.gz -R index.html*,*_from_genomic.fna.gz -e robots=off $ftp_path"
    eval $cmd
fi

# Check if wget was successful
if [ $? -ne 0 ]; then
    echo "Error downloading files for species: $species_lower" >&2
    exit 1
fi


fasta_files=($(find . -type f -name "*.fa.gz"))
if [ ${#fasta_files[@]} -gt 1 ]; then
    gunzip -c *.fa.gz | gzip > "${species}.dna.primary_assembly.fa.gz"

    rm *.primary_assembly.*.fa.gz
fi

## unzip the fasta file
gunzip *.fa.gz