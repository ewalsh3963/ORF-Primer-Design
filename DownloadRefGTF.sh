#!/bin/bash

## Inputs from command line
outdir=$1
species=$2
species_lower=${species,,}

## GTF file download
ftp_root="https://ftp.ensembl.org/pub"
ftp_path="${ftp_root}/current_gtf/${species_lower}/"
cmd="wget -r -nd -np -A *.gtf.gz -R index.html*,*chr.gtf.gz,*abinitio.gtf.gz,*scaff.gtf.gz -e robots=off $ftp_path"
eval $cmd

# Check if wget was successful
if [ $? -ne 0 ]; then
    echo "Error downloading files for species: $species_lower" >&2

    exit 1
fi


