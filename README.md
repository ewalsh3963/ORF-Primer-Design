# ORF-Primer-Design

This repository contains tools and pipelines for automated primer design targeting open reading frames (ORFs) based on input coding sequences (CDS) and genome reference files. It is intended for use in high-throughput cloning and expression studies.

## 🎯 Purpose

The project streamlines the process of:

- Parsing CDS and genome files
- Designing primers flanking ORFs
- Exporting primer sequences for synthesis or downstream molecular biology workflows

## ⚙️ Features

- Accepts FASTA/GenBank-formatted genome and CDS files
- Automatically identifies ORF regions and generates compatible primer pairs
- Accepts Homo sapiens gene symbols and can convert them to other species to support multi-species expression assays
- Batch processing of multiple genes
- Outputs primer sequences and metadata in CSV or tab-delimited formats

## 🧰 Technologies Used

- Python 3.x
- `Biopython` for sequence handling
- `pandas` for data manipulation
- `argparse` for command-line interfaces

## 🚀 Getting Started

1. **Clone the repository:**

   ```bash
   git clone https://github.com/ewalsh3963/ORF-Primer-Design.git
   cd ORF-Primer-Design

2. **Install dependencies:**
```
pip install -r requirements.txt
```

3. **Run the primer design script:**

Example usage:
```
nextflow run main.nf -bg -process.echo --species Mus_musculus \
 --outDir /path/to/outdir \
--features /path/to/geneSymbols/inCSV/format
```
