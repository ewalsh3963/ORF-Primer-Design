#!/usr/bin/env nextflow

include { GetFasta } from './modules/GetFasta.nf'
include { GetFasta as GetFastaRedesign } from './modules/GetFasta.nf'
include { PrimerDesign } from './modules/GetFasta.nf'
include { PrimerDesign as PrimerRedesign } from './modules/GetFasta.nf'

process DownloadRefGenome {
    input:
        val species
        val outdir
    
    output:
        path "*.dna.primary_assembly.fa" // Specify the expected output files    

    publishDir params.outDir, mode: 'copy'
    script:
    
    """    

    /home/ewalsh/PrimerDesign/DownloadRefGenome.sh ${outdir} ${species}
    
    """
}

process DownloadRefGTF {
    input:
        val species 
        val outdir
        
    output:
        path "*.gtf.gz" // Specify the expected output files

    publishDir params.outDir, mode: 'copy'
    script:
    
    """    
    
    /home/ewalsh/PrimerDesign/DownloadRefGTF.sh $outdir $species

    """
}

process CreateBedFile {
    input:
        path gtffile // Expecting a file input from the previous process
        val features
        val species

    output:
    	path("features.bed")
    publishDir params.outDir, mode: 'copy'

    script:
    
    """

    echo "Processing file: ${gtffile}"
    python3 /home/ewalsh/PrimerDesign/FilterGTF.py -s ${species} -gtf ${gtffile} -f ${features} > features.bed

    """
}

process CopyBedFile {
    input: 
        val bedFile
    
    output: 
        path("features.bed")

    publishDir params.outDir, mode: 'copy'
    
    script:

    """
        cp ${bedFile} features_oi.bed
    """
}

process BedEdit {
    input:
    path gtfFiles 
    path orf_fasta
    path orf_primers
    path bedFile

    output:
    	path("features_redesign.bed")
    // publishDir params.outDir, mode: 'copy'

    script:
    
    """
    echo $gtfFiles

    python3 /home/ewalsh/PrimerDesign/primerReDesign.py -gtf ${gtfFiles} -fa ${orf_fasta} -bed ${bedFile} -pd ${orf_primers} > features_redesign.bed

    """

}

/*
 * Define the workflow
 */
workflow {
    if (params.ref_fasta == null) {
        println "Downloading genome reference fasta file..."
        fastaFile = DownloadRefGenome(params.species, params.outDir)
    } else {
        fastaFile = file(params.ref_fasta)
    }

    if (params.bed != null) {
        CopyBedFile(params.bed)
    } else if (params.gtf == null) {
        println "Downloading genome reference GTF file..."

        // Call the DownloadRefGTF process and capture its output
        gtfFiles = DownloadRefGTF(params.species, params.outDir)

        // View the files captured in gtfFiles
        gtfFiles.view { file -> "Found file: ${file}" }
        gtfFiles.count().view { count -> "Number of elements in gtfFiles: ${count}" }

        println "+++++++++++++++++++++++++++++++++++++++++"
        println "Creating bed file from ${params.features}"
        
        // Pass the output of DownloadRefGTF directly to CreateBedFile
        bedFile = CreateBedFile(gtfFiles, params.features, params.species)
    }

    orf_fasta = GetFasta(fastaFile, bedFile, "PrimaryDesign")
    orf_primers = PrimerDesign(orf_fasta, bedFile, "PrimaryDesign")

    // Change the search region for UTR with no hits in the first pass 
    bed_redesign = BedEdit(gtfFiles, orf_fasta, orf_primers, bedFile)
    orf_fasta_redesign = GetFastaRedesign(fastaFile, bed_redesign, "Redesign")
    orf_primers = PrimerRedesign(orf_fasta_redesign, bed_redesign, "Redesign")

    // Merge design and redesign files 
    // Channel
    // .fromPath("${params.outDir}/ORF*.fa")
    // .collectFile(name: "${params.outDir}/ORF.fa")

    Channel
    .fromPath("${params.outDir}/ORF_primers*.csv")
    .collectFile(name: "${params.outDir}/ORF_primers.csv")

}


// nextflow run main.nf  -process.echo --species Mus_musculus \
// --outDir /home/ewalsh/scratch \
// --features /home/ewalsh/PrimerDesign/GeneSymbols.csv