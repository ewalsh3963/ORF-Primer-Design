process GetFasta {
    input:
    path fasta
    path bedFile // Expecting a file input from the previous process
    val suffix
    output:
        path("ORF_${suffix}.fa")
    
    publishDir params.outDir, mode: 'copy'

    script:
    
    """
    echo ${fasta}
    echo ${bedFile}
    
    /usr/bin/bedtools getfasta -s -name -fi ${fasta} -bed ${bedFile} -fo ORF_${suffix}.fa

    """

}

process PrimerDesign {
    input:
    path orf_fasta
    path bedFile
    val suffix

    output:
    path("ORF_primers_${suffix}.csv")

    publishDir params.outDir, mode: 'copy'

    script:
    """
    if [ "${suffix}" == "PrimaryDesign" ]; then
        python3 /home/ewalsh/PrimerDesign/primerDesign.py -title -fa "${orf_fasta}" -bed "${bedFile}" > ORF_primers_${suffix}.csv
    else
        python3 /home/ewalsh/PrimerDesign/primerDesign.py -fa "${orf_fasta}" -bed "${bedFile}" > ORF_primers_${suffix}.csv
    fi
    """
}


// /usr/bin/bedtools getfasta -s -name -fi /home/ewalsh/PrimerDesign/work/dd/57460bad61dcbd2fd8fc00e7443cbd/Mus_musculus.GRCm39.dna.primary_assembly.fa -bed /home/ewalsh/PrimerDesign/work/a2/905b56b667859d31d83897376738d7/features_redesign.bed

// grep "Ccer2_five_prime_utr" /home/ewalsh/PrimerDesign/work/a2/905b56b667859d31d83897376738d7/features_redesign.bed


// Ccl27al_three_prime_utr