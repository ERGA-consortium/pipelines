process EXTRACT_PROTEOME {
    label 'process_single'
    label 'biopython'

    input:
    path(ch_fasta)
    path(ch_gff)
    
    output:
    path("predicted_proteome.fa"), emit: proteome

    script:
    """
    python3 ${projectDir}/bin/gff2fasta_biopython.py --gff ${ch_gff} --seq ${ch_fasta} --concat --out predicted_proteome.fa --tr --gc ${params.genetic_code}
    """
}
