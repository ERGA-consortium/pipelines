process EXTRACT_TRANSCRIPTOME {
    label 'process_single'
    label 'process_medium_memory'
    label 'biopython'

    input:
    path(ch_gff)
    path(ch_fasta)
     
    output:
    path("predicted_transcriptome.fa"), emit: transcriptome

    script:
    """
    python3 ${projectDir}/bin/gff2fasta_biopython.py --gff ${ch_gff} --seq ${ch_fasta} --concat --out predicted_transcriptome.fa
    """
}
