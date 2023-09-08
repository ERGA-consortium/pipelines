
rule fastqc:
    """RNA-seq quality check using FastQC"""

    input:
        RNA_dir = config['RNA_dir']
    output:
        outdir = directory(os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/out_fastQC"))
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/align_RNA/fastqc.log')
    conda:
        '../envs/fastqc.yaml'
    threads: 40
    shell:
        """
        mkdir {output.outdir}
        cd {input.RNA_dir}
        for file in $(ls {input.RNA_dir})
        do
            SAMPLE=$(basename $file)
            fastqc -t {threads} $SAMPLE -o {output.outdir}
        done    
        """


rule trimm:
    """Trim reads using trim_galore"""

    input:
        RNA_dir = config['RNA_dir']
    output:
        outTrimm = directory(os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/out_trim_galore"))
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/align_RNA/trimm.log')
    conda:
        '../envs/trimm.yaml'
    threads: 40
    params:
        trimlog = os.path.join(config['snakemake_dir_path'], 'results/2_braker/align_RNA/trimlog.log'),
    shell:
        """
        mkdir results/2_braker/align_RNA/out_trim_galore
        for filepath in {input.RNA_dir}/*_1.fastq.gz
        do
            i=$(basename "$filepath" _1.fastq.gz)
            echo "$i"
            trim_galore --paired -j {threads} --gzip -o results/2_braker/align_RNA/out_trim_galore/ {input.RNA_dir}/$i\_1.fastq.gz {input.RNA_dir}/$i\_2.fastq.gz
        done
        """



rule fastqc_trimmed:
    """Create a Database for RepeatModeler and execute RepeatModeler"""
    input:
        outTrimm = directory(os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/out_trim_galore"))
    output:
        outdir = directory(os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/out_fastQC_trimmed"))
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/align_RNA/fastqc_trimmed.log')
    conda:
        '../envs/fastqc.yaml'
    threads: 40
    shell:
        """
        mkdir {output.outdir}
        cd {input.outTrimm}
        for file in $(ls {input.outTrimm})
        do
            SAMPLE=$(basename $file)
            fastqc -t {threads} $SAMPLE -o {output.outdir}
        done
        """


rule hisat2:
    """ (Be careful, if genome fasta has long header it will cause problems afterward, so cut header with cut -d ' ' -f1 your_file.fa > new_file.fa) """
    input:
        out_trimm = directory(os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/out_trim_galore")),
        asm = config['asm']
    output:
        aln_out = directory(os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/hisat2")),
        aln_summary = os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/hisat2/splicesite.txt")
    threads: 40
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/align_RNA/hisat2/hisat2.log')
    params:
        index_dir = os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/index"),
    conda:
        '../envs/hisat.yaml'
    shell:
        """
        mkdir {params.index_dir}
        (hisat2-build {input.asm} {params.index_dir}/genome_index) 2> {log}
        for filepath in {input.out_trimm}/*_R1_paired.fq.gz
        do
             i=$(basename "$filepath" _R1_paired.fq.gz)
             echo "$i"
             hisat2 --phred33 --new-summary --novel-splicesite-outfile {output.aln_summary} -p {threads} -x {params.index_dir}/genome_index -1 {input.out_trimm}/$i\_R1_paired.fq.gz -2 {input.out_trimm}/$i\_R2_paired.fq.gz -S {output.aln_out}/$i\_accepted_hits.sam
        done
        """



rule to_bam:
    """ sam into bam file and sort + index for Apollo """
    input:
        aln_out = directory(os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/hisat2"))
    output:
        aln_sam = os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/hisat2/test.txt")
    threads: 4
    params:
        index_dir = os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/index"),
    conda:
        '../envs/samtools.yaml'
    shell:
        """
        cd {input.aln_out}
        for file in *.sam
        do
            SAMPLE=${{file%%_*}}
            samtools view -bS -o ${{SAMPLE}}_accepted_hits.bam $file
            samtools sort -o ${{SAMPLE}}_accepted_hits.sorted.bam ${{SAMPLE}}_accepted_hits.bam
            samtools index ${{SAMPLE}}_accepted_hits.sorted.bam 
        done
        echo "ok" > {output.aln_sam}
        """



# I should do like a output function to have proper output file not sheeting as I do















