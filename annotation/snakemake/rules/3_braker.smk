

rule braker:
    input:
        asm_masked = os.path.join(config['snakemake_dir_path'], 'results/1_MaskRepeat/RepeatMasker/', os.path.basename(config['asm'])) + '.masked',
        aln_bam = os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/hisat2/merge.sorted.bam")
    output:
        os.path.join(config['snakemake_dir_path'],"results/2_braker/out_braker/braker/braker.aa") 
    threads: config['max_threads']
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/out_braker/out_braker.log')
    params:
        protDB = os.path.join(config['snakemake_dir_path'], 'files/uniprot_sprot.fasta'),
        bams_dir = os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/hisat2"),
        out_dir = directory(os.path.join(config['snakemake_dir_path'],"results/2_braker/out_braker"))
    singularity:
        docker://teambraker/braker3:v.1.0.4'
    shell:
        """
        sed -i '/^>/ s/ .*//' {input.asm_masked}
        cd {params.out_dir}
        (braker.pl --genome={input.asm_masked} --prot_seq={params.protDB} --bam={input.aln_bam} --softmasking --threads {threads} --gff3) 2> {log}
        """



rule eval:
    input: 
        braker_aa = os.path.join(config['snakemake_dir_path'],"results/2_braker/out_braker/braker/braker.aa")
    output: 
        outdir = directory(os.path.join(config['snakemake_dir_path'], 'results/2_braker/braker_busco'))
    log:
        os.path.join(config['snakemake_dir_path'],'logs/2_braker/busco/busco.log')
    singularity:
        'docker://ezlabgva/busco:v5.5.0_cv1'
    threads: config['max_threads']
    params:
        out_path =  os.path.join(config['snakemake_dir_path'], 'results/2_braker'),
        out_name = 'braker_busco',
        lineage = config['busco_phylum'],
    shell:
        """
        cd {params.out_path}
        (busco -f -m prot -i {input.braker_aa} -o {params.out_name} --out_path {params.out_path} -l {params.lineage} -c {threads}) 2> {log}
        rm -r busco_downloads/
        """

