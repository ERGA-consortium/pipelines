rule RepeatModeler:
    """Create a Database for RepeatModeler and execute RepeatModeler"""

    input:
        asm = config['asm']
    output:
        os.path.join(config['snakemake_dir_path'],"results/1_MaskRepeat/RepeatModeler/repeatMM_summary.txt"),
        os.path.join(config['snakemake_dir_path'], 'results/1_MaskRepeat/RepeatModeler', config['name']+"-families.fa")
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/1_MaskRepeat/RepeatModeler/RepeatModeler.log')
    singularity:
        '/srv/public/users/brown/pipelines/annotation/tetools_latest.sif'
    threads: 40
    params:
        out_db = directory(os.path.join(config['snakemake_dir_path'], 'results/1_MaskRepeat/RepeatModeler')),
        name = config['name']
    shell:
        """
        cd {params.out_db}
        (BuildDatabase -name {params.name} {input.asm}) 2> {log}
        (RepeatModeler -database {params.name} -threads {threads} -LTRStruct) 2>> {log}
        cat > repeatMM_summary.txt <<EOF
        'total_nb_of_consensus_repeats:' $(grep -c '>' $name-families.fa)
        'nbr_repeats_per_cat:' $(grep '>' $name-families.fa | sed -r 's/.+#//' | sed -r 's/\s+.+//' | sort | uniq -c)
        EOF
        """


PROT_DIR=os.path.dirname(os.path.join(config['snakemake_dir_path'],"files/uniprot_sprot.fasta.gz")) + "/"
PROT_FASTA_GZ=os.path.basename(os.path.join(config['snakemake_dir_path'],"files/uniprot_sprot.fasta.gz"))
PROT_FASTA, PROT_EXT_GZ=os.path.splitext(PROT_FASTA_GZ)
PROT_NAME, PROT_EXT=os.path.splitext(PROT_FASTA)


rule split_uniprot:
    """Split the UniProt/Swissprot protein database into chunks for transposonPSI"""
    input: 
        prot = os.path.join(config['snakemake_dir_path'],"files/uniprot_sprot.fasta.gz")# curated proteins from swissprot/uniprot
    output:
        chunks = temp(expand(PROT_DIR + "split_result/" + PROT_NAME + "_chunk{nr}.fa", nr=range(1, 101)))
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/1_MaskRepeat/split_uniprot/split_uniprot.log')
    params:
        dir = PROT_DIR,
        prot_f = PROT_FASTA,
        prot = PROT_NAME,
        prot_path = PROT_DIR + PROT_NAME
    conda: "../envs/gaas.yaml"
    shell:
        """
        cd {params.dir}
        gunzip -k {input.prot}
        gaas_fasta_splitter.pl -f {params.prot_f} --nb_chunks 100 -o tmp/
        mv tmp/*fa split_result/
        rm -r tmp/
        """

rule transposonPSI:
    """Identify transposons in the UniProt/Swissprot protein dataset"""
    input:
        chunk = PROT_DIR + "split_result/" + PROT_NAME + "_chunk{nr}.fa"
    output:
        allHits = temp(PROT_DIR + "split_result/" + PROT_NAME + "_chunk{nr}.fa.TPSI.allHits"),
        topHits = temp(PROT_DIR + "split_result/" + PROT_NAME + "_chunk{nr}.fa.TPSI.topHits")
    params:
        dir = PROT_DIR + "split_result/",
    conda: "../envs/tePSI.yaml"
    shell:
        """
        cd {params.dir}
        transposonPSI.pl {input.chunk} prot
        """


rule list_tePSI_hits:
    input:
        topHits = expand(PROT_DIR + "split_result/" + PROT_NAME + "_chunk{nr}.fa.TPSI.topHits", nr=range(1, 101))
    output:
        allTopHits = PROT_DIR + PROT_NAME + ".TPSI.topHits",
        prot_list = PROT_DIR + PROT_NAME + ".TPSI.topHits.accessions.txt"
    shell:
        """
        cat {input.topHits} > {output.allTopHits} &&
        awk '{{if($0 ~ /^[^\/\/.*]/) print $5}}' {output.allTopHits} | sort -u > {output.prot_list}
        """


rule filter_uniprot_fasta:
    """Remove transposons from the UniProt/Swissprot protein dataset"""
    input:
        prot = os.path.join(config['snakemake_dir_path'],"files/uniprot_sprot.fasta.gz"),
        prot_list = PROT_DIR + PROT_NAME + ".TPSI.topHits.accessions.txt"
    output:
        prot_filtered = PROT_DIR + PROT_NAME + ".noTEs.fa"
    log: os.path.join(config['snakemake_dir_path'], 'logs/1_MaskRepeat/filter_uniprot_fasta/filter_uniprot_fasta.log')
    params:
        dir = PROT_DIR,
        prot_path = PROT_DIR + PROT_FASTA
    conda: "../envs/gaas.yaml"
    shell:
        """
        echo 'etape 2.0' > {log}
        cd {params.dir}
        echo 'etape 2.1' >> {log}
        gaas_fasta_removeSeqFromIDlist.pl -f {params.prot_path} -l {input.prot_list} -o {output.prot_filtered}
        """


rule filtered_blast_db:
    """Generate BLAST database from filtered UniProt/Swissprot protein dataset"""
    input: 
        prot_filtered = PROT_DIR + PROT_NAME + ".noTEs.fa"
    output:
        phr = PROT_DIR + PROT_NAME + ".noTEs.fa.phr",
        pin = PROT_DIR + PROT_NAME + ".noTEs.fa.pin",
        psq = PROT_DIR + PROT_NAME + ".noTEs.fa.psq"
    params:
        dir = PROT_DIR
    conda: "../envs/blast.yaml"
    shell:
        """
        cd {params.dir}
        makeblastdb -in {input.prot_filtered} -dbtype prot
        """

rule symbolic_links:
    """Create symbolic links to repeat libraries in output directory"""
    input:
        repmo_raw = os.path.join(config['snakemake_dir_path'], 'results/1_MaskRepeat/RepeatModeler', config['name']+"-families.fa")
    output:
        sym_link = os.path.join(config['snakemake_dir_path'], 'results/1_MaskRepeat/filter_TEprot', config['name']+"-families.fa", config['name']+"-families.fa")
    shell:
        """
        ln -s {input.repmo_raw} {output.sym_link}
        """


rule blast_repeat_library:
    """Blastx repeat library to filtered Uniprot/Swissprot database"""
    input:
        repmo_raw = os.path.join(config['snakemake_dir_path'], 'results/1_MaskRepeat/filter_TEprot', config['name']+"-families.fa", config['name']+"-families.fa"), 
        blast_db_idx = rules.filtered_blast_db.output,
        blast_db = PROT_DIR + PROT_NAME + ".noTEs.fa"
    output:
        blast = os.path.join(config['snakemake_dir_path'], 'results/1_MaskRepeat/filter_TEprot', config['name']+"-families.fa", config['name']+"-families.fa.blastx.out")
    params:
        dir = os.path.join(config['snakemake_dir_path'], 'results/1_MaskRepeat/filter_TEprot', config['name']+"-families.fa")
    threads: 8
    conda: "../envs/blast.yaml"
    shell:
        """
        cd {params.dir}
        blastx -num_threads {threads} -db {input.blast_db} -query {input.repmo_raw} -out {output.blast}
        """

rule protexcluder:
    """Remove blast hits from repeat library"""
    input:
        repmo_raw = os.path.join(config['snakemake_dir_path'], 'results/1_MaskRepeat/filter_TEprot', config['name']+"-families.fa", config['name']+"-families.fa"),
        blast = os.path.join(config['snakemake_dir_path'], 'results/1_MaskRepeat/filter_TEprot', config['name']+"-families.fa", config['name']+"-families.fa.blastx.out")
    output:
        repmo_fil = os.path.join(config['snakemake_dir_path'], 'results/1_MaskRepeat/filter_TEprot', config['name']+"-families.fa", config['name']+"-families.fanoProtFinal") 
    params:
        dir = os.path.join(config['snakemake_dir_path'], 'results/1_MaskRepeat/filter_TEprot', config['name']+"-families.fa") 
    conda: "../envs/protexcluder.yaml"
    shell:
        """
        cd {params.dir}
        ProtExcluder.pl {input.blast} {input.repmo_raw}
        """

rule mask:
    input:
        repmo_fil = os.path.join(config['snakemake_dir_path'], 'results/1_MaskRepeat/filter_TEprot', config['name']+"-families.fa", config['name']+"-families.fanoProtFinal"),
        asm = config['asm']
    output:
        out_mask = directory(os.path.join(config['snakemake_dir_path'], 'results/1_MaskRepeat/RepeatMasker')),
        asm_masked = os.path.join(config['snakemake_dir_path'], 'results/1_MaskRepeat/RepeatMasker/asm_decontaminated.fasta.masked')
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/1_MaskRepeat/RepeatMasker/RepeatMasker.log')
    singularity:
        '/srv/public/users/brown/pipelines/annotation/tetools_latest.sif'
    threads: 40
    shell:
        """
        (RepeatMasker -pa {threads} -gff -xsmall -dir {output.out_mask} -lib {input.repmo_fil} {input.asm}) 2> {log}
        """



