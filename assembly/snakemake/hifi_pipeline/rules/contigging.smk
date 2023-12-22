rule hifiasm:
	"""Execute Hifiasm to generate primary and alternative assemblies."""

	input:
		hifi_fastq = expand("{sample}", sample=samples['hifi_data']),
		maximum_depth = os.path.join(config['results'], prefix, "genome_profiling", "maximum_depth.txt")
	output:
		asm_p_fa = os.path.join(config['results'], prefix, "contigging/hifiasm", prefix + '.asm.p_ctg.gfa'),
		asm_a_fa = os.path.join(config['results'], prefix, "contigging/hifiasm", prefix + '.asm.a_ctg.gfa')
	log:
		os.path.join(config['results'], "logs", prefix, "hifiasm.log")
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/hifiasm.yml')
	threads:
		resource['hifiasm']['threads']
	resources:
		mem_mb=resource['hifiasm']['mem_mb'],
		time=resource['hifiasm']['time']
	params:
		hifiasm_prefix = os.path.join(config['results'], prefix, "contigging/hifiasm", prefix + '.asm'),
		purgel=purgel
	shell:
		"""
		MAX_DEPTH=$(cat {input.maximum_depth})
		(hifiasm -o {params.hifiasm_prefix}  --purge-max $MAX_DEPTH --primary -t {threads} -l {params.purgel} {input.hifi_fastq}) 2> {log}
		"""

rule hifiasm_hic:
	"""Execute Hifiasm to generate primary and alternative assemblies in phasing mode. HiC data required"""

	input:
		hifi_fastq = expand("{sample}", sample=samples['hifi_data']),
		hic1 = hic_R1,
		hic2 = hic_R2,
		maximum_depth = os.path.join(config['results'], prefix, "genome_profiling", "maximum_depth.txt")
	output:
		asm_p_fa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic", prefix + '.asm.hic.p_ctg.gfa'),
		asm_hap1_fa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic", prefix + '.asm.hic.hap1.p_ctg.gfa'),
		asm_hap2_fa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic", prefix + '.asm.hic.hap2.p_ctg.gfa'),
		asm_a_fa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic", prefix + '.asm.hic.a_ctg.gfa')
	params:
		hifiasm_prefix = os.path.join(config['results'], prefix, "contigging/hifiasm_hic", prefix + '.asm'),
		purgel = purgel
	log:
		os.path.join(config['results'], "logs", prefix, 'hifiasm_hic.log')
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/hifiasm.yml')
	threads:
		resource['hifiasm_hic']['threads']
	resources:
		mem_mb=resource['hifiasm_hic']['mem_mb'],
		time=resource['hifiasm_hic']['time']

	shell:
		"""
		MAX_DEPTH=$(cat {input.maximum_depth})
		(hifiasm -o {params.hifiasm_prefix} --purge-max $MAX_DEPTH --primary -t {threads} --h1 {input.hic1} --h2 {input.hic2} -l {params.purgel} {input.hifi_fastq}) 2> {log}
		"""

rule hifiasm_hic_purging:
	"""Execute Hifiasm with phasing mode to generate primary and alternative assemblies. HiC data required. Multiple purging levels are tested: 1, 2 and 3"""
	
	input:
		hifi_fastq = expand("{sample}", sample=samples['hifi_data']),
		hic1 = hic_R1,
		hic2 = hic_R2,
		maximum_depth = os.path.join(config['results'], prefix, "genome_profiling", "maximum_depth.txt")
	output:
		l1_hap1_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l1/", prefix + '.asm.hic.hap1.p_ctg.gfa'),
		l1_hap2_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l1/", prefix + '.asm.hic.hap2.p_ctg.gfa'),
		l1_p_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l1/", prefix + '.asm.hic.p_ctg.gfa'),
		l1_a_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l1/", prefix + '.asm.hic.a_ctg.gfa'),
		l2_hap1_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l2/", prefix + '.asm.hic.hap1.p_ctg.gfa'),
		l2_hap2_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l2/", prefix + '.asm.hic.hap2.p_ctg.gfa'),
		l2_p_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l2/", prefix + '.asm.hic.p_ctg.gfa'),
		l2_a_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l2/", prefix + '.asm.hic.a_ctg.gfa'),
		l3_hap1_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l3/", prefix + '.asm.hic.hap1.p_ctg.gfa'),
		l3_hap2_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l3/", prefix + '.asm.hic.hap2.p_ctg.gfa'),
		l3_p_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l3/", prefix + '.asm.hic.p_ctg.gfa'),
		l3_a_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l3/", prefix + '.asm.hic.a_ctg.gfa')
	params:
		hifiasm_path = os.path.join(config['results'], prefix, "contigging/hifiasm_hic"),
		hifiasm_prefix_l1 = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l1/", prefix + '.asm'),
		hifiasm_prefix_l2 = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l2/", prefix + '.asm'),
		hifiasm_prefix_l3 = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l3/", prefix + '.asm'),
	log:
		l1_log = os.path.join(config['results'], "logs", prefix, 'hifiasm_hic_purging_l1.log'),
		l2_log = os.path.join(config['results'], "logs", prefix, 'hifiasm_hic_purging_l2.log'),
		l3_log = os.path.join(config['results'], "logs", prefix, 'hifiasm_hic_purging_l3.log')
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/hifiasm.yml')
	threads:
		resource['hifiasm_hic_purging']['threads']
	resources:
		mem_mb=resource['hifiasm_hic_purging']['mem_mb'],
		time=resource['hifiasm_hic_purging']['time']
	shell:
		"""
		MAX_DEPTH=$(cat {input.maximum_depth})

		(hifiasm -o {params.hifiasm_prefix_l1} -t {threads} --primary -l1 --purge-max $MAX_DEPTH --h1 {input.hic1} --h2 {input.hic2} {input.hifi_fastq}) 2> {log.l1_log}

		mv {params.hifiasm_prefix_l1}.ovlp.source.bin {params.hifiasm_prefix_l2}.ovlp.source.bin
		mv {params.hifiasm_prefix_l1}.ovlp.reverse.bin {params.hifiasm_prefix_l2}.ovlp.reverse.bin
		mv {params.hifiasm_prefix_l1}.ec.bin {params.hifiasm_prefix_l2}.ec.bin

		(hifiasm -o {params.hifiasm_prefix_l2} -t {threads} --primary -l2 --purge-max $MAX_DEPTH --h1 {input.hic1} --h2 {input.hic2} {input.hifi_fastq}) 2> {log.l2_log}

		mv {params.hifiasm_prefix_l2}.ovlp.source.bin {params.hifiasm_prefix_l3}.ovlp.source.bin
		mv {params.hifiasm_prefix_l2}.ovlp.reverse.bin {params.hifiasm_prefix_l3}.ovlp.reverse.bin
		mv {params.hifiasm_prefix_l2}.ec.bin {params.hifiasm_prefix_l3}.ec.bin

		(hifiasm -o {params.hifiasm_prefix_l3} -t {threads} --primary -l3 --purge-max $MAX_DEPTH --h1 {input.hic1} --h2 {input.hic2} {input.hifi_fastq}) 2> {log.l3_log}	
		"""


rule hicanu:
	"""Execute HiCanu to generate primary and alternative assemblies."""

	input:
		hifi_fastq = expand("{sample}", sample=samples['hifi_data']),
		estimated_genome_size = os.path.join(config['results'], prefix, "genome_profiling", "estimated_genome_size.txt")
	output:
		outdir = directory(os.path.join(config['results'], prefix, 'contigging/hicanu')),
		outasm = os.path.join(config['results'], prefix, 'contigging/hicanu/asm.contigs.fasta')
	log:
		os.path.join(config['results'], "logs", prefix, 'hicanu.log')
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/hicanu.yml')
	threads:
		resource['hicanu']['threads']
	resources:
		mem_mb=resource['hicanu']['mem_mb'],
		time=resource['hicanu']['time']
	shell:
		"""
		#mem_gb=$(({resources.mem_mb} / 1000))
		estimated_genome_size=$(cat {input.estimated_genome_size})
		(canu -p asm -d {output.outdir} genomeSize=$estimated_genome_size -pacbio-hifi {input.hifi_fastq} useGrid=false) 2> {log}
		"""

rule flye:
	"""Execute Flye to generate primary and alternative assemblies."""

	input:
		hifi_fastq = expand("{sample}", sample=samples['hifi_data']),
		estimated_genome_size = os.path.join(config['results'], prefix, "genome_profiling", "estimated_genome_size.txt")
	output:
		outdir = directory(os.path.join(config['results'], prefix, 'contigging/flye')),
		asm_flye = os.path.join(config['results'],  prefix, 'contigging/flye/assembly.fasta')
	params:
		estimated_genome_size = os.path.join(config['results'], prefix, "genome_profiling", "estimated_genome_size.txt"),
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/flye.yml')
	threads:
		resource['flye']['threads']
	resources:
		mem_mb=resource['flye']['mem_mb'],
		time=resource['flye']['time']
	log:
		os.path.join(config['results'], "logs", prefix, 'flye.log')
	shell:
		"""
		estimated_genome_size=$(cat {input.estimated_genome_size})
		(flye --threads {threads} --out-dir {output.outdir} --genome-size $estimated_genome_size --pacbio-hifi {input.hifi_fastq}) 2> {log}
		"""

rule hifiasm_gfa_to_fasta:
	"""Convert the GFA output from Hifiasm into a FASTA format."""
	input:
		asm_p_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm", prefix + '.asm.p_ctg.gfa'),
		asm_a_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm", prefix + '.asm.a_ctg.gfa'),
	output:
		asm_p_fa = os.path.join(config['results'], prefix, "contigging/hifiasm", prefix + '.asm.p_ctg.fasta'),
		asm_a_fa = os.path.join(config['results'], prefix, "contigging/hifiasm", prefix + '.asm.a_ctg.fasta')
	threads:
		resource['hifiasm_gfa_to_fasta']['threads']
	resources:
		mem_mb=resource['hifiasm_gfa_to_fasta']['mem_mb'],
		time=resource['hifiasm_gfa_to_fasta']['time']
	log:
		os.path.join(config['results'], "logs", prefix, 'hifiasm_gfa_to_fa.log')
	shell:
		"""awk '/^S/{{print ">"$2;print $3}}' {input.asm_p_gfa} > {output.asm_p_fa} 2> {log};
		awk '/^S/{{print ">"$2;print $3}}' {input.asm_a_gfa} > {output.asm_a_fa} 2>> {log}
		"""

rule hifiasm_hic_gfa_to_fasta:
	"""Convert the GFA output from Hifiasm into a FASTA format."""
	input:
		asm_p_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic", prefix + '.asm.hic.p_ctg.gfa'),
		asm_hap1_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic", prefix + '.asm.hic.hap1.p_ctg.gfa'),
		asm_hap2_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic", prefix + '.asm.hic.hap2.p_ctg.gfa'),
		asm_a_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic", prefix + '.asm.hic.a_ctg.gfa')
	output:
		asm_p_fa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic", prefix + '.asm.hic.p_ctg.fasta'),
		asm_hap1_fa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic", prefix + '.asm.hic.hap1.p_ctg.fasta'),
		asm_hap2_fa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic", prefix + '.asm.hic.hap2.p_ctg.fasta'),
		asm_a_fa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic", prefix + '.asm.hic.a_ctg.fasta')
	threads:
		resource['hifiasm_hic_gfa_to_fasta']['threads']
	resources: 
		mem_mb=resource['hifiasm_hic_gfa_to_fasta']['mem_mb'],
		time=resource['hifiasm_hic_gfa_to_fasta']['time']
	log:
		os.path.join(config['results'], "logs", prefix, 'hifiasm_hic_gfa_to_fa.log')
	shell:
		"""awk '/^S/{{print ">"$2;print $3}}' {input.asm_p_gfa} > {output.asm_p_fa} 2> {log};
		awk '/^S/{{print ">"$2;print $3}}' {input.asm_a_gfa} > {output.asm_a_fa} 2>> {log};
		awk '/^S/{{print ">"$2;print $3}}' {input.asm_hap1_gfa} > {output.asm_hap1_fa} 2>> {log};
		awk '/^S/{{print ">"$2;print $3}}' {input.asm_hap2_gfa} > {output.asm_hap2_fa} 2>> {log};
		"""


rule hifiasm_hic_purging_gfa_to_fasta:
	"""Convert the GFA output from Hifiasm into a FASTA format. """
	input:
		asm_hap1_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l{l}", prefix + '.asm.hic.hap1.p_ctg.gfa'),
		asm_hap2_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l{l}", prefix + '.asm.hic.hap2.p_ctg.gfa'),
		asm_p_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l{l}", prefix + '.asm.hic.p_ctg.gfa'),
		asm_a_gfa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l{l}", prefix + '.asm.hic.a_ctg.gfa')
	output:
		asm_hap1_fa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l{l}", prefix + '.asm.hic.hap1.p_ctg.fa'),
		asm_hap2_fa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l{l}", prefix + '.asm.hic.hap2.p_ctg.fa'),
		asm_p_fa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l{l}", prefix + '.asm.hic.p_ctg.fa'),
		asm_a_fa = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l{l}", prefix + '.asm.hic.a_ctg.fa'),	
	threads:
		resource['hifiasm_hic_gfa_to_fasta']['threads']
	resources: 
		mem_mb=resource['hifiasm_hic_gfa_to_fasta']['mem_mb'],
		time=resource['hifiasm_hic_gfa_to_fasta']['time']
	log:
		os.path.join(config['results'], "logs", prefix, 'hifiasm_hic_purging_gfa_to_fa_{l}.log')
	shell:
		"""awk '/^S/{{print ">"$2;print $3}}' {input.asm_hap1_gfa} > {output.asm_hap1_fa} 2> {log};
		awk '/^S/{{print ">"$2;print $3}}' {input.asm_hap2_gfa} > {output.asm_hap2_fa} 2> {log};
		awk '/^S/{{print ">"$2;print $3}}' {input.asm_p_gfa} > {output.asm_p_fa} 2> {log};
		awk '/^S/{{print ">"$2;print $3}}' {input.asm_a_gfa} > {output.asm_a_fa} 2> {log};
		"""
