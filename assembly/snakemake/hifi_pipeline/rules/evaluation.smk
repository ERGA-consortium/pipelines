rule merqury:
	""" Execute Merqury following Hifiasm to evaluate the quality and accuracy of the assembled genome."""

	input:
		asm_p = contig_output_primary(contig_assembler, phasing_mode),
		meryl_db = os.path.join(config['results'] , prefix, "genome_profiling", (prefix +".meryl"))
	output:
		out_qv = os.path.join(evaluation_folder(contig_assembler, phasing_mode), 'merqury', prefix + '.qv'),
		out_completeness = os.path.join(evaluation_output_folder, 'merqury', prefix + '.completeness.stats')
	params:
		merqury_prefix = prefix, 
		out_folder = evaluation_folder(contig_assembler, phasing_mode),
		asm_a = os.path.join(config['results'], prefix, "contigging/hifiasm", prefix + '.asm.a_ctg.fasta') if (contig_assembler == "hifiasm" and not phasing_mode) else "" 
	conda:
		os.path.join(workflow.basedir, 'envs/evaluation.yml')
	threads:
		resource['merqury']['threads']
	resources:
		mem_mb=resource['merqury']['mem_mb'],
		time=resource['merqury']['time']
	log:
		os.path.join(config['results'], "logs", prefix, 'merqury.log')
	shell:
		"""
		cd {params.out_folder}
		mkdir -p merqury
		cd merqury
		export OMP_NUM_THREADS={threads}
		(merqury.sh {input.meryl_db} {input.asm_p} {params.asm_a} {params.merqury_prefix}) &> {log}
		"""


rule gfastats:
	""" run gfastats on assembly"""

	input:
		asm_p = contig_output_primary(contig_assembler, phasing_mode)
	output:
		out_gfastats = os.path.join(evaluation_output_folder, "stats", 'gfastats.txt')
	conda:
		os.path.join(workflow.basedir, 'envs/evaluation.yml')
	threads:
		resource['gfastats']['threads']
	resources:
		mem_mb=resource['gfastats']['mem_mb'],
		time=resource['gfastats']['time']
	log:
		os.path.join(config['results'], "logs", prefix, 'gfastats.log')
	shell:
		"(gfastats --nstar-report {input.asm_p} > {output}) 2> {log};"


rule merqury_purging:
	""" Execute Merqury following Hifiasm to evaluate the quality and accuracy of the assembled genome."""
	input:
		asm_p = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l{l}", prefix + '.asm.hic.p_ctg.fa'),
		meryl_db = os.path.join(config['results'] , prefix, "genome_profiling", (prefix +".meryl"))
	output:
		out_qv = os.path.join(config['results'], prefix, "assembly_evaluation/hifiasm_hic_l{l}/merqury", prefix + ".qv"),
		out_completeness = os.path.join(config['results'], prefix, "assembly_evaluation/hifiasm_hic_l{l}/merqury", prefix + ".completeness.stats")
	params:
		merqury_prefix = prefix,
		out_folder = os.path.join(config['results'], prefix, "assembly_evaluation/hifiasm_hic_l{l}"),
		asm_a = os.path.join(config['results'], prefix, "contigging/hifiasm", prefix + '.asm.a_ctg.fasta') if (contig_assembler == "hifiasm" and not phasing_mode) else "" 
	conda:
		os.path.join(workflow.basedir, 'envs/evaluation.yml')
	threads:
		resource['merqury']['threads']
	resources:
		mem_mb=resource['merqury']['mem_mb'],
		time=resource['merqury']['time']
	log:
		os.path.join(config['results'], "logs", prefix, 'merqury_hifiasm_l{l}.log')
	shell:
		"""
		cd {params.out_folder}
		mkdir -p merqury
		cd merqury
		export OMP_NUM_THREADS={threads}
		(merqury.sh {input.meryl_db} {input.asm_p} {params.asm_a} {params.merqury_prefix}) 2> {log}
		"""

rule gfastats_purging:
	""" run gfastats on assembly"""

	input:
		asm_p = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l{l}", prefix + '.asm.hic.p_ctg.fa'),
	output:
		out_gfastats = os.path.join(config['results'], prefix, "assembly_evaluation/hifiasm_hic_l{l}", "stats", 'gfastats.txt')
	conda:
		os.path.join(workflow.basedir, 'envs/evaluation.yml')
	threads:
		resource['gfastats']['threads']
	resources:
		mem_mb=resource['gfastats']['mem_mb'],
		time=resource['gfastats']['time']
	log:
		os.path.join(config['results'], "logs", prefix, 'gfastats_purging_l{l}.log')
	shell:
		"(gfastats --nstar-report {input.asm_p} > {output}) 2> {log};"

