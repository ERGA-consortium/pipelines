localrules:
	get_genomescope_stats

### meryl

rule count_kmers:
	input:
		reads=lambda wildcards: samples[samples["sample_nr"] == wildcards.sample_nr]["hifi_data"]
	output:
		temp(directory(os.path.join(config['results'], prefix, "genome_profiling/{sample_nr}.meryl")))
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/genome_profiling.yml")
	threads:
		resource['meryl_count']['threads']
	resources:
		mem_mb=resource['meryl_count']['mem_mb'],
		time=resource['meryl_count']['time']
	shell:
		"meryl count k={kmer} memory={resources.mem_mb} threads={threads} {input.reads} output {output}"

rule merge_kmers:
	input:
		expand(os.path.join(config['results'], prefix, "genome_profiling/{sample_nr}.meryl"), sample_nr=samples['sample_nr']),  
	output:
		directory(os.path.join(config['results'] , prefix, "genome_profiling", (prefix +".meryl")))
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/genome_profiling.yml")
	threads:
		resource['meryl_merge']['threads']
	resources:
		mem_mb=resource['meryl_merge']['mem_mb'],
		time=resource['meryl_merge']['time']
	shell:
		"meryl union-sum memory={resources.mem_mb} threads={threads} output {output} {input} "

rule create_hist:
	input:
		os.path.join(config['results'], prefix, "genome_profiling", (prefix + ".meryl"))
	output:
		os.path.join(config['results'], prefix, "genome_profiling", (prefix + ".hist"))
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/genome_profiling.yml")
	threads:
		resource['meryl_histogram']['threads']
	resources:
		mem_mb=resource['meryl_histogram']['mem_mb'],
		time=resource['meryl_histogram']['time']
	shell:
		"meryl histogram {input} > {output}"

#### smudgeplot

rule run_smudgeplot:
	input:
		hist=os.path.join(config['results'], prefix, "genome_profiling", (prefix + ".hist"))
	output:
		plot1=os.path.join(config['results'], prefix, "genome_profiling", "smudgeplot", "results_smudgeplot.png"),
		plot2=os.path.join(config['results'], prefix, "genome_profiling", "smudgeplot", "results_smudgeplot_log10.png"),
		summary_table=os.path.join(config['results'], prefix, "genome_profiling", "smudgeplot", "results_summary_table.tsv")
	params:
		output_path=os.path.join(config['results'], prefix, "genome_profiling/smudgeplot"),
		meryl_path=os.path.join(config['results'], prefix, "genome_profiling", (prefix + ".meryl")),
		kmer=kmer
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/genome_profiling.yml")
	threads:
		resource['run_smudgeplot']['threads']
	resources:
		mem_mb=resource['run_smudgeplot']['mem_mb'],
		time=resource['run_smudgeplot']['time']
	log:
		os.path.join(config['results'], "logs",  prefix, "smudgeplot.log")
	shell:
		"""
		L=$(smudgeplot.py cutoff {input.hist} L)
		U=$(smudgeplot.py cutoff {input.hist} U)

		(meryl print less-than ${{U}} greater-than ${{L}} threads={threads} memory={resources.mem_mb} {params.meryl_path} | sort > {params.output_path}/meryl_L${{L}}_U${{U}}.dump) &> {log}

		(smudgeplot.py hetkmers -o {params.output_path}/meryl_L${{L}}_U${{U}} --middle {params.output_path}/meryl_L${{L}}_U${{U}}.dump) &> {log}

		(smudgeplot.py plot -o {params.output_path}/results {params.output_path}/meryl_L${{L}}_U${{U}}_coverages.tsv -k {params.kmer} ) &> {log}
		"""

#### genomescope

rule run_genomescope:
	input:
		os.path.join(config['results'], prefix, "genome_profiling", (prefix + ".hist"))
	output:
		directory(os.path.join(config['results'], prefix, "genome_profiling", "genomescope"))
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/genome_profiling.yml")
	threads:
		resource['run_genomescope']['threads']
	resources:
		mem_mb=resource['run_genomescope']['mem_mb'],
		time=resource['run_genomescope']['time']
	shell:
		"genomescope2 -i {input} -p {ploidy} -k {kmer} -o {output}"


rule get_genomescope_stats:
	input:
		os.path.join(config['results'], prefix, "genome_profiling", "genomescope")
	output:
		estimated_genome_size = os.path.join(config['results'], prefix, "genome_profiling", "estimated_genome_size.txt"),
		maximum_depth = os.path.join(config['results'], prefix, "genome_profiling", "maximum_depth.txt")
	threads:
		1
	shell:
		"""cd {input}
		VAR="$(grep -n "Genome Haploid Length" summary.txt | cut -f1 -d:)"
		sed -n $VAR\p summary.txt | sed -e 's/  \+/\\t/g' | cut -f3 | sed -e 's/,//g' | sed -e 's/ bp//g' > {output.estimated_genome_size}
		VAR="$(grep -n "kmercov " model.txt | cut -f1 -d:)"
		KCOV="$(printf "%.2f\\n" $(sed -n $VAR\p model.txt | sed -e 's/ \+/\\t/g' | cut -f2))"
		printf "%.0f\\n" $(echo "$KCOV * 1.5" | bc) > transition_parameter
		printf "%.0f\\n" $(echo ""$(cat transition_parameter)" * 3" | bc) > {output.maximum_depth}
		"""

