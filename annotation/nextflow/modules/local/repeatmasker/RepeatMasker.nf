process REPEATMASKER {
	label 'process_medium_high'
	label 'repeatmasker'

	input:
	path(familyFa)
	path(cleanfasta)
	
	output:
	path("genome.fasta.masked.clean"), emit: masked_genome
	path("genome.fasta.masked")      , emit: raw_masked_genome
	path("genome.fasta.tbl")         , emit: summary
	path("genome.fasta.out")         , emit: annotation
	path("genome.fasta.out.gff")     , emit: annotation_gff
	path("genome.fasta.cat.gz")      , emit: concatenate
	path "versions.yml"              , emit: versions

	script:
	"""
	#Changing the name in case no masking
	mv ${cleanfasta} genome.fasta

	#Masking with RepeatMasker
	RepeatMasker -gff -xsmall -s -pa ${task.cpus} -lib ${familyFa} genome.fasta

	#Rename fasta header to only contain the first element
	awk '/^>/ {sub(/ .*/, ""); printf("%s\\n", \$0); next;} {print;}' genome.fasta.masked > genome.fasta.masked.clean

	cat <<-VERSIONS > versions.yml
	"${task.process}":
	    RepeatMasker:  \$( RepeatMasker -v | sed -e "s/RepeatMasker version //g")
	VERSIONS
	"""
}
