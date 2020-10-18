rule samtools_index_fasta:
	input:
		"{fasta}"
	output:
		 "{fasta}.fai"
	conda:
		"../envs/samtools.yaml"
	shell:
		"samtools faidx {input}"

localrules: samtools_index_fasta

rule rtg_format:
	input:
		"{prefix}.fa"
	output:
		directory("{prefix}.sdf")
	conda:
		"../envs/rtg.yaml"
	shell:
		"rtg format {input} -o {output}"

localrules: rtg_format
