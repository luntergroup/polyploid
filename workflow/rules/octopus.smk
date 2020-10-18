rule octopus_call:
	input:
		reference="data/references/{reference}.fa",
		bam="data/reads/mapped/{sample}.{library}.{depth}x.{reference}.{mapper}.bam",
		bai="data/reads/mapped/{sample}.{library}.{depth}x.{reference}.{mapper}.bam.bai",
		bed="data/references/{reference}.chromosomes.bed"
	output:
		vcf="results/calls/{sample}.{library}.{depth}x.{reference}.{mapper}.Octopus.vcf.gz",
		vcf_index="results/calls/{sample}.{library}.{depth}x.{reference}.{mapper}.Octopus.vcf.gz.tbi"
	params:
		ploidy=config["sample_ploidy"]*len(config["samples"]),
		other=config["caller_options"]["Octopus"]
	log:
		"logs/octopus/{sample}.{library}.{depth}x.{reference}.{mapper}.log"
	benchmark:
		"results/benchmarks/octopus/{sample}.{library}.{depth}x.{reference}.{mapper}.tsv"
	threads: 20
	shell:
		"(octopus \
		 -R {input.reference} \
		 -I {input.bam} \
		 -t {input.bed} \
		 -o {output.vcf} \
		 -P {params.ploidy} \
		 --threads {threads}) \
		 {params.other} \
		 2> {log}"
