rule freebayes_call:
	input:
		reference="data/references/{reference}.fa",
		bam="data/reads/mapped/{sample}.{library}.{depth}x.{reference}.{mapper}.bam",
		bai="data/reads/mapped/{sample}.{library}.{depth}x.{reference}.{mapper}.bam.bai",
		bed="data/references/{reference}.autosomes.bed"
	output:
		vcf="results/calls/{sample}.{library}.{depth}x.{reference}.{mapper}.FreeBayes.vcf.gz",
		vcf_index="results/calls/{sample}.{library}.{depth}x.{reference}.{mapper}.FreeBayes.vcf.gz.tbi"
	params:
		ploidy=config["sample_ploidy"]*len(config["samples"]),
		chunk_size=50000000,
		other=config["caller_options"]["FreeBayes"] if "FreeBayes" in config["caller_options"] else ""
	log:
		"logs/freebayes/{sample}.{library}.{depth}x.{reference}.{mapper}.log"
	benchmark:
		"results/benchmarks/freebayes/{sample}.{library}.{depth}x.{reference}.{mapper}.tsv"
	threads: 20
	conda:
		"../envs/freebayes.yaml"
	shell:
		"(freebayes-parallel \
			 <(bedtools makewindows -b {input.bed} -w {params.chunk_size} | awk '{{print $1\":\"$2\"-\"$3}}') {threads} \
			 -f {input.reference} \
			 -b {input.bam} \
			 --ploidy {params.ploidy} \
			 -= \
			 {params.other} | \
		bcftools filter \
			 -i 'QUAL > 1 & GQ > 1 & SAF > 0 & SAR > 0' \
			 -s FAIL \
			 -Oz -o {output.vcf} && \
		tabix {output.vcf}) \
			 2> {log}"
