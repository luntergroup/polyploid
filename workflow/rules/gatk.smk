rule gatk_create_sequence_dictionary:
	input:
		"data/references/{reference}.fa"
	output:
		"data/references/{reference}.dict"
	conda:
		"../envs/gatk4.yaml"
	shell:
		"gatk CreateSequenceDictionary -R {input}"

localrules: gatk_create_sequence_dictionary

rule gatk_mark_duplicates:
	input:
		bam="data/reads/mapped/{sample}.{library}.{depth}x.{reference}.{mapper}.bam",
		bai="data/reads/mapped/{sample}.{library}.{depth}x.{reference}.{mapper}.bam.bai"
	output:
		bam="data/reads/mapped/{sample}.{library}.{depth}x.{reference}.{mapper}.dedup.bam",
		bai="data/reads/mapped/{sample}.{library}.{depth}x.{reference}.{mapper}.dedup.bai",
		metrics=temp("data/reads/mapped/{sample}.{library}.{depth}x.{reference}.{mapper}.dedup.bam.metric.txt")
	log:
		"logs/gatk/{sample}.{library}.{depth}x.{reference}.{mapper}.dedup.log"
	benchmark:
		"benchmarks/gatk/{sample}.{library}.{depth}x.{reference}.{mapper}.dedup.tsv"
	conda:
		"../envs/gatk4.yaml"
	shell:
		"(gatk \
		 --java-options '-Xmx4g -XX:ParallelGCThreads=1' \
		 MarkDuplicates \
		 -I {input.bam} \
		 -M {output.metrics} \
		 -O {output.bam} \
		 --CREATE_INDEX) \
		 2> {log}"

rule gatk_call:
	input:
		bam=rules.gatk_mark_duplicates.output.bam,
		bai=rules.gatk_mark_duplicates.output.bai,
		reference="data/references/{reference}.fa",
		reference_dict="data/references/{reference}.dict",
		bed="data/references/{reference}.chromosomes.bed"
	output:
		vcf="results/calls/{sample}.{library}.{depth}x.{reference}.{mapper}.GATK4.raw.vcf.gz",
		vcf_index="results/calls/{sample}.{library}.{depth}x.{reference}.{mapper}.GATK4.raw.vcf.gz.tbi"
	params:
		ploidy=config["sample_ploidy"]*len(config["samples"]),
		other=config["caller_options"]["GATK4"]
	benchmark:
		"results/benchmarks/gatk/{sample}.{library}.{depth}x.{reference}.{mapper}.tsv"
	threads: 20
	conda:
		"../envs/gatk4.yaml"
	shell:
		"workflow/scripts/gatk_parallel.py \
			-R {input.reference} \
			-I {input.bam} \
			-L {input.bed} \
			-O {output.vcf} \
			--sample-ploidy {params.ploidy} \
			--threads {threads} \
			{params.other}"

rule gatk_filter:
	input:
		vcf=rules.gatk_call.output.vcf,
		vcf_index=rules.gatk_call.output.vcf_index
	output:
		vcf="results/calls/{sample}.{library}.{depth}x.{reference}.{mapper}.GATK4.vcf.gz",
		vcf_index="results/calls/{sample}.{library}.{depth}x.{reference}.{mapper}.GATK4.vcf.gz.tbi"
	log:
		"logs/gatk/{sample}.{library}.{depth}x.{reference}.{mapper}.filter.log"
	benchmark:
		"results/benchmarks/gatk/{sample}.{library}.{depth}x.{reference}.{mapper}.filter.tsv"
	conda:
		"../envs/gatk4.yaml"
	shell:
		"(gatk \
		 --java-options '-Xmx4g -XX:ParallelGCThreads=1' \
		 VariantFiltration \
		 -V {input.vcf} \
		 -O {output.vcf}) \
		 -filter 'QD < 2.0' --filter-name 'QD2' \
		 -filter 'QUAL < 50' --filter-name 'Q50' \
		 -filter 'GQ < 5' --filter-name 'GQ5' \
		 -filter 'FS > 60.0' --filter-name 'FS60' \
		 -filter 'SOR > 3.0' --filter-name 'SOR3' \
		 -filter 'MQ < 40.0' --filter-name 'MQ40' \
		 -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' \
		 -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' \
		 2> {log}"
