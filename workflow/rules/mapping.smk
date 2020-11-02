rule bwa_index:
	input:
		fa="data/references/{reference}.fa",
		fai="data/references/{reference}.fa.fai"
	output:
		"data/references/{reference}.fa.amb",
		"data/references/{reference}.fa.ann",
		"data/references/{reference}.fa.bwt",
		"data/references/{reference}.fa.pac",
		"data/references/{reference}.fa.sa"
	conda:
		"../envs/bwa.yaml"
	shell:
		"bwa index {input.fa}"

localrules: bwa_index

rule bwa_map:
	input:
		rules.bwa_index.output,
		fa="data/references/{reference}.fa",
		fai="data/references/{reference}.fa.fai",
		fq1="data/reads/raw/{sample}.{library}.{depth}x.R1.fastq.gz",
		fq2="data/reads/raw/{sample}.{library}.{depth}x.R2.fastq.gz"
	output:
		"data/reads/mapped/{sample}.{library}.{depth}x.{reference}.bwa.bam"
	params:
		rg=r"@RG\tID:{sample}\tSM:{sample}\tLB:{library}\tPU:Illumina",
		sort_memory_per_thread="4G"
	log:
		"logs/bwa/{sample}.{library}.{depth}x.{reference}.log"
	threads: 20
	conda:
		"../envs/bwa.yaml"
	shell:
		"(bwa mem -t {threads} -R '{params.rg}' {input.fa} {input.fq1} {input.fq2} | \
		  samtools view -bh | \
		  samtools sort -@ {threads} -m {params.sort_memory_per_thread} -o {output}) \
		 2> {log}"

rule pbmm2_map:
	input:
		"data/references/{reference}.fa",
		"data/reads/raw/{sample}.{library}.{depth}x.fastq.gz"
	output:
		"data/reads/mapped/{sample}.{library}.{depth}x.{reference}.pbmm2.bam"
	params:
		rg=r"@RG\tID:{sample}\tSM:{sample}\tLB:{library}\tPU:PacBio"
	log:
		"logs/pbmm2/{sample}.{library}.{depth}x.{reference}.log"
	threads: 20
	conda:
		"../envs/pbmm2.yaml"
	shell:
		"(pbmm2 align {input} -o {output} \
		--preset HIFI \
		-j {threads} -J {threads} \
		-rg {params.rg} \
		--sort) \
		 2> {log}"

rule samtools_index_bam:
	input:
		"data/reads/mapped/{sample}.{library}.{depth}x.{reference}.{mapper}.bam"
	output:
		"data/reads/mapped/{sample}.{library}.{depth}x.{reference}.{mapper}.bam.bai"
	conda:
		"../envs/samtools.yaml"
	shell:
		"samtools index {input}"

localrules: samtools_index_bam
