rule vcfeval:
	input:
		reference="data/references/{reference}.sdf",
		baseline_vcf="data/truth/{sample}.{reference}.vcf.gz",
		baseline_vcf_index="data/truth/{sample}.{reference}.vcf.gz.tbi",
		evaluation_regions="data/truth/{sample}.{reference}.bed",
		calls_vcf="results/calls/{sample}.{library}.{depth}x.{reference}.{mapper}.{caller}.vcf.gz",
		calls_vcf_index="results/calls/{sample}.{library}.{depth}x.{reference}.{mapper}.{caller}.vcf.gz.tbi"
	output:
		directory("results/eval/{sample}.{library}.{depth}x.{reference}.{mapper}.{caller}.{filter}.{match}.vcfeval")
	params:
		score_field=lambda wildcards: config["score_fields"][wildcards.caller],
		ploidy=config["sample_ploidy"]*len(config["samples"]),
		all_records=lambda wildcards: "--all-records" if wildcards.filter=="raw" else "",
		squash_ploidy=lambda wildcards: "--squash-ploidy" if wildcards.match=="AL" else "",
		decompose=lambda wildcards: "--decompose" if wildcards.match=="AL" else "",
		output_mode="split",
		memory="40g"
	threads: 20
	conda:
		"../envs/rtg.yaml"
	shell:
		"workflow/scripts/vcfeval.py \
			-t {input.reference} \
			-b {input.baseline_vcf} \
			--evaluation-regions {input.evaluation_regions} \
			-c {input.calls_vcf} \
			-o {output} \
			-f {params.score_field} \
			--ploidy {params.ploidy} \
			--ref-overlap \
			--output-mode {params.output_mode} \
			--threads {threads} \
			--memory {params.memory} \
			{params.all_records} \
			{params.squash_ploidy} \
			{params.decompose}"

rule install_starfish:
	output:
		"workflow/scripts/starfish.py"
	shell:
		"""
		rm -rf starfish || 1
		git clone https://github.com/dancooke/starfish
		mv starfish/starfish.py {output}
		rm -rf starfish
		"""

rule starfish:
	input:
		starfish="workflow/scripts/starfish.py",
		reference="data/references/{reference}.sdf",
		calls=expand("results/calls/{{sample}}.{{library}}.{{depth}}x.{{reference}}.{{mapper}}.{caller}.vcf.gz", caller=config["callers"])
	output:
		"results/eval/{sample}.{library}.{depth}x.{reference}.{mapper}.{filter}.{match}.isec"
	params:
		ploidy=config["sample_ploidy"]*len(config["samples"]),
		all_records=lambda wildcards: "--all-records" if wildcards.filter=="raw" else "",
		squash_ploidy=lambda wildcards: "--squash-ploidy" if wildcards.match=="AL" else "",
		decompose=lambda wildcards: "--decompose" if wildcards.match=="AL" else ""
	conda:
		"../envs/rtg.yaml"
	shell:
		"{input.starfish} \
		 -t {input.reference} \
		 -V {input.calls} \
		 -O {output} \
		 --ref-overlap \
		 --ploidy {params.ploidy} \
		 --threads {threads} \
		 {params.all_records} \
		 {params.squash_ploidy} \
		 {params.decompose}"