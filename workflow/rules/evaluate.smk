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
	log:
		"logs/eval/{sample}.{library}.{depth}x.{reference}.{mapper}.{caller}.{filter}.{match}.vcfeval.log"
	params:
		score_field=lambda wildcards: config["score_fields"][wildcards.caller],
		ploidy=config["sample_ploidy"]*len(config["samples"]),
		all_records=lambda wildcards: "--all-records" if wildcards.filter=="raw" else "",
		squash_ploidy=lambda wildcards: "--squash-ploidy" if wildcards.match in ["AL", "CB"] else "",
		decompose=lambda wildcards: "--decompose" if wildcards.match=="AL" else "",
		output_mode=lambda wildcards: "combine" if wildcards.match=="CB" else "split",
		memory="40g"
	threads: 20
	resources:
		mem_gb=40
	conda:
		"../envs/rtg.yaml"
	shell:
		"(workflow/scripts/vcfeval.py \
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
			--memory {resources.mem_gb}g \
			{params.all_records} \
			{params.squash_ploidy} \
			{params.decompose} \
			)2> {log}"

rule install_starfish:
	output:
		"workflow/scripts/starfish.py"
	params:
		url="https://github.com/dancooke/starfish/raw/master/starfish.py"
	shell:
		"curl -L -o {output} {params.url}"
localrules: install_starfish

def _get_starfish_calls(wildcards):
	res = []
	for library, depth in config["runs"].items():
		for caller in config["callers"]:
			res.append("results/calls/" + wildcards.sample + "." + library + "." + str(depth) + "x." + wildcards.reference + "." + wildcards.mapper + "." + caller + ".vcf.gz")
	return res

def _get_starfish_labels(wildcards):
	res = []
	for library in config["runs"].keys():
		for caller in config["callers"]:
			res.append(caller + "." + library)
	return res

rule starfish:
	input:
		starfish="workflow/scripts/starfish.py",
		reference="data/references/{reference}.sdf",
		calls=_get_starfish_calls
	output:
		directory("results/eval/{sample}.{reference}.{mapper}.{filter}.{match}.isec")
	params:
		ploidy=config["sample_ploidy"]*len(config["samples"]),
		all_records=lambda wildcards: "--all-records" if wildcards.filter=="raw" else "",
		squash_ploidy=lambda wildcards: "--squash-ploidy" if wildcards.match=="AL" else "",
		decompose=lambda wildcards: "--decompose" if wildcards.match=="AL" else "",
		names=_get_starfish_labels
	log:
		"logs/eval/{sample}.{reference}.{mapper}.{filter}.{match}.isec.log"
	threads: 20
	resources:
		mem_gb=40
	conda:
		"../envs/rtg.yaml"
	shell:
		"(python {input.starfish} \
		 -t {input.reference} \
		 -V {input.calls} \
		 -O {output} \
		 --ref-overlap \
		 --ploidy {params.ploidy} \
		 --names {params.names} \
		 --threads {threads} \
		 --memory {resources.mem_gb}g \
		 {params.all_records} \
		 {params.squash_ploidy} \
		 {params.decompose} \
		 --verbose \
		 )2> {log}"
