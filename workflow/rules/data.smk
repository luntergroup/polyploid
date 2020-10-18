rule download_hs38DH:
	output:
		"data/references/hs38DH.fa",
		"data/references/hs38DH.fa.alt"
	conda:
		"../envs/bwakit.yaml"
	shell:
		"""
		run-gen-ref hs38DH
		mv hs38DH.fa* $(dirname {output[0]})
		"""

rule download_GRCh38:
	output:
		"data/references/GRCh38.fa"
	params:
		url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
	shell:
		"curl {params.url} | gzip -d > {output}"

rule generate_hs38DH_chromosomes_bed:
	input:
		"data/references/hs38DH.fa.fai"
	output:
		"data/references/hs38DH.chromosomes.bed"
	shell:
		"head -23 {input} | awk -v OFS='\t' '{{print $1,0,$2}}' > {output}"

rule generate_GRCh38_chromosomes_bed:
	input:
		"data/references/GRCh38.fa.fai"
	output:
		"data/references/GRCh38.chromosomes.bed"
	shell:
		"head -23 {input} | awk -v OFS='\t' '{{print $1,0,$2}}' > {output}"

localrules: download_hs38DH, download_GRCh38, generate_hs38DH_chromosomes_bed, generate_GRCh38_chromosomes_bed

import re

rule download_giab:
	output:
		vcf="data/truth/{sample}.{reference}.vcf.gz",
		vcf_index="data/truth/{sample}.{reference}.vcf.gz.tbi",
		bed="data/truth/{sample}.{reference}.bed"
	wildcard_constraints:
		sample='|'.join([re.escape(s) for s in config["samples"]])
	params:
		url_prefix="ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020"
	shell:
		"""
		curl -o {output.vcf} {params.url_prefix}/{wildcards.sample}_GRCh38_1_22_v4.2_benchmark.vcf.gz
		curl -o {output.vcf_index} {params.url_prefix}/{wildcards.sample}_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi
		curl -o {output.bed} {params.url_prefix}/{wildcards.sample}_GRCh38_1_22_v4.2_benchmark.bed
		"""
		
localrules: download_giab

MAX_DEPTH = str(int(config["sample_depth"]) * len(config["samples"]))
MIXED_SAMPLE = '+'.join(config["samples"])

rule mix_paired_reads:
	input:
		expand("data/reads/raw/{sample}.{{library}}." + str(config["sample_depth"]) + "x.{{strand}}.fastq.gz", sample=config["samples"])
	output:
		"data/reads/raw/" + MIXED_SAMPLE + ".{library}." + MAX_DEPTH + "x.{strand}.fastq.gz"
	shell:
		"cat {input} > {output}"

localrules: mix_paired_reads

rule mix_reads:
	input:
		expand("data/reads/raw/{sample}.{{library}}." + str(config["sample_depth"]) + "x.fastq.gz", sample=config["samples"])
	output:
		"data/reads/raw/" + MIXED_SAMPLE + ".{library}." + MAX_DEPTH + "x.fastq.gz"
	shell:
		"cat {input} > {output}"

localrules: mix_reads

rule downsample_paired_fastq:
	input:
		rules.mix_paired_reads.output
	output:
		"data/reads/raw/" + MIXED_SAMPLE + ".{library}.{depth}x.{strand}.fastq.gz"
	params:
		max_depth=MAX_DEPTH
	conda:
		"../envs/seqtk.yaml"
	shell:
		"bc<<<'scale=10; {wildcards.depth}/{params.max_depth}' | \
		 xargs -I{{}} seqtk sample {input} {{}} \
		 | gzip > {output}"

rule downsample_fastq:
	input:
		rules.mix_reads.output
	output:
		"data/reads/raw/" + MIXED_SAMPLE + ".{library}.{depth}x.fastq.gz"
	params:
		max_depth=MAX_DEPTH
	conda:
		"../envs/seqtk.yaml"
	shell:
		"seqtk sample {input} \
		 <(bc<<<'scale=10; {wildcards.depth}/{params.max_depth}') \
		 | gzip > {output}"
		 
def read_sample_name(vcf_filename):
	return ps.VariantFile(vcf_filename).header.samples[0]

rule make_truth_vcf:
	input:
		vcfs=expand("data/truth/{sample}.{{reference}}.vcf.gz", sample=config["samples"]),
		vcf_indices=expand("data/truth/{sample}.{{reference}}.vcf.gz.tbi", sample=config["samples"])
	output:
		vcf="data/truth/" + MIXED_SAMPLE + ".{reference}.vcf.gz",
		vcf_index="data/truth/" + MIXED_SAMPLE + ".{reference}.vcf.gz.tbi"
	params:
		gt = ' + "|" + '.join(["SAMPLES[" + str(i) + "].GT" for i in range(len(config["samples"]))])
	conda:
		"../envs/rtg.yaml"
	shell:
		"bcftools merge -0 --force-samples {input.vcfs} | \
		rtg vcffilter -i - -o - \
		 --javascript 'function record() {{SAMPLES[0].GT = {params.gt}}}' | \
		rtg vcfannotate -i - \
		 --relabel <(cat <(bcftools view -h {input.vcfs[0]} | tail -1 | awk '{{print $NF}}') <(echo TRUTH) | tr '\n' ' ') \
		 -o - | \
		rtg vcfsubset -i - --keep-sample TRUTH -o {output.vcf}"

localrules: make_truth_vcf

rule make_truth_bed:
	input:
		expand("data/truth/{sample}.{{reference}}.bed", sample=config["samples"])
	output:
		"data/truth/" + MIXED_SAMPLE + ".{reference}.bed"
	conda:
		"../envs/bedtools.yaml"
	shell:
		"bedtools multiinter -i {input} | cut -f 1,2,3 > {output}"

localrules: make_truth_bed