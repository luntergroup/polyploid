from pathlib import Path
import shutil
import wget

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
		url = "ftp://ftp-trace.ncbi.nlm.nih.gov//genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
	shell:
		"curl {params.url} | gzip -d > {output}"

rule download_musa_acuminata_fasta:
	output:
		"data/references/musa_acuminata.fa"
	params:
		url = "https://banana-genome-hub.southgreen.fr/filebrowser/download/42948"
	shell:
		"curl -o {output} {params.url}"

localrules: download_musa_acuminata_fasta

rule generate_hs38DH_autosomes_bed:
	input:
		"data/references/hs38DH.fa.fai"
	output:
		"data/references/hs38DH.autosomes.bed"
	shell:
		"head -22 {input} | awk -v OFS='\t' '{{print $1,0,$2}}' > {output}"
		
rule generate_hs38DH_chromosomes_bed:
	input:
		"data/references/hs38DH.fa.fai"
	output:
		"data/references/hs38DH.chromosomes.bed"
	shell:
		"head -25 {input} | awk -v OFS='\t' '{{print $1,0,$2}}' > {output}"

rule generate_GRCh38_autosomes_bed:
	input:
		"data/references/GRCh38.fa.fai"
	output:
		"data/references/GRCh38.autosomes.bed"
	shell:
		"head -22 {input} | awk -v OFS='\t' '{{print $1,0,$2}}' > {output}"
		
rule generate_GRCh38_chromosomes_bed:
	input:
		"data/references/GRCh38.fa.fai"
	output:
		"data/references/GRCh38.chromosomes.bed"
	shell:
		"head -25 {input} | awk -v OFS='\t' '{{print $1,0,$2}}' > {output}"

localrules: download_hs38DH, download_GRCh38, generate_hs38DH_autosomes_bed, generate_hs38DH_chromosomes_bed, generate_GRCh38_autosomes_bed, generate_GRCh38_chromosomes_bed

rule generate_musa_acuminata_chromosomes_bed:
	input:
		"data/references/musa_acuminata.fa.fai"
	output:
		"data/references/musa_acuminata.chromosomes.bed"
	shell:
		"head -11 {input} | awk -v OFS='\t' '{{print $1,0,$2}}' > {output}"

rule generate_musa_acuminata_autosome_bed:
	input:
		"data/references/musa_acuminata.fa.fai"
	output:
		"data/references/musa_acuminata.autosomes.bed"
	shell:
		"head -11 {input} | awk -v OFS='\t' '{{print $1,0,$2}}' > {output}"

localrules: generate_musa_acuminata_chromosomes_bed, generate_musa_acuminata_autosome_bed

import re

def _get_ajtrio_giab_prefix(sample, reference):
	base_url = "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio"
	subdirs = {
		"HG002": "HG002_NA24385_son",
		"HG003": "HG003_NA24149_father",
		"HG004": "HG004_NA24143_mother",
	}
	giab_version = "NISTv4.2.1"
	if reference == "hs38DH":
		reference = "GRCh38"
	return f"{base_url}/{subdirs[sample]}/{giab_version}/{reference}/{sample}_{reference}_1_22_v4.2.1_benchmark"

def _get_ajtrio_giab_vcf_url(wildcards):
	return _get_ajtrio_giab_prefix(wildcards.sample, wildcards.reference) + ".vcf.gz"

def _get_ajtrio_giab_bed_url(wildcards):
	return _get_ajtrio_giab_prefix(wildcards.sample, wildcards.reference) + "_noinconsistent.bed"

rule download_giab:
	output:
		vcf="data/truth/{sample}.{reference}.vcf.gz",
		vcf_index="data/truth/{sample}.{reference}.vcf.gz.tbi",
		bed="data/truth/{sample}.{reference}.bed"
	wildcard_constraints:
		sample='|'.join([re.escape(s) for s in ["HG002", "HG003", "HG004"]])
	params:
		vcf_url=_get_ajtrio_giab_vcf_url,
		bed_url=_get_ajtrio_giab_bed_url
	shell:
		"""
		curl -o {output.vcf} {params.vcf_url}
		curl -o {output.vcf_index} {params.vcf_url}.tbi
		curl -o {output.bed} {params.bed_url}
		"""

localrules: download_giab

def try_get_link(config, sample, strand=None):
	try:
		if strand is None:
			return config["links"][sample]
		else:
			return config["links"][sample][strand]
	except:
		return None

rule download_novaseq_precision_fda_reads:
	output:
		"data/reads/raw/{sample}.NovaSeq.35x.{strand}.fastq.gz"
	wildcard_constraints:
		sample='|'.join([re.escape(s) for s in ["HG002", "HG003", "HG004"]])
	params:
		url=lambda wildcards: try_get_link(config, wildcards.sample, wildcards.strand)
	shell:
		"curl -O {output} {params.url}"

localrules: download_novaseq_precision_fda_reads

rule download_pacbio_precision_fda_reads:
	output:
		"data/reads/raw/{sample}.PacBioHiFi.35x.fastq.gz"
	wildcard_constraints:
		sample='|'.join([re.escape(s) for s in ["HG002", "HG003", "HG004"]])
	params:
		url=lambda wildcards: try_get_link(config, wildcards.sample)
	shell:
		"curl -O {output} {params.url}"

localrules: download_pacbio_precision_fda_reads

def download_ftp(link, save_path):
	wget.download(link, out=str(save_path))

def cat_files(source_filenames, destination_filename):
	with destination_filename.open('wb') as destination:
		for source_filename in source_filenames:
			shutil.copyfileobj(source_filename.open('rb'), destination)

def download_ena(bucket_url, run_accessions, strand, output):
	tmp_fastqs = []
	for run_accession, d in run_accessions:
		link_prefix = bucket_url + "/" + d + "/" + run_accession + '/' + run_accession + "_"
		temp_fastq = output.parent / (run_accession + ".R" + strand + ".fastq.gz")
		download_ftp(link_prefix + strand + ".fastq.gz", temp_fastq)
		tmp_fastqs.append(temp_fastq)
	cat_files(tmp_fastqs, output)
 	for fq in tmp_fastqs:
		fq.unlink()

rule download_banana_hiseq_reads:
	output:
		"data/reads/raw/banana.HiSeq.55x.{strand}.fastq.gz"
	params:
		bucket_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR341",
		run_accessions = [("ERR3412983", "003"), ("ERR3412984", "004")],
		strand = lambda wildcards: wildcards.strand[-1:]
	run:
		download_ena(params.bucket_url, params.run_accessions, params.strand, Path(output[0]))

rule download_banana_nextseq_reads:
	output:
		"data/reads/raw/banana.NextSeq.65x.{strand}.fastq.gz"
	params:
		bucket_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR341",
		run_accessions = [("ERR3413471", "001"), ("ERR3413472", "002"), ("ERR3413473", "003"), ("ERR3413474", "004")],
		strand = lambda wildcards: wildcards.strand[-1:]
	run:
		download_ena(params.bucket_url, params.run_accessions, params.strand, Path(output[0]))

localrules: download_banana_hiseq_reads, download_banana_nextseq_reads

if config["samples"][0] != "banana":
	SAMPLE_DEPTH = int(config["sample_depth"])
	MAX_DEPTH = str(SAMPLE_DEPTH * len(config["samples"]))
else:
	SAMPLE_DEPTH = 0
	MAX_DEPTH = 0
MIXED_SAMPLE = '+'.join(config["samples"])

rule mix_paired_reads:
	input:
		expand("data/reads/raw/{sample}.{{library}}." + str(SAMPLE_DEPTH) + "x.{{strand}}.fastq.gz", sample=config["samples"])
	output:
		"data/reads/raw/" + MIXED_SAMPLE + ".{library}." + str(MAX_DEPTH) + "x.{strand}.fastq.gz"
	shell:
		"cat {input} > {output}"

localrules: mix_paired_reads

rule mix_reads:
	input:
		expand("data/reads/raw/{sample}.{{library}}." + str(SAMPLE_DEPTH) + "x.fastq.gz", sample=config["samples"])
	output:
		"data/reads/raw/" + MIXED_SAMPLE + ".{library}." + str(MAX_DEPTH) + "x.fastq.gz"
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
		"bc<<<'scale=10; {wildcards.depth}/{params.max_depth}' | \
		 xargs -I{{}} seqtk sample {input} {{}} \
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
	params:
		num_samples=str(len(config["samples"]))
	conda:
		"../envs/bedtools.yaml"
	shell:
		"bedtools multiinter -i {input} | awk '{{if ($4=={params.num_samples}) print}}' | cut -f 1,2,3 > {output}"

localrules: make_truth_bed