rule download_germline_forest:
    output:
        "workflow/resources/octopus/forests/polyploid.forest"
    shell:
        "curl https://storage.googleapis.com/luntergroup/octopus/forests/germline.v0.7.1.forest.gz | gunzip > {output}"
localrules: download_germline_forest

rule download_polyploid_forest:
    output:
        "workflow/resources/octopus/forests/polyploid.forest"
    shell:
        "curl https://storage.googleapis.com/luntergroup/octopus/forests/polyploid.v0.7.1.forest.gz | gunzip > {output}"
localrules: download_polyploid_forest

rule octopus_call:
    input:
        reference="data/references/{reference}.fa",
        bam="data/reads/mapped/{sample}.{library}.{depth}x.{reference}.{mapper}.bam",
        bai="data/reads/mapped/{sample}.{library}.{depth}x.{reference}.{mapper}.bam.bai",
        bed="data/references/{reference}.autosomes.bed"
        # forest=rules.download_polyploid_forest.output if config["sample_ploidy"]*len(config["samples"])>2 else rules.download_germline_forest.output
    output:
        vcf="results/calls/{sample}.{library}.{depth}x.{reference}.{mapper}.Octopus.vcf.gz",
        vcf_index="results/calls/{sample}.{library}.{depth}x.{reference}.{mapper}.Octopus.vcf.gz.tbi"
    params:
        ploidy=config["sample_ploidy"]*len(config["samples"]),
        other=config["caller_options"]["Octopus"] if "Octopus" in config["caller_options"] else ""
    log:
        "logs/octopus/{sample}.{library}.{depth}x.{reference}.{mapper}.log"
    benchmark:
        "results/benchmarks/octopus/{sample}.{library}.{depth}x.{reference}.{mapper}.tsv"
    threads: 20
    container:
        "docker://dancooke/octopus"
    shell:
        "(octopus \
         -R {input.reference} \
         -I {input.bam} \
         -t {input.bed} \
         -o {output} \
         -P {params.ploidy} \
         --threads {threads} \
         {params.other}) \
         2> {log}"

rule octopus_call_longhaps:
    input:
        reference="data/references/{reference}.fa",
        bam="data/reads/mapped/{sample}.{library}.{depth}x.{reference}.{mapper}.bam",
        bai="data/reads/mapped/{sample}.{library}.{depth}x.{reference}.{mapper}.bam.bai",
        bed="data/references/{reference}.chromosomes.bed",
        source_vcf=rules.octopus_call.output.vcf
    output:
        vcf="results/calls/{sample}.{library}.{depth}x.{reference}.{mapper}.Octopus_longhaps.vcf.gz",
        vcf_index="results/calls/{sample}.{library}.{depth}x.{reference}.{mapper}.Octopus_longhaps.vcf.gz.tbi"
    params:
        ploidy=config["sample_ploidy"]*len(config["samples"]),
        other=config["caller_options"]["Octopus_longhaps"] if "Octopus_longhaps" in config["caller_options"] else ""
    log:
        "logs/octopus/longhaps/{sample}.{library}.{depth}x.{reference}.{mapper}.longhaps.log"
    benchmark:
        "results/benchmarks/octopus/longhaps/{sample}.{library}.{depth}x.{reference}.{mapper}.tsv"
    threads: 20
    container:
        "docker://dancooke/octopus:develop"
    shell:
        "(octopus \
         -R {input.reference} \
         -I {input.bam} \
         -t {input.bed} \
         -o {output} \
         -P {params.ploidy} \
         --threads {threads} \
         --disable-denovo-variant-discovery \
         --source-candidates {input.source_vcf} \
         --use-filtered-source-candidates \
         -x 400 \
         --lagging-level OPTIMISTIC \
         --backtrack-level AGGRESSIVE \
         {params.other}) \
         2> {log}"

rule octopus_recall_banana:
    input:
        reference="data/references/{reference}.fa",
        bam="data/reads/mapped/{sample}.{library}.{depth}x.{reference}.{mapper}.bam",
        bai="data/reads/mapped/{sample}.{library}.{depth}x.{reference}.{mapper}.bam.bai",
        bed="data/references/{reference}.chromosomes.bed",
        source_vcfs=expand("results/calls/banana.{library}.musa_acuminata.{caller}.vcf.gz", \
                           library=["HiSeq-1500.55x", "NextSeq-500.65x"], \
                           caller=["Octopus", "GATK4"])
    output:
        vcf="results/calls/{sample}.{library}.{depth}x.{reference}.{mapper}.Octopus_recall.vcf.gz",
        vcf_index="results/calls/{sample}.{library}.{depth}x.{reference}.{mapper}.Octopus_recall.vcf.gz.tbi"
    params:
        ploidy=config["sample_ploidy"]*len(config["samples"]),
        other=config["caller_options"]["Octopus_longhaps"] if "Octopus_recall" in config["caller_options"] else ""
    log:
        "logs/octopus/longhaps/{sample}.{library}.{depth}x.{reference}.{mapper}.longhaps.log"
    benchmark:
        "results/benchmarks/octopus/longhaps/{sample}.{library}.{depth}x.{reference}.{mapper}.tsv"
    threads: 20
    container:
        "docker://dancooke/octopus"
    shell:
        "(octopus \
         -R {input.reference} \
         -I {input.bam} \
         -t {input.bed} \
         -o {output} \
         -P {params.ploidy} \
         --threads {threads} \
         --disable-denovo-variant-discovery \
         --source-candidates {input.source_vcf} \
         --use-filtered-source-candidates \
         {params.other}) \
         2> {log}"
