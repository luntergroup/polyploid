#!/usr/bin/env python3

import argparse
import csv
import concurrent.futures
import subprocess as sp
from pathlib import Path

def read_bed_regions(bed_filename):
    with bed_filename.open() as bed:
        bedreader = csv.reader(bed, delimiter='\t')
        res = []
        for row in bedreader:
            res.append(row[0] + ':' + str(int(row[1]) + 1) + '-' + row[2])
        return res

def remove_vcf(vcf_filename, remove_index=True):
    vcf_filename.unlink()
    if remove_index:
        vcf_index_filename = vcf_filename.with_suffix(vcf_filename.suffix + '.tbi')
        if vcf_index_filename:
            vcf_index_filename.unlink()

def gatk_call_helper(reference, bam, region, output, ploidy):
    sp.call([
        "gatk", \
        "--java-options", "-Xmx4g -XX:ParallelGCThreads=1", \
        "HaplotypeCaller", \
        "-R", reference, \
        "-I", bam, \
        "-L", region, \
        "-O", output, \
        "-ploidy", str(ploidy), \
        "--native-pair-hmm-threads", "1", \
        "-stand-call-conf", "10"
    ])

def main(args):
    regions = read_bed_regions(args.regions)
    tmp_dir = args.output.with_suffix('.tmp')
    tmp_dir.mkdir(exist_ok=True)
    tmp_vcfs = [tmp_dir / (region.replace(':', '_') + '.vcf.gz') for region in regions]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        for region, tmp_vcf in zip(regions, tmp_vcfs):
            futures.append(executor.submit(gatk_call_helper, reference=args.reference, bam=args.bam, region=region, output=tmp_vcf, ploidy=args.sample_ploidy))
        for future in concurrent.futures.as_completed(futures):
            print(future.result())
    sp.call(['bcftools', 'concat', '-Oz', '-o', args.output] + tmp_vcfs)
    sp.call(['tabix', args.output])
    for tmp_vcf in tmp_vcfs:
        remove_vcf(tmp_vcf)
    tmp_dir.unlink()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-R', '--reference',
                        type=Path,
                        required=True,
                        help='Reference FASTA')
    parser.add_argument('-I', '--bam',
                        type=Path,
                        required=True,
                        help='Input BAM')
    parser.add_argument('-L', '--regions',
                        type=Path,
                        required=True,
                        help='Regions BED')
    parser.add_argument('-O', '--output',
                        type=Path,
                        required=True,
                        help='Output VCF')
    parser.add_argument('--sample-ploidy',
                        type=int,
                        default=2,
                        help='Sample ploidy')
    parser.add_argument('--threads',
                        type=int,
                        default=1,
                        help='Threads')
    parsed, unparsed = parser.parse_known_args()
    main(parsed)
