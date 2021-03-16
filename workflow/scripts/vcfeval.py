#!/usr/bin/env python3

import argparse
import csv
import pysam as ps
import subprocess as sp
from pathlib import Path
import random

rtg_bin = Path("rtg") 

def run_vcfeval(ref_sdf, 
                baseline_vcf, 
                call_vcf, 
                out_dir,
                bed_regions=None, 
                evaluation_regions=None,
                all_records=False,
                ref_overlap=False, 
                squash_ploidy=False,
                sample=None,
                score_field=None,
                decompose=False,
                ploidy=None,
                output_mode=None,
                flag_alternates=False,
                threads=None,
                memory=None):
	cmd = [rtg_bin]
	if memory is not None:
		cmd.append('RTG_MEM=' + memory)
	cmd += ['vcfeval', '-t', ref_sdf, '-b', baseline_vcf, '-c', call_vcf, '-o', out_dir]
	if bed_regions is not None:
		cmd += ['--bed-regions', bed_regions]
	if evaluation_regions is not None:
		cmd += ['--evaluation-regions', evaluation_regions]
	if score_field is not None:
		cmd += ['-f', score_field]
	if all_records:
		cmd.append('--all-records')
	if ref_overlap:
		cmd.append('--ref-overlap')
	if squash_ploidy:
		cmd.append('--squash-ploidy')
	if sample is not None:
		cmd += ['--sample', sample]
	if decompose:
		cmd.append('--decompose')
	if threads is not None:
		cmd += ['--threads', str(threads)]
	if ploidy is not None:
		cmd += ['--sample-ploidy', str(ploidy)]
	if output_mode is not None:
		cmd += ['--output-mode', output_mode]
	if flag_alternates:
		cmd.append("--XXcom.rtg.vcf.eval.flag-alternates=true")
	sp.call(cmd)
	
def index_vcf(vcf_filename, overwrite=True):
	if overwrite:
		sp.call(['tabix', '-f', vcf_filename])
	else:
		sp.call(['tabix', vcf_filename])

def decompose_multiallelic(source_vcf, dest_vcf):
	cmd = ['bcftools', 'norm', '-m', '-', '--force', '-Oz', '-o', dest_vcf, source_vcf]
	sp.call(cmd)
	index_vcf(dest_vcf)

def split_baseline_annotated(baseline_vcf_filename, tp_vcf_filename, fn_vcf_filename):
    baseline_vcf = ps.VariantFile(baseline_vcf_filename)
    tp_vcf = ps.VariantFile(tp_vcf_filename, 'wz', header=baseline_vcf.header)
    fn_vcf = ps.VariantFile(fn_vcf_filename, 'wz', header=baseline_vcf.header)
    for rec in baseline_vcf:
        if rec.info["BASE"] == "TP":
            tp_vcf.write(rec)
        elif rec.info["BASE"] == "FN":
            if "BASE_ALTERNATE" in rec.info:
                tp_vcf.write(rec)
            else:
                fn_vcf.write(rec)
    tp_vcf.close()
    index_vcf(tp_vcf_filename)
    fn_vcf.close()
    index_vcf(fn_vcf_filename)

def split_calls_annotated(calls_vcf_filename, tp_vcf_filename, fp_vcf_filename):
    calls_vcf = ps.VariantFile(calls_vcf_filename)
    tp_vcf = ps.VariantFile(tp_vcf_filename, 'wz', header=calls_vcf.header)
    fp_vcf = ps.VariantFile(fp_vcf_filename, 'wz', header=calls_vcf.header)
    for rec in calls_vcf:
        if rec.info["CALL"] == "TP":
            tp_vcf.write(rec)
        elif rec.info["CALL"] == "FP":
            if "CALL_ALTERNATE" in rec.info:
                tp_vcf.write(rec)
            else:
                fp_vcf.write(rec)
    tp_vcf.close()
    index_vcf(tp_vcf_filename)
    fp_vcf.close()
    index_vcf(fp_vcf_filename)

def count_records(vcf_filename):
    vcf = ps.VariantFile(vcf_filename, 'r')
    return sum(1 for rec in vcf)

def vcf_index_exists(vcf_filename):
    vcf_index_filaneme = vcf_filename.with_suffix(vcf_filename.suffix + '.tbi')
    return vcf_index_filaneme.exists()

def remove_vcf_index(vcf_filename):
    vcf_index_filaneme = vcf_filename.with_suffix(vcf_filename.suffix + '.tbi')
    vcf_index_filaneme.unlink()

def remove_vcf(vcf_filename, index=True):
    vcf_filename.unlink()
    if index and vcf_index_exists(vcf_filename):
        remove_vcf_index(vcf_filename)

def vcfeval_alleles_helper(ref, baseline, calls, out, bed_regions=None, evaluation_regions=None, sample=None, ref_overlap=False, ploidy=None, all_records=False, decompose=False, threads=None, memory=None):
	normed_calls = calls.with_suffix('.norm.tmp' + str(random.randint(0, 1e5)) + '.vcf.gz')
	decompose_multiallelic(calls, normed_calls)
	run_vcfeval(ref, baseline, normed_calls, out, evaluation_regions=evaluation_regions, sample=sample, ref_overlap=ref_overlap, ploidy=ploidy, all_records=all_records, decompose=decompose, squash_ploidy=True, output_mode='annotate', flag_alternates=True, threads=threads, memory=memory)
	split_baseline_annotated(out / 'baseline.vcf.gz', out / "tp-baseline.vcf.gz", out / "fn.vcf.gz")
	split_calls_annotated(out / 'calls.vcf.gz', out / "tp.vcf.gz", out / "fp.vcf.gz")
	remove_vcf(normed_calls)
	tp_baseline = count_records(out / "tp-baseline.vcf.gz")
	tp = count_records(out / "tp.vcf.gz")
	fn = count_records(out / "fn.vcf.gz")
	fp = count_records(out / "fp.vcf.gz")
	sensitivity = tp_baseline / (tp_baseline + fn)
	precision = tp_baseline / (tp_baseline + fp)
	f_measure = 2 * sensitivity * precision / (sensitivity + precision)
	with open(out / "new_summary.txt", 'w') as summary:
		summarywriter = csv.writer(summary, delimiter='\t')
		summarywriter.writerow(["Threshold", "True-pos-baseline", "True-pos-call", "False-pos", "False-neg", "Precision", "Sensitivity", "F-measure"])
		summarywriter.writerow(["None", str(tp_baseline), str(tp), str(fp), str(fn), str(precision), str(sensitivity), str(f_measure)])

def main(args):
    global rtg_bin
    rtg_bin = args.rtg
    if args.squash_ploidy and args.output_mode == "split":
        vcfeval_alleles_helper(args.sdf, args.baseline, args.calls, args.output, \
                    bed_regions=args.bed_regions, \
                    evaluation_regions=args.evaluation_regions, \
                    sample=args.sample, \
                    all_records=args.all_records, \
                    ref_overlap=args.ref_overlap, \
                    ploidy=args.ploidy, \
                    decompose=args.decompose, \
                    threads=args.threads, \
                    memory=args.memory)
    else:
        run_vcfeval(args.sdf, args.baseline, args.calls, args.output, \
                    bed_regions=args.bed_regions, \
                    evaluation_regions=args.evaluation_regions, \
                    sample=args.sample, \
                    all_records=args.all_records, \
                    ref_overlap=args.ref_overlap, \
                    ploidy=args.ploidy, \
                    decompose=args.decompose, \
                    output_mode=args.output_mode, \
                    threads=args.threads, \
                    memory=args.memory)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--sdf',
                        type=Path,
                        required=True,
                        help='RTG Tools SDF reference index')
    parser.add_argument('-b', '--baseline',
                        type=Path,
                        required=True,
                        help='Baseline VCF')
    parser.add_argument('-c', '--calls',
                        type=Path,
                        required=True,
                        help='Calls VCF')
    parser.add_argument('-o', '--output',
                        type=Path,
                        required=True,
                        help='Output directory')
    parser.add_argument('--bed-regions',
                        type=Path,
                        required=False,
                        help='regions in BED format to perform intersection')
    parser.add_argument('--evaluation-regions',
                        type=Path,
                        required=False,
                        help='regions in BED format to perform intersection')
    parser.add_argument('--squash-ploidy',
                        default=False,
                        action='store_true',
                        help='Perform haplploid matching - ignore genotype mismatches')
    parser.add_argument('--sample',
                        type=str,
                        required=False,
                        help='Sample to compare (if multiple samples in VCFs) or ALT to ignore genotypes')
    parser.add_argument('--all-records',
                        default=False,
                        action='store_true',
                        help='Intersect all records')
    parser.add_argument('--output-mode',
                        type=str,
                        required=False,
                        help='Set output mode')
    parser.add_argument('--ploidy', 
                        type=int,
                        required=False,
                        help='Set the default ploidy for comparison')
    parser.add_argument('--threads',
                        type=int,
                        required=False,
                        help='Maximum number of threads to use (default is all cores)')
    parser.add_argument('--memory',
                        type=str,
                        required=False,
                        help='Set maximum memory that can be used')
    parser.add_argument('--ref-overlap',
                        default=False,
                        action='store_true',
                        help='Call RTG vcfeval with "ref-overlap" option')
    parser.add_argument('--decompose',
                        default=False,
                        action='store_true',
                        help='Decompose multi-allelic and complex alleles')
    parser.add_argument('--rtg',
                        type=Path,
                        default=Path("rtg"),
                        help="RTG Tools binary")
    parsed, unparsed = parser.parse_known_args()
    main(parsed)
