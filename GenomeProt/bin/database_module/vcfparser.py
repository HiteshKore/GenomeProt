#!/usr/bin/env python3

# This script parses the VCF file (v4.2) and segregates homozygous and heterozygous variants into two files.
# It can also handle a merged VCF file generated using bcftools merge or gatk joint variant calling. It simplifies the VCF by
# taking allele with maximum depth in case of mutiple variant alleles. It modifies the allele attributes accordingly

import os, re, sys
from vcfpy import Reader, Substitution
from collections import Counter

def remove_file_if_exists(filename):
    if os.path.exists(filename):
        try:
            os.remove(filename)
        except:
            pass

# Convert a VCF record to a plain VCF line
def record_to_vcf_line(record): #returns atributes for first sample in case of multiple samples
    # Extract basic fields
    chrom = record.CHROM
    pos = record.POS
    ref = record.REF

    alts = ','.join([alt.value for alt in record.ALT])          # Extract alternative alleles
    qual = record.QUAL if record.QUAL is not None else '.'      # Handle missing quality
    filters = ','.join(record.FILTER) if record.FILTER else '.'
    info = '.'                                                  # INFO field (empty in this case)

    # Join all parts to form the VCF line
    vcf_line = '\t'.join([chrom, str(pos), '.', ref, alts, qual, filters, info, "GT", "0/1"])   # Simplified VCF compatible with vcftools consensus
    return vcf_line

# Get the variant with maximum depth in case there are multiple variants at a genomic locus
def getMostFrequentVariant(record, most_likely_genotype, sum_depths, sum_ref_allele_depths):
    varmap = dict.fromkeys([alt.value for alt in record.ALT], 0)    # k: variant, v: total depth

    for call in record.calls:
        ad = call.data.get("AD", 0)
        if isinstance(ad, int) or len(ad) == 0 or ad is None:
            continue

        ad.pop(0)
        alt_alleles = list(varmap)
        for i in range(len(ad)):
            depth = ad[i]
            if depth is not None:
                varmap[alt_alleles[i]] += int(depth)

    snv_vars = {}   # k: variant, v: total depth; stores SNVs
    for k in varmap:
        if len(k) == 1 and k != '*':
            snv_vars[k] = varmap[k]

    if len(snv_vars) == 0:
        return ''

    alt_allele = max(snv_vars, key = snv_vars.get)  # Variant with maximum depth
    alt_allele_depth = snv_vars[alt_allele]         # Alt allele depth
    alt_allele_depth_new = [sum_ref_allele_depths, alt_allele_depth]
    record.ALT = [Substitution(type_ = "SNV", value = alt_allele)]

    for call in record.calls:
        call.data["AD"] = alt_allele_depth_new  # Update alternate depth
        call.data["DP"] = sum_depths            # Update total depth
        call.data["GT"] = most_likely_genotype  # Update overall genotype

    return record_to_vcf_line(record)

def main():
    args = sys.argv
    if len(args) != 3:
        print("Usage: python vcfparser.py <vcf_file> <outdir>")
        sys.exit(1)

    # Store the arguments
    [arg_vcf_path, arg_outdir] = args[1:]

    # Check the arguments
    error_message = ""
    if not os.path.isfile(arg_vcf_path):
        error_message = f"The VCF file '{arg_vcf_path}' either does not exist or is not a file."
    elif os.path.exists(arg_outdir) and not os.path.isdir(arg_outdir):
        error_message = f"'{arg_outdir}' exists but it is not a directory."

    if error_message:
        print(error_message)
        sys.exit(1)

    # If the output directory does not exist, try to create it
    if not os.path.exists(arg_outdir):
        try:
            os.mkdir(arg_outdir)
        except Exception as e:
            print(f"Failed to create the output directory '{arg_outdir}'. Error message: '{e}'.")
            sys.exit(1)

    vcf_filename = os.path.basename(arg_vcf_path)
    vcf_homo_output   = os.path.join(arg_outdir, vcf_filename.replace(".vcf",   "_homozygous.vcf"))
    vcf_hetero_output = os.path.join(arg_outdir, vcf_filename.replace(".vcf", "_heterozygous.vcf"))
    remove_file_if_exists(vcf_homo_output)
    remove_file_if_exists(vcf_hetero_output)

    with open(arg_vcf_path, 'r') as f:
        lines_to_write = ''
        for raw_line in f:
            line = raw_line.strip()
            if not line.startswith('#'):
                break
            line_to_write = line + '\n'
            if line.startswith("#CHROM"):
                cols = line.split('\t')
                if "FORMAT" not in cols:
                    continue
                modified_header = '\t'.join(cols[:cols.index("FORMAT") + 2])
                line_to_write = modified_header + '\n'
            lines_to_write += line_to_write
        with open(vcf_homo_output, 'a') as g:
            _ = g.write(lines_to_write)
        with open(vcf_hetero_output, 'a') as g:
            _ = g.write(lines_to_write)

    # Read the records in the VCF file
    lines_for_homo = ''
    lines_for_hetero = ''
    vcf_reader = Reader.from_path(arg_vcf_path)
    for record in vcf_reader:
        if len(record.REF) != 1:        # Check that the reference allele is a SNV
            continue

        num_samples = len(record.calls) # Number of samples in the analysis
        if num_samples == 0:
            continue

        sum_depths = 0
        gts = []
        sum_ref_allele_depths = 0

        for call in record.calls:
            dp = call.data.get("DP", 0)
            if dp is not None:
                sum_depths += dp

            gt = call.data.get("GT", 0)
            if '.' not in gt:
                gts.append(gt)

            ad = call.data.get("AD", 0)
            if not isinstance(ad, int) and len(ad) > 0 and ad[0] is not None:   # Handles format change in joint VCF produced by GATk and bcftool merge
                sum_ref_allele_depths += int(ad[0])

        # Calculate the fraction of samples the genotype is detected in
        per_samples_with_genotype = (len(gts) / num_samples) * 100
        if per_samples_with_genotype < 60 or sum_depths < 10:
            continue

        unique_gts_map = dict(Counter(gts))
        most_likely_genotype = max(unique_gts_map, key = unique_gts_map.get)

        splitted = re.split(r"[\/|]", most_likely_genotype)[:2] # Homozygous records '|' phased genotype: when maternal and paternal alleles are known
        if int(splitted[0]) > 0 and int(splitted[1]) > 0 and sum_ref_allele_depths < 10:    # There might be some misalignments if the reference allele depth is less than 10' mutation will still be considered homozygous
            vcf_rec = getMostFrequentVariant(record, most_likely_genotype, sum_depths, sum_ref_allele_depths)
            if vcf_rec.strip():
                lines_for_homo += vcf_rec + '\n'
        elif sum_ref_allele_depths >= 10:   # Check if reference allele depth is not zero
            vcf_rec = getMostFrequentVariant(record, most_likely_genotype, sum_depths, sum_ref_allele_depths)
            if vcf_rec.strip():
                lines_for_hetero += vcf_rec + '\n'
        else:
            continue

    if lines_for_homo:
        with open(vcf_homo_output, 'a') as f:
            _ = f.write(lines_for_homo)

    if lines_for_hetero:
        with open(vcf_hetero_output, 'a') as f:
            _ = f.write(lines_for_hetero)

if __name__ == "__main__":
    main()