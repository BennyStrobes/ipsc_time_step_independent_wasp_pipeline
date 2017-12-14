#!/bin/bash
#SBATCH --time=7:00:00 --mem=10GB --partition=broadwl

INDIVIDUAL="$1"
time_step="$2"
genotype_dir="$3"
raw_allelic_counts_dir="$4"
chrom_info_file="$5"
parameter_string="$6"
target_regions_dir="$7"
cht_input_file_dir="$8"


ALL_SAMPLES_FILE=$genotype_dir"all_genotyped_samples.txt"
target_region_file=$target_regions_dir"target_regions_"$parameter_string"_merged.txt"

date
#
# create CHT input file for this individual
#
python extract_haplotype_read_counts.py \
    --chrom $chrom_info_file \
    --snp_index $genotype_dir"snp_index.h5" \
    --snp_tab $genotype_dir"snp_tab.h5" \
    --geno_prob $genotype_dir"geno_probs.h5" \
    --haplotype $genotype_dir"haps.h5" \
    --samples $ALL_SAMPLES_FILE \
    --individual $INDIVIDUAL \
    --ref_as_counts $raw_allelic_counts_dir"ref_as_counts."$INDIVIDUAL"_"$time_step".h5" \
    --alt_as_counts $raw_allelic_counts_dir"alt_as_counts."$INDIVIDUAL"_"$time_step".h5" \
    --other_as_counts $raw_allelic_counts_dir"other_as_counts."$INDIVIDUAL"_"$time_step".h5" \
    --read_counts $raw_allelic_counts_dir"read_counts."$INDIVIDUAL"_"$time_step".h5" \
    $target_region_file \
    | gzip > $cht_input_file_dir"haplotype_read_counts_"$parameter_string"."$INDIVIDUAL"_"$time_step".txt.gz"


date