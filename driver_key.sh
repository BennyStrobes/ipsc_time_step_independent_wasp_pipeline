
######################################
# TIME STEP INDEPENDENT WASP ANALYSIS
######################################


# This scripts assumes you have run the ipsc_preproccess_pipeline first (https://github.com/BennyStrobes/ipsc_preprocess_pipeline)
# And have save the results here:
preprocess_dir="/project2/gilad/bstrober/ipsc_differentiation/preprocess/"



##########################################################################
# Input Files
##########################################################################
# Gencode hg19 gene annotation file
gencode_gene_annotation_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/gencode.v19.annotation.gtf.gz"

# Genotype directory created in preprocess pipeline (contains genotype data in h5 format)
genotype_dir=$preprocess_dir"genotype/"

# Directory created in preprocess pipeline that contains allelic count data in h5 format
raw_allelic_counts_dir=$preprocess_dir"raw_allelic_counts/"

# Downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz on 10/20/17
# Required by WASP
chrom_info_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/chromInfo.txt"

# Heterozygous probability genotype file for all cell_lines in our analysis
# These heterozygous probabilities come from impute2
het_prob_genotype_file=$preprocess_dir"genotype/YRI_het_prob_genotype.vcf"

# Dosage genotype for all cell lines in our analysis
dosage_genotype_file=$preprocess_dir"genotype/YRI_genotype.vcf"

# PCA loadings for our samples (done independently in each time step)
#  Therefor, there is one file for each time step
pca_loading_file_stem=$preprocess_dir"covariates/pca_loadings_time_step_independent_"

# Processed total expression data
# This is the quantile normalized expression data
quantile_normalized_expression=$preprocess_dir"processed_total_expression/quantile_normalized.txt"

#  Processed total expression data that was processed independently at each time point
# This is the quantile normalized data
quantile_normalized_time_step_independent_expression=$preprocess_dir"processed_total_expression/time_step_independent_quantile_normalized.txt"

# Processed and corrected total expression data
# This the quantile normalized expression data with the SVA latent factors regressed out
corrected_quantile_normalized_expression=$preprocess_dir"processed_total_expression/quant_expr_sva_corrected.txt"




##########################################################################
# Output Directories (Assumes these directories exist before script starts)
##########################################################################
# Root directory for ipsc qtl pipelines
wasp_qtl_root="/project2/gilad/bstrober/ipsc_differentiation/time_step_independent_qtl_pipelines/wasp/"

# Directory containing target regions for WASP Test
target_regions_dir=$wasp_qtl_root"target_regions/"

# Directory containing CHT input files
cht_input_file_dir=$wasp_qtl_root"cht_input_files/"

# Directory containing CHT output files 
cht_output_dir=$wasp_qtl_root"cht_output/"



##########################################################################
# Model Parameters
##########################################################################
# Max distance between gene TSS and variant
cis_distance=50000

# Test variant must MAF >= maf_cutoff in each time step
maf_cutoff=.1

#  Must be >= $min_read_count (summed across all cell lines in 1 time step) reads at heterozygous variant!
min_read_count=60

# Must be >= $min_as_read_count (summed across all cell lines in 1 time step) reads on less popular (min transformation) allele
min_as_read_count=10

# Must be $min_het_count heterozygous variant in exon of gene of interest
min_het_count=2

# Keep track of model parameters
parameter_string="cis_distance_"$cis_distance"_maf_cutoff_0"$maf_cutoff"_min_reads_"$min_read_count"_min_as_reads_"$min_as_read_count"_min_het_counts_"$min_het_count







##################################################################
### Prepare input files to Combined Haplotype Test (CHT)
##################################################################



# Run the following 4 steps in series
##################################################################

###################################################################
### Step 1: get_wasp_target_regions.sh (run in parallel for each time step)
### This script makes a list of variant gene pairs (to be tested with CHT) for each time step that pass the following filters:
######## A. Distance between variant and gene TSS is <= $cis_distance
######## B. Test variant must have MAF >= $maf_cutoff
######## C. The gene must have >= $min_het_count heterozygous variants (summed across cell lines)
######## D. The number of reads overlapping heterozygous variants (summed across cell lines) must be >= $min_read_count
######## E. The number of reads overlapping heterozygous variants that map to less popular allele (min transformation) (summed across cell lines) must be >= $min_as_read_count
### Takes about 7 hours to run per time step
### This code was initially get_target_regions.py (in WASP repo), but I made some minor edits to it
if false; then
for time_step in $(seq 0 15); do
    sbatch get_wasp_target_regions.sh $time_step $target_regions_dir $corrected_quantile_normalized_expression $genotype_dir $raw_allelic_counts_dir $chrom_info_file $dosage_genotype_file $cis_distance $maf_cutoff $min_read_count $min_as_read_count $min_het_count $gencode_gene_annotation_file
done
fi


###################################################################
### Step 2: merge_wasp_target_regions.sh
### This script makes a list of variant gene pairs (to be tested with CHT) that pass step 1's filters in ALL time steps
### Basically just take the variant-gene pairs that are found in ALL time steps
### Takes about 15 minutes to run
if false; then
sbatch merge_target_regions.sh $target_regions_dir $cis_distance $maf_cutoff $min_read_count $min_as_read_count $min_het_count
fi

###################################################################
### Step 3: wasp_haplotype_shell.sh (run in parrallel for each sample)
### This script makes a CHT input file for each sample
### This script calls extract_haplotype_read_counts.py (WASP script)
### Takes about 4 hours per sample
if false; then
for time_step in $(seq 0 15); do
    sample_file=$target_regions_dir"rna_seq_samples_"$time_step".txt"
    while read cell_line; do
        sbatch wasp_haplotype_shell.sh $cell_line $time_step $genotype_dir $raw_allelic_counts_dir $chrom_info_file $parameter_string $target_regions_dir $cht_input_file_dir
    done<$sample_file
done
fi

###################################################################
### Step 4: fit_overdispersion_parameter.sh (run in parrallel for each time step)
### This script consists of two parts:
#### A. fit_as_coefficients.py: Estimate overdispersion parameters for allele-specific test (beta binomial)
#### B. fit_bnb_coefficients.py: Estimate overdispersion parameters for association test (beta-negative binomial)
#### C. get_PCs.R: Extract PCs based on haplotype read count data
# Takes about 3 hours per time step
if false; then
for time_step in $(seq 0 15); do
    sbatch fit_overdispersion_parameters.sh $time_step $parameter_string $cht_input_file_dir $target_regions_dir
done
fi





##################################################################
### Optimize number of PCs for Combined Haplotype Test (CHT)
##################################################################

# 2-6 (inclusive are done)
for pc_num in $(seq 0 5); do
    for time_step in $(seq 0 15); do 
        for chrom_num in $(seq 1 1); do 
            sbatch submit_chrom_parallel_cht_test.sh $time_step $chrom_num $pc_num $parameter_string $cht_input_file_dir $cht_output_dir
        done
    done
done


# No Need to do this in parallel
if false; then
for pc_num in $(seq 3 3); do
    for time_step in $(seq 0 15); do 
        sbatch organize_and_visualize_wasp_cht_test_optimize_number_pcs_results.sh $parameter_string $cht_input_file_dir $cht_output_pc_opti_dir $target_regions_dir $dosage_genotype_file $gencode_gene_annotation_file $corrected_quantile_normalized_expression $pc_num $time_step
    done
done
fi
if false; then
sh organize_and_visualize_wasp_cht_test_optimize_number_pcs_results.sh $parameter_string $cht_input_file_dir $cht_output_pc_opti_dir $target_regions_dir $dosage_genotype_file $gencode_gene_annotation_file $corrected_quantile_normalized_expression
fi
