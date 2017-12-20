#!/bin/bash
#SBATCH --time=4:00:00 --mem=10GB --partition=broadwl


parameter_string="$1"
cht_input_file_dir="$2"
cht_output_dir="$3"
target_regions_dir="$4"
dosage_genotype_file="$5"
gencode_gene_annotation_file="$6"
corrected_quantile_normalized_expression="$7"
time_step="$8"
eqtl_data_set_file="$9"
mvalue_file="${10}"
cht_enrichment_dir="${11}"
cis_distance="${12}"





# Some parameters to this function
qval_thresh=".1"
target_region_file=$target_regions_dir"target_regions_"$parameter_string"_merged.txt"


# Do these for varying number of PCs
for pc_num in $(seq 0 5); do
    echo $pc_num"_"$time_step
    
    # Concatenate all 22 chromosomes of real data
    # Also, modify/subset some of the columns from WASP output to make the output easier to handle
    wasp_results_stem=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_"
    python organize_wasp_output.py $time_step $cht_output_dir $dosage_genotype_file $corrected_quantile_normalized_expression $gencode_gene_annotation_file $target_region_file $wasp_results_stem
    
    # Concatenate all 22 chromosomes of permuted data
    # Also, modify/subset some of the columns from WASP output to make the output easier to handle    
    wasp_results_stem=$cht_output_dir"cht_perm_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_"
    python organize_wasp_output.py $time_step $cht_output_dir $dosage_genotype_file $corrected_quantile_normalized_expression $gencode_gene_annotation_file $target_region_file $wasp_results_stem
    
    # Compute qvalues on null data
    null_file=$cht_output_dir"cht_perm_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_eqtl_results.txt"
    real_file=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_eqtl_results.txt"
    qvalue_file=$cht_output_dir"null_qvalue_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step".txt"
    Rscript compute_qvalues.R $null_file $qvalue_file

    # Assess genome wide significance of actual data based on the qvalue threshold (FDR <= .1) in null data
    significant_results=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_qval_"$qval_thresh"_significant.txt"
    significant_gene_results=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_qval_"$qval_thresh"_significant_egenes.txt"
    python assess_wasp_significance.py $null_file $real_file $qvalue_file $significant_results $significant_gene_results $qval_thresh

    # Using reference eqtl data sets, extract pvalues belonging to reference variant-gene pairs from OUR DATA
    eqtl_file=$cht_output_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_eqtl_results.txt"
    enrichment_output_stem=$cht_enrichment_dir"enrichment_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step
    python get_snp_gene_pairs_from_other_data_and_match_for_background.py $eqtl_file $enrichment_output_stem $eqtl_data_set_file $cis_distance $mvalue_file
done