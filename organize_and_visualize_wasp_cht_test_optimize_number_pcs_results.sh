#!/bin/bash
#SBATCH --time=1:00:00 --mem=10GB --partition=gilad


parameter_string="$1"
cht_input_file_dir="$2"
cht_output_pc_opti_dir="$3"
target_regions_dir="$4"
dosage_genotype_file="$5"
gencode_gene_annotation_file="$6"
corrected_quantile_normalized_expression="$7"


# Some parameters to this function
qval_thresh=".1"
target_region_file=$target_regions_dir"target_regions_"$parameter_string"_merged.txt"
start_chrom="2"
end_chrom="6"

if false; then
# NEED TO LOOP #PCS 0-2 (eventually 5)
# NEED TO LOOP TIME STEPS 0-15
    echo $pc_num"_"$time_step
    # Run for real data
    wasp_results_stem=$cht_output_pc_opti_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_"
    python organize_wasp_output.py $time_step $cht_output_pc_opti_dir $dosage_genotype_file $corrected_quantile_normalized_expression $gencode_gene_annotation_file $target_region_file $wasp_results_stem $start_chrom $end_chrom
    
    wasp_results_stem=$cht_output_pc_opti_dir"cht_perm_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_"
    #python organize_wasp_output.py $time_step $cht_output_pc_opti_dir $dosage_genotype_file $corrected_quantile_normalized_expression $gencode_gene_annotation_file $target_region_file $wasp_results_stem $start_chrom $end_chrom
    
    # Compute qvalues on null data
    null_file=$cht_output_pc_opti_dir"cht_perm_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_eqtl_results.txt"
    real_file=$cht_output_pc_opti_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_eqtl_results.txt"
    qvalue_file=$cht_output_pc_opti_dir"null_qvalue_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step".txt"
    # Rscript compute_qvalues.R $null_file $qvalue_file

    # Assess genome wide significance
    significant_results=$cht_output_pc_opti_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_qval_"$qval_thresh"_significant.txt"
    significant_gene_results=$cht_output_pc_opti_dir"cht_results_"$parameter_string"_num_pc_"$pc_num"_time_"$time_step"_qval_"$qval_thresh"_significant_egenes.txt"
    # python assess_wasp_significance.py $null_file $real_file $qvalue_file $significant_results $significant_gene_results $qval_thresh
fi

Rscript visualize_pc_optimization.R $cht_output_pc_opti_dir $parameter_string