#!/bin/bash
#SBATCH --time=8:00:00 --partition=broadwl --mem=5GB


parameter_string="$1"
cht_output_dir="$2"
cht_visualization_dir="$3"
pc_num="$4"
chrom_hmm_input_dir="$5"
cht_enrichment_dir="$6"
num_permutations="$7"
efdr="$8"
time_step="$9"

echo $pc_num"_"$time_step


#########################
marker_type="enhancer"
#########################

cell_line_version="ipsc_only_cell_lines"
python organize_chrom_hmm_enrichment.py $parameter_string $cht_output_dir $cht_visualization_dir $pc_num $chrom_hmm_input_dir $cht_enrichment_dir $num_permutations $time_step $marker_type $cell_line_version $efdr


cell_line_version="heart_only_cell_lines"
python organize_chrom_hmm_enrichment.py $parameter_string $cht_output_dir $cht_visualization_dir $pc_num $chrom_hmm_input_dir $cht_enrichment_dir $num_permutations $time_step $marker_type $cell_line_version $efdr

cell_line_version="heart_and_ipsc_cell_lines"
python organize_chrom_hmm_enrichment.py $parameter_string $cht_output_dir $cht_visualization_dir $pc_num $chrom_hmm_input_dir $cht_enrichment_dir $num_permutations $time_step $marker_type $cell_line_version $efdr


cell_line_version="ipsc_cell_lines"
python organize_chrom_hmm_enrichment.py $parameter_string $cht_output_dir $cht_visualization_dir $pc_num $chrom_hmm_input_dir $cht_enrichment_dir $num_permutations $time_step $marker_type $cell_line_version $efdr


cell_line_version="heart_cell_lines"
python organize_chrom_hmm_enrichment.py $parameter_string $cht_output_dir $cht_visualization_dir $pc_num $chrom_hmm_input_dir $cht_enrichment_dir $num_permutations $time_step $marker_type $cell_line_version $efdr

cell_line_version="all_cell_lines"
python organize_chrom_hmm_enrichment.py $parameter_string $cht_output_dir $cht_visualization_dir $pc_num $chrom_hmm_input_dir $cht_enrichment_dir $num_permutations $time_step $marker_type $cell_line_version $efdr



#########################
marker_type="promotor"
#########################

cell_line_version="ipsc_only_cell_lines"
python organize_chrom_hmm_enrichment.py $parameter_string $cht_output_dir $cht_visualization_dir $pc_num $chrom_hmm_input_dir $cht_enrichment_dir $num_permutations $time_step $marker_type $cell_line_version $efdr


cell_line_version="heart_only_cell_lines"
python organize_chrom_hmm_enrichment.py $parameter_string $cht_output_dir $cht_visualization_dir $pc_num $chrom_hmm_input_dir $cht_enrichment_dir $num_permutations $time_step $marker_type $cell_line_version $efdr

cell_line_version="heart_and_ipsc_cell_lines"
python organize_chrom_hmm_enrichment.py $parameter_string $cht_output_dir $cht_visualization_dir $pc_num $chrom_hmm_input_dir $cht_enrichment_dir $num_permutations $time_step $marker_type $cell_line_version $efdr


cell_line_version="ipsc_cell_lines"
python organize_chrom_hmm_enrichment.py $parameter_string $cht_output_dir $cht_visualization_dir $pc_num $chrom_hmm_input_dir $cht_enrichment_dir $num_permutations $time_step $marker_type $cell_line_version $efdr


cell_line_version="heart_cell_lines"
python organize_chrom_hmm_enrichment.py $parameter_string $cht_output_dir $cht_visualization_dir $pc_num $chrom_hmm_input_dir $cht_enrichment_dir $num_permutations $time_step $marker_type $cell_line_version $efdr

cell_line_version="all_cell_lines"
python organize_chrom_hmm_enrichment.py $parameter_string $cht_output_dir $cht_visualization_dir $pc_num $chrom_hmm_input_dir $cht_enrichment_dir $num_permutations $time_step $marker_type $cell_line_version $efdr
