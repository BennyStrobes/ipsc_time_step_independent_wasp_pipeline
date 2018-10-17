#!/bin/bash
#SBATCH --time=8:00:00 --partition=broadwl --mem=5GB


cht_output_dir="$1"
chrom_hmm_input_dir="$2"
num_permutations="$3"
ipsc_egene_file="$4"
cm_egene_file="$5"
num_permutations="$6"
ipsc_eqtl_file="$7"
cm_eqtl_file="$8"
parameter_string="$9"
cht_enrichment_dir="${10}"
cht_visualization_dir="${11}"


if false; then

#########################
marker_type="enhancer"
#########################
cell_line_version="all_cell_lines"
python banovich_ipsc_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $ipsc_egene_file $ipsc_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir


cell_line_version="heart_and_ipsc_cell_lines"
python banovich_ipsc_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $ipsc_egene_file $ipsc_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir


cell_line_version="heart_only_cell_lines"
python banovich_ipsc_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $ipsc_egene_file $ipsc_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir


cell_line_version="ipsc_only_cell_lines"
python banovich_ipsc_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $ipsc_egene_file $ipsc_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir


cell_line_version="heart_cell_lines"
python banovich_ipsc_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $ipsc_egene_file $ipsc_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir


cell_line_version="ipsc_cell_lines"
python banovich_ipsc_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $ipsc_egene_file $ipsc_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir


#########################
marker_type="promotor"
#########################
cell_line_version="all_cell_lines"
python banovich_ipsc_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $ipsc_egene_file $ipsc_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir


cell_line_version="heart_and_ipsc_cell_lines"
python banovich_ipsc_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $ipsc_egene_file $ipsc_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir


cell_line_version="heart_only_cell_lines"
python banovich_ipsc_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $ipsc_egene_file $ipsc_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir


cell_line_version="ipsc_only_cell_lines"
python banovich_ipsc_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $ipsc_egene_file $ipsc_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir


cell_line_version="heart_cell_lines"
python banovich_ipsc_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $ipsc_egene_file $ipsc_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir


cell_line_version="ipsc_cell_lines"
python banovich_ipsc_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $ipsc_egene_file $ipsc_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir


#########################
marker_type="enhancer"
#########################
cell_line_version="all_cell_lines"
python banovich_cm_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $cm_egene_file $cm_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir

cell_line_version="heart_and_ipsc_cell_lines"
python banovich_cm_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $cm_egene_file $cm_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir

cell_line_version="heart_only_cell_lines"
python banovich_cm_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $cm_egene_file $cm_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir

cell_line_version="ipsc_only_cell_lines"
python banovich_cm_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $cm_egene_file $cm_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir

cell_line_version="heart_cell_lines"
python banovich_cm_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $cm_egene_file $cm_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir

cell_line_version="ipsc_cell_lines"
python banovich_cm_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $cm_egene_file $cm_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir

#########################
marker_type="promotor"
#########################


cell_line_version="all_cell_lines"
python banovich_cm_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $cm_egene_file $cm_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir

cell_line_version="heart_and_ipsc_cell_lines"
python banovich_cm_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $cm_egene_file $cm_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir

cell_line_version="heart_only_cell_lines"
python banovich_cm_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $cm_egene_file $cm_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir

cell_line_version="ipsc_only_cell_lines"
python banovich_cm_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $cm_egene_file $cm_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir

cell_line_version="heart_cell_lines"
python banovich_cm_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $cm_egene_file $cm_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir

cell_line_version="ipsc_cell_lines"
python banovich_cm_organize_chrom_hmm_enrichment.py $cht_output_dir $chrom_hmm_input_dir $num_permutations $cm_egene_file $cm_eqtl_file $marker_type $cell_line_version $parameter_string $cht_enrichment_dir
fi



Rscript visualize_banovich_chrom_hmm_enrichment_analysis.R $cht_enrichment_dir $cht_visualization_dir $num_permutations

