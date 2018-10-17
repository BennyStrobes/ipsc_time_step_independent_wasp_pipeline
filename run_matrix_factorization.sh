#!/bin/bash
#SBATCH --time=8:00:00 --partition=broadwl --mem=5GB


parameter_string="$1"
cht_output_dir="$2"
cht_visualization_dir="$3"
target_regions_dir="$4"
gencode_gene_annotation_file="$5"
num_genes="$6"
alpha="$7"
fdr="$8"
pc_num="$9"

module load Anaconda2


python run_matrix_factorization.py $parameter_string $cht_output_dir $cht_visualization_dir $target_regions_dir $gencode_gene_annotation_file $num_genes $alpha $fdr $pc_num