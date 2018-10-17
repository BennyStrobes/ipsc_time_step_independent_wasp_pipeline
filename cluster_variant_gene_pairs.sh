#!/bin/bash
#SBATCH --time=5:00:00 --mem=10GB --partition=broadwl


parameter_string="$1"
cht_output_dir="$2"
cht_visualization_dir="$3"

module load Anaconda2

pc_num="3"


python cluster_variant_gene_pairs.py $parameter_string $cht_output_dir $cht_visualization_dir $pc_num