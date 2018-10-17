#!/bin/bash
#SBATCH --time=18:00:00 --partition=broadwl --mem=12GB



parameter_string="$1"
cht_output_dir="$2"
pc_num="$3"
used_gtex_tissues_file="$4"
snp_id_to_rs_id_file="$5"
dosage_genotype_file="$6"
cm_eqtl_file="$7"
cm_egene_file="$8"
ipsc_eqtl_file="$9"
ipsc_egene_file="${10}"

python organize_tests_across_studies.py $parameter_string $cht_output_dir $pc_num $used_gtex_tissues_file $snp_id_to_rs_id_file $dosage_genotype_file $cm_eqtl_file $cm_egene_file $ipsc_eqtl_file $ipsc_egene_file
