import numpy as np 
import os
import sys
import gzip
import pdb



def convert_from_ensamble_to_gene_symbol(ensamble_hits, gencode_file):
    f = gzip.open(gencode_file)
    gene_symbol_hits = []
    used = {}
    for line in f:
        line = line.decode('utf-8').rstrip()
        data = line.split()
        if line.startswith('#'):
            continue
        line_ensamble_id = data[9].split('"')[1].split('.')[0]
        line_gene_symbol = data[17].split('"')[1]
        if line_ensamble_id in ensamble_hits:
            gene_symbol_hits.append(line_gene_symbol)
    return np.unique(gene_symbol_hits)

def sort_gsea(save_file, new_save_file):
    f = open(save_file)
    t = open(new_save_file,'w')
    pairs = []
    for i,line in enumerate(f):
        line = line.rstrip()
        data = line.split()
        if i < 4:
            continue
        pvalue = float(data[6])
        pairs.append((pvalue, line))
    sorted_pairs = sorted(pairs, key=lambda x: x[0])
    for pair in sorted_pairs:
        liner = pair[1]
        t.write(liner + '\n')
    t.close()

def print_array(file_name, array):
    t = open(file_name,'w')
    for ele in array:
        t.write(ele + '\n')
    t.close()

def gene_set_enrichment(ensamble_hits, ensamble_background, output_root, gencode_file):
    gene_symbol_hits = convert_from_ensamble_to_gene_symbol(ensamble_hits, gencode_file)
    gene_symbol_background = convert_from_ensamble_to_gene_symbol(ensamble_background, gencode_file)
    hits_file = output_root + 'temp_hit_genes.txt'
    background_file = output_root + 'temp_background_genes.txt'
    print_array(hits_file, gene_symbol_hits)
    print_array(background_file, gene_symbol_background)

    gene_sets = {'kegg':'c2.cp.kegg.v5.1.symbols.gmt.txt', 'biocarta':'c2.cp.biocarta.v5.1.symbols.gmt.txt', 'reactome':'c2.cp.reactome.v5.1.symbols.gmt.txt', 'hallmark':'h.all.v5.1.symbols.gmt.txt'}

    for gene_set in gene_sets.keys():
        save_file = output_root + gene_set + '_gsea_output.txt'
        geneset_file = '/project2/gilad/bstrober/tools/tools/gsea/data/' + gene_sets[gene_set]
        os.system('gsea ' + hits_file + ' ' + background_file + ' ' + geneset_file + ' ' + save_file)
        new_save_file = output_root + gene_set + '_gsea_sorted_output.txt'
        sort_gsea(save_file, new_save_file)
        os.system('rm ' + save_file)
    #os.system('rm ' + hits_file)
    #os.system('rm ' + background_file)


def extract_genes(egene_file):
	f = open(egene_file)
	head_count = 0
	egenes = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble_id = data[1]
		egenes[ensamble_id] = 1
	return egenes

def extract_ipsc_egenes(egene_file):
    f = open(egene_file)
    head_count = 0
    egenes = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[0].split('.')[0]
        egenes[ensamble_id] = 0
    return egenes

def extract_top_n_ipsc_egenes(egene_file, num_genes):
    f = open(egene_file)
    head_count = 0
    egenes = {}
    arr = []
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[0].split('.')[0]
        pvalue = float(data[2])
        arr.append((ensamble_id, pvalue))
    arr.sort(key=lambda x: x[1])
    for i, tupler in enumerate(arr):
        if i > (num_genes-1):
            continue
        egenes[tupler[0]] = 0
    return egenes

def extract_ipsc_all_genes(ipsc_eqtl_file):
    f = open(ipsc_eqtl_file)
    egenes = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        ensamble_id = data[0].split('.')[0]
        egenes[ensamble_id] = 0
    return egenes

def extract_cm_egenes(egene_file):
    f = open(egene_file)
    head_count = 0
    egenes = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        egenes[data[0]] = 1
    return egenes

def extract_cm_all_genes(cm_eqtl_file):
    f = open(cm_eqtl_file)
    head_count = 0
    genes = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        genes[data[0]] = 1
    return genes

def extract_top_n_cm_egenes(eqtl_file, num_genes):
    f = open(cm_eqtl_file)
    head_count = 0
    genes = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[0]
        pvalue = float(data[2])
        if ensamble_id not in genes:
            genes[ensamble_id] = pvalue
        else:
            old_pvalue = genes[ensamble_id]
            genes[ensamble_id] = min(pvalue, old_pvalue)
    gene_arr = []
    egenes = {}
    for gene in genes.keys():
        gene_arr.append((gene,genes[gene]))
    gene_arr.sort(key=lambda x: x[1])
    for i, tupler in enumerate(gene_arr):
        if i > (num_genes-1):
            continue
        egenes[tupler[0]] = 0
    return egenes

def extract_gtex_tissue_significant_genes(egene_file):
    f = open(egene_file)
    head_count = 0
    all_genes = {}
    egenes = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[0].split('.')[0]
        rs_id = data[18]
        qvalue = float(data[28])
        pvalue = float(data[23])
        all_genes[ensamble_id] = 1
        if qvalue <= .05:
            egenes[ensamble_id] = 1
    return egenes, all_genes

def extract_gtex_tissue_top_genes(egene_file):
    f = open(egene_file)
    head_count = 0
    all_genes = {}
    egenes = {}
    genes = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[0].split('.')[0]
        rs_id = data[18]
        qvalue = float(data[28])
        pvalue = float(data[23])
        all_genes[ensamble_id] = 1
        if ensamble_id not in genes:
            genes[ensamble_id] = pvalue
        else:
            genes[ensamble_id] = min(genes[ensamble_id], pvalue)
    f.close()
    gene_arr = []
    for geney in genes.keys():
        gene_arr.append((geney,genes[geney]))
    gene_arr.sort(key=lambda x: x[1])
    for i, tupler in enumerate(gene_arr):
        if i > (num_genes-1):
            continue
        egenes[tupler[0]] = 0
    return egenes, all_genes


parameter_string = sys.argv[1]
cht_output_dir = sys.argv[2]
cht_visualization_dir = sys.argv[3]
target_regions_dir = sys.argv[4]
gencode_gene_annotation_file = sys.argv[5]
ipsc_egene_file = sys.argv[6]
ipsc_eqtl_file = sys.argv[7]
cm_egene_file = sys.argv[8]
cm_eqtl_file = sys.argv[9]
eqtl_data_set_file = sys.argv[10]

#############################################################
# Gene set enrichment for GTEx v7
#############################################################
# Loop through desired gtex tissues
f = open(eqtl_data_set_file)
for line in f:
    line = line.rstrip()
    data = line.split()
    if data[2] != 'gtex_egene':
        continue
    egene_file = data[1]
    tissue_name = egene_file.split('/')[-1].split('.')[0]
    print(tissue_name)
    # Do gene set enrichment based on egenes (fdr <= .05 in this tissue.. done with qvalue)
    tissue_egenes, background_genes = extract_gtex_tissue_significant_genes(egene_file)
    output_root = cht_visualization_dir + 'gtex_v7_gene_set_enrichment_significant_genes_05_fdr_' + tissue_name + '_'
    gene_set_enrichment(tissue_egenes, background_genes, output_root, gencode_gene_annotation_file)
    num_genes_arr = [50, 100, 200, 300, 400, 500]
    for num_genes in num_genes_arr:
        tissue_top_genes, background_genes = extract_gtex_tissue_top_genes(egene_file)
        output_root = cht_visualization_dir + 'gtex_v7_gene_set_enrichment_top_' + str(num_genes) + '_genes_' + tissue_name + '_'
        gene_set_enrichment(tissue_top_genes, background_genes, output_root, gencode_gene_annotation_file)


f.close()
'''



#############################################################
# Gene set enrichment for Banovich ipsc and ipsc-cm eqtls
#############################################################
# ipsc-cm eqtls
cm_all_genes = extract_cm_all_genes(cm_eqtl_file)
num_genes_arr = [50, 100, 200, 300, 400, 500]
for num_genes in num_genes_arr:
    print(num_genes)
    top_n_cm_egenes = extract_top_n_cm_egenes(cm_eqtl_file, num_genes)
    output_root = cht_visualization_dir + 'banovich_ipsc_cm_gene_set_enrichment_top_' + str(num_genes) + '_'
    gene_set_enrichment(top_n_cm_egenes, cm_all_genes, output_root, gencode_gene_annotation_file)
    




# ipsc eqtls
ipsc_egenes = extract_ipsc_egenes(ipsc_egene_file)
ipsc_all_genes = extract_ipsc_all_genes(ipsc_eqtl_file)
output_root = cht_visualization_dir + 'banovich_ipsc_gene_set_enrichment_all_sig_'
gene_set_enrichment(ipsc_egenes, ipsc_all_genes, output_root, gencode_gene_annotation_file)
# Do ipsc eqtls (only selecting top n genes)
num_genes_arr = [100, 200, 300, 400, 500]
for num_genes in num_genes_arr:
    print(num_genes)
    ipsc_egenes_top_n = extract_top_n_ipsc_egenes(ipsc_egene_file, num_genes)
    output_root = cht_visualization_dir + 'banovich_ipsc_gene_set_enrichment_top_' + str(num_genes) + '_'
    gene_set_enrichment(ipsc_egenes_top_n, ipsc_all_genes, output_root, gencode_gene_annotation_file)



#############################################################
# Gene set enrichment for per-time step eqtls
#############################################################
num_pc = '3'
efdr = '.05'


for time_step in range(16):
	print(time_step)
	egene_file = cht_output_dir + 'cht_results_' + parameter_string + '_num_pc_' + num_pc + '_time_' + str(time_step) + '_efdr_thresh_' + efdr + '_significant_egenes.txt'
	time_step_egenes = extract_genes(egene_file)
	all_tests_file = cht_output_dir + 'cht_results_' + parameter_string + '_num_pc_' + num_pc + '_time_' + str(time_step) + '_eqtl_results.txt'
	background_genes = extract_genes(all_tests_file)


	output_root = cht_visualization_dir + parameter_string + '_gene_set_enrichment_per_time_' + str(time_step) + '_'

	gene_set_enrichment(time_step_egenes, background_genes, output_root, gencode_gene_annotation_file)
'''
