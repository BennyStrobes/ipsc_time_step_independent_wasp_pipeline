import numpy as np 
import os
import sys
import pdb
import math
import gzip
import random
import scipy.stats



# Extract list of cell line ids used for this cell_line_version
def get_cell_line_ids(cell_line_version, chrom_hmm_input_dir):
    cell_line_ids = []
    if cell_line_version == 'heart_cell_lines':
        cell_line_ids.append('E095')
        #cell_line_ids.append('E104')
        cell_line_ids.append('E105')
        cell_line_ids.append('E083')
    elif cell_line_version == 'ipsc_cell_lines':
        cell_line_ids.append('E018')
        cell_line_ids.append('E019')
        cell_line_ids.append('E020')
        cell_line_ids.append('E021')
        cell_line_ids.append('E022')

    elif cell_line_version == 'all_cell_lines':
        for file_name in os.listdir(chrom_hmm_input_dir):
            if file_name.endswith('mnemonics.bed.gz') == False:
                continue
            cell_line_id = file_name.split('_')[0]
            cell_line_ids.append(cell_line_id)
    return np.asarray(cell_line_ids)



# Create Mapping from variant-gene pair to quartuple (chrom_num, variant_position, distToTss, MAF)
def extract_variant_gene_pair_info(time_step_independent_file):
    head_count = 0
    # Open time_step_independent file that contains distance to tss knowledge for each variant gene pair
    f = open(time_step_independent_file)
    # Dictionary for each variant-gene pair
    variant_gene_pair_info = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        # Extract info
        ensamble_id = data[1]
        rs_id = data[3]
        pair_name = rs_id + '_' + ensamble_id
        # compute distance to tss
        gene_tss_pos = float(data[2])
        rs_pos = float(data[4])
        disty = abs(rs_pos - gene_tss_pos)
        maf = float(data[5])
        # Add relevent info to small dictionary
        mini_dicti = {}
        mini_dicti['chrom_num'] = int(data[0])
        mini_dicti['dist_to_tss'] = disty
        mini_dicti['maf'] = maf
        mini_dicti['variant_position'] = rs_pos
        # Add variant-gene pair to dictionary
        if pair_name in variant_gene_pair_info:
            print('assumption errror')
            pdb.set_trace()
        variant_gene_pair_info[pair_name] = mini_dicti
    return variant_gene_pair_info


def extract_significant_variant_gene_pairs(egenes_file):
    f = open(egenes_file)
    dicti = {}  # dictionary to keep variant gene pairs
    head_count = 0  # Used to skip header
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # Skip header
            head_count = head_count + 1
            continue
        #  A significant variant-gene pairs
        rs_id = data[3]
        ensamble_id = data[1]
        # Simple check
        if rs_id + '_' + ensamble_id in dicti:
            print('fundamental assumption error')
            pdb.set_trace()
        # Add variant gene pair to dictionary
        dicti[rs_id + '_' + ensamble_id] = 1
    return dicti

# Return the bin number corresponding to this distance
def get_distance_bin(distance, distance_bin_size):
    return int(math.floor(distance/distance_bin_size))


# Return the bin number corresponding to this distance
def get_maf_bin(maf, maf_bin_size):
    return int(math.floor(maf/maf_bin_size))

def make_background_object(time_step_independent_file):
    distance_bin_size = 10000
    maf_bin_size = .05
    eqtl_distance = 50000
    ####################
    # Initialize object
    ####################
    background_qtls = []
    # number of bins needed for maf and distance
    num_distance_bins = int(math.ceil(eqtl_distance/distance_bin_size + 1))
    num_maf_bins = int(math.ceil(.5/maf_bin_size + 1))
    # Add each possible bin
    for distance_bin in range(num_distance_bins):
        background_qtls.append([])
        for maf_bin in range(num_maf_bins):
            background_qtls[distance_bin].append([])
    f = open(time_step_independent_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[1]
        rs_id = data[3]
        test_name = rs_id + '_' + ensamble_id
        distance = abs(float(data[2]) - float(data[4]))
        maf = float(data[5])
        # Return the bin number corresponding to this distance
        distance_bin = get_distance_bin(distance, distance_bin_size)
        # Return the bin number corresponding to this distance
        maf_bin = get_maf_bin(maf, maf_bin_size)
        background_qtls[distance_bin][maf_bin].append(test_name)
    f.close()
    return background_qtls


def get_valid_markers(marker_type):
    valid_markers = {}
    if marker_type == 'promotor':
        valid_markers['1_TssA'] = 1
        valid_markers['2_TssAFlnk'] = 1
        valid_markers['10_TssBiv'] = 1
        valid_markers['11_BivFlnk'] = 1
    elif marker_type == 'enhancer':
        valid_markers['7_Enh'] = 1
        valid_markers['6_EnhG'] = 1
        valid_markers['12_EnhBiv'] = 1
        valid_markers['11_BivFlnk'] = 1
    return valid_markers

# Count number of variants that overlap a marker on this chromosome
def count_variant_overlap(chrom_num, chromosome, sig_variant_gene_pairs, variant_gene_pair_info):
    county = 0
    for variant_gene in sig_variant_gene_pairs.keys():
        variant_gene_pair_dict = variant_gene_pair_info[variant_gene]
        if variant_gene_pair_dict['chrom_num'] != chrom_num:
            continue
        variant_position = int(variant_gene_pair_dict['variant_position'])
        county = county + chromosome[variant_position]
    return county

def count_variant_overlap_specificity(chrom_num, heart_chromosome, ipsc_chromosome, sig_variant_gene_pairs, variant_gene_pair_info, cell_line_version):
    county = 0
    for variant_gene in sig_variant_gene_pairs.keys():
        variant_gene_pair_dict = variant_gene_pair_info[variant_gene]
        if variant_gene_pair_dict['chrom_num'] != chrom_num:
            continue
        variant_position = int(variant_gene_pair_dict['variant_position'])
        if cell_line_version == 'heart_and_ipsc_cell_lines':
            if heart_chromosome[variant_position] == 1.0 and ipsc_chromosome[variant_position] == 1.0:
                county = county + 1
        elif cell_line_version == 'heart_only_cell_lines':
            if heart_chromosome[variant_position] == 1.0 and ipsc_chromosome[variant_position] == 0.0:
                county = county + 1
        elif cell_line_version == 'ipsc_only_cell_lines':
            if heart_chromosome[variant_position] == 0.0 and ipsc_chromosome[variant_position] == 1.0:
                county = county + 1
    return county



# Make binary array length of a chromosome. If array == 0, no marker there. If array == 1, there is a marker there
def make_binary_chromosome(chrom_num, chrom_hmm_input_dir, cell_line_ids, marker_type):
    chrom_num_string = 'chr' + str(chrom_num)
    chromosome = np.zeros(259250621)
    valid_markers = get_valid_markers(marker_type)
    for cell_line_id in cell_line_ids:
        chrom_hmm_file = chrom_hmm_input_dir + cell_line_id + '_15_coreMarks_mnemonics.bed.gz'
        f = gzip.open(chrom_hmm_file)
        for line in f:
            line = line.rstrip()
            data = line.split()
            chromer = data[0]
            # Only consider marks on this chromsome
            if chromer != chrom_num_string:
                continue
            line_marker = data[3]
            if line_marker in valid_markers:
                start = int(data[1])
                end = int(data[2])
                if start > end:
                    print('assumption errororo')
                chromosome[start:end] = np.ones(end-start)
    return chromosome

def sample_background_variant_gene_pairs(background_object, variant_gene_pair_info, sig_variant_gene_pairs):
    distance_bin_size = 10000
    maf_bin_size = .05
    eqtl_distance = 50000
    background_variant_gene_pairs = {}
    for test_name in sig_variant_gene_pairs.keys():
        dist_to_tss = variant_gene_pair_info[test_name]['dist_to_tss']
        maf = variant_gene_pair_info[test_name]['maf']
        # Return the bin number corresponding to this distance
        distance_bin = get_distance_bin(dist_to_tss, distance_bin_size)
        # Return the bin number corresponding to this distance
        maf_bin = get_maf_bin(maf, maf_bin_size)
        converged = False
        while converged == False:
            if len(background_object[distance_bin][maf_bin]) < 10:
                print('small backgroudn: should investigate')
                pdb.set_trace()
            randomly_selected_pair = random.choice(background_object[distance_bin][maf_bin])
            if randomly_selected_pair not in background_variant_gene_pairs:
                background_variant_gene_pairs[randomly_selected_pair] = 1
                converged = True
    if len(background_variant_gene_pairs) != len(sig_variant_gene_pairs):
        print('ASSUMPTION EROROR')
    return background_variant_gene_pairs

# Return list of length num_permutations where each element of the dictionary list of variant-gene pairs (of len(sig_variant_gene_pairs)) matched for dist_to_tss and maf
def extract_perm_variant_gene_pairs(sig_variant_gene_pairs, time_step_independent_file, num_permutations, variant_gene_pair_info):
    wrapper_list = []
    background_object = make_background_object(time_step_independent_file)
    for perm_num in range(num_permutations):
        background_variant_gene_pairs = sample_background_variant_gene_pairs(background_object, variant_gene_pair_info, sig_variant_gene_pairs)
        wrapper_list.append(background_variant_gene_pairs)
    return wrapper_list


parameter_string = sys.argv[1]
cht_output_dir = sys.argv[2]
cht_visualization_dir = sys.argv[3]
pc_num = sys.argv[4]
chrom_hmm_input_dir = sys.argv[5]
cht_enrichment_dir = sys.argv[6]
num_permutations = int(sys.argv[7])
time_step = sys.argv[8]
marker_type = sys.argv[9]
cell_line_version = sys.argv[10]
efdr = sys.argv[11]


# Extract list of cell line ids used for this cell_line_version
cell_line_ids = get_cell_line_ids(cell_line_version, chrom_hmm_input_dir)

# file containing all variant gene pairs for this time step, and their corresponding pvalue
all_hits_file = cht_output_dir + 'cht_results_' + parameter_string + '_num_pc_' + pc_num + '_time_' + time_step + '_eqtl_results.txt'

# file containing all significant genes (eFDR <= $efdr) for this time step, and their strongest associated variant
egenes_file = cht_output_dir + 'cht_results_' + parameter_string + '_num_pc_' + pc_num + '_time_' + time_step + '_efdr_thresh_' + efdr + '_significant_egenes.txt'


# Create Mapping from variant-gene pair to quartuple (chrom_num, variant_position, distToTss, MAF)
variant_gene_pair_info = extract_variant_gene_pair_info(all_hits_file)

# First create dictionary list of the significant variant gene pairs where each key is of form $variantID"_"$geneID
sig_variant_gene_pairs = extract_significant_variant_gene_pairs(egenes_file)
print(len(sig_variant_gene_pairs))

# Return list of length num_permutations where each element of the dictionary list of variant-gene pairs (of len(sig_variant_gene_pairs)) matched for dist_to_tss and maf
perm_variant_gene_pairs = extract_perm_variant_gene_pairs(sig_variant_gene_pairs, all_hits_file, num_permutations, variant_gene_pair_info)

#############################
# Count up number of times a variant overlaps one of our markers
###############################
real_overlaps = 0  # Initialze real counts
perm_overlaps = np.zeros(num_permutations) # initialize vector of perm counts

# Do seperately for each chromosome
for chrom_num in range(1,23):
    print(chrom_num)
    if cell_line_version == 'heart_cell_lines' or cell_line_version == 'ipsc_cell_lines' or cell_line_version == 'all_cell_lines':
        # Make binary array length of a chromosome. If array == 0, no marker there. If array == 1, there is a marker there
        chromosome = make_binary_chromosome(chrom_num, chrom_hmm_input_dir, cell_line_ids, marker_type)

        # Count number of variants that overlap a marker on this chromosome
        real_overlaps = real_overlaps + count_variant_overlap(chrom_num, chromosome, sig_variant_gene_pairs, variant_gene_pair_info)
     
        # Count number of variants that overlap a marker on this chromosome
        for perm_num in range(num_permutations):
            perm_overlaps[perm_num] = perm_overlaps[perm_num] + count_variant_overlap(chrom_num, chromosome, perm_variant_gene_pairs[perm_num], variant_gene_pair_info)
    elif cell_line_version == 'heart_and_ipsc_cell_lines' or cell_line_version == 'heart_only_cell_lines' or cell_line_version == 'ipsc_only_cell_lines':
        heart_cell_line_ids = get_cell_line_ids('heart_cell_lines', chrom_hmm_input_dir)
        ipsc_cell_line_ids = get_cell_line_ids('ipsc_cell_lines', chrom_hmm_input_dir)
        # Make binary array length of a chromosome. If array == 0, no marker there. If array == 1, there is a marker there
        heart_chromosome = make_binary_chromosome(chrom_num, chrom_hmm_input_dir, heart_cell_line_ids, marker_type)
        ipsc_chromosome = make_binary_chromosome(chrom_num, chrom_hmm_input_dir, ipsc_cell_line_ids, marker_type)


        # Count number of variants that overlap a marker on this chromosome
        real_overlaps = real_overlaps + count_variant_overlap_specificity(chrom_num, heart_chromosome, ipsc_chromosome, sig_variant_gene_pairs, variant_gene_pair_info, cell_line_version)

        # Count number of variants that overlap a marker on this chromosome
        for perm_num in range(num_permutations):
            perm_overlaps[perm_num] = perm_overlaps[perm_num] + count_variant_overlap_specificity(chrom_num, heart_chromosome, ipsc_chromosome, perm_variant_gene_pairs[perm_num], variant_gene_pair_info, cell_line_version)

# For each permutation run, compute fischers exact test and print to output file

output_file = cht_enrichment_dir + parameter_string + '_num_pc_' + pc_num + '_time_' + time_step + '_efdr_thresh_' + efdr + '_' + cell_line_version + '_' + marker_type + '_' + str(num_permutations) + '_enrich.txt'
t = open(output_file, 'w')
num_samp = len(sig_variant_gene_pairs)
t.write('real_overlaps\treal_misses\tperm_overlaps\tperm_misses\todds_ratio\tpvalue\n')
for perm_num in range(num_permutations):
    aa = real_overlaps
    bb = num_samp - real_overlaps
    cc = perm_overlaps[perm_num]
    dd = num_samp - cc
    odds_ratio, pvalue = scipy.stats.fisher_exact([[aa,bb],[cc,dd]])
    t.write(str(aa) + '\t' + str(bb) + '\t' + str(cc) + '\t' + str(dd) + '\t' + str(odds_ratio) + '\t' + str(pvalue) + '\n')
t.close()

