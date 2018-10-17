import numpy as np
import os
import sys
import pdb
import random
import math

# Return the bin number corresponding to this distance
def get_distance_bin(distance, distance_bin_size):
    return int(math.floor(distance/distance_bin_size))


# Return the bin number corresponding to this distance
def get_maf_bin(maf, maf_bin_size):
    return int(math.floor(maf/maf_bin_size))

# Make object that input you can input distance and maf into, and it will return an array of all the pvalues contained in that bin
def make_background_qtls_object(our_eqtl_file, distance_bin_size, maf_bin_size, eqtl_distance):
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

    ####################
    # Add pvalues from our_eqtl_file to the appropriate bin
    ####################
    f = open(our_eqtl_file)
    head_count = 0  # Used to skip header
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1  # Skip Header
            continue
        # Extract relevent fields from column delimited file
        gene_position = float(data[2])
        variant_position = float(data[4])
        distance = abs(gene_position-variant_position)
        maf = float(data[5])
        pvalue = float(data[8])
        rs_id = data[3]
        if maf > .5 or maf < .1:
            print('maf errro')
            pdb.set_trace()

        # Check to make sure we are only dealing with labeled rs_ids
        # This should have been taken care of in `prepare_independent_time_step_matrix_eqtl_files.py`. But just checking here
        if rs_id == '.':
            print('ERROR')
            continue

        # Return the bin number corresponding to this distance
        distance_bin = get_distance_bin(distance, distance_bin_size)
        # Return the bin number corresponding to this distance
        maf_bin = get_maf_bin(maf, maf_bin_size)

        # Add pvalue of this variant gene pair to the appropriate bin
        background_qtls[distance_bin][maf_bin].append(pvalue)
    return background_qtls

# Extract dictionary list of variant_geneID pairs that are found to be eqtls in reference_eqtl_file
def extract_reference_eqtls(reference_eqtl_file,version):
    eqtls = {}
    # Stream file
    f = open(reference_eqtl_file)
    head_count = 0  # used to skip header
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # Skip header
            head_count = head_count + 1
            continue
        # Extract gene id and variant id in different ways depending on where the data comes from
        if version == 'ipsc':
            # Extract relevent fields from line
            gene_id = data[0].split('.')[0]
            rsid = data[1].split('.')[0]
        elif version == 'gtex_egene':
            gene_id = data[0].split('.')[0]
            rsid = data[18]
            qvalue = float(data[-2])
            if qvalue > .05:
                continue
        elif version == 'ipsc_cardiac':
            gene_id = data[0]
            rsid = data[1]
        elif version == 'gtex_egene_beta_filter':
            gene_id = data[0].split('.')[0]
            rsid = data[18]
            beta = float(data[24])
            qvalue = float(data[-2])
            if abs(beta) < .4 or qvalue > .05:
                continue
        # Concatenate geneID and rsID to get name of test
        test_name = gene_id + '_' + rsid
        # Add eqtl to dictionary
        eqtls[test_name] = 1
    return eqtls

# Stream our_eqtl_file and find those that are in our_eqtl_file
# For each one of those hits, randomly select background pvalue
def extract_and_print_pvalues(background_qtls, reference_eqtls, our_eqtl_file, output_file):
    # Initialize output handle
    t = open(output_file, 'w')
    t.write('real_pvalue\tmatched_pvalue\n')

    # Stream our_eqtl_file
    head_count = 0
    f = open(our_eqtl_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1  # Skip Header
            continue
        # Extract relevent fields
        pvalue = float(data[8])
        gene_id = data[1]
        rs_id = data[3]
        test_id = gene_id + '_' + rs_id
        gene_position = float(data[2])
        variant_position = float(data[4])
        distance = abs(gene_position-variant_position)
        maf = float(data[5])
        # Check to see if test_id in reference_eqtls
        if test_id not in reference_eqtls:
            continue
        # Now randomly select a background variant-gene pair's pvalue for reference
        # Return the bin number corresponding to this distance
        distance_bin = get_distance_bin(distance, distance_bin_size)
        # Return the bin number corresponding to this distance
        maf_bin = get_maf_bin(maf, maf_bin_size)
        
        # Extract list of pvalue corresponding to this distance bin and maf_bin
        pvalues = background_qtls[distance_bin][maf_bin]

        if len(pvalues) < 5:
            print('small bin error')
            print(len(pvalues))
            print(distance_bin)
            print(maf_bin)
            print(distance)
            print(maf)
        # Randomly select one element from the list
        random_pvalue = random.choice(pvalues)

        # Write to output file
        t.write(str(pvalue) + '\t' + str(random_pvalue) + '\n')

def remove_NA(vector):
    new_vector = []
    for ele in vector:
        if ele == 'NA':
            new_vector.append('NaN')
        else:
            new_vector.append(ele)
    return new_vector

def extract_mvalue_tissue_specificities(mvalue_file):
    tissues_of_interest = [('breast',19),('hlv',28),('stomach',42),('skin_not_sun',38)]
    f = open(mvalue_file)
    head_count = 0
    eqtl_tissue_specificities = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1 
            continue
        dicti = {}

        mvalz = np.asarray(remove_NA(data[-48:])).astype(float)
        gene_id = data[0].split(',')[1]
        snp_id = data[0].split(',')[0]
        for tissue_of_interest in tissues_of_interest:
            tissue_name = tissue_of_interest[0]
            tissue_index = tissue_of_interest[1]
            # Check to see if the snp-gene pair is tissue specific for tissue_of_interest
            if mvalz[tissue_index] > .8 and len(np.where(mvalz > .8)[0]) < 20:
                dicti[tissue_name] = 1
        eqtl_tissue_specificities[gene_id + ',' + snp_id] = dicti
    return eqtl_tissue_specificities
###############################################################################################
# Objective
###############################################################################################
# We are given significant eqtls from another eqtl analysis.
# We extract both:
#  1. the p-values of those eqtls in our data.
#  2. p-values of matched snp-gene pairs (matched for distance to tss and maf)


###############################################################################################
# Input Data
###############################################################################################
# eqtl results from our analysis
our_eqtl_file = sys.argv[1]
# Output file stem
output_stem = sys.argv[2]
# File where each line is a eqtl data set (with information on that data set)
eqtl_data_set_file = sys.argv[3]
# Max distance a variant and gene can be from one another
eqtl_distance = float(sys.argv[4])
# File containing mvalue data for gtex v7
mvalue_file = sys.argv[5]


# eqtl_tissue_specificities = extract_mvalue_tissue_specificities(mvalue_file)

# Size of bins used for background matching
distance_bin_size = 14000
maf_bin_size = .08

# Create file handle for our data set file
f = open(eqtl_data_set_file)
# Loop through data set file
# Ie. each line is a data set
for line in f:
    # Extract information regarding the data set
    line = line.rstrip()
    data = line.split()
    data_set_name = data[0]
    reference_eqtl_file = data[1]
    version = data[2]

    # Name of output file for this data set
    output_file = output_stem + '_' + data_set_name + '_real_v_matched_controls.txt'
    
    # Make object that input you can input distance and maf into, and it will return an array of all the pvalues contained in that bin
    background_qtls = make_background_qtls_object(our_eqtl_file, distance_bin_size, maf_bin_size, eqtl_distance)


    # Extract dictionary list of geneID_variantID pairs that are found to be eqtls in reference_eqtl_file
    reference_eqtls = extract_reference_eqtls(reference_eqtl_file,version)


    # Stream our_eqtl_file and find those that are in our_eqtl_file
    # For each one of those hits, randomly select background pvalue
    extract_and_print_pvalues(background_qtls, reference_eqtls, our_eqtl_file, output_file)