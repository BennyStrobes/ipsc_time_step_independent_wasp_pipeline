import numpy as np 
import os
import sys
import pdb
import math
import gzip
import random
import scipy.stats
from sklearn.cluster import KMeans
from sklearn import mixture




def extract_variant_gene_pair_time_step_info(time_step_independent_stem):
    # Dictionary to keep track of mapping from variant-gene pairs to vector of length time-steps taht contains pvalues
    dicti = {}
    dicti_pval = {}
    # loop through all time steps
    for time_step in range(16):
        file_name = time_step_independent_stem + str(time_step) + '_eqtl_results.txt'
        head_count = 0  # Used to skip header
        f = open(file_name)
        for line in f:
            line = line.rstrip()
            data = line.split()
            if head_count == 0:
                head_count = head_count + 1
                continue
            # Extract relevent info
            rs_id = data[3]
            ensamble_id = data[1]
            pvalue = float(data[-1])
            alpha = float(data[-3])
            beta = float(data[-2])
            allelic_fraction = alpha/(alpha+beta)
            test_name = rs_id + '_' + ensamble_id
            # If we've never seen this test before
            if test_name not in dicti:
                if time_step != 0:
                    print('assumpriton erroro')
                    pdb.set_trace()
                dicti[test_name] = np.zeros(16)
            dicti[test_name][time_step] = allelic_fraction
            # If we've never seen this test before
            if test_name not in dicti_pval:
                if time_step != 0:
                    print('assumpriton erroro')
                    pdb.set_trace()
                dicti_pval[test_name] = np.zeros(16)
            dicti_pval[test_name][time_step] = pvalue
        f.close()
    return dicti, dicti_pval


def in_same_direction(estimated_t0, estimated_t15, error_bound):
    # Both going in same direction
    if estimated_t0 <= .5 and estimated_t15 <= .5:
        return True
    elif estimated_t0 > .5 and estimated_t15 > .5:
        return True 
    # Different directions
    diff_t0 = abs(estimated_t0 - .5)
    diff_t15 = abs(estimated_t15 - .5)
    if diff_t0 >= diff_t15:
        if diff_t15 <= error_bound:
            return True 
    elif diff_t15 > diff_t0:
        if diff_t0 <= error_bound:
            return True
    return False


def linear_regression_allelic_fraction_classification(dicti, error_bound, output_file):
    t = open(output_file, 'w')
    t.write('variant_id\tgene_id\tcluster\n')
    for test_name in dicti.keys():
        rs_id = test_name.split('_')[0]
        gene_id = test_name.split('_')[1]
        allelic_vec = dicti[test_name]
        time_steps = np.arange(16)
        result = scipy.stats.linregress(time_steps,allelic_vec)
        slope = result[0]
        intercept = result[1]
        estimated_t0 = intercept
        estimated_t15 = intercept + (slope*15.0)
        #estimated_t0 = allelic_vec[0]
        #estimated_t15 = allelic_vec[15]
        if abs(estimated_t15 - .5) > abs(estimated_t0 - .5) and in_same_direction(estimated_t0, estimated_t15, error_bound):
            cluster = 'late_time_step_hits'
        elif abs(estimated_t15 - .5) <= abs(estimated_t0 - .5) and in_same_direction(estimated_t0, estimated_t15, error_bound):
            cluster = 'early_time_step_hits'
        else:
            cluster = 'change_in_sign_hits'
        t.write(rs_id + '\t' + gene_id + '\t' + cluster + '\n')
    t.close()


def convert_from_dictionary_to_matrix(allelic_fraction_dict,pval_dict):
    test_names = []
    matrix = []
    for test_name in allelic_fraction_dict.keys():
        allelic_vec = allelic_fraction_dict[test_name]
        pval_vec = pval_dict[test_name]
        #if min(pval_vec) > .05:
        #    continue
        allelic_vec = np.abs(allelic_vec - .5)
        allelic_vec = allelic_vec - np.mean(allelic_vec)
        matrix.append(allelic_vec)
        test_names.append(test_name)
    return np.asarray(test_names), np.asmatrix(matrix)

def convert_from_dictionary_to_matrix_pvalue(allelic_fraction_dict,pval_dict, pvalue_thresh):
    test_names = []
    matrix = []
    for test_name in allelic_fraction_dict.keys():
        allelic_vec = allelic_fraction_dict[test_name]
        pval_vec = pval_dict[test_name]
        if min(pval_vec) > pvalue_thresh:
            continue
        allelic_vec = np.abs(allelic_vec - .5)
        allelic_vec = allelic_vec - np.mean(allelic_vec)
        matrix.append(allelic_vec)
        test_names.append(test_name)
    return np.asarray(test_names), np.asmatrix(matrix)

def convert_from_dictionary_to_matrix_standardized(allelic_fraction_dict,pval_dict):
    test_names = []
    matrix = []
    for test_name in allelic_fraction_dict.keys():
        allelic_vec = allelic_fraction_dict[test_name]
        pval_vec = pval_dict[test_name]
        #if min(pval_vec) > .05:
        #    continue
        allelic_vec = np.abs(allelic_vec - .5)
        allelic_vec = (allelic_vec - np.mean(allelic_vec))/np.std(allelic_vec)
        matrix.append(allelic_vec)
        test_names.append(test_name)
    return np.asarray(test_names), np.asmatrix(matrix)


def convert_from_dictionary_to_matrix_standardized_pvalue(allelic_fraction_dict,pval_dict,pvalue_thresh):
    test_names = []
    matrix = []
    for test_name in allelic_fraction_dict.keys():
        allelic_vec = allelic_fraction_dict[test_name]
        pval_vec = pval_dict[test_name]
        if min(pval_vec) > pvalue_thresh:
            continue
        allelic_vec = np.abs(allelic_vec - .5)
        allelic_vec = (allelic_vec - np.mean(allelic_vec))/np.std(allelic_vec)
        matrix.append(allelic_vec)
        test_names.append(test_name)
    return np.asarray(test_names), np.asmatrix(matrix)

def unsupervised_clustering_kmeans(allelic_fraction_dict, pval_dict, output_file):
    test_names, matrix = convert_from_dictionary_to_matrix(allelic_fraction_dict, pval_dict)
    kmeans = KMeans(n_clusters=3).fit(matrix)
    cluster1_center = kmeans.cluster_centers_[0,:]
    cluster2_center = kmeans.cluster_centers_[1,:]
    cluster3_center = kmeans.cluster_centers_[2,:]
    time_steps = np.arange(16)
    
    cluster1_slope = scipy.stats.linregress(time_steps,cluster1_center)[0]
    cluster2_slope = scipy.stats.linregress(time_steps,cluster2_center)[0]
    cluster3_slope = scipy.stats.linregress(time_steps,cluster3_center)[0]

    # To convert from cluster number to cluster name
    converter = {}
    slopes = [cluster1_slope, cluster2_slope, cluster3_slope]
    # Get index of cluster that decreases over time (ie early time step qtls)
    index = np.where(slopes == min(slopes))[0][0]
    converter[index] = 'early_time_step_hits'
    # Get index of cluster that increases over time (ie late time step qtls)
    index = np.where(slopes == max(slopes))[0][0]
    converter[index] = 'late_time_step_hits'
    # Other
    index = np.where(slopes == np.median(slopes))[0][0]
    converter[index] = 'change_in_sign_hits'

    cluster_assignments = kmeans.labels_
    t = open(output_file, 'w')
    t.write('variant_id\tgene_id\tcluster\n')
    for i,test_name in enumerate(test_names):
        rs_id = test_name.split('_')[0]
        gene_id = test_name.split('_')[1]
        cluster_name = converter[cluster_assignments[i]]
        t.write(rs_id + '\t' + gene_id + '\t' + cluster_name + '\n')
    t.close()

def unsupervised_clustering_kmeans_pvalue_thresh(allelic_fraction_dict, pval_dict, output_file,pvalue_thresh):
    test_names, matrix = convert_from_dictionary_to_matrix_pvalue(allelic_fraction_dict, pval_dict,pvalue_thresh)
    kmeans = KMeans(n_clusters=3).fit(matrix)
    cluster1_center = kmeans.cluster_centers_[0,:]
    cluster2_center = kmeans.cluster_centers_[1,:]
    cluster3_center = kmeans.cluster_centers_[2,:]
    time_steps = np.arange(16)
    
    cluster1_slope = scipy.stats.linregress(time_steps,cluster1_center)[0]
    cluster2_slope = scipy.stats.linregress(time_steps,cluster2_center)[0]
    cluster3_slope = scipy.stats.linregress(time_steps,cluster3_center)[0]

    # To convert from cluster number to cluster name
    converter = {}
    slopes = [cluster1_slope, cluster2_slope, cluster3_slope]
    # Get index of cluster that decreases over time (ie early time step qtls)
    index = np.where(slopes == min(slopes))[0][0]
    converter[index] = 'early_time_step_hits'
    # Get index of cluster that increases over time (ie late time step qtls)
    index = np.where(slopes == max(slopes))[0][0]
    converter[index] = 'late_time_step_hits'
    # Other
    index = np.where(slopes == np.median(slopes))[0][0]
    converter[index] = 'change_in_sign_hits'

    cluster_assignments = kmeans.labels_
    t = open(output_file, 'w')
    t.write('variant_id\tgene_id\tcluster\n')
    used = {}
    for i,test_name in enumerate(test_names):
        rs_id = test_name.split('_')[0]
        used[test_name] = 1
        gene_id = test_name.split('_')[1]
        cluster_name = converter[cluster_assignments[i]]
        t.write(rs_id + '\t' + gene_id + '\t' + cluster_name + '\n')
    for test_name in pval_dict.keys():
        if test_name not in used:
            t.write(test_name.split('_')[0] + '\t' + test_name.split('_')[1] + '\t' + 'change_in_sign_hits' + '\n')
    t.close()

def unsupervised_clustering_full_gmm(allelic_fraction_dict, pval_dict, output_file):
    test_names, matrix = convert_from_dictionary_to_matrix_standardized(allelic_fraction_dict, pval_dict)
  
    g = mixture.GMM(n_components=3,covariance_type='full')

    fit_gmm = g.fit(matrix)

    cluster1_center = fit_gmm.means_[0,:]
    cluster2_center = fit_gmm.means_[1,:]
    cluster3_center = fit_gmm.means_[2,:]
    time_steps = np.arange(16)
    
    cluster1_slope = scipy.stats.linregress(time_steps,cluster1_center)[0]
    cluster2_slope = scipy.stats.linregress(time_steps,cluster2_center)[0]
    cluster3_slope = scipy.stats.linregress(time_steps,cluster3_center)[0]

    # To convert from cluster number to cluster name
    converter = {}
    slopes = [cluster1_slope, cluster2_slope, cluster3_slope]
    # Get index of cluster that decreases over time (ie early time step qtls)
    index = np.where(slopes == min(slopes))[0][0]
    converter[index] = 'early_time_step_hits'
    # Get index of cluster that increases over time (ie late time step qtls)
    index = np.where(slopes == max(slopes))[0][0]
    converter[index] = 'late_time_step_hits'
    # Other
    index = np.where(slopes == np.median(slopes))[0][0]
    converter[index] = 'change_in_sign_hits'

    cluster_assignments = fit_gmm.predict(matrix)
    t = open(output_file, 'w')
    t.write('variant_id\tgene_id\tcluster\n')
    for i,test_name in enumerate(test_names):
        rs_id = test_name.split('_')[0]
        gene_id = test_name.split('_')[1]
        cluster_name = converter[cluster_assignments[i]]
        t.write(rs_id + '\t' + gene_id + '\t' + cluster_name + '\n')
    t.close()

def unsupervised_clustering_full_gmm_pvalue_thresh(allelic_fraction_dict, pval_dict, output_file, pvalue_thresh):
    test_names, matrix = convert_from_dictionary_to_matrix_standardized_pvalue(allelic_fraction_dict, pval_dict, pvalue_thresh)
  
    g = mixture.GMM(n_components=3,covariance_type='full')

    fit_gmm = g.fit(matrix)

    cluster1_center = fit_gmm.means_[0,:]
    cluster2_center = fit_gmm.means_[1,:]
    cluster3_center = fit_gmm.means_[2,:]
    time_steps = np.arange(16)
    
    cluster1_slope = scipy.stats.linregress(time_steps,cluster1_center)[0]
    cluster2_slope = scipy.stats.linregress(time_steps,cluster2_center)[0]
    cluster3_slope = scipy.stats.linregress(time_steps,cluster3_center)[0]

    # To convert from cluster number to cluster name
    converter = {}
    slopes = [cluster1_slope, cluster2_slope, cluster3_slope]
    # Get index of cluster that decreases over time (ie early time step qtls)
    index = np.where(slopes == min(slopes))[0][0]
    converter[index] = 'early_time_step_hits'
    # Get index of cluster that increases over time (ie late time step qtls)
    index = np.where(slopes == max(slopes))[0][0]
    converter[index] = 'late_time_step_hits'
    # Other
    index = np.where(slopes == np.median(slopes))[0][0]
    converter[index] = 'change_in_sign_hits'

    cluster_assignments = fit_gmm.predict(matrix)
    t = open(output_file, 'w')
    t.write('variant_id\tgene_id\tcluster\n')
    used = {}
    for i,test_name in enumerate(test_names):
        rs_id = test_name.split('_')[0]
        used[test_name] = 1
        gene_id = test_name.split('_')[1]
        cluster_name = converter[cluster_assignments[i]]
        t.write(rs_id + '\t' + gene_id + '\t' + cluster_name + '\n')
    for test_name in pval_dict.keys():
        if test_name not in used:
            t.write(test_name.split('_')[0] + '\t' + test_name.split('_')[1] + '\t' + 'change_in_sign_hits' + '\n')
    t.close()


parameter_string = sys.argv[1]
cht_output_dir = sys.argv[2]
cht_visualization_dir = sys.argv[3]
pc_num = sys.argv[4]


input_root = cht_output_dir + 'cht_results_' + parameter_string + '_num_pc_' + pc_num + '_time_'

allelic_fraction_dict, pval_dict = extract_variant_gene_pair_time_step_info(input_root)

############################################################
# Classify based on linear regression of allelic fraction
############################################################
error_bound = 0
output_file = cht_output_dir + parameter_string + '_num_pc_' + pc_num + '_variant_gene_pair_clusters_af_lr_' + str(error_bound) + '.txt'
linear_regression_allelic_fraction_classification(allelic_fraction_dict, error_bound, output_file)

error_bound = .005
output_file = cht_output_dir + parameter_string + '_num_pc_' + pc_num + '_variant_gene_pair_clusters_af_lr_' + str(error_bound) + '.txt'
linear_regression_allelic_fraction_classification(allelic_fraction_dict, error_bound, output_file)

error_bound = .01
output_file = cht_output_dir + parameter_string + '_num_pc_' + pc_num + '_variant_gene_pair_clusters_af_lr_' + str(error_bound) + '.txt'
linear_regression_allelic_fraction_classification(allelic_fraction_dict, error_bound, output_file)


############################################################
# Classify based unsupervised clustering
############################################################
output_file = cht_output_dir + parameter_string + '_num_pc_' + pc_num + '_variant_gene_pair_clusters_kmeans.txt'
# unsupervised_clustering_kmeans(allelic_fraction_dict, pval_dict, output_file)

output_file = cht_output_dir + parameter_string + '_num_pc_' + pc_num + '_variant_gene_pair_clusters_full_gmm.txt'
# unsupervised_clustering_full_gmm(allelic_fraction_dict, pval_dict, output_file)

pvalue_thresh = .01
output_file = cht_output_dir + parameter_string + '_num_pc_' + pc_num + '_variant_gene_pair_clusters_full_gmm_' + str(pvalue_thresh) + '.txt'
# unsupervised_clustering_full_gmm_pvalue_thresh(allelic_fraction_dict, pval_dict, output_file, pvalue_thresh)

pvalue_thresh = .05
output_file = cht_output_dir + parameter_string + '_num_pc_' + pc_num + '_variant_gene_pair_clusters_full_gmm_' + str(pvalue_thresh) + '.txt'
# unsupervised_clustering_full_gmm_pvalue_thresh(allelic_fraction_dict, pval_dict, output_file, pvalue_thresh)

pvalue_thresh = .01
output_file = cht_output_dir + parameter_string + '_num_pc_' + pc_num + '_variant_gene_pair_clusters_kmeans_' + str(pvalue_thresh) + '.txt'
#unsupervised_clustering_kmeans_pvalue_thresh(allelic_fraction_dict, pval_dict, output_file, pvalue_thresh)

pvalue_thresh = .05
output_file = cht_output_dir + parameter_string + '_num_pc_' + pc_num + '_variant_gene_pair_clusters_kmeans_' + str(pvalue_thresh) + '.txt'
#unsupervised_clustering_kmeans_pvalue_thresh(allelic_fraction_dict, pval_dict, output_file, pvalue_thresh)
