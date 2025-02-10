#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 09:24:31 2021

@author: zhenzhang
"""
# this version  is based on the modification of the original 3/ leving out the dynamic preserving part
# finalize version
#from scipy.stats import fisher_exact
import scipy.stats as stats
from scipy.stats import chisquare
import sys
import numpy as np
#import re
import pandas as pd
#from pandas.core.frame import DataFrame
#from random import sample, choice
from Bio import Align
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pysam
#from scipy import stats
import biotite.sequence as seqbio
import biotite.sequence.align as align
import itertools as it
#import math
import multiprocessing as mp
import networkx as nx
import copy
#import seaborn as sns
import itertools
import matplotlib.pyplot as plt
#from fitter import Fitter, get_common_distributions, get_distributions
#from distfit import distfit
#from scipy.special import comb, perm
#from random import sample
from sklearn.cluster import SpectralClustering,DBSCAN,OPTICS,Birch
#from mpl_toolkits import mplot3d
#from scipy.special import psi
import scipy
from scipy.optimize import linear_sum_assignment
from random import sample
from communities.algorithms import louvain_method, girvan_newman, hierarchical_clustering
from networkx.algorithms.community import greedy_modularity_communities, naive_greedy_modularity_communities,\
    lukes_partitioning, asyn_lpa_communities, label_propagation_communities, asyn_fluidc
#sys.path.append('/Users/zhenzhang/Desktop/codehome/code_for_SNV_inference/') 
#import snv_inference_pac as sip
from sklearn.decomposition import PCA
#import concurrent.futures
#from mpire import WorkerPool
#from functools import partial
import time
from scipy.optimize import fsolve
import os
import collections as col





Standardbase = ['A', 'T', 'C', 'G', '-', 'N']

#main_url_save = '/lustre/project/Stat/s1155133512/Zhen/gene/bamdata/Simulation_2_AESM3/Data/'
#length_dir = '/lustre/project/Stat/s1155133512/Zhen/Data/npy_legnth/'
#length_dir = '/home/zhen/Data/npy_legnth/'
#fasta_file_path ='/lustre/project/Stat/s1155133512/Zhen/Data/hg38/ncbi_dataset/data/GCF_000001405.39/chr8.fna'
#fasta_file_path ='/home/zhen/Data/hg38/ncbi_dataset/data/GCF_000001405.39/chr8.fna'

def extimation_p(p_0,*import_data):
    corretion_eps_from_zero = 0.000000000001
    c_last, num_data = import_data
    pp_0 = p_0*np.log(p_0+corretion_eps_from_zero)*c_last+num_data*(1-p_0)
    return(pp_0)

def bam_2_msa_ins(start,seq,cigar):
    
    #start = Result_table.iloc[2,2]
    #seq = Result_table.iloc[2,1]
    #cigar = Result_table.iloc[2,5]
    df_tempt_ins2_xun_dic= {}
    df_tempt_ins2_xun_dic_bone = {}
    #df_tempt_bone_dic = {}
    
    back_bone_start = start
    present_seq = seq
    present_cigar = cigar
    ### insertion part

    bone_index_relative = 0
    seq_index_relative = 0
    #updated_str = ''
    for cig_xun in present_cigar:
        if cig_xun[0] in [0,7,8]:
            #updated_str= updated_str + present_seq[bone_index_relative:int(seq_index_relative+cig_xun[1])]
            for ll_per in range(cig_xun[1]):
                df_tempt_ins2_xun_dic_bone[back_bone_start+bone_index_relative+ll_per] = present_seq[seq_index_relative+ll_per]
            
            seq_index_relative = seq_index_relative + cig_xun[1]
            bone_index_relative = bone_index_relative + cig_xun[1]
            
        elif cig_xun[0] in [2]:
            #updated_str = updated_str + 'n'*cig_xun[1]
            for ll_per in range(cig_xun[1]):
                df_tempt_ins2_xun_dic_bone[back_bone_start+bone_index_relative+ll_per] = '-'
            bone_index_relative = bone_index_relative + cig_xun[1]
        elif cig_xun[0] in [1]:
            #updated_str = updated_str + present_seq[bone_index_relative:int(seq_index_relative+cig_xun[1])]
            df_tempt_ins2_xun_dic[back_bone_start+bone_index_relative-1] = present_seq[seq_index_relative:int(seq_index_relative+cig_xun[1])]
            seq_index_relative = seq_index_relative + cig_xun[1]
        elif cig_xun[0] in [4]:
            seq_index_relative = seq_index_relative + cig_xun[1]
    return([back_bone_start+bone_index_relative-1,df_tempt_ins2_xun_dic,df_tempt_ins2_xun_dic_bone]) 




def bam_2_msa_ins2(start,seq,cigar,region_start):
    
    #start = Result_table.iloc[2,2]
    #seq = Result_table.iloc[2,1]
    #cigar = Result_table.iloc[2,5]
    df_tempt_ins2_xun_dic= {}
    #df_tempt_bone_dic = {}
    
    back_bone_start = start - region_start
    present_seq = seq
    present_cigar = cigar
    ### insertion part

    bone_index_relative = 0
    seq_index_relative = 0
    #updated_str = ''
    for cig_xun in present_cigar:
        if cig_xun[0] in [0,7,8]:
            #updated_str= updated_str + present_seq[bone_index_relative:int(seq_index_relative+cig_xun[1])]
            seq_index_relative = seq_index_relative + cig_xun[1]
            bone_index_relative = bone_index_relative + cig_xun[1]
        elif cig_xun[0] in [2]:
            #updated_str = updated_str + 'n'*cig_xun[1]
            bone_index_relative = bone_index_relative + cig_xun[1]
        elif cig_xun[0] in [1]:
            #updated_str = updated_str + present_seq[bone_index_relative:int(seq_index_relative+cig_xun[1])]
            df_tempt_ins2_xun_dic[back_bone_start+bone_index_relative-1] = present_seq[seq_index_relative:int(seq_index_relative+cig_xun[1])]
            seq_index_relative = seq_index_relative + cig_xun[1]
        elif cig_xun[0] in [4]:
            seq_index_relative = seq_index_relative + cig_xun[1]
    return([back_bone_start+bone_index_relative-1,df_tempt_ins2_xun_dic]) 


       
def refine_the_fasta(ins_num_dic,reference_fsatq):
    finalized_ref = copy.deepcopy(reference_fsatq)
    n_c = len(reference_fsatq)
    cal_bone_index = 0
    tem_seq_bone = ''
    seg_ment_bone_seq = {}
    while cal_bone_index <=  n_c-1:
        if cal_bone_index not in ins_num_dic:
            tem_seq_bone = tem_seq_bone + finalized_ref[int(cal_bone_index)]
        else:
            tem_seq_bone = tem_seq_bone + finalized_ref[int(cal_bone_index)]
            seg_ment_bone_seq[cal_bone_index] = tem_seq_bone
            tem_seq_bone = ''
            
            #tem_seq_bone = ''
        
        cal_bone_index = cal_bone_index + 1    
    if len(tem_seq_bone)>0:
        seg_ment_bone_seq[-1] = tem_seq_bone
    final_ref = ''
    final_ref = ''
    for kk_key in sorted(ins_num_dic):
        #print(kk_key,ins_inserted_dic[kk_key])
        final_ref = final_ref + seg_ment_bone_seq[kk_key] + ins_num_dic[kk_key]*'N'
    if -1 in seg_ment_bone_seq:
      final_ref_standard = final_ref + seg_ment_bone_seq[-1]
    else:
      final_ref_standard = final_ref
    
    return(final_ref_standard)


def string_2_code(string_in):
    out_code = []
    for ss in range(len(string_in)):
        out_code.append(Standardbase.index(string_in[ss]))
    return(out_code)


def finding_sub_region(seq, cigar, strat_binding, standar_start, standar_end):
    #Notes_ = []
    cigar_tuples_all = cigar
    # read_sequence_name = read.query_name#reads name
    read_start = strat_binding
    read_sequence_all = seq
    collection_N_B = []
    start_pos = standar_start+1  # just for simulation
    end_pos = standar_end+1
    len_BestRefSeq = end_pos - start_pos + 1
    first_cigar_index = 0
    last_cigar_index = 0

    for first_index in range(len(cigar_tuples_all)):
        if cigar_tuples_all[first_index][0] in [3, 9]:
            collection_N_B.append(first_index)

    if len(collection_N_B) < 1:
        if read_start+1 <= end_pos:
            read_sequence = read_sequence_all
            cigar_tuples = cigar_tuples_all

            if read_start+1 >= start_pos:
                start_read = 0  # read conditional start
                start_con_pos = read_start-start_pos+1  # read conditional end
            else:
                start_read = start_pos - read_start - 1  # 0 base
                start_con_pos = 0

            # begining of the getting the end area
            count_reads_base = 0
            loops_real_tempt = 0
            loops_real = 0
            # get the precise start point of the read sequence
            if start_read == 0:
                for count_index in range(len(cigar_tuples)):
                    if cigar_tuples[count_index][0] == 4 or cigar_tuples[count_index][0] == 5 or cigar_tuples[count_index][0] == 6:
                        if cigar_tuples[count_index][0] == 4:
                            count_reads_base = count_reads_base + \
                                cigar_tuples[count_index][1]
                    else:
                        break
                count_reads_base_initial = count_reads_base  # 0 base, next is present
                remain_trure_base_gap = 0  # before the binds
                start_tupe_g2 = count_index
                loops_real = start_con_pos
                first_cigar_index = start_tupe_g2
                first_cigar_index_left_base = cigar_tuples[first_cigar_index][1]
            else:

                for count_index_test in range(len(cigar_tuples)):
                    if cigar_tuples[count_index_test][0] == 0 or cigar_tuples[count_index_test][0] == 7:
                        break  # getting the starting binding

                # Then count the previous reads index
                for count_index_reads in range(count_index_test):
                    if cigar_tuples[count_index_reads][0] == 4 or cigar_tuples[count_index_reads][0] == 5 or\
                            cigar_tuples[count_index_reads][0] == 6 or cigar_tuples[count_index_reads][0] == 1 or cigar_tuples[count_index_reads][0] == 8:
                        if cigar_tuples[count_index_reads][0] == 4:
                            count_reads_base = count_reads_base + \
                                cigar_tuples[count_index_reads][1]
                # after binding
                for count_index in range(count_index_test, len(cigar_tuples)):
                    if (cigar_tuples[count_index][0] == 0 or cigar_tuples[count_index][0] == 2 or cigar_tuples[count_index][0] == 3 or
                            cigar_tuples[count_index][0] == 7 or cigar_tuples[count_index][0] == 8):
                        loops_real_tempt = loops_real_tempt + \
                            cigar_tuples[count_index][1]
                        if loops_real_tempt - 1 >= (start_read):  # change
                            break
                    if cigar_tuples[count_index][0] == 0 or cigar_tuples[count_index][0] == 1 or cigar_tuples[count_index][0] == 4 or\
                            cigar_tuples[count_index][0] == 7 or cigar_tuples[count_index][0] == 8:
                        count_reads_base = count_reads_base + \
                            cigar_tuples[count_index][1]

                loops_real_tempt = (loops_real_tempt - 1 -
                                    cigar_tuples[count_index][1])  # 1base
                remain_trure_base_gap = start_read - loops_real_tempt
                if cigar_tuples[count_index][0] in [0, 7, 8, 1, 4]:

                    count_reads_base_initial = (
                        count_reads_base-1) + remain_trure_base_gap
                    count_reads_base = count_reads_base_initial  # change

                else:
                    count_reads_base = count_reads_base  # change it to 1 index
                    count_reads_base_initial = count_reads_base + 1

                loops_real = start_con_pos  # 0 base
                start_tupe_g2 = count_index

            if remain_trure_base_gap == 0:
                existing_flag = 0
                for remaining_tupe in range(start_tupe_g2, len(cigar_tuples)):
                    # for per_base in range(cigar_tuples[remaining_tupe][1]):
                    if cigar_tuples[remaining_tupe][0] == 0 or cigar_tuples[remaining_tupe][0] == 7:
                        loops_real = loops_real + \
                            cigar_tuples[remaining_tupe][1]
                        count_reads_base = count_reads_base + \
                            cigar_tuples[remaining_tupe][1]
                        if loops_real > len_BestRefSeq-1:

                            difference_gap = loops_real-(len_BestRefSeq)+1
                            loops_real = loops_real - difference_gap
                            count_reads_base = count_reads_base - difference_gap
                            last_cigar_index = remaining_tupe
                            last_cigar_index_left_base = cigar_tuples[last_cigar_index][1] - \
                                difference_gap+1

                            existing_flag = 1
                            break

                    # deletion and skip does not consume the base of the reads
                    elif cigar_tuples[remaining_tupe][0] == 2:
                        loops_real = loops_real + \
                            cigar_tuples[remaining_tupe][1]
                        if loops_real > len_BestRefSeq-1:

                            difference_gap = loops_real-(len_BestRefSeq)+1
                            loops_real = loops_real - difference_gap
                            last_cigar_index = remaining_tupe
                            last_cigar_index_left_base = cigar_tuples[last_cigar_index][1] - \
                                difference_gap+1
                            existing_flag = 1
                            break

                    elif cigar_tuples[remaining_tupe][0] == 8:
                        loops_real = loops_real + \
                            cigar_tuples[remaining_tupe][1]
                        count_reads_base = count_reads_base + \
                            cigar_tuples[remaining_tupe][1]
                        if loops_real > len_BestRefSeq-1:

                            difference_gap = loops_real-(len_BestRefSeq)+1
                            loops_real = loops_real - difference_gap
                            count_reads_base = count_reads_base - difference_gap
                            last_cigar_index = remaining_tupe
                            last_cigar_index_left_base = cigar_tuples[last_cigar_index][1] - \
                                difference_gap+1
                            existing_flag = 1
                            break

                    elif cigar_tuples[remaining_tupe][0] == 4 and remaining_tupe != (len(cigar_tuples)-1):
                        # Notes_.append(str(loops_real+1) + ':' + '#'+':4')# represnt the soft copy, did not use the sequence infomation
                        count_reads_base = count_reads_base + \
                            cigar_tuples[remaining_tupe][1]

                    elif cigar_tuples[remaining_tupe][0] == 1:
                        count_reads_base = count_reads_base + \
                            cigar_tuples[remaining_tupe][1]

                # here to set if the cigar is normally finished
                if not existing_flag:
                    # the last is 1 or other complete in === of 2,1,7,8
                    if cigar_tuples[remaining_tupe][0] != 4:
                        last_cigar_index = remaining_tupe
                        last_cigar_index_left_base = cigar_tuples[last_cigar_index][1]
                    else:  # the last is 4
                        if remaining_tupe - first_cigar_index > 1:
                            last_cigar_index = remaining_tupe - 1
                            last_cigar_index_left_base = cigar_tuples[last_cigar_index][1]
                        else:
                            last_cigar_index = remaining_tupe - 1
                            last_cigar_index_left_base = first_cigar_index_left_base

            else:
                # after mapping the begining, start get the recording
                pre_base = (cigar_tuples[start_tupe_g2]
                            [1]) - (remain_trure_base_gap)

                first_cigar_index = start_tupe_g2
                first_cigar_index_left_base = pre_base + 1

                existing_flag = 0
                if cigar_tuples[start_tupe_g2][0] == 0 or cigar_tuples[start_tupe_g2][0] == 7:
                    loops_real = loops_real + pre_base
                    count_reads_base = count_reads_base + pre_base
                    if loops_real > len_BestRefSeq-1:

                        difference_gap = loops_real-(len_BestRefSeq)+1
                        loops_real = loops_real - difference_gap
                        count_reads_base = count_reads_base - difference_gap
                        last_cigar_index = start_tupe_g2
                        last_cigar_index_left_base = cigar_tuples[last_cigar_index][1] - \
                            difference_gap+1
                        existing_flag = 1
                        # break

                # deletion and skip does not consume the base of the reads
                elif cigar_tuples[start_tupe_g2][0] == 2:
                    loops_real = loops_real + pre_base
                    if loops_real > len_BestRefSeq-1:

                        difference_gap = loops_real-(len_BestRefSeq)+1
                        loops_real = loops_real - difference_gap
                        last_cigar_index = start_tupe_g2
                        last_cigar_index_left_base = cigar_tuples[last_cigar_index][1] - \
                            difference_gap+1
                        existing_flag = 1
                        # break

                elif cigar_tuples[start_tupe_g2][0] == 8:
                    loops_real = loops_real + pre_base
                    count_reads_base = count_reads_base + pre_base
                    if loops_real > len_BestRefSeq-1:

                        difference_gap = loops_real-(len_BestRefSeq)+1
                        loops_real = loops_real - difference_gap
                        count_reads_base = count_reads_base - difference_gap
                        last_cigar_index = start_tupe_g2
                        last_cigar_index_left_base = cigar_tuples[last_cigar_index][1] - \
                            difference_gap+1
                        existing_flag = 1
                        # break

                elif cigar_tuples[start_tupe_g2][0] == 4 and start_tupe_g2 != (len(cigar_tuples)):
                    # Notes_.append(str(loops_real+1) + ':' + '#'+':4')# represnt the soft copy, did not use the sequence infomation
                    count_reads_base = count_reads_base + pre_base

                elif cigar_tuples[start_tupe_g2][0] == 1:
                    count_reads_base = count_reads_base + pre_base

    # other complete tupes
                if loops_real < len_BestRefSeq-1 and start_tupe_g2+1 < len(cigar_tuples):
                    existing_flag = 0

                    for remaining_tupe in range(start_tupe_g2+1, len(cigar_tuples)):
                        if cigar_tuples[remaining_tupe][0] == 0 or cigar_tuples[remaining_tupe][0] == 7:
                            loops_real = loops_real + \
                                cigar_tuples[remaining_tupe][1]
                            count_reads_base = count_reads_base + \
                                cigar_tuples[remaining_tupe][1]
                            if loops_real > len_BestRefSeq-1:

                                difference_gap = loops_real-(len_BestRefSeq)+1
                                loops_real = loops_real - difference_gap
                                count_reads_base = count_reads_base - difference_gap
                                last_cigar_index = remaining_tupe
                                last_cigar_index_left_base = cigar_tuples[last_cigar_index][1] - \
                                    difference_gap+1
                                existing_flag = 1
                                break

                        # deletion and skip does not consume the base of the reads
                        elif cigar_tuples[remaining_tupe][0] == 2:
                            loops_real = loops_real + \
                                cigar_tuples[remaining_tupe][1]
                            if loops_real > len_BestRefSeq-1:

                                difference_gap = loops_real-(len_BestRefSeq)+1
                                loops_real = loops_real - difference_gap
                                last_cigar_index = remaining_tupe
                                last_cigar_index_left_base = cigar_tuples[last_cigar_index][1] - \
                                    difference_gap+1
                                existing_flag = 1
                                break

                        elif cigar_tuples[remaining_tupe][0] == 8:
                            loops_real = loops_real + \
                                cigar_tuples[remaining_tupe][1]
                            count_reads_base = count_reads_base + \
                                cigar_tuples[remaining_tupe][1]
                            if loops_real > len_BestRefSeq-1:

                                difference_gap = loops_real-(len_BestRefSeq)+1
                                loops_real = loops_real - difference_gap
                                count_reads_base = count_reads_base - difference_gap
                                last_cigar_index = remaining_tupe
                                last_cigar_index_left_base = cigar_tuples[last_cigar_index][1] - \
                                    difference_gap+1
                                existing_flag = 1
                                break

                        elif cigar_tuples[remaining_tupe][0] == 4 and remaining_tupe != (len(cigar_tuples)-1):
                            # Notes_.append(str(loops_real+1) + ':' + '#'+':4')# represnt the soft copy, did not use the sequence infomation

                            count_reads_base = count_reads_base + \
                                cigar_tuples[remaining_tupe][1]
                            existing_flag = 0

                        elif cigar_tuples[remaining_tupe][0] == 1:
                            count_reads_base = count_reads_base + \
                                cigar_tuples[remaining_tupe][1]
                            existing_flag = 0

                    # here to set if the cigar is normally finished
                    if not existing_flag:
                        # the last is 1 or other complete in === of 2,1,7,8
                        if cigar_tuples[remaining_tupe][0] != 4:
                            last_cigar_index = remaining_tupe
                            last_cigar_index_left_base = cigar_tuples[last_cigar_index][1]
                        else:  # the last is 4
                            if remaining_tupe - first_cigar_index > 1:
                                last_cigar_index = remaining_tupe - 1
                                last_cigar_index_left_base = cigar_tuples[last_cigar_index][1]
                            else:
                                last_cigar_index = remaining_tupe - 1
                                last_cigar_index_left_base = first_cigar_index_left_base

                        # Omit another situation due to no possibility
                else:
                    if not existing_flag:
                        # the last is 1 or other complete in === of 2,1,7,8
                        if cigar_tuples[start_tupe_g2][0] != 4:
                            last_cigar_index = start_tupe_g2
                            last_cigar_index_left_base = first_cigar_index_left_base

            #start_bam = read.reference_start
            if start_read == 0 and loops_real < len_BestRefSeq-1:  # make sure every condition is 0 base end
                sequence_start = start_con_pos+start_pos-1
                sequence_end = sequence_start + (loops_real-start_con_pos)-1
                relative_start = start_con_pos
                relative_end = loops_real-1
                true_start_insequence = count_reads_base_initial
                true_end_insequence = count_reads_base-1
            else:
                sequence_start = start_con_pos+start_pos-1
                sequence_end = sequence_start + (loops_real-start_con_pos)
                relative_start = start_con_pos
                relative_end = loops_real
                true_start_insequence = count_reads_base_initial
                true_end_insequence = count_reads_base

            # Cigar information revised
            if last_cigar_index-first_cigar_index > 0:
                revised_cigar = [(cigar_tuples_all[first_cigar_index][0], first_cigar_index_left_base)]+cigar_tuples_all[(first_cigar_index+1):last_cigar_index] +\
                    [(cigar_tuples_all[last_cigar_index]
                      [0], last_cigar_index_left_base)]
            else:
                revised_cigar = [
                    (cigar_tuples_all[last_cigar_index][0], last_cigar_index_left_base)]

            # revised sequence
            read_temp = read_sequence[int(
                true_start_insequence):int(true_end_insequence+1)]
            #df_tempt.loc[loopl] = [base,read_sequence_name,start_bam,sequence_start,sequence_end,true_start_insequence,true_end_insequence,relative_start,relative_end,Seq_read,revised_cigar,0]
            #loopl = loopl + 1

    else:
        collection_N_B.insert(0, -1)
        collection_N_B.insert(len(collection_N_B), len(cigar_tuples_all))

        temp_consume_read = 0
        read_sequence = read_sequence_all
        #start_bam = read.reference_start
        for second_index in range(len(collection_N_B)-1):
            #Notes_ = []
            cigar_tuples = cigar_tuples_all[collection_N_B[second_index] +
                                            1:collection_N_B[second_index+1]]
            seg = 0
            ##previous_code##
            # if not second_index:###initial does start with 0 base
            if read_start+1 <= end_pos:
                if read_start+1 >= start_pos:
                    start_read = 0  # read conditional start
                    start_con_pos = read_start-start_pos+1  # read conditional end
                else:
                    start_read = start_pos - read_start - 1  # 0 base
                    start_con_pos = 0

                # begining of the getting the end area
                count_reads_base = 0
                loops_real_tempt = 0
                loops_real = 0
                # get the precise start point of the read sequence
                if start_read == 0:
                    for count_index in range(len(cigar_tuples)):
                        if cigar_tuples[count_index][0] == 4 or cigar_tuples[count_index][0] == 5 or cigar_tuples[count_index][0] == 6:
                            if cigar_tuples[count_index][0] == 4:
                                count_reads_base = count_reads_base + \
                                    cigar_tuples[count_index][1]
                        else:
                            break
                    count_reads_base_initial = count_reads_base  # 0 base, next is present
                    remain_trure_base_gap = 0  # before the binds
                    start_tupe_g2 = count_index
                    loops_real = start_con_pos
                    first_cigar_index = start_tupe_g2
                    first_cigar_index_left_base = cigar_tuples[first_cigar_index][1]
                    forward_flag = 1
                else:
                    for count_index_test in range(len(cigar_tuples)):
                        if cigar_tuples[count_index_test][0] == 0 or cigar_tuples[count_index_test][0] == 7:
                            break  # getting the starting binding

                    # Then count the previous reads index
                    for count_index_reads in range(count_index_test):
                        if cigar_tuples[count_index_reads][0] == 4 or cigar_tuples[count_index_reads][0] == 5 or\
                                cigar_tuples[count_index_reads][0] == 6 or cigar_tuples[count_index_reads][0] == 1 or cigar_tuples[count_index_reads][0] == 8:
                            if cigar_tuples[count_index_reads][0] == 4:
                                count_reads_base = count_reads_base + \
                                    cigar_tuples[count_index_reads][1]
                    # after binding
                    for count_index in range(count_index_test, len(cigar_tuples)):
                        if (cigar_tuples[count_index][0] == 0 or cigar_tuples[count_index][0] == 2 or cigar_tuples[count_index][0] == 3 or
                                cigar_tuples[count_index][0] == 7 or cigar_tuples[count_index][0] == 8):
                            loops_real_tempt = loops_real_tempt + \
                                cigar_tuples[count_index][1]
                            if loops_real_tempt - 1 >= (start_read):  # change
                                break
                        if cigar_tuples[count_index][0] == 0 or cigar_tuples[count_index][0] == 1 or cigar_tuples[count_index][0] == 4 or\
                                cigar_tuples[count_index][0] == 7 or cigar_tuples[count_index][0] == 8:
                            count_reads_base = count_reads_base + \
                                cigar_tuples[count_index][1]

                    if loops_real_tempt - 1 >= start_read:
                        loops_real_tempt = (
                            loops_real_tempt - 1 - cigar_tuples[count_index][1])  # 1base
                        remain_trure_base_gap = start_read - loops_real_tempt
                        if cigar_tuples[count_index][0] in [0, 7, 8, 1, 4]:
                            count_reads_base_initial = count_reads_base-1 + remain_trure_base_gap
                            count_reads_base = count_reads_base_initial  # change

                        else:
                            count_reads_base = count_reads_base  # change it to 1 index
                            count_reads_base_initial = count_reads_base + 1

                        loops_real = start_con_pos  # 0 base
                        start_tupe_g2 = count_index
                        forward_flag = 1

                    else:
                        forward_flag = 0
                        sequence_start = read_start
                        sequence_end = read_start
                        relative_start = 0
                        relative_end = 0
                        true_start_insequence = temp_consume_read
                        true_end_insequence = temp_consume_read + count_reads_base-1
                if forward_flag:
                    if remain_trure_base_gap == 0:
                        existing_flag = 0
                        for remaining_tupe in range(start_tupe_g2, len(cigar_tuples)):
                            # for per_base in range(cigar_tuples[remaining_tupe][1]):
                            if cigar_tuples[remaining_tupe][0] == 0 or cigar_tuples[remaining_tupe][0] == 7:
                                loops_real = loops_real + \
                                    cigar_tuples[remaining_tupe][1]
                                count_reads_base = count_reads_base + \
                                    cigar_tuples[remaining_tupe][1]
                                if loops_real > len_BestRefSeq-1:

                                    difference_gap = loops_real - \
                                        (len_BestRefSeq)+1
                                    loops_real = loops_real - difference_gap
                                    count_reads_base = count_reads_base - difference_gap
                                    last_cigar_index = remaining_tupe
                                    last_cigar_index_left_base = cigar_tuples[
                                        last_cigar_index][1]-difference_gap+1
                                    existing_flag = 1
                                    break

                            # deletion and skip does not consume the base of the reads
                            elif cigar_tuples[remaining_tupe][0] == 2:
                                loops_real = loops_real + \
                                    cigar_tuples[remaining_tupe][1]
                                if loops_real > len_BestRefSeq-1:

                                    difference_gap = loops_real - \
                                        (len_BestRefSeq)+1
                                    loops_real = loops_real - difference_gap
                                    last_cigar_index = remaining_tupe
                                    last_cigar_index_left_base = cigar_tuples[
                                        last_cigar_index][1]-difference_gap+1
                                    existing_flag = 1
                                    break

                            elif cigar_tuples[remaining_tupe][0] == 8:
                                loops_real = loops_real + \
                                    cigar_tuples[remaining_tupe][1]
                                count_reads_base = count_reads_base + \
                                    cigar_tuples[remaining_tupe][1]
                                if loops_real > len_BestRefSeq-1:

                                    difference_gap = loops_real - \
                                        (len_BestRefSeq)+1
                                    loops_real = loops_real - difference_gap
                                    count_reads_base = count_reads_base - difference_gap
                                    last_cigar_index = remaining_tupe
                                    last_cigar_index_left_base = cigar_tuples[
                                        last_cigar_index][1]-difference_gap+1
                                    existing_flag = 1
                                    break

                            elif cigar_tuples[remaining_tupe][0] == 4 and remaining_tupe != (len(cigar_tuples)-1):
                                # Notes_.append(str(loops_real+1) + ':' + '#'+':4')# represnt the soft copy, did not use the sequence infomation
                                count_reads_base = count_reads_base + \
                                    cigar_tuples[remaining_tupe][1]

                            elif cigar_tuples[remaining_tupe][0] == 1:
                                count_reads_base = count_reads_base + \
                                    cigar_tuples[remaining_tupe][1]

                        # here to set if the cigar is normally finished
                        if not existing_flag:
                            # the last is 1 or other complete in === of 2,1,7,8
                            if cigar_tuples[remaining_tupe][0] != 4:
                                last_cigar_index = remaining_tupe
                                last_cigar_index_left_base = cigar_tuples[last_cigar_index][1]
                            else:  # the last is 4
                                if remaining_tupe - first_cigar_index > 1:
                                    last_cigar_index = remaining_tupe - 1
                                    last_cigar_index_left_base = cigar_tuples[last_cigar_index][1]
                                else:
                                    last_cigar_index = remaining_tupe - 1
                                    last_cigar_index_left_base = first_cigar_index_left_base

                    else:
                        # after mapping the begining, start get the recording
                        pre_base = (
                            cigar_tuples[start_tupe_g2][1]) - (remain_trure_base_gap)

                        first_cigar_index = start_tupe_g2
                        first_cigar_index_left_base = pre_base + 1

                        existing_flag = 0
                        if cigar_tuples[start_tupe_g2][0] == 0 or cigar_tuples[start_tupe_g2][0] == 7:
                            loops_real = loops_real + pre_base
                            count_reads_base = count_reads_base + pre_base
                            if loops_real > len_BestRefSeq-1:

                                difference_gap = loops_real-(len_BestRefSeq)+1
                                loops_real = loops_real - difference_gap
                                count_reads_base = count_reads_base - difference_gap
                                last_cigar_index = start_tupe_g2
                                last_cigar_index_left_base = cigar_tuples[last_cigar_index][1] - \
                                    difference_gap+1
                                existing_flag = 1
                                break
                        # deletion and skip does not consume the base of the reads
                        elif cigar_tuples[start_tupe_g2][0] == 2:
                            loops_real = loops_real + pre_base
                            if loops_real > len_BestRefSeq-1:

                                difference_gap = loops_real-(len_BestRefSeq)+1
                                loops_real = loops_real - difference_gap
                                last_cigar_index = start_tupe_g2
                                last_cigar_index_left_base = cigar_tuples[last_cigar_index][1] - \
                                    difference_gap+1
                                existing_flag = 1
                                break
                        elif cigar_tuples[start_tupe_g2][0] == 8:
                            loops_real = loops_real + pre_base
                            count_reads_base = count_reads_base + pre_base
                            if loops_real > len_BestRefSeq-1:

                                difference_gap = loops_real-(len_BestRefSeq)+1
                                loops_real = loops_real - difference_gap
                                count_reads_base = count_reads_base - difference_gap
                                last_cigar_index = start_tupe_g2
                                last_cigar_index_left_base = cigar_tuples[last_cigar_index][1] - \
                                    difference_gap+1
                                existing_flag = 1
                                break

                        elif cigar_tuples[start_tupe_g2][0] == 4 and start_tupe_g2 != (len(cigar_tuples)):
                            # Notes_.append(str(loops_real+1) + ':' + '#'+':4')# represnt the soft copy, did not use the sequence infomation
                            count_reads_base = count_reads_base + pre_base

                        elif cigar_tuples[start_tupe_g2][0] == 1:
                            count_reads_base = count_reads_base + pre_base

            # other complete tupes
                        if loops_real < len_BestRefSeq-1 and start_tupe_g2+1 < len(cigar_tuples):
                            existing_flag = 0
                            for remaining_tupe in range(start_tupe_g2+1, len(cigar_tuples)):
                                if cigar_tuples[remaining_tupe][0] == 0 or cigar_tuples[remaining_tupe][0] == 7:
                                    loops_real = loops_real + \
                                        cigar_tuples[remaining_tupe][1]
                                    count_reads_base = count_reads_base + \
                                        cigar_tuples[remaining_tupe][1]
                                    if loops_real > len_BestRefSeq-1:

                                        difference_gap = loops_real - \
                                            (len_BestRefSeq)+1
                                        loops_real = loops_real - difference_gap
                                        count_reads_base = count_reads_base - difference_gap
                                        last_cigar_index = remaining_tupe
                                        last_cigar_index_left_base = cigar_tuples[
                                            last_cigar_index][1]-difference_gap+1
                                        existing_flag = 1
                                        break

                                # deletion and skip does not consume the base of the reads
                                elif cigar_tuples[remaining_tupe][0] == 2:
                                    loops_real = loops_real + \
                                        cigar_tuples[remaining_tupe][1]
                                    if loops_real > len_BestRefSeq-1:

                                        difference_gap = loops_real - \
                                            (len_BestRefSeq)+1
                                        loops_real = loops_real - difference_gap
                                        last_cigar_index = remaining_tupe
                                        last_cigar_index_left_base = cigar_tuples[
                                            last_cigar_index][1]-difference_gap+1
                                        existing_flag = 1
                                        break

                                elif cigar_tuples[remaining_tupe][0] == 8:
                                    loops_real = loops_real + \
                                        cigar_tuples[remaining_tupe][1]
                                    count_reads_base = count_reads_base + \
                                        cigar_tuples[remaining_tupe][1]
                                    if loops_real > len_BestRefSeq-1:

                                        difference_gap = loops_real - \
                                            (len_BestRefSeq)+1
                                        loops_real = loops_real - difference_gap
                                        count_reads_base = count_reads_base - difference_gap
                                        last_cigar_index = remaining_tupe
                                        last_cigar_index_left_base = cigar_tuples[
                                            last_cigar_index][1]-difference_gap+1
                                        existing_flag = 1
                                        break

                                elif cigar_tuples[remaining_tupe][0] == 4 and remaining_tupe != (len(cigar_tuples)-1):
                                    # Notes_.append(str(loops_real+1) + ':' + '#'+':4')# represnt the soft copy, did not use the sequence infomation
                                    count_reads_base = count_reads_base + \
                                        cigar_tuples[remaining_tupe][1]
                                    existing_flag = 0

                                elif cigar_tuples[remaining_tupe][0] == 1:
                                    count_reads_base = count_reads_base + \
                                        cigar_tuples[remaining_tupe][1]
                                    existing_flag = 0

                            # here to set if the cigar is normally finished
                            if not existing_flag:
                                # the last is 1 or other complete in === of 2,1,7,8
                                if cigar_tuples[remaining_tupe][0] != 4:
                                    last_cigar_index = remaining_tupe
                                    last_cigar_index_left_base = cigar_tuples[last_cigar_index][1]
                                else:  # the last is 4
                                    if remaining_tupe - first_cigar_index > 1:
                                        last_cigar_index = remaining_tupe - 1
                                        last_cigar_index_left_base = cigar_tuples[last_cigar_index][1]
                            else:
                                last_cigar_index = remaining_tupe - 1
                                last_cigar_index_left_base = first_cigar_index_left_base
                        else:
                            if not existing_flag:
                                # the last is 1 or other complete in === of 2,1,7,8
                                if cigar_tuples[start_tupe_g2][0] != 4:
                                    last_cigar_index = start_tupe_g2
                                    last_cigar_index_left_base = first_cigar_index_left_base

                    #start_bam = read.reference_start
                    if start_read == 0 and loops_real < len_BestRefSeq-1:  # make sure every condition is 0 base end
                        sequence_start = start_con_pos+start_pos-1
                        sequence_end = sequence_start + \
                            (loops_real-start_con_pos)-1
                        relative_start = start_con_pos
                        relative_end = loops_real-1
                        true_start_insequence = temp_consume_read+count_reads_base_initial
                        true_end_insequence = temp_consume_read+count_reads_base-1
                    else:
                        sequence_start = start_con_pos+start_pos-1
                        sequence_end = sequence_start + \
                            (loops_real-start_con_pos)
                        relative_start = start_con_pos
                        relative_end = loops_real
                        true_start_insequence = temp_consume_read+count_reads_base_initial
                        true_end_insequence = temp_consume_read+count_reads_base

                    # Cigar information revised
                    if last_cigar_index-first_cigar_index > 0:
                        revised_cigar = [(cigar_tuples_all[first_cigar_index][0], first_cigar_index_left_base)]+cigar_tuples_all[first_cigar_index+1:last_cigar_index] +\
                            [(cigar_tuples_all[last_cigar_index]
                              [0], last_cigar_index_left_base)]
                    else:
                        revised_cigar = [
                            (cigar_tuples_all[last_cigar_index][0], last_cigar_index_left_base)]

                    # revised sequence
                    read_temp = read_sequence[int(
                        true_start_insequence):int(true_end_insequence+1)]
                    #df_tempt.loc[loopl] = [base,read_sequence_name,start_bam,sequence_start,sequence_end,true_start_insequence,true_end_insequence,relative_start,relative_end,Seq_read,revised_cigar,0]
                    #loopl = loopl + 1

                # else:
                #    df_tempt.loc[loopl] = [base,read_sequence_name,start_bam,sequence_start,sequence_end,true_start_insequence,true_end_insequence,relative_start,relative_end,[],[],0]
                #    loopl = loopl + 1
                    #print('not enough',sample,info,read_sequence_name,flush=True)
                    ###end of the previous code ##
            # else:
                    #print('over the end',sample,info,read_sequence_name,read_start,end_pos,flush=True)
            #        df_tempt.loc[loopl] = [base,read_sequence_name,start_bam,read_start,read_start,0,0,0,0,[],[],0]
            #        loopl = loopl + 1
            #        continue
            for third_index in range(len(cigar_tuples)):
                if cigar_tuples[third_index][0] in [0, 2, 7, 8]:
                    seg = seg + cigar_tuples[third_index][1]
            if collection_N_B[second_index+1] < len(cigar_tuples_all):
                if cigar_tuples_all[collection_N_B[second_index+1]][0] in [3, 9]:
                    if cigar_tuples_all[collection_N_B[second_index+1]][0] == 3:
                        read_start = read_start + seg + \
                            cigar_tuples_all[collection_N_B[second_index+1]][1]
                    else:
                        read_start = read_start + seg - \
                            cigar_tuples_all[collection_N_B[second_index+1]][1]

            temp_consume_read = true_end_insequence + 1

# eps
    return([read_temp, revised_cigar, relative_start, relative_end])


def coding_mismatch_all_insertion_error2_per(ii):
    #seq = Re
    global out_put_df_list
    start_pos = out_put_df_list[ii][0]
    seq_read = out_put_df_list[ii][1]
    end_pos = start_pos + len(seq_read) - 1
    seq_read_code = string_2_code(seq_read)
    refined_count_table_matrix = np.zeros([len(seq_read_code),7],dtype=int)
    tem_dic = {}
    tem_final_code_ = []
    for ii_base in range(len(seq_read_code)):
        tem_dic[start_pos+ii_base] = seq_read_code[ii_base]
        tem_final_code_.append(str(seq_read_code[ii_base]+6*(start_pos+ii_base)))
        refined_count_table_matrix[ii_base,1+int(seq_read_code[ii_base])] = 1
    refined_count_table_matrix[:,0] = list(range(start_pos,end_pos+1))   
    return([start_pos,end_pos,refined_count_table_matrix,tem_dic,tem_final_code_])


def coding_mismatch_all_insertion_error2(in_df, refied_reference):
    #seq = Re
    reference_code = string_2_code(refied_reference)
    whole_length = len(reference_code)
    refined_count_table = pd.DataFrame(data=np.zeros([whole_length, 10]), columns=[
                                       'ID', 'pos', 'A', 'T', 'C', 'G', '-', 'N', 'Final', 'flag'])
    final_indicator_table = pd.DataFrame(
        columns=['ID', 'start', 'end', 'collection_wrong', 'collection_correct'])
    for read_index in range(len(in_df)):
        start_pos = in_df.iloc[read_index, 0]
        seq_read = in_df.iloc[read_index, 1]
        seq_read_code = string_2_code(seq_read)
        tem_collection_wrong = []
        tem_collection_correct = []
        for ii in range(len(seq_read)):
            tem_code = (start_pos+ii)*6+seq_read_code[ii]
            if seq_read_code[ii] != reference_code[ii+int(start_pos)]:

                tem_collection_wrong.append(tem_code)
            else:
                tem_collection_correct.append(tem_code)
            refined_count_table.iloc[int(start_pos+ii), 2+int(seq_read_code[ii])
                                     ] = refined_count_table.iloc[int(start_pos+ii), 2+int(seq_read_code[ii])] + 1
        final_indicator_table.loc[read_index] = [
            read_index, start_pos, start_pos+len(seq_read)-1, tem_collection_wrong, tem_collection_correct]
    return([refined_count_table, final_indicator_table])


def adding_ins_to_bone(Result_table, ins_inserted_dic, dic_multiple_ali_ins, final_ref_code):
    #Result_table = Result_table_xun
    #ins_inserted_dic = ins_inserted_dic_nei
    #dic_multiple_ali_ins = dic_multiple_ali_ins_nei
    #final_ref_code = final_ref_code_standard_nei
    # whole_ins_dic=ins_inserted_dic
    key_ = sorted(ins_inserted_dic.keys())
    collection_of_sites = []
    #collection_of_sites_correct = []
    #new_data_frame = pd.DataFrame(columns=['start','refined_seq'])
    out_put_df = pd.DataFrame(columns=['start', 'refined_seq_'])
    for ii in range(len(Result_table)):
        seq = Result_table.iloc[ii, 1]
        cigar_info = Result_table.iloc[ii, 5]# 5 for remove the predic start and end
        start_pos = Result_table.iloc[ii, 2]
        end_pos = Result_table.iloc[ii, 3]
        ## since we do not care the clipping information/we cut the cigar
        pre_clip = 0
        after_clip = 0
        l_clip_flag = 0
        r_clip_flag = 0
        if cigar_info[0][0] in [4]:
            pre_clip = cigar_info[0][1]
            l_clip_flag = 1
        if cigar_info[-1][0] in [4]:
            after_clip = cigar_info[-1][1]
            r_clip_flag = 1
        seq_refined = seq[pre_clip:len(seq)-after_clip]
        refined_cigar = cigar_info[l_clip_flag:len(cigar_info)-r_clip_flag]
        
        # returned_dic={}
        tem_right_shift = 0
        if len(key_)>0:
            if start_pos > min(key_):
    
                st_index = 0
    
                while start_pos > key_[st_index]:
                    max_len = ins_inserted_dic[key_[st_index]]
                    tem_right_shift = tem_right_shift + max_len
                    st_index = st_index + 1
                    if st_index == len(key_):
                        break
        register_strat = start_pos + tem_right_shift
        key_ = np.array(key_)
        remaining_key_ = key_[np.where((key_ >= start_pos) & (key_ < end_pos))]
        if len(remaining_key_) > 0:
            # catch the bone
            tem_seq = ''  # only cosider the bone
            seq_index_o = 0
            index_bone = start_pos
            #segment_insertion = []
            #bound_index = []
            seq_specific_dic = {}
            for cigar_index in range(len(refined_cigar)):
                if refined_cigar[cigar_index][0] in [2]:
                    tem_seq = tem_seq + '-'*refined_cigar[cigar_index][1]
                    index_bone = index_bone + refined_cigar[cigar_index][1]
                elif refined_cigar[cigar_index][0] in [0, 7, 8]:

                    tem_seq = tem_seq + \
                        seq_refined[seq_index_o:seq_index_o+refined_cigar[cigar_index][1]]
                    index_bone = index_bone + refined_cigar[cigar_index][1]
                    seq_index_o = seq_index_o + refined_cigar[cigar_index][1]
                elif refined_cigar[cigar_index][0] in [1]:
                    # bound_index.append([seq_index_o,seq_index_o+ cigar_info[cigar_index][1]]) # the form is [),we can not get the largest value
                    seq_specific_dic[index_bone -
                                     1] = seq_refined[seq_index_o:seq_index_o+refined_cigar[cigar_index][1]]
                    seq_index_o = seq_index_o + refined_cigar[cigar_index][1]
                #elif  cigar_info[cigar_index][0] in [4]:
                    
                    #seq_index_o = seq_index_o + cigar_info[cigar_index][1]

            # cut the bone index according to overall situation
            cal_bone_index = start_pos
            tem_seq_bone = ''
            seg_ment_bone_seq_xun = {}
            while cal_bone_index <= end_pos:
                if cal_bone_index not in remaining_key_:
                    tem_seq_bone = tem_seq_bone + \
                        tem_seq[int(cal_bone_index-start_pos)]
                else:
                    tem_seq_bone = tem_seq_bone + \
                        tem_seq[int(cal_bone_index-start_pos)]
                    seg_ment_bone_seq_xun[cal_bone_index] = tem_seq_bone
                    tem_seq_bone = ''

                    #tem_seq_bone = ''

                cal_bone_index = cal_bone_index + 1
            if len(tem_seq_bone) > 0:
                seg_ment_bone_seq_xun[-1] = tem_seq_bone

            # makeup insertion segment
            finalized_ins_dic = {}
            for kkey in remaining_key_:
                if ins_inserted_dic[kkey] < 2:  # only one insertion
                    if kkey in seq_specific_dic:
                        finalized_ins_dic[kkey] = seq_specific_dic[kkey]
                    else:
                        finalized_ins_dic[kkey] = 'N'
                        #back_ground_ins = max(whole_ins_dic[ss_key_],key=len,default = '')
                else:
                    if kkey in seq_specific_dic:
                        present_ins = seq_specific_dic[kkey]
                        tem_data_frame = pd.DataFrame(
                            dic_multiple_ali_ins[kkey])
                        refined_ins = (
                            tem_data_frame[tem_data_frame[0] == present_ins]).iloc[0, 1]
                        finalized_ins_dic[kkey] = refined_ins
                    else:
                        finalized_ins_dic[kkey] = 'N'*ins_inserted_dic[kkey]

            # finished line for insertion and begin for make up the whole sequence
            final_string = ''
            for key_key in sorted(finalized_ins_dic):
                final_string = final_string + \
                    seg_ment_bone_seq_xun[key_key] + finalized_ins_dic[key_key]

            final_string = final_string + seg_ment_bone_seq_xun[-1]

        else:
            final_string = seq

        ###
        # begin mapping the mismatch and extracting insertion information
        final_string_code = string_2_code(final_string)
        ref_start = register_strat
        #final_dic = {}
        #cigar_o = []
        for ss in range(len(final_string_code)):
            if final_string_code[ss] != final_ref_code[int(ss+ref_start)] and (ss+ref_start) not in collection_of_sites:
                collection_of_sites.append(ss+ref_start)

        #new_data_frame.loc[ii] = [register_strat,final_string_code]
        # based on final_string to infer

        out_put_df.loc[ii] = [register_strat,
                              final_string]

    collection_of_sites_array = np.array(sorted(collection_of_sites))
    return([out_put_df, collection_of_sites_array])



def adding_ins_to_bone_per(ii):
    #Result_table = Result_table_xun
    #ins_inserted_dic = ins_inserted_dic_nei
    #dic_multiple_ali_ins = dic_multiple_ali_ins_nei
    #final_ref_code = final_ref_code_standard_nei
    # whole_ins_dic=ins_inserted_dic
    global Result_table,ins_inserted_dic, dic_multiple_ali_ins,final_ref_code
    key_ = sorted(ins_inserted_dic.keys())
    #collection_of_sites = []
    #collection_of_sites_correct = []
    #new_data_frame = pd.DataFrame(columns=['start','refined_seq'])
    #out_put_df = pd.DataFrame(columns=['start', 'refined_seq_'])
    #for ii in range(len(Result_table)):
    seq = Result_table.iloc[ii, 1]
    cigar_info = Result_table.iloc[ii, 5]# 5 for remove the predic start and end
    start_pos = Result_table.iloc[ii, 2]
    end_pos = Result_table.iloc[ii, 3]
    ## since we do not care the clipping information/we cut the cigar
    pre_clip = 0
    after_clip = 0
    l_clip_flag = 0
    r_clip_flag = 0
    if cigar_info[0][0] in [4]:
        pre_clip = cigar_info[0][1]
        l_clip_flag = 1
    if cigar_info[-1][0] in [4]:
        after_clip = cigar_info[-1][1]
        r_clip_flag = 1
    seq_refined = seq[pre_clip:len(seq)-after_clip]
    refined_cigar = cigar_info[l_clip_flag:len(cigar_info)-r_clip_flag]
    
    # returned_dic={}
    tem_right_shift = 0
    if len(key_)>0:
        if start_pos > min(key_):
    
            st_index = 0
    
            while start_pos > key_[st_index]:
                max_len = ins_inserted_dic[key_[st_index]]
                tem_right_shift = tem_right_shift + max_len
                st_index = st_index + 1
                if st_index == len(key_):
                    break
    register_strat = start_pos + tem_right_shift
    key_ = np.array(key_)
    remaining_key_ = key_[np.where((key_ >= start_pos) & (key_ < end_pos))]
    if len(remaining_key_) > 0:
        # catch the bone
        tem_seq = ''  # only cosider the bone
        seq_index_o = 0
        index_bone = start_pos
        #segment_insertion = []
        #bound_index = []
        seq_specific_dic = {}
        for cigar_index in range(len(refined_cigar)):
            if refined_cigar[cigar_index][0] in [2]:
                tem_seq = tem_seq + '-'*refined_cigar[cigar_index][1]
                index_bone = index_bone + refined_cigar[cigar_index][1]
            elif refined_cigar[cigar_index][0] in [0, 7, 8]:

                tem_seq = tem_seq + \
                    seq_refined[seq_index_o:seq_index_o+refined_cigar[cigar_index][1]]
                index_bone = index_bone + refined_cigar[cigar_index][1]
                seq_index_o = seq_index_o + refined_cigar[cigar_index][1]
            elif refined_cigar[cigar_index][0] in [1]:
                # bound_index.append([seq_index_o,seq_index_o+ cigar_info[cigar_index][1]]) # the form is [),we can not get the largest value
                seq_specific_dic[index_bone -
                                 1] = seq_refined[seq_index_o:seq_index_o+refined_cigar[cigar_index][1]]
                seq_index_o = seq_index_o + refined_cigar[cigar_index][1]
            #elif  cigar_info[cigar_index][0] in [4]:
                
                #seq_index_o = seq_index_o + cigar_info[cigar_index][1]

        # cut the bone index according to overall situation
        cal_bone_index = start_pos
        tem_seq_bone = ''
        seg_ment_bone_seq_xun = {}
        while cal_bone_index <= end_pos:
            if cal_bone_index not in remaining_key_:
                tem_seq_bone = tem_seq_bone + \
                    tem_seq[int(cal_bone_index-start_pos)]
            else:
                tem_seq_bone = tem_seq_bone + \
                    tem_seq[int(cal_bone_index-start_pos)]
                seg_ment_bone_seq_xun[cal_bone_index] = tem_seq_bone
                tem_seq_bone = ''

                #tem_seq_bone = ''

            cal_bone_index = cal_bone_index + 1
        if len(tem_seq_bone) > 0:
            seg_ment_bone_seq_xun[-1] = tem_seq_bone

        # makeup insertion segment
        finalized_ins_dic = {}
        for kkey in remaining_key_:
            if ins_inserted_dic[kkey] < 2:  # only one insertion
                if kkey in seq_specific_dic:
                    finalized_ins_dic[kkey] = seq_specific_dic[kkey]
                else:
                    finalized_ins_dic[kkey] = 'N'
                    #back_ground_ins = max(whole_ins_dic[ss_key_],key=len,default = '')
            else:
                if kkey in seq_specific_dic:
                    present_ins = seq_specific_dic[kkey]
                    tem_data_frame = pd.DataFrame(
                        dic_multiple_ali_ins[kkey])
                    refined_ins = (
                        tem_data_frame[tem_data_frame[0] == present_ins]).iloc[0, 1]
                    finalized_ins_dic[kkey] = refined_ins
                else:
                    finalized_ins_dic[kkey] = 'N'*ins_inserted_dic[kkey]

        # finished line for insertion and begin for make up the whole sequence
        final_string = ''
        for key_key in sorted(finalized_ins_dic):
            final_string = final_string + \
                seg_ment_bone_seq_xun[key_key] + finalized_ins_dic[key_key]

        final_string = final_string + seg_ment_bone_seq_xun[-1]

    else:
        final_string = seq_refined

    ###
    # begin mapping the mismatch and extracting insertion information
    final_string_code = string_2_code(final_string)
    ref_start = register_strat
    #final_dic = {}
    #cigar_o = []
    #for ss in range(len(final_string_code)):
    #    if final_string_code[ss] != final_ref_code[int(ss+ref_start)] and (ss+ref_start) not in collection_of_sites:
    #        collection_of_sites.append(ss+ref_start)

    #new_data_frame.loc[ii] = [register_strat,final_string_code]
    # based on final_string to infer

    #out_put_df.loc[ii] = [register_strat,
    #                      final_string]

    #collection_of_sites_array = np.array(sorted(collection_of_sites))
    return([register_strat, final_string])




def extract_s_e(s, e):
    return([s, e])


def Jaccard_pair8(read_index1, read_index2):
    # here we should refine the condition
    #threshold = 0
    #read_index1 = 2
    #read_index2 = 3
    #ref_code = finalized_ref_code_coll
    #coverage_df = reads_coverage
    #theta_r = 0.9
    #proba_list = labels_initial0_matrix
    #reads_infer_clu = reads_infer_clu
    global final_indicator_table
    start1 = final_full_result[read_index1][0]#final_indicator_table.iloc[read_index1, 1]
    end1 = final_full_result[read_index1][1]#final_indicator_table.iloc[read_index1, 2]
    start2 = final_full_result[read_index2][0]#final_indicator_table.iloc[read_index2, 1]
    end2 =  final_full_result[read_index2][1]#final_indicator_table.iloc[read_index2, 2]
    # read_flag1=reads_infer_clu.iloc[read_index1,1]
    # read_flag2=reads_infer_clu.iloc[read_index2,1]
    #l_min =threshold*(end1-start1+1)

    l_min = 30
    max_s = max([start1, start2])
    min_e = min([end1, end2])
    return_coe = 0

    if min_e - max_s > l_min:
        code_s = max_s
        code_e = (min_e+1)
        #mutation_list_select = mutation_list[np.where(
        #    (mutation_list > code_s) & ((mutation_list < code_e)))]
        #mutation_list_select_id = {
        #    int(xx/6): xx % 6 for xx in mutation_list_select}
        #mutation_list_select_id_unicque = set(mutation_list_select_id.keys())
        #if len(mutation_list_select) > 0:
        #if len(mutation_list_select_id_unicque) > 0:
            #cover_s = max_s
            #cover_e = min_e
            #sub_cover = coverage_df.iloc[int(cover_s):int(cover_e+1),:]
            #sub_whole_cover = sum(sub_cover.iloc[:,1])
            #mismatched_col_whole1 = final_indicator_table.iloc[read_index1, 3]
            #matched_col_whole1 = final_indicator_table.iloc[read_index1, 4]
            #mismatched_col_whole2 = final_indicator_table.iloc[read_index2, 3]
            #matched_col_whole2 = final_indicator_table.iloc[read_index2, 4]
            #whole1 = np.array(sorted(mismatched_col_whole1+matched_col_whole1))
            #whole2 = np.array(sorted(mismatched_col_whole2+matched_col_whole2))

            #selected_base1 = whole1[np.where((whole1>code_s)&((whole1<code_e)))]
        selected_base1_dic = final_full_result[read_index1][3]
        #selected_base2 = whole2[np.where((whole2>code_s)&((whole2<code_e)))]
        selected_base2_dic = final_full_result[read_index2][3]
        # for read_flag1 in range(len(proba_list_1)):
        #ref1 = ref_code[int(read_flag1)]

        Jaccard_coe = 0

        fenmu_ = 0
        for base_r in range(code_s,code_e):
        #for base_r in list(mutation_list_select_id_unicque):
            # if selected_base2_dic[base_r]!=5 and selected_base1_dic[base_r]!=5:
            fenmu_ = fenmu_ + 1
            if selected_base2_dic[base_r] == selected_base1_dic[base_r]:
                Jaccard_coe = Jaccard_coe + 1
                #Jaccard_coe = Jaccard_coe + 1/(cover_e-cover_s+1)

        # return_coe = np.exp((cover_1_c+cover_2_c)/(2*(min_e-max_s+1))*np.log(p_correct)+\
        #    (cover_1_w+cover_2_w)/(2*(min_e-max_s+1))*np.log(p_wrong+adjust_r))*Jaccard_coe
        if fenmu_ > 5:
            return_coe = Jaccard_coe/fenmu_
    return ([read_index1, read_index2, return_coe])



def Jaccard_pair9(read_index1, read_index2):
    # here we should refine the condition
    #threshold = 0.5
    #read_index1 = 0
    #read_index2 = 2
    #ref_code = finalized_ref_code_coll
    #coverage_df = reads_coverage
    #theta_r = 0.9
    #proba_list = labels_initial0_matrix
    #reads_infer_clu = reads_infer_clu
    global mutation_list, final_full_result
    start1 = final_full_result[read_index1][0]
    end1 = final_full_result[read_index1][1]
    start2 = final_full_result[read_index2][0]
    end2 = final_full_result[read_index2][1]
    # read_flag1=reads_infer_clu.iloc[read_index1,1]
    # read_flag2=reads_infer_clu.iloc[read_index2,1]
    #l_min =threshold*min((end1-start1+1),(end2-start2+1))
    #NN= 5
    l_min = 0
    max_s = max([start1, start2])
    min_e = min([end1, end2])
    num_base = 6
    return_coe = 0
    #adjust_r = pow(0.1,10)
    #p_correct = (theta_r*gamma_r+(1-gamma_r)*(1-theta_r)/(NN-1))
    #p_wrong = (1/((NN-1)**2)*(theta_r+gamma_r+NN*(1-theta_r*gamma_r)-2))
   # mutation_list = copy.deepcopy(new_keys)
    mutation_list = sorted(mutation_list)
    mutation_list = np.array(mutation_list)
    if min_e - max_s > l_min:
        code_s = max_s
        code_e = (min_e+1)
        mutation_list_select = mutation_list[np.where(
            (mutation_list > code_s) & ((mutation_list < code_e)))]
        #mutation_list_select_id = {
        #    int(xx/6): xx % 6 for xx in mutation_list_select}
        #mutation_list_select_id_unicque = set(mutation_list_select_id.keys())
        if len(mutation_list_select) > 0:
        #if len(mutation_list_select_id_unicque) > 0:
            #cover_s = max_s
            #cover_e = min_e
            #sub_cover = coverage_df.iloc[int(cover_s):int(cover_e+1),:]
            #sub_whole_cover = sum(sub_cover.iloc[:,1])
            #mismatched_col_whole1 = final_indicator_table.iloc[read_index1, 3]
            #matched_col_whole1 = final_indicator_table.iloc[read_index1, 4]
            #mismatched_col_whole2 = final_indicator_table.iloc[read_index2, 3]
            #matched_col_whole2 = final_indicator_table.iloc[read_index2, 4]
            #whole1 = np.array(sorted(mismatched_col_whole1+matched_col_whole1))
            #whole2 = np.array(sorted(mismatched_col_whole2+matched_col_whole2))

            #selected_base1 = whole1[np.where((whole1>code_s)&((whole1<code_e)))]
            selected_base1_dic = final_full_result[read_index1][3]
            #selected_base2 = whole2[np.where((whole2>code_s)&((whole2<code_e)))]
            selected_base2_dic = final_full_result[read_index2][3]
            # for read_flag1 in range(len(proba_list_1)):
            #ref1 = ref_code[int(read_flag1)]

            Jaccard_coe = 0

            fenmu_ = 0
            for base_r in list(mutation_list_select):
            #for base_r in list(mutation_list_select_id_unicque):
                # if selected_base2_dic[base_r]!=5 and selected_base1_dic[base_r]!=5:
                fenmu_ = fenmu_ + 1
                if selected_base2_dic[base_r] == selected_base1_dic[base_r]:
                    Jaccard_coe = Jaccard_coe + 1
                    #Jaccard_coe = Jaccard_coe + 1/(cover_e-cover_s+1)

            # return_coe = np.exp((cover_1_c+cover_2_c)/(2*(min_e-max_s+1))*np.log(p_correct)+\
            #    (cover_1_w+cover_2_w)/(2*(min_e-max_s+1))*np.log(p_wrong+adjust_r))*Jaccard_coe
            if fenmu_ > 5:
                return_coe = Jaccard_coe/fenmu_
    return ([read_index1, read_index2, return_coe])


def change_label_matrix(label_array, K_):
    #K_ = len(set(label_array))---
    #label_array = labels_initial0
    #K_ = 1
    new_label = pd.DataFrame(np.zeros([len(label_array), K_+1]))
    new_label.iloc[:, 0] = range(len(label_array))
    for ll in range(len(label_array)):
        new_label.iloc[ll, int(label_array[ll])+1] = 1
    return(new_label)


def converting_code_2_seq(code, seq_o):
    tem_seq = ''
    for code_unit in code:
        if code_unit == -1:
            tem_seq = tem_seq + 'N'
        else:
            tem_seq = tem_seq + seq_o[int(code_unit)]
    return(tem_seq)



def multiple_seq_alignment2_1(seq_collection_key):
    seq_collection = whole_insertion[seq_collection_key]
    if len(seq_collection)>1:
        recollection_fatsa = []
        recollection_fatsa = (SeqRecord(Seq(seq_collection[s]),id=str(s)) for s in range(len(seq_collection)))
        present_pid_string_o = str(mp.current_process())
        #print(present_pid_string)
        start_pos_string  = present_pid_string_o.index('ForkPoolWorker')
        end_pos_string =  present_pid_string_o.index('parent')
        present_pid_string = present_pid_string_o[start_pos_string:end_pos_string-2]
        #present_pid=str(os.getpid())
        SeqIO.write(recollection_fatsa, present_pid_string+"_msa.fasta", "fasta")
        os.system('mafft --auto ' + present_pid_string+"_msa.fasta > " + present_pid_string+  '_msa_result.fasta')
        #read_recollection_fatsa = open(present_pid+  '_msa.txt')
        os.system('samtools faidx ' + present_pid_string+"_msa_result.fasta")
        fasta_original_file = pysam.FastaFile(present_pid_string+"_msa_result.fasta")
        result_mapping = []
        #print(fasta_original_file.references)
        for line_per in range(len(seq_collection)):
            tem_seq_msa = fasta_original_file.fetch(str(line_per))
            line2=(tem_seq_msa.upper()).replace('-','N')
            result_mapping.append([seq_collection[line_per],line2])
                
    else:
         result_mapping = [seq_collection+seq_collection]
    
    return(result_mapping)




def multiple_seq_alignment2(seq_collection):
    #seq_collection = whole_insertion[6528]
    if len(seq_collection)>1:
        recollection_fatsa = []
        recollection_fatsa = (SeqRecord(Seq(seq_collection[s]),id=str(s)) for s in range(len(seq_collection)))
        present_pid_string = str(mp.current_process())
        print(present_pid_string)
        #start_pos_string  = present_pid_string.index('-')
        #end_pos_string = present_pid_string[start_pos_string:].index("\'") 
        present_pid=str(os.getpid())
        SeqIO.write(recollection_fatsa, present_pid+"_msa.fasta", "fasta")
        os.system('mafft --auto ' + present_pid+"_msa.fasta > " + present_pid+  '_msa_result.fasta')
        #read_recollection_fatsa = open(present_pid+  '_msa.txt')
        os.system('samtools faidx ' + present_pid+"_msa_result.fasta")
        fasta_original_file = pysam.FastaFile(present_pid+"_msa_result.fasta")
        result_mapping = []
        print(fasta_original_file.references)
        for line_per in range(len(seq_collection)):
            tem_seq_msa = fasta_original_file.fetch(str(line_per))
            line2=(tem_seq_msa.upper()).replace('-','N')
            result_mapping.append([seq_collection[line_per],line2])
                
    else:
         result_mapping = [seq_collection+seq_collection]
    
    return(result_mapping)


def multiple_seq_alignment(seq_collection):
    #seq_collection = whole_insertion[6528]
    recollection = []
    for seq_index in range(len(seq_collection)):
        tem_seq = seqbio.NucleotideSequence(seq_collection[seq_index])
        recollection.append(tem_seq)

    try:
        alignment, order, guide_tree, distance_matrix = align.align_multiple(
            recollection,
            matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
            gap_penalty=(-10, -0.5),
            terminal_penalty=False
        )
        alignment_trace = alignment.trace
        length_compare, num_col = np.shape(alignment_trace)
        result_mapping = []
        for col_index in range(num_col):
            #tem_result_mapping = [recollection[col_index]]
            coverted_seq = converting_code_2_seq(
                alignment_trace[:, col_index], recollection[col_index])
            # tem_result_mapping.append(coverted_seq)
            result_mapping.append(
                [seq_collection[col_index], coverted_seq])
    except (ValueError, ZeroDivisionError):
        max_len = len(seq_collection[0])
        for col_index in range(len(seq_collection)):
            if len(seq_collection[col_index]) > max_len:
                max_len = len(seq_collection[col_index])
        result_mapping = []
        for col_index in range(len(seq_collection)):
            if len(seq_collection[col_index]) < max_len:
                refined_ins_ = seq_collection[col_index] + \
                    'N'*(max_len-len(seq_collection[col_index]))

            else:
                refined_ins_ = seq_collection[col_index]
            result_mapping.append([seq_collection[col_index], refined_ins_])

    return(result_mapping)

# flag is develop to discriminate original and refined
def select_susscipious_sites(selected_dic, final_ref_code, collection_of_sites):

    #selected_dic = select_dic2
    key_collection = pd.DataFrame(selected_dic.keys())
    collection_sites = []
    k_mer = len(key_collection.iloc[0, :])-1  # excluding the first site
    collection_of_sites = sorted(collection_of_sites)
    for index_key in range(len(key_collection)):
        #index_key= 1
        pos = key_collection.iloc[index_key, 0]
        index_pos = collection_of_sites.index(pos)
        for neighbour_index in range(k_mer):
            #neighbour_index = 0
            tem_pos = int(collection_of_sites[int(index_pos)+neighbour_index])
            # print(tem_pos)
            if key_collection.iloc[index_key, 1+neighbour_index] != final_ref_code[tem_pos]:
                if tem_pos not in collection_sites:
                    collection_sites.append(tem_pos)

    return(collection_sites)


def true_variation_label(ins_inserted_dic, true_sites_list):

    true_sites_var = {}
    keys_array = np.array(sorted(ins_inserted_dic.keys()))
    for sites in true_sites_list:
        #sites = 34
        index_ = np.where(keys_array < sites)
        selected_sites = keys_array[index_]
        if keys_array[index_].size:
            tem_right_shift = 0
            for ii in range(len(selected_sites)):
                tem_right_shift = tem_right_shift + \
                    ins_inserted_dic[selected_sites[ii]]
            new_pos = sites + tem_right_shift
            true_sites_var[sites] = new_pos
        else:
            true_sites_var[sites] = sites

    return(true_sites_var)

def trans_(ins_inserted_dic,true_sites_list):
    transferred_dic = {}
    ins_inserted_dic_keys = np.array(sorted(list(ins_inserted_dic.keys())))
    for site in true_sites_list:
        tem_sum = 0
        partial_list = ins_inserted_dic_keys[ins_inserted_dic_keys<site]
        for ins_id in partial_list:
            tem_sum = tem_sum + ins_inserted_dic[ins_id]
        transferred_dic[site]= site+  tem_sum
    return(transferred_dic)

def reverse_trans(ins_inserted_dic, true_sites_list):
    original_sites_var = {}
    #ins_inserted_dic= {0:2,3:2}
    #true_sites_list = list(range(9))
    keys_list = list(sorted(ins_inserted_dic.keys()))

    for sites in true_sites_list:
        #sites = 34
        #print(sites)
        ins_num = 0
        #sites = 6
        if len(keys_list)>0:
            if sites <= keys_list[0]:
                original_sites_var[sites] = [sites, -1]
            else:
                tem_distance_sum = keys_list[0]
                if tem_distance_sum+ins_inserted_dic[keys_list[0]] >= sites:
                    small_ins = sites-tem_distance_sum-1  # 0 index
                    original_sites_var[sites] = [keys_list[0], small_ins]
                else:
    
                    tem_distance_sum = tem_distance_sum + \
                        ins_inserted_dic[keys_list[0]]
                    ins_num = ins_inserted_dic[keys_list[0]]
                    loopll = 1
                    while tem_distance_sum < sites:
                        if loopll < len(keys_list):
                            if tem_distance_sum+keys_list[loopll]-keys_list[loopll-1] >= sites:
    
                                small_bone = sites - ins_num
                                original_sites_var[sites] = [small_bone, -1]
                                tem_distance_sum = sites
                                end_flag = 0
                                break
                            else:
                                tem_distance_sum = tem_distance_sum + \
                                    keys_list[loopll]-keys_list[loopll-1]
                                if tem_distance_sum+ins_inserted_dic[keys_list[loopll]] >= sites:
                                    small_ins = sites-tem_distance_sum-1  # 0 index
                                    original_sites_var[sites] = [
                                        keys_list[loopll], small_ins]
                                    end_flag = 0
                                    break
                                else:
                                    tem_distance_sum = tem_distance_sum + \
                                        ins_inserted_dic[keys_list[loopll]]
                                    ins_num = ins_num + \
                                        ins_inserted_dic[keys_list[loopll]]
                                    loopll = loopll + 1
    
                        else:
                            end_flag = 1
                            break
                    if end_flag:
                        small_bone = sites - ins_num
                        original_sites_var[sites] = [small_bone, -1]
        else:
            original_sites_var[sites] = [sites, -1]
    return(original_sites_var)


def reverse_trans2(ins_inserted_dic_df, true_sites_list,full_flag):
    original_sites_var = []
    #ins_inserted_dic= {0:2,3:2}
    #true_sites_list = list(range(9))
    #keys_list = list(sorted(ins_inserted_dic.keys()))
    if len(ins_inserted_dic_df)>0:
        if full_flag:
            tem_infer_var_site = copy.deepcopy(full_var_site)
        else:
            tem_infer_var_site = copy.deepcopy(infer_var_site)
        for sites in true_sites_list:
            ins_inserted_dic_df_sub = ins_inserted_dic_df[ins_inserted_dic_df.iloc[:,2]<=sites]
            tem_reduce_num = 0
            if len(ins_inserted_dic_df_sub)>0:
                if sites==ins_inserted_dic_df_sub.iloc[-1,2]:
                    #tem_reduce_num= sum(ins_inserted_dic_df_sub.iloc[:,1]) -1
                    reverse_original_site = ins_inserted_dic_df_sub.iloc[-1,0]
                    
                    if sites in tem_infer_var_site:
                        original_sites_var.append([sites,reverse_original_site,ins_inserted_dic_df_sub.iloc[-1,1]-1,tem_infer_var_site.index(sites)])
                    else:
                        original_sites_var.append([sites,reverse_original_site,ins_inserted_dic_df_sub.iloc[-1,1]-1,-2])
                else:
                    if len(ins_inserted_dic_df_sub)<len(ins_inserted_dic_df):## this is we have the last inserted pos
                        next_id_hang = len(ins_inserted_dic_df_sub)
                        tem_next_end_bone = ins_inserted_dic_df.iloc[next_id_hang,2] - ins_inserted_dic_df.iloc[next_id_hang,1]
                        if tem_next_end_bone>=sites:
                            reverse_original_site = ins_inserted_dic_df_sub.iloc[-1,0] + sites - ins_inserted_dic_df_sub.iloc[-1,2]
                            if sites in tem_infer_var_site:
                                original_sites_var.append([sites,reverse_original_site,-1,tem_infer_var_site.index(sites)])
                            else:
                                original_sites_var.append([sites,reverse_original_site,-1,-2])
                        else:
                            reverse_original_site = ins_inserted_dic_df.iloc[next_id_hang,0]
                            if sites in tem_infer_var_site:
                                original_sites_var.append([sites,reverse_original_site,sites-tem_next_end_bone-1,tem_infer_var_site.index(sites)])
                            else:
                                original_sites_var.append([sites,reverse_original_site,sites-tem_next_end_bone-1,-2])
                    else:
                        reverse_original_site = sites - ins_inserted_dic_df.iloc[-1,2] + ins_inserted_dic_df.iloc[-1,0]
                        if sites in tem_infer_var_site:
                            original_sites_var.append( [sites,reverse_original_site,-1,tem_infer_var_site.index(sites)])
                        else:
                            original_sites_var.append([sites,reverse_original_site,-1,-2])
                       
            else:
                if sites in tem_infer_var_site:
                    original_sites_var.append([sites,sites,-1,tem_infer_var_site.index(sites)])
                else:
                    original_sites_var.append([sites,sites,-1,-2])
    else:
        
        for sites in true_sites_list:
            original_sites_var.append([sites,sites,-1,full_var_site.index(sites)])

    return(original_sites_var)


def reverse_trans2_per(sites):
    global ins_inserted_dic_df,infer_var_site
    original_sites_var = []
    #ins_inserted_dic= {0:2,3:2}
    #true_sites_list = list(range(9))
    #keys_list = list(sorted(ins_inserted_dic.keys()))
    #sites = 647949
    #for sites in true_sites_list:
    if len(ins_inserted_dic_df)>0:
        ins_inserted_dic_df_sub = ins_inserted_dic_df[ins_inserted_dic_df.iloc[:,2]<=sites]
        #tem_reduce_num = 0
        if len(ins_inserted_dic_df_sub)>0:
            if sites==ins_inserted_dic_df_sub.iloc[-1,2]:
                #tem_reduce_num= sum(ins_inserted_dic_df_sub.iloc[:,1]) -1
                reverse_original_site = ins_inserted_dic_df_sub.iloc[-1,0]
                
                if sites in infer_var_site:
                    original_sites_var=[sites,reverse_original_site,ins_inserted_dic_df_sub.iloc[-1,1]-1,infer_var_site.index(sites)]
                else:
                    original_sites_var=[sites,reverse_original_site,ins_inserted_dic_df_sub.iloc[-1,1]-1,-2]
            else:
                if len(ins_inserted_dic_df_sub)<len(ins_inserted_dic_df):## this is we have the last inserted pos
                    next_id_hang = len(ins_inserted_dic_df_sub)
                    tem_next_end_bone = ins_inserted_dic_df.iloc[next_id_hang,2] - ins_inserted_dic_df.iloc[next_id_hang,1]
                    if tem_next_end_bone>=sites:
                        reverse_original_site = ins_inserted_dic_df_sub.iloc[-1,0] + sites - ins_inserted_dic_df_sub.iloc[-1,2]
                        if sites in infer_var_site:
                            original_sites_var=[sites,reverse_original_site,-1,infer_var_site.index(sites)]
                        else:
                            original_sites_var=[sites,reverse_original_site,-1,-2]
                    else:
                        reverse_original_site = ins_inserted_dic_df.iloc[next_id_hang,0]
                        if sites in infer_var_site:
                            original_sites_var=[sites,reverse_original_site,sites-tem_next_end_bone-1,infer_var_site.index(sites)]
                        else:
                            original_sites_var=[sites,reverse_original_site,sites-tem_next_end_bone-1,-2]
                else:
                    reverse_original_site = sites - ins_inserted_dic_df.iloc[-1,2] + ins_inserted_dic_df.iloc[-1,0]
                    if sites in infer_var_site:
                        original_sites_var = [sites,reverse_original_site,-1,infer_var_site.index(sites)]
                    else:
                        original_sites_var = [sites,reverse_original_site,-1,-2]
            
                
                   
        else:
            if sites in infer_var_site:
                reverse_original_site = ins_inserted_dic_df.iloc[0,0]
                original_sites_var=[sites,reverse_original_site,sites-reverse_original_site-1,infer_var_site.index(sites)]
            else:
                original_sites_var=[sites,sites,-1,-2]
    else:
        original_sites_var=[sites,sites,-1,-2]
    return(original_sites_var)



def trans_to_dic(final_indicator_table, specific_var):
    len_read = len(final_indicator_table)
    feed_back = []
    for rr in range(len_read):
        tem_dic = {}
        start_ = final_indicator_table.iloc[rr, 1]
        end_ = final_indicator_table.iloc[rr, 2]
        correct_one = final_indicator_table.iloc[rr, 3]
        wrong_one = final_indicator_table.iloc[rr, 4]
        target_var = correct_one+wrong_one
        if len(target_var) > 0:

            full_list = correct_one + wrong_one
            tem_dic = {int(xx/6): (xx % 6) for xx in full_list}

            feed_back.append([start_, end_, tem_dic])
        else:
            feed_back.append([0, 0, tem_dic])

    return(feed_back)


def select_read_blocks(start_block, end_block):
    #start_block = list_block_nei[0][0]
    #end_block = list_block[0][1]
    global ele_list,rselct_result
    tem_return_res = []
    #reads_info = rselct_result
    #ele_list_tem = ele_list
    for read_label_xun in ele_list:
        read_start = rselct_result[read_label_xun][0]
        read_end = rselct_result[read_label_xun][1]

        if not (start_block > read_end) and not (read_start >= end_block):
            tem_dic = rselct_result[read_label_xun][2]
            select_dic = {k: v for k, v in tem_dic.items() if (
                (k >= start_block) and (k <= end_block))}
            tem_return_res.append([read_start, read_end, select_dic])
    return(tem_return_res)

def select_read_blocks2(start_block_id, end_block_id):
    #global ele_list,rselct_result,infer_var_site,final_full_result
    start_block = infer_var_site[int(start_block_id)]
    end_block = infer_var_site[int(end_block_id)]
    #start_block = list_block_nei[0][0]
    #end_block = list_block[0][1]
    #global ele_list,rselct_result
    tem_return_res = []
    #reads_info = rselct_result
    #ele_list_tem = ele_list
    for read_label_xun in ele_list:## abs collection
        read_start = final_full_result[read_label_xun][0]
        read_end = final_full_result[read_label_xun][1]

        if not (start_block > read_end) and not (read_start >= end_block):
            tem_dic = final_full_result[read_label_xun][3]
            select_dic = {k: v for k, v in tem_dic.items() if (
                (k >= start_block) and (k <= end_block))}
            tem_return_res.append([read_start, read_end, select_dic])
    return(tem_return_res)


def select_ref_block(start_block, end_block, kk):
    #tem_return_res = []
    #global Best_ref
    select_dic = {k: v for k, v in Best_ref[kk].items() if (
        (k >= start_block) and (k <= end_block))}
    return(select_dic)

def select_ref_block2(start_block_id, end_block_id, kk):
    #tem_return_res = []
    #start_block_id = 10413
    #end_block_id = 10416
    #global Best_ref,infer_var_site,Best_ref_ins_full
#for block_fullid in range(len(remaining_dependent)):
#    start_block_id = remaining_dependent.iloc[block_fullid,0]
#    end_block_id = remaining_dependent.iloc[block_fullid,1]
    start_block_id = int(start_block_id)
    end_block_id = int(end_block_id)
    select_dic = {}
    for ll_tem in range(start_block_id,end_block_id+1):
        original_infor = reverse_trans2_per(infer_var_site[ll_tem])
        original_infor_pos = original_infor[1]
        original_infor_state = original_infor[2]
        if original_infor_state <0:
            inferred_base = Best_ref[kk][original_infor_pos]
        else:
            full_infer_ins_string = Best_ref_ins_full[kk][original_infor_pos]
            inferred_base = full_infer_ins_string[original_infor_state]
        select_dic[infer_var_site[ll_tem]] = int(inferred_base)
    #start_block = infer_var_site[start_block_id]
    #end_block = infer_var_site[end_block_id]
    #select_dic = {k: v for k, v in Best_ref[kk].items() if (
    #    (k >= start_block) and (k <= end_block))}
    return(select_dic)


def pre_max_(site_id, base_number):
    
    #global remaining_independent,ref_code_standard
#for index in range(len(remaining_independent)):
#    site_id = index#remaining_independent[index]
#    base_number = present_clu_dis2[kk][(index), :]
    base_number_list = list(base_number)
    max_ele = max(base_number_list)
    if max_ele<1:
        standard_bone_original_list = reverse_trans2_per(remaining_independent[site_id])
        standard_bone_original_pos = standard_bone_original_list[1]
        standard_bone_original_state = standard_bone_original_list[2]
        if standard_bone_original_state<0:
            max_ele_id = ref_code_standard[standard_bone_original_pos]
        else:
            max_ele_id = 5
    else:
        max_ele_id = base_number_list.index(max_ele)
    return([remaining_independent[site_id], max_ele_id])


def signle_sites_infer_max2(site_id, base_number, standard_ref_code_site,kk):## this is a problem
    #kk = 0
    #site_id = 70
    #base_number = present_clu_disdf[kk][70,:]
    #base_number = present_clu_dis2cdf.iloc[9810,:]
    #base_number = np.array([0,1,0,3,0])
    #standard_ref_code_site = standard_ref_code[9810]
    #kk = 0
    #gamma_ = [gamma1[kk], gamma2[kk], gamma3[kk], gamma4[kk]]
    #theta_ = [pcorrect, p_mis_loc, pins1_loc,
    #          pdel1_loc, lamuna_ins, lamuna_del]
    #gamma_ = gamma[kk]
    #for  ll in remaining_independent:
        #base_number = present_clu_dis2[kk][ll,:]
        #standard_ref_code_site = standard_ref_code[ll]
        #gamma_ = [0.98838986,0.01145738,0.00015277,0]
        #theta_ = copy.deepcopy(parameters_f[4:])
        global gamma1, gamma2, gamma3, gamma4,theta_,del_beta,del_p,ins_beta,ins_p
        gamma_c =  gamma1[kk]
        gamma_mis = gamma2[kk]
        gamma_del = gamma3[kk]
        gamma_ins = gamma4[kk]
        gamma_del_p = del_p[kk]
        gamma_del_beta = del_beta[kk]
        gamma_ins_p = ins_p[kk]
        gamma_ins_beta = ins_beta[kk]
        pcorrect = theta_[0]
        p_mis_loc = theta_[1]
        pins1_loc = theta_[2]
        pdel1_loc = theta_[3]
        lamuna_ins = theta_[4]
        lamuna_del = theta_[5]
        #potential_probability = np.zeros(5)
        main_ref = standard_ref_code_site
        corretion_eps_from_zero = 0.000001
        #base_number = [0,0,0,2,0,0]
        #main_ref = 2
        if main_ref > 4:  # ins in the bone
            
            base_number_ins = sum(base_number[:4])
            if base_number_ins>1:## ins only consider the maximum one
                max_num_bei= max(base_number[:4])
                max_num_bei_id = int(list(base_number[:4]).index(max_num_bei))
                prb = np.zeros(2)
                #for main_base_xun in range(4):
                same_base = base_number[max_num_bei_id]
                not_same_base = sum(base_number[:4]) - same_base
                del_num = base_number[5]
                prb[0] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero) +\
                    (same_base+not_same_base)*np.log(1-pdel1_loc+corretion_eps_from_zero) + del_num*np.log(pdel1_loc+corretion_eps_from_zero) +\
                    del_num*(possion_log_probability(0,lamuna_ins)) + \
                    np.log(gamma_ins+corretion_eps_from_zero)+expotential_log_log_probability(0,gamma_ins_p,gamma_ins_beta)
                prb[1] = (del_num)*np.log(1-pins1_loc+corretion_eps_from_zero) +\
                    sum(base_number[:4])*(np.log(pins1_loc+corretion_eps_from_zero) +
                                          possion_log_probability(0,lamuna_ins))+np.log(1-gamma_ins+corretion_eps_from_zero)
                #tem_probability = copy.deepcopy(potential_probability/tem_pri)
                #tem_probability = tem_probability - max(tem_probability)
                #tem_probability_main = np.exp(tem_probability)
                #tem_probability_main = tem_probability_main/sum(tem_probability_main)
                #post_sampling = np.random.multinomial(n= 1,pvals = tem_probability_main)
                #max_ID = list(post_sampling).index(1)
                #    n=1, pvals=[6/8, 1/8, 1/16, 1/32, 1/32])
                #l_i = list(l_i_vec).index(1)+1
                max_ID = list(prb).index(max(prb))
                if max_ID == 1:
                    max_ID = 5
                else:
                    max_ID = max_num_bei_id
            else:
                max_ID = 5
        else:
            if sum(base_number) > 0:
                ### ref_case
                del_num = base_number[4]
                if del_num>0:
                    prb = np.zeros(3)
                    same_base = base_number[main_ref]
                    not_same_base = sum(base_number[:4]) - same_base
                    
                    prb[0] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero) +\
                        (same_base+not_same_base)*np.log(1-pdel1_loc+corretion_eps_from_zero) + del_num*np.log(pdel1_loc+corretion_eps_from_zero) + np.log(1-gamma_del+corretion_eps_from_zero)\
                        +del_num*possion_log_probability(0,lamuna_del)+np.log(gamma_c-corretion_eps_from_zero)
                    prb[1] = (del_num)*np.log(1-pins1_loc+corretion_eps_from_zero) +\
                        sum(base_number[:4])*(np.log(pins1_loc+corretion_eps_from_zero) +
                                             possion_log_probability(0,lamuna_ins))+np.log(gamma_del+corretion_eps_from_zero)+ expotential_log_log_probability(0,gamma_del_p,gamma_del_beta)
                    remain_base = list(base_number[:4])
                    del remain_base[main_ref]
                    base_index = list(range(4))
                    del base_index[main_ref]
                    
                    if max(remain_base)>1:
                        id_max_bei = base_index[int(list(remain_base).index(max(remain_base)))]
                        same_base = base_number[int(id_max_bei)]
                        not_same_base = sum(base_number[:4]) - same_base
                        
                        prb[2] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero) +\
                            (same_base+not_same_base)*np.log(1-pdel1_loc+corretion_eps_from_zero) + del_num*np.log(pdel1_loc+corretion_eps_from_zero)  + np.log(1-gamma_del+corretion_eps_from_zero)\
                            + del_num*possion_log_probability(0,lamuna_del)+np.log(gamma_mis+corretion_eps_from_zero)
                    else:
                        prb[2] = float('-inf')
                    max_ID = list(prb).index(max(prb))
                    if max_ID==0:
                        max_ID = main_ref
                    elif max_ID==1:
                        max_ID = 4
                    else:
                        max_ID = id_max_bei
                else:
                    prb = np.zeros(2)
                    same_base = base_number[main_ref]
                    not_same_base = sum(base_number[:4]) - same_base
                    
                    prb[0] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero)
                       
                   
                    remain_base = list(base_number[:4])
                    del remain_base[main_ref]
                    base_index = list(range(4))
                    del base_index[main_ref]
                    
                    if max(remain_base)>1:
                        id_max_bei = base_index[int(list(remain_base).index(max(remain_base)))]
                        same_base = base_number[int(id_max_bei)]
                        not_same_base = sum(base_number[:4]) - same_base
                        
                        prb[1] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero) 
                            
                    else:
                        prb[1] = float('-inf')
                    max_ID = list(prb).index(max(prb))
                    if max_ID==0:
                        max_ID = main_ref

                    else:
                        max_ID = id_max_bei
    
            else:
                max_ID = main_ref
        #max_ID = list(potential_probability).index(max(potential_probability))
        return([site_id, max_ID])



def signle_sites_infer_max3(site_id, base_number, standard_ref_code_site,kk):## this is a problem
    #kk = 0
    #site_id = 70
    #base_number = present_clu_disdf[kk][70,:]
    #base_number = present_clu_dis2cdf.iloc[9810,:]
    #base_number = np.array([0,1,0,3,0])
    #standard_ref_code_site = standard_ref_code[9810]
    #kk = 0
    #gamma_ = [gamma1[kk], gamma2[kk], gamma3[kk], gamma4[kk]]
    #theta_ = [pcorrect, p_mis_loc, pins1_loc,
    #          pdel1_loc, lamuna_ins, lamuna_del]
    #gamma_ = gamma[kk]
    #for  ll in remaining_independent:
    #    base_number = present_clu_dis2[kk][ll,:]
    #    standard_ref_code_site = standard_ref_code[infer_var_site[ll]]
        #gamma_ = [0.98838986,0.01145738,0.00015277,0]
        #theta_ = copy.deepcopy(parameters_f[4:])
        #global gamma1, gamma2, gamma3, gamma4,theta_,\
        #    del_beta,del_p,ins_beta,ins_p,infer_var_site
        gamma_c =  gamma1[kk]
        gamma_mis = gamma2[kk]
        gamma_del = gamma3[kk]
        gamma_ins = gamma4[kk]
        gamma_del_p = del_p[kk]
        gamma_del_beta = del_beta[kk]
        gamma_ins_p = ins_p[kk]
        gamma_ins_beta = ins_beta[kk]
        pcorrect = theta_[0]
        p_mis_loc = theta_[1]
        pins1_loc = theta_[2]
        pdel1_loc = theta_[3]
        lamuna_ins = theta_[4]
        lamuna_del = theta_[5]
        #potential_probability = np.zeros(5)
        main_ref = standard_ref_code_site
        corretion_eps_from_zero = 0.000001
        #base_number = [0,0,0,2,0,0]
        #main_ref = 2
        if main_ref > 4:  # ins in the bone
            
            base_number_ins = sum(base_number[:4])
            if base_number_ins>1:## ins only consider the maximum one
                max_num_bei= max(base_number[:4])
                max_num_bei_id = int(list(base_number[:4]).index(max_num_bei))
                prb = np.zeros(2)
                #for main_base_xun in range(4):
                same_base = base_number[max_num_bei_id]
                not_same_base = sum(base_number[:4]) - same_base
                del_num = base_number[5]
                prb[0] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero) +\
                    (same_base+not_same_base)*np.log(1-pdel1_loc+corretion_eps_from_zero) + del_num*np.log(pdel1_loc+corretion_eps_from_zero) +\
                    del_num*(possion_log_probability(0,lamuna_ins)) + \
                    np.log(gamma_ins+corretion_eps_from_zero)+expotential_log_log_probability(0,gamma_ins_p,gamma_ins_beta)
                prb[1] = (del_num)*np.log(1-pins1_loc+corretion_eps_from_zero) +\
                    sum(base_number[:4])*(np.log(pins1_loc+corretion_eps_from_zero) +
                                          possion_log_probability(0,lamuna_ins))+np.log(1-gamma_ins+corretion_eps_from_zero)
                #tem_probability = copy.deepcopy(potential_probability/tem_pri)
                #tem_probability = tem_probability - max(tem_probability)
                #tem_probability_main = np.exp(tem_probability)
                #tem_probability_main = tem_probability_main/sum(tem_probability_main)
                #post_sampling = np.random.multinomial(n= 1,pvals = tem_probability_main)
                #max_ID = list(post_sampling).index(1)
                #    n=1, pvals=[6/8, 1/8, 1/16, 1/32, 1/32])
                #l_i = list(l_i_vec).index(1)+1
                max_ID = list(prb).index(max(prb))
                if max_ID == 1:
                    max_ID = 5
                else:
                    max_ID = max_num_bei_id
            else:
                max_ID = 5
        else:
            if sum(base_number) > 0:
                ### ref_case
                del_num = base_number[4]
                if del_num>0:
                    prb = np.zeros(3)
                    same_base = base_number[main_ref]
                    not_same_base = sum(base_number[:4]) - same_base
                    
                    prb[0] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero) +\
                        (same_base+not_same_base)*np.log(1-pdel1_loc+corretion_eps_from_zero) + del_num*np.log(pdel1_loc+corretion_eps_from_zero) + np.log(1-gamma_del+corretion_eps_from_zero)\
                        +del_num*possion_log_probability(0,lamuna_del)+np.log(gamma_c-corretion_eps_from_zero)
                    prb[1] = (del_num)*np.log(1-pins1_loc+corretion_eps_from_zero) +\
                        sum(base_number[:4])*(np.log(pins1_loc+corretion_eps_from_zero) +
                                             possion_log_probability(0,lamuna_ins))+np.log(gamma_del+corretion_eps_from_zero)+ expotential_log_log_probability(0,gamma_del_p,gamma_del_beta)
                    remain_base = list(base_number[:4])
                    del remain_base[main_ref]
                    base_index = list(range(4))
                    del base_index[main_ref]
                    
                    if max(remain_base)>1:
                        id_max_bei = base_index[int(list(remain_base).index(max(remain_base)))]
                        same_base = base_number[int(id_max_bei)]
                        not_same_base = sum(base_number[:4]) - same_base
                        
                        prb[2] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero) +\
                            (same_base+not_same_base)*np.log(1-pdel1_loc+corretion_eps_from_zero) + del_num*np.log(pdel1_loc+corretion_eps_from_zero)  + np.log(1-gamma_del+corretion_eps_from_zero)\
                            + del_num*possion_log_probability(0,lamuna_del)+np.log(gamma_mis+corretion_eps_from_zero)
                    else:
                        prb[2] = float('-inf')
                    max_ID = list(prb).index(max(prb))
                    if max_ID==0:
                        max_ID = main_ref
                    elif max_ID==1:
                        max_ID = 4
                    else:
                        max_ID = id_max_bei
                else:
                    prb = np.zeros(2)
                    same_base = base_number[main_ref]
                    not_same_base = sum(base_number[:4]) - same_base
                    
                    prb[0] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero)
                       
                   
                    remain_base = list(base_number[:4])
                    del remain_base[main_ref]
                    base_index = list(range(4))
                    del base_index[main_ref]
                    
                    if max(remain_base)>1:
                        id_max_bei = base_index[int(list(remain_base).index(max(remain_base)))]
                        same_base = base_number[int(id_max_bei)]
                        not_same_base = sum(base_number[:4]) - same_base
                        
                        prb[1] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero) 
                            
                    else:
                        prb[1] = float('-inf')
                    max_ID = list(prb).index(max(prb))
                    if max_ID==0:
                        max_ID = main_ref

                    else:
                        max_ID = id_max_bei
    
            else:
                max_ID = main_ref
        #max_ID = list(potential_probability).index(max(potential_probability))
        return([infer_var_site[site_id], max_ID])


def signle_sites_infer_max4(site_id, base_number, standard_ref_code_site,kk):## this is a problem
    #kk = 0
    #site_id = 70
    #base_number = present_clu_disdf[kk][70,:]
    #base_number = present_clu_dis2cdf.iloc[9810,:]
    #base_number = np.array([0,1,0,3,0])
    #standard_ref_code_site = standard_ref_code[9810]
    #kk = 0
    #gamma_ = [gamma1[kk], gamma2[kk], gamma3[kk], gamma4[kk]]
    #theta_ = [pcorrect, p_mis_loc, pins1_loc,
    #          pdel1_loc, lamuna_ins, lamuna_del]
    #gamma_ = gamma[kk]
    #for  ll in remaining_independent:
        #base_number = present_clu_dis2[kk][ll,:]
        #standard_ref_code_site = standard_ref_code[infer_var_site[ll]]
        #gamma_ = [0.98838986,0.01145738,0.00015277,0]
        #theta_ = copy.deepcopy(parameters_f[4:])
        #global gamma1, gamma2, gamma3, gamma4,theta_,\
        #    del_beta,del_p,ins_beta,ins_p,infer_var_site
        gamma_c =  gamma1[kk]
        gamma_mis = gamma2[kk]
        gamma_del = gamma3[kk]
        gamma_ins = gamma4[kk]
        gamma_del_p = del_p[kk]
        gamma_del_beta = del_beta[kk]
        gamma_ins_p = ins_p[kk]
        gamma_ins_beta = ins_beta[kk]
        pcorrect = theta_[0]
        p_mis_loc = theta_[1]
        pins1_loc = theta_[2]
        pdel1_loc = theta_[3]
        lamuna_ins = theta_[4]
        lamuna_del = theta_[5]
        #potential_probability = np.zeros(5)
        main_ref = standard_ref_code_site
        corretion_eps_from_zero = 0.000001
        #base_number = [0,0,0,2,0,0]
        #main_ref = 2
        if main_ref > 4:  # ins in the bone
            
            base_number_ins = sum(base_number[:4])
            if base_number_ins>1:## ins only consider the maximum one
                max_num_bei= max(base_number[:4])
                max_num_bei_id = int(list(base_number[:4]).index(max_num_bei))
                prb = np.zeros(2)
                #for main_base_xun in range(4):
                same_base = base_number[max_num_bei_id]
                not_same_base = sum(base_number[:4]) - same_base
                del_num = base_number[5]
                prb[0] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero) +\
                    (same_base+not_same_base)*np.log(1-pdel1_loc+corretion_eps_from_zero) + del_num*np.log(pdel1_loc+corretion_eps_from_zero) +\
                    del_num*(possion_log_probability(0,lamuna_ins)) + \
                    np.log(gamma_ins+corretion_eps_from_zero)+expotential_log_log_probability(0,gamma_ins_p,gamma_ins_beta)
                prb[1] = (del_num)*np.log(1-pins1_loc+corretion_eps_from_zero) +\
                    sum(base_number[:4])*(np.log(pins1_loc+corretion_eps_from_zero) +
                                          possion_log_probability(0,lamuna_ins))+np.log(1-gamma_ins+corretion_eps_from_zero)
                #tem_probability = copy.deepcopy(potential_probability/tem_pri)
                #tem_probability = tem_probability - max(tem_probability)
                #tem_probability_main = np.exp(tem_probability)
                #tem_probability_main = tem_probability_main/sum(tem_probability_main)
                #post_sampling = np.random.multinomial(n= 1,pvals = tem_probability_main)
                #max_ID = list(post_sampling).index(1)
                #    n=1, pvals=[6/8, 1/8, 1/16, 1/32, 1/32])
                #l_i = list(l_i_vec).index(1)+1
                max_ID = list(prb).index(max(prb))
                if max_ID == 1:
                    max_ID = 5
                else:
                    max_ID = max_num_bei_id
            else:
                max_ID = 5
        else:
            if sum(base_number) > 0:
                ### ref_case
                del_num = base_number[4]
                if del_num>0:
                    prb = np.zeros(3)
                    same_base = base_number[main_ref]
                    not_same_base = sum(base_number[:4]) - same_base
                    
                    prb[0] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero) +\
                        (same_base+not_same_base)*np.log(1-pdel1_loc+corretion_eps_from_zero) + del_num*np.log(pdel1_loc+corretion_eps_from_zero) + np.log(1-gamma_del+corretion_eps_from_zero)\
                        +del_num*possion_log_probability(0,lamuna_del)+np.log(gamma_c-corretion_eps_from_zero)
                    prb[1] = (del_num)*np.log(1-pins1_loc+corretion_eps_from_zero) +\
                        sum(base_number[:4])*(np.log(pins1_loc+corretion_eps_from_zero) +
                                             possion_log_probability(0,lamuna_ins))+np.log(gamma_del+corretion_eps_from_zero)+ expotential_log_log_probability(0,gamma_del_p,gamma_del_beta)
                    remain_base = list(base_number[:4])
                    del remain_base[main_ref]
                    base_index = list(range(4))
                    del base_index[main_ref]
                    
                    if max(remain_base)>1:
                        id_max_bei = base_index[int(list(remain_base).index(max(remain_base)))]
                        same_base = base_number[int(id_max_bei)]
                        not_same_base = sum(base_number[:4]) - same_base
                        
                        prb[2] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero) +\
                            (same_base+not_same_base)*np.log(1-pdel1_loc+corretion_eps_from_zero) + del_num*np.log(pdel1_loc+corretion_eps_from_zero)  + np.log(1-gamma_del+corretion_eps_from_zero)\
                            + del_num*possion_log_probability(0,lamuna_del)+np.log(gamma_mis+corretion_eps_from_zero)
                    else:
                        prb[2] = float('-inf')
                    max_ID = list(prb).index(max(prb))
                    if max_ID==0:
                        max_ID = main_ref
                    elif max_ID==1:
                        max_ID = 4
                    else:
                        max_ID = id_max_bei
                else:
                    prb = np.zeros(2)
                    same_base = base_number[main_ref]
                    not_same_base = sum(base_number[:4]) - same_base
                    
                    prb[0] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero)
                       
                   
                    remain_base = list(base_number[:4])
                    del remain_base[main_ref]
                    base_index = list(range(4))
                    del base_index[main_ref]
                    
                    if max(remain_base)>1:
                        id_max_bei = base_index[int(list(remain_base).index(max(remain_base)))]
                        same_base = base_number[int(id_max_bei)]
                        not_same_base = sum(base_number[:4]) - same_base
                        
                        prb[1] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero) 
                            
                    else:
                        prb[1] = float('-inf')
                    max_ID = list(prb).index(max(prb))
                    if max_ID==0:
                        max_ID = main_ref

                    else:
                        max_ID = id_max_bei
    
            else:
                max_ID = main_ref
        #max_ID = list(potential_probability).index(max(potential_probability))
        return([infer_var_site[site_id], max_ID,0])




def block_sites_infer_ECM2(start_block, end_block, blobk_base_number, standard_ref_code_block, updated_segment, segment_reads,kk):

#for  ll in range(len(list_block)):
#    start_block = list_block[ll][0]
#    end_block = list_block[ll][1]
#    blobk_base_number = present_clu_dis2[kk][start_block:end_block+1,:]
#    standard_ref_code_block = final_ref_code_standard[start_block:end_block+1]
    #gamma_ = [0.99786129,0.00152765,0.00061106,0]
    # theta_ =
    #ll = 9
#    updated_segment = copy.deepcopy(Ref_dic_block[ll])
#    segment_reads = Dep_dic[ll]
    #start_block = 1327
    #end_block = 1328
    #blobk_base_number = copy.deepcopy(present_clu_dis2[kk][1327:1329,:])
    #dep_s = sip.select_ref_block(1327,1328,Best_ref_collection[8][2])
    #updated_segment = copy.deepcopy(dep_s)
    #segment_reads = copy.deepcopy(dep)
    read_number_tem = len(segment_reads)
    #read_id_tem = list(segment_reads.keys())
    #global gamma1, gamma2, gamma3, gamma4,theta_,del_beta,del_p,ins_beta,ins_p
    gamma_c =  gamma1[kk]
    gamma_mis = gamma2[kk]
    gamma_del = gamma3[kk]
    gamma_ins = gamma4[kk]
    gamma_del_p = del_p[kk]
    gamma_del_beta = del_beta[kk]
    gamma_ins_p = ins_p[kk]
    gamma_ins_beta = ins_beta[kk]
    pcorrect = theta_[0]
    p_mis_loc = theta_[1]
    pins1_loc = theta_[2]
    pdel1_loc = theta_[3]
    lamuna_ins = theta_[4]
    lamuna_del = theta_[5]
    corretion_eps_from_zero = 0.0000000001
    return_result = []
    # Indel_bound
    # sreening
    state_mapping = []
    original_true_state_ref = []
    for basss_index in range(start_block, end_block+1):
        if blobk_base_number[int(basss_index-start_block), 5] > 0:
            original_true_state_ref.append(-1)
            state_mapping.append(1)
        else:
            state_mapping.append(2)
            original_true_state_ref.append(0)
    #state_mapping = [1,1,2,2,2,1,1,2,2]
    # linking
    state_dep_dic = {}
    for basss_index in range(start_block, end_block+1):
        # left-ward direction
        state_now = state_mapping[basss_index-start_block]
        if basss_index == start_block:
            left_collection = []
        else:
            left_collection = []
            if state_now == 2:
                shif_l = 1
                while basss_index-shif_l >= start_block:
                    if state_mapping[int(basss_index-shif_l-start_block)] == 2:
                        #left_collection.append(int(basss_index-shif_l))
                        break
                    else:
                        left_collection.append(int(basss_index-shif_l))
                        shif_l = shif_l + 1

            else:
                shif_l = 1
                while basss_index-shif_l >= start_block:
                    if state_mapping[int(basss_index-shif_l-start_block)] == 1:
                        shif_l = shif_l + 1
                    else:
                        left_collection.append(int(basss_index-shif_l))
                        break

        # right-ward direction
        if basss_index == end_block:
            right_collection = []
        else:
            right_collection = []
            if state_now == 2:
                shif_l = 1
                while basss_index+shif_l <= end_block:
                    if state_mapping[int(basss_index+shif_l-start_block)] == 2:
                        break
                    else:
                        right_collection.append(int(basss_index+shif_l))
                        shif_l = shif_l + 1

            else:
                shif_l = 1
                while basss_index+shif_l <= end_block:
                    if state_mapping[int(basss_index+shif_l-start_block)] == 1:
                        shif_l = shif_l + 1
                    else:
                        right_collection.append(int(basss_index+shif_l))
                        break
        final_dep_collec = left_collection+right_collection
        if len(final_dep_collec) > 0:
            state_dep_dic[basss_index] = final_dep_collec
    updated_state_mapping = copy.deepcopy(state_mapping)
    for basss_index in range(start_block, end_block+1):
        # print(updated_segment)
        #basss_index = 899
        potential_probability = np.zeros(5)
        updated_segment_cp = copy.deepcopy(updated_segment)
        #updated_state_mapping_cp = copy.deepcopy()
        # if

        if blobk_base_number[int(basss_index-start_block), 5] > 0:  # ins in the bone
            if basss_index in state_dep_dic:
                tem_vei = state_dep_dic[basss_index]  # neighbour
                ll_xun_ = 0
                case_flag = 1
                while ll_xun_ < len(tem_vei):
                    if updated_segment_cp[tem_vei[ll_xun_]] == 4:
                        case_flag = 0
                        break
                    else:
                        ll_xun_ = ll_xun_ + 1
            else:
                case_flag = 1

            if case_flag == 1:

                for main_base_xun in [0, 1, 2, 3, 5]:
                    if blobk_base_number[int(basss_index-start_block),main_base_xun]>0:
                        updated_segment_cp[basss_index] = main_base_xun
                        full_full_prob = 0
                        for read_xun_tem in range(read_number_tem):
                            tem_dic = segment_reads[int(read_xun_tem)][2]
                            begin_id = max(
                                segment_reads[int(read_xun_tem)][0], start_block)
                            end_id = min(
                                segment_reads[int(read_xun_tem)][1], end_block)
    
                            num_correct = 0
                            num_mis = 0
                            num_del = 0
                            num_ins = 0
                            len_del = {}
                            len_ins = {}
                            start_del_flag = 0
                            end_del_flag = 0
    
                            base_xun = int(begin_id)
                            while base_xun <= end_id:
                                if updated_segment_cp[base_xun] < 4:
                                    #bone_flag = 1
                                    if tem_dic[base_xun] == updated_segment_cp[base_xun]:
    
                                        if num_del > 0:
                                            if num_del not in len_del:
                                                len_del[num_del] = 1
                                            else:
                                                len_del[num_del] = len_del[num_del] + 1
                                            num_del = 0
    
                                        if num_ins > 0:
                                            if num_ins not in len_ins:
                                                len_ins[num_ins] = 1
                                            else:
                                                len_ins[num_ins] = len_ins[num_ins] + 1
                                            num_ins = 0
                                        num_correct = num_correct + 1
    
                                    elif tem_dic[base_xun] < 4 and tem_dic[base_xun] != updated_segment_cp[base_xun]:
                                        if num_del > 0:
                                            if num_del not in len_del:
                                                len_del[num_del] = 1
                                            else:
                                                len_del[num_del] = len_del[num_del] + 1
                                            num_del = 0
                                        if num_ins > 0:
                                            if num_ins not in len_ins:
                                                len_ins[num_ins] = 1
                                            else:
                                                len_ins[num_ins] = len_ins[num_ins] + 1
                                            num_ins = 0
                                        num_mis = num_mis + 1
                                    else:
                                        if num_ins > 0:
                                            if num_ins not in len_ins:
                                                len_ins[num_ins] = 1
                                            else:
                                                len_ins[num_ins] = len_ins[num_ins] + 1
                                            num_ins = 0
    
                                        num_del = num_del + 1
                                        if base_xun == begin_id:
                                            start_del_flag = 1
                                        if base_xun == end_id:
                                            end_del_flag = 1
                                elif updated_segment_cp[base_xun] == 4:
                                    if tem_dic[base_xun] != 4:
                                        if tem_dic[base_xun] < 4:
                                            num_ins = num_ins + 1
                                        else:
                                            print('error1')
                                else:
                                    if tem_dic[base_xun] == 4:
                                        print('error2')
                                    elif tem_dic[base_xun] < 4:
                                        num_ins = num_ins + 1
                                base_xun = base_xun + 1
    
                            
                            if num_ins>0:
                                if num_ins in len_ins:
                                    len_ins[num_ins] = len_ins[num_ins] + 1
                                else:
                                    len_ins[num_ins] = 1
                            if num_del>0:
                                if num_del in len_del:
                                    len_del[num_del] = len_del[num_del] + 1
                                else:
                                    len_del[num_del] = 1
                            tem_probabity_correct = num_correct*np.log(pcorrect+corretion_eps_from_zero)
                            tem_probabity_mis = num_mis*np.log(p_mis_loc+corretion_eps_from_zero)
                            if len(len_ins) > 0:
                                values_ = len_ins.values()
                                num_ins_loc = sum(np.array(list(values_)))
                                tem_probabity_ins = num_ins_loc*np.log(pins1_loc+corretion_eps_from_zero)
                                for ii in len_ins:
                                    # np.log(stats.expon.pdf(ii,lamuna_ins))
                                    tem_probabity_ins = tem_probabity_ins + \
                                        len_ins[ii] * \
                                        possion_log_probability(ii-1,lamuna_ins)
                                    # print(tem_probabity_ins)
                                    # tem_probabity_ins = tem_probabity_ins + len_ins[ii]*np.log(stats.nbinom.pmf(ii,n_ins,p_ins))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                    # tem_probabity_ins = tem_probabity_ins + len_ins[ii]*np.log(1/(ii*sig_ins*np.sqrt(2*math.pi))*np.exp(-(np.log(ii)-mu_ins)**2/(2*sig_ins**2)))#np.log(stats.expon.pdf(ii,lamuna_ins))
                            else:
                                tem_probabity_ins = 0
                                num_ins_loc = 0
    
                            del_whole = 0
                            if len(len_del) > 0:
                                values_ = len_del.values()
                                num_del_loc = sum(np.array(list(values_)))
                                tem_probabity_del = num_del_loc*np.log(pdel1_loc+corretion_eps_from_zero)
                                for dd in len_del:
                                    # np.log(stats.expon.pdf(ii,lamuna_ins))
                                    tem_probabity_del = tem_probabity_del + \
                                        len_del[dd] * \
                                         possion_log_probability(dd-1,lamuna_del)
                                    # tem_probabity_del = tem_probabity_del + len_del[dd]*np.log(stats.nbinom.pmf(dd,n_del,p_del))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                    # tem_probabity_del = tem_probabity_del + len_del[dd]*np.log(1/(dd*sig_del*np.sqrt(2*math.pi))*np.exp(-(np.log(dd)-mu_del)**2/(2*sig_del**2)))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                    del_whole = del_whole + dd*len_del[dd]
    
                            else:
                                tem_probabity_del = 0
                                num_del_loc = 0
    
                            other_gaps = num_correct + num_mis - 1-num_del_loc+start_del_flag + end_del_flag
                            gap_not_ins_pro = other_gaps*np.log(1-pins1_loc+corretion_eps_from_zero)
                            not_del_num = num_correct + num_mis
                            gap_not_del_pro = not_del_num*np.log(1-pdel1_loc+corretion_eps_from_zero)
                            full_prob = gap_not_del_pro+tem_probabity_correct + tem_probabity_mis +\
                                tem_probabity_ins+tem_probabity_del+gap_not_ins_pro
    
                            full_full_prob = full_full_prob + full_prob
                        
                        ### evo comparision
                        

                        ID_tem = [0, 1, 2, 3, 5].index(main_base_xun)
                        ### here we should judge the influence of the variations
                        
                        
                        
                        
                        Tem_ins_evo = {}
                        Tem_del_evo = {}
                        ins_len_evo = 0
                        del_len_evo = 0
                        num_mis_evo = 0
                        num_cor_evo = 0 
                        for ll_check in range(start_block, end_block+1):
                            relative_id =  ll_check - start_block
                            if original_true_state_ref[relative_id]<0:## meaning this should ins evo fromthe refence
                                
                                if del_len_evo>0:
                                    if del_len_evo in Tem_del_evo:
                                        Tem_del_evo[del_len_evo] = Tem_del_evo[del_len_evo] + 1
                                    else:
                                        Tem_del_evo[del_len_evo] = 1
                                    del_len_evo = 0
                                if updated_state_mapping[relative_id]>0:
                                    ins_len_evo = ins_len_evo + 1
                                
                            else:
                                if ins_len_evo>0:
                                    if ins_len_evo in Tem_ins_evo:
                                        Tem_ins_evo[ins_len_evo] = Tem_ins_evo[ins_len_evo] + 1
                                    else:
                                        Tem_ins_evo[ins_len_evo] =  1
                                    ins_len_evo = 0
                                if updated_state_mapping[relative_id]>0:
                                    del_len_evo = del_len_evo + 1
                                else:
                                    if ll_check == basss_index:
                                        if main_base_xun == standard_ref_code_block[relative_id]:
                                            num_cor_evo = num_cor_evo + 1
                                        else:
                                            num_mis_evo = num_mis_evo + 1
                                    else:
                                        if updated_segment_cp[ll_check] ==  standard_ref_code_block[relative_id]:
                                            num_cor_evo = num_cor_evo + 1
                                        else:
                                            num_mis_evo = num_mis_evo + 1
                        if del_len_evo>0:
                            if del_len_evo in Tem_del_evo:
                                Tem_del_evo[del_len_evo] = Tem_del_evo[del_len_evo] + 1
                            else:
                                Tem_del_evo[del_len_evo] = 1
                        if ins_len_evo>0:
                            if ins_len_evo in Tem_ins_evo:
                                Tem_ins_evo[ins_len_evo] = Tem_ins_evo[ins_len_evo] + 1
                            else:
                                Tem_ins_evo[ins_len_evo] =  1
                        ### infer the prob
                        tem_evo = 0
                        
                        for del_per in Tem_del_evo:
                            tem_evo = tem_evo + Tem_del_evo[del_per]*(expotential_log_log_probability(del_per-1,del_p[kk],del_beta[kk]) + np.log(gamma_del+corretion_eps_from_zero) )
                        for ins_per in Tem_ins_evo:
                            tem_evo = tem_evo + Tem_ins_evo[ins_per]*(expotential_log_log_probability(ins_per-1,ins_p[kk],ins_beta[kk]) + np.log(gamma_ins+corretion_eps_from_zero) )
                        tem_evo = tem_evo + (num_cor_evo+ num_mis_evo)*np.log((1-gamma_del)) + (num_cor_evo)*np.log(gamma_c) + (num_mis_evo)*np.log(gamma_mis+corretion_eps_from_zero)              
                       
                        
                        potential_probability[int(ID_tem)] = full_full_prob + tem_evo
                    else:
                        ID_tem = [0, 1, 2, 3, 5].index(main_base_xun)
                        potential_probability[int(ID_tem)] = float('-inf')
                #tem_probability = copy.deepcopy(potential_probability/tem_pri)
                #tem_probability = tem_probability - max(tem_probability)
                #tem_probability_main = np.exp(tem_probability)
                #tem_probability_main = tem_probability_main/sum(tem_probability_main)
                #post_sampling = np.random.multinomial(n= 1,pvals = tem_probability_main)
                #max_ID = list(post_sampling).index(1)
                #    n=1, pvals=[6/8, 1/8, 1/16, 1/32, 1/32])
                #l_i = list(l_i_vec).index(1)+1
                max_ID = list(potential_probability).index(max(potential_probability))
                
                if max_ID >3:
                    updated_state_mapping[basss_index- start_block] = -1
                else:
                    updated_state_mapping[basss_index- start_block] = 1
                    
                if max_ID == 4:
                    updated_segment[basss_index] = 5
                    
            else:
                updated_segment[basss_index] = 5
                updated_state_mapping[basss_index- start_block] = -1
        else:  # del condition
            if basss_index in state_dep_dic:
                tem_vei = state_dep_dic[basss_index]
                ll_xun_ = 0
                case_flag = 1
                while ll_xun_ < len(tem_vei):
                    if updated_segment_cp[tem_vei[ll_xun_]] < 4:
                        case_flag = 0
                        break
                    else:
                        ll_xun_ = ll_xun_ + 1
            else:
                case_flag = 1

            if case_flag:
                for main_base_xun in [0, 1, 2, 3, 4]:
                    if blobk_base_number[int(basss_index-start_block),main_base_xun]>0:
                        updated_segment_cp[basss_index] = main_base_xun
                        full_full_prob = 0
                        for read_xun_tem in range(read_number_tem):
                            tem_dic = segment_reads[int(read_xun_tem)][2]
                            begin_id = max(
                                segment_reads[int(read_xun_tem)][0], start_block)
                            end_id = min(
                                segment_reads[int(read_xun_tem)][1], end_block)
    
                            num_correct = 0
                            num_mis = 0
                            num_del = 0
                            num_ins = 0
                            len_del = {}
                            len_ins = {}
                            start_del_flag = 0
                            end_del_flag = 0
    
                            base_xun = int(begin_id)
                            while base_xun <= end_id:
                                if updated_segment_cp[base_xun] < 4:
                                    #bone_flag = 1
                                    if tem_dic[base_xun] == updated_segment_cp[base_xun]:
    
                                        if num_del > 0:
                                            if num_del not in len_del:
                                                len_del[num_del] = 1
                                            else:
                                                len_del[num_del] = len_del[num_del] + 1
                                            num_del = 0
    
                                        if num_ins > 0:
                                            if num_ins not in len_ins:
                                                len_ins[num_ins] = 1
                                            else:
                                                len_ins[num_ins] = len_ins[num_ins] + 1
                                            num_ins = 0
                                        num_correct = num_correct + 1
    
                                    elif tem_dic[base_xun] < 4 and tem_dic[base_xun] != updated_segment_cp[base_xun]:
                                        if num_del > 0:
                                            if num_del not in len_del:
                                                len_del[num_del] = 1
                                            else:
                                                len_del[num_del] = len_del[num_del] + 1
                                            num_del = 0
                                        if num_ins > 0:
                                            if num_ins not in len_ins:
                                                len_ins[num_ins] = 1
                                            else:
                                                len_ins[num_ins] = len_ins[num_ins] + 1
                                            num_ins = 0
                                        num_mis = num_mis + 1
                                    else:
                                        if num_ins > 0:
                                            if num_ins not in len_ins:
                                                len_ins[num_ins] = 1
                                            else:
                                                len_ins[num_ins] = len_ins[num_ins] + 1
                                            num_ins = 0
    
                                        num_del = num_del + 1
                                        if base_xun == begin_id:
                                            start_del_flag = 1
                                        if base_xun == end_id:
                                            end_del_flag = 1
                                elif updated_segment_cp[base_xun] == 4:
                                    if tem_dic[base_xun] != 4:
                                        if tem_dic[base_xun] < 4:
                                            num_ins = num_ins + 1
                                        else:
                                            print('error1')
                                            print(basss_index)
                                            #break
                                else:
                                    if tem_dic[base_xun] == 4:
                                        print('error2')
                                    elif tem_dic[base_xun] < 4:
                                        num_ins = num_ins + 1
                                base_xun = base_xun + 1
                            
                            if num_ins>0:
                                if num_ins in len_ins:
                                    len_ins[num_ins] = len_ins[num_ins] + 1
                                else:
                                    len_ins[num_ins] = 1
                            if num_del>0:
                                if num_del in len_del:
                                    len_del[num_del] = len_del[num_del] + 1
                                else:
                                    len_del[num_del] = 1
                            tem_probabity_correct = num_correct*np.log(pcorrect+corretion_eps_from_zero)
                            tem_probabity_mis = num_mis*np.log(p_mis_loc+corretion_eps_from_zero)
                            if len(len_ins) > 0:
                                values_ = len_ins.values()
                                num_ins_loc = sum(np.array(list(values_)))
                                tem_probabity_ins = num_ins_loc*np.log(pins1_loc+corretion_eps_from_zero)
                                for ii in len_ins:
                                    # np.log(stats.expon.pdf(ii,lamuna_ins))
                                    tem_probabity_ins = tem_probabity_ins + \
                                        len_ins[ii] * \
                                         possion_log_probability(ii-1,lamuna_ins)
                                    # print(tem_probabity_ins)
                                    # tem_probabity_ins = tem_probabity_ins + len_ins[ii]*np.log(stats.nbinom.pmf(ii,n_ins,p_ins))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                    # tem_probabity_ins = tem_probabity_ins + len_ins[ii]*np.log(1/(ii*sig_ins*np.sqrt(2*math.pi))*np.exp(-(np.log(ii)-mu_ins)**2/(2*sig_ins**2)))#np.log(stats.expon.pdf(ii,lamuna_ins))
                            else:
                                tem_probabity_ins = 0
                                num_ins_loc = 0
    
                            del_whole = 0
                            if len(len_del) > 0:
                                values_ = len_del.values()
                                num_del_loc = sum(np.array(list(values_)))
                                tem_probabity_del = num_del_loc*np.log(pdel1_loc+corretion_eps_from_zero)
                                for dd in len_del:
                                    # np.log(stats.expon.pdf(ii,lamuna_ins))
                                    tem_probabity_del = tem_probabity_del + \
                                        len_del[dd] * \
                                         possion_log_probability(dd-1,lamuna_del)                                    # tem_probabity_del = tem_probabity_del + len_del[dd]*np.log(stats.nbinom.pmf(dd,n_del,p_del))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                    # tem_probabity_del = tem_probabity_del + len_del[dd]*np.log(1/(dd*sig_del*np.sqrt(2*math.pi))*np.exp(-(np.log(dd)-mu_del)**2/(2*sig_del**2)))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                    del_whole = del_whole + dd*len_del[dd]
    
                            else:
                                tem_probabity_del = 0
                                num_del_loc = 0
    
                            other_gaps = num_correct + num_mis - 1-num_del_loc+start_del_flag + end_del_flag
                            gap_not_ins_pro = other_gaps*np.log(1-pins1_loc+corretion_eps_from_zero)
                            not_del_num = num_correct + num_mis
                            gap_not_del_pro = not_del_num*np.log(1-pdel1_loc+corretion_eps_from_zero)
                            full_prob = gap_not_del_pro+tem_probabity_correct + tem_probabity_mis +\
                                tem_probabity_ins+tem_probabity_del+gap_not_ins_pro
    
                            full_full_prob = full_full_prob + full_prob
    
                        ID_tem = [0, 1, 2, 3, 4].index(main_base_xun)
                        
                        
                        Tem_ins_evo = {}
                        Tem_del_evo = {}
                        ins_len_evo = 0
                        del_len_evo = 0
                        num_mis_evo = 0
                        num_cor_evo = 0 
                        for ll_check in range(start_block, end_block+1):
                            relative_id =  ll_check - start_block
                            if original_true_state_ref[relative_id]<0:## meaning this should ins evo fromthe refence
                                
                                if del_len_evo>0:
                                    if del_len_evo in Tem_del_evo:
                                        Tem_del_evo[del_len_evo] = Tem_del_evo[del_len_evo] + 1
                                    else:
                                        Tem_del_evo[del_len_evo] = 1
                                    del_len_evo = 0
                                if updated_state_mapping[relative_id]>0:
                                    ins_len_evo = ins_len_evo + 1
                                
                            else:
                                if ins_len_evo>0:
                                    if ins_len_evo in Tem_ins_evo:
                                        Tem_ins_evo[ins_len_evo] = Tem_ins_evo[ins_len_evo] + 1
                                    else:
                                        Tem_ins_evo[ins_len_evo] =  1
                                    ins_len_evo = 0
                                if updated_state_mapping[relative_id]>0:
                                    del_len_evo = del_len_evo + 1
                                else:
                                    if ll_check == basss_index:
                                        if main_base_xun == standard_ref_code_block[relative_id]:
                                            num_cor_evo = num_cor_evo + 1
                                        else:
                                            num_mis_evo = num_mis_evo + 1
                                    else:
                                        if  updated_segment_cp[ll_check] ==  standard_ref_code_block[relative_id]:
                                            num_cor_evo = num_cor_evo + 1
                                        else:
                                            num_mis_evo = num_mis_evo + 1
                        
                        
                        if del_len_evo>0:
                            if del_len_evo in Tem_del_evo:
                                Tem_del_evo[del_len_evo] = Tem_del_evo[del_len_evo] + 1
                            else:
                                Tem_del_evo[del_len_evo] = 1
                        
                        if ins_len_evo>0:
                            if ins_len_evo in Tem_ins_evo:
                                Tem_ins_evo[ins_len_evo] = Tem_ins_evo[ins_len_evo] + 1
                            else:
                                Tem_ins_evo[ins_len_evo] =  1
                        ### infer the prob
                        tem_evo = 0
                        
                        for del_per in Tem_del_evo:
                            tem_evo = tem_evo + Tem_del_evo[del_per]*(expotential_log_log_probability(del_per-1,del_p[kk],del_beta[kk]) + np.log(gamma_del+corretion_eps_from_zero) )
                        
                        for ins_per in Tem_ins_evo:
                            tem_evo = tem_evo + Tem_ins_evo[ins_per]*(expotential_log_log_probability(ins_per-1,ins_p[kk],ins_beta[kk]) + np.log(gamma_ins+corretion_eps_from_zero) )
                        tem_evo = tem_evo + (num_cor_evo+ num_mis_evo)*np.log((1-gamma_del)) + (num_cor_evo)*np.log(gamma_c) + (num_mis_evo)*np.log(gamma_mis+corretion_eps_from_zero)     
                        
                        potential_probability[int(ID_tem)] = full_full_prob + tem_evo
                    else:
                        ID_tem = [0, 1, 2, 3, 4].index(main_base_xun)
                        potential_probability[int(ID_tem)] = float('-inf')
                #tem_probability = copy.deepcopy(potential_probability/tem_pri)
                #tem_probability = tem_probability - max(tem_probability)
                #tem_probability_main = np.exp(tem_probability)
                #tem_probability_main = tem_probability_main/sum(tem_probability_main)
                #post_sampling = np.random.multinomial(n= 1,pvals = tem_probability_main)
                #max_ID = list(post_sampling).index(1)
                #    n=1, pvals=[6/8, 1/8, 1/16, 1/32, 1/32])
                #l_i = list(l_i_vec).index(1)+1
                max_ID = list(potential_probability).index(max(potential_probability))
                
                if max_ID >3:
                    updated_state_mapping[basss_index- start_block] = 2
                else:
                    updated_state_mapping[basss_index- start_block] = 0
                

            else:
                for main_base_xun in [0, 1, 2, 3]:
                    if blobk_base_number[int(basss_index-start_block),main_base_xun]>0:
                        updated_segment_cp[basss_index] = main_base_xun
                        full_full_prob = 0
                        for read_xun_tem in range(read_number_tem):
                            tem_dic = segment_reads[int(read_xun_tem)][2]
                            begin_id = max(
                                segment_reads[int(read_xun_tem)][0], start_block)
                            end_id = min(
                                segment_reads[int(read_xun_tem)][1], end_block)
    
                            num_correct = 0
                            num_mis = 0
                            num_del = 0
                            num_ins = 0
                            len_del = {}
                            len_ins = {}
                            start_del_flag = 0
                            end_del_flag = 0
    
                            base_xun = int(begin_id)
                            while base_xun <= end_id:
                                if updated_segment_cp[base_xun] < 4:
                                    #bone_flag = 1
                                    if tem_dic[base_xun] == updated_segment_cp[base_xun]:
    
                                        if num_del > 0:
                                            if num_del not in len_del:
                                                len_del[num_del] = 1
                                            else:
                                                len_del[num_del] = len_del[num_del] + 1
                                            num_del = 0
    
                                        if num_ins > 0:
                                            if num_ins not in len_ins:
                                                len_ins[num_ins] = 1
                                            else:
                                                len_ins[num_ins] = len_ins[num_ins] + 1
                                            num_ins = 0
                                        num_correct = num_correct + 1
    
                                    elif tem_dic[base_xun] < 4 and tem_dic[base_xun] != updated_segment_cp[base_xun]:
                                        if num_del > 0:
                                            if num_del not in len_del:
                                                len_del[num_del] = 1
                                            else:
                                                len_del[num_del] = len_del[num_del] + 1
                                            num_del = 0
                                        if num_ins > 0:
                                            if num_ins not in len_ins:
                                                len_ins[num_ins] = 1
                                            else:
                                                len_ins[num_ins] = len_ins[num_ins] + 1
                                            num_ins = 0
                                        num_mis = num_mis + 1
                                    else:
                                        if num_ins > 0:
                                            if num_ins not in len_ins:
                                                len_ins[num_ins] = 1
                                            else:
                                                len_ins[num_ins] = len_ins[num_ins] + 1
                                            num_ins = 0
    
                                        num_del = num_del + 1
                                        if base_xun == begin_id:
                                            start_del_flag = 1
                                        if base_xun == end_id:
                                            end_del_flag = 1
                                elif updated_segment_cp[base_xun] == 4:
                                    if tem_dic[base_xun] != 4:
                                        if tem_dic[base_xun] < 4:
                                            num_ins = num_ins + 1
                                        else:
                                            print('error1')
                                else:
                                    if tem_dic[base_xun] == 4:
                                        print('error2')
                                    elif tem_dic[base_xun] < 4:
                                        num_ins = num_ins + 1
                                base_xun = base_xun + 1
    
                            
                            if num_ins>0:
                                if num_ins in len_ins:
                                    len_ins[num_ins] = len_ins[num_ins] + 1
                                else:
                                    len_ins[num_ins] = 1
                            if num_del>0:
                                if num_del in len_del:
                                    len_del[num_del] = len_del[num_del] + 1
                                else:
                                    len_del[num_del] = 1
                            tem_probabity_correct = num_correct*np.log(pcorrect+corretion_eps_from_zero)
                            tem_probabity_mis = num_mis*np.log(p_mis_loc+corretion_eps_from_zero)
                            if len(len_ins) > 0:
                                values_ = len_ins.values()
                                num_ins_loc = sum(np.array(list(values_)))
                                tem_probabity_ins = num_ins_loc*np.log(pins1_loc+corretion_eps_from_zero)
                                for ii in len_ins:
                                    # np.log(stats.expon.pdf(ii,lamuna_ins))
                                    tem_probabity_ins = tem_probabity_ins + \
                                        len_ins[ii] * \
                                         possion_log_probability(ii-1,lamuna_ins)
                                    # print(tem_probabity_ins)
                                    # tem_probabity_ins = tem_probabity_ins + len_ins[ii]*np.log(stats.nbinom.pmf(ii,n_ins,p_ins))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                    # tem_probabity_ins = tem_probabity_ins + len_ins[ii]*np.log(1/(ii*sig_ins*np.sqrt(2*math.pi))*np.exp(-(np.log(ii)-mu_ins)**2/(2*sig_ins**2)))#np.log(stats.expon.pdf(ii,lamuna_ins))
                            else:
                                tem_probabity_ins = 0
                                num_ins_loc = 0
    
                            del_whole = 0
                            if len(len_del) > 0:
                                values_ = len_del.values()
                                num_del_loc = sum(np.array(list(values_)))
                                tem_probabity_del = num_del_loc*np.log(pdel1_loc+corretion_eps_from_zero)
                                for dd in len_del:
                                    # np.log(stats.expon.pdf(ii,lamuna_ins))
                                    tem_probabity_del = tem_probabity_del + \
                                        len_del[dd] * \
                                         possion_log_probability(dd-1,lamuna_del)
                                    # tem_probabity_del = tem_probabity_del + len_del[dd]*np.log(stats.nbinom.pmf(dd,n_del,p_del))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                    # tem_probabity_del = tem_probabity_del + len_del[dd]*np.log(1/(dd*sig_del*np.sqrt(2*math.pi))*np.exp(-(np.log(dd)-mu_del)**2/(2*sig_del**2)))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                    del_whole = del_whole + dd*len_del[dd]
    
                            else:
                                tem_probabity_del = 0
                                num_del_loc = 0
    
                            other_gaps = num_correct + num_mis - 1-num_del_loc+start_del_flag + end_del_flag
                            gap_not_ins_pro = other_gaps*np.log(1-pins1_loc+corretion_eps_from_zero)
                            not_del_num = num_correct + num_mis
                            gap_not_del_pro = not_del_num*np.log(1-pdel1_loc+corretion_eps_from_zero)
                            full_prob = gap_not_del_pro+tem_probabity_correct + tem_probabity_mis +\
                                tem_probabity_ins+tem_probabity_del+gap_not_ins_pro
    
                            full_full_prob = full_full_prob + full_prob
    
                        ID_tem = [0, 1, 2, 3].index(main_base_xun)
                 
                        updated_state_mapping[basss_index-start_block] = 0
                        Tem_ins_evo = {}
                        
                        Tem_del_evo = {}
                        ins_len_evo = 0
                        del_len_evo = 0
                        num_mis_evo = 0
                        num_cor_evo = 0 
                        for ll_check in range(start_block, end_block+1):
                            relative_id =  ll_check - start_block
                            if original_true_state_ref[relative_id]<0:## meaning this should ins evo fromthe refence
                                
                                if del_len_evo>0:
                                    if del_len_evo in Tem_del_evo:
                                        Tem_del_evo[del_len_evo] = Tem_del_evo[del_len_evo] + 1
                                    else:
                                        Tem_del_evo[del_len_evo] = 1
                                    del_len_evo = 0
                                if updated_state_mapping[relative_id]>0:
                                    ins_len_evo = ins_len_evo + 1
                                
                            else:
                                if ins_len_evo>0:
                                    if ins_len_evo in Tem_ins_evo:
                                        Tem_ins_evo[ins_len_evo] = Tem_ins_evo[ins_len_evo] + 1
                                    else:
                                        Tem_ins_evo[ins_len_evo] =  1
                                    ins_len_evo = 0
                                if updated_state_mapping[relative_id]>0:
                                    del_len_evo = del_len_evo + 1
                                else:
                                    if ll_check == basss_index:
                                        if main_base_xun == standard_ref_code_block[relative_id]:
                                            num_cor_evo = num_cor_evo + 1
                                        else:
                                            num_mis_evo = num_mis_evo + 1
                                    else:
                                        if  updated_segment_cp[ll_check] ==  standard_ref_code_block[relative_id]:
                                            num_cor_evo = num_cor_evo + 1
                                        else:
                                            num_mis_evo = num_mis_evo + 1
                        
                        
                        if del_len_evo>0:
                            if del_len_evo in Tem_del_evo:
                                Tem_del_evo[del_len_evo] = Tem_del_evo[del_len_evo] + 1
                            else:
                                Tem_del_evo[del_len_evo] = 1
                        if ins_len_evo>0:
                            if ins_len_evo in Tem_ins_evo:
                                Tem_ins_evo[ins_len_evo] = Tem_ins_evo[ins_len_evo] + 1
                            else:
                                Tem_ins_evo[ins_len_evo] =  1
                        ### infer the prob
                        tem_evo = 0
                        
                        for del_per in Tem_del_evo:
                            tem_evo = tem_evo + Tem_del_evo[del_per]*(expotential_log_log_probability(del_per-1,del_p[kk],del_beta[kk]) + np.log(gamma_del+corretion_eps_from_zero) )
                        for ins_per in Tem_ins_evo:
                            tem_evo = tem_evo + Tem_ins_evo[ins_per]*(expotential_log_log_probability(ins_per-1,ins_p[kk],ins_beta[kk]) + np.log(gamma_ins+corretion_eps_from_zero) )
                        tem_evo = tem_evo + (num_cor_evo+ num_mis_evo)*np.log((1-gamma_del)) + (num_cor_evo)*np.log(gamma_c) + (num_mis_evo)*np.log(gamma_mis+corretion_eps_from_zero)     
                        
                        potential_probability[int(ID_tem)] = full_full_prob + tem_evo
                    
                    else:
                        ID_tem = [0, 1, 2, 3].index(main_base_xun)
                        potential_probability[int(ID_tem)] = float('-inf')
                #tem_probability = copy.deepcopy(potential_probability/tem_pri)
                #tem_probability = tem_probability - max(tem_probability)
                #tem_probability_main = np.exp(tem_probability)
                #tem_probability_main = tem_probability_main/sum(tem_probability_main)
                #post_sampling = np.random.multinomial(n= 1,pvals = tem_probability_main)
                #max_ID = list(post_sampling).index(1)
                #    n=1, pvals=[6/8, 1/8, 1/16, 1/32, 1/32])
                #l_i = list(l_i_vec).index(1)+1
                max_ID = list(potential_probability).index(max(potential_probability))
                updated_state_mapping[basss_index - start_block] = 0

            updated_segment[basss_index] = max_ID

        return_result.append([basss_index, updated_segment[basss_index]])

    
    return(return_result)



def block_sites_infer_ECM3(start_block_id, end_block_id, blobk_base_number, standard_ref_code_block, updated_segment,\
                           segment_reads,common_del_per,common_ins_per,kk):
    #for ll_test in range(len(remaining_dependent)):
    #    print(ll_test)
        #ll_test = 26
        global infer_var_site,common_del_block_df,common_ins_block_dic,infer_var_ori_bone,infer_var_ori_state
    #    start_block_id = int(remaining_dependent.iloc[ll_test,0])
    #    end_block_id = int(remaining_dependent.iloc[ll_test,1])
        start_block = infer_var_site[start_block_id]
        end_block = infer_var_site[end_block_id]
    #    blobk_base_number = present_clu_dis2[kk][int(start_block_id):int(end_block_id)+1, :]
    #    standard_ref_code_block = []
    #    for ll_b_id in range(start_block_id,end_block_id+1):
    #        standard_ref_code_block.append(final_ref_code_standard[infer_var_site[ll_b_id]])
    #    updated_segment = copy.deepcopy(Ref_dic_block[ll_test])
    #    segment_reads = Dep_dic[ll_test]
    #    common_del_per = dep_spe_infor_del_col_full[ll_test]
    #    common_ins_per = dep_spe_infor_ins_col_full[ll_test]
       
        #common_ins_per_df = pd.DataFrame(common_ins_per)
        common_del_per_df = pd.DataFrame(common_del_per)
        
        ### common del part we present the neighbour condition
        common_del_per_df_list = []
        if len(common_del_per_df)>0:
            for ll_del_id in range(len(common_del_per_df)):
                #print(common_del_per_df)
                loc_common_del_df = common_del_per_df.iloc[ll_del_id,2]## common del df  row
                tem_start_del = common_del_block_df.iloc[loc_common_del_df,0]
                tem_end_del = common_del_block_df.iloc[loc_common_del_df,1]
                if [tem_start_del,tem_end_del] not in common_del_per_df_list:
                    common_del_per_df_list.append([tem_start_del,tem_end_del])
        #common_del_per_df_list = set(common_del_per_df_list)
        common_del_per_df_list_df = pd.DataFrame(common_del_per_df_list)
                #if common_del_per_df.iloc[ll_del_id,1] not in common_del_per_df_dic:
                #    common_del_per_df_dic[common_del_per_df.iloc[ll_del_id,1]] = [common_del_per_df.iloc[ll_del_id,0],common_del_per_df.iloc[ll_del_id,2]]
                #else:
                #    common_del_per_df_dic[common_del_per_df.iloc[ll_del_id,1]].append([common_del_per_df.iloc[ll_del_id,0],common_del_per_df.iloc[ll_del_id,2]])
       
        ### ins part we know the common ins position
        common_ins_per_df_list = []
        if len(common_ins_per)>0:
            for ll_ins_id in common_ins_per:
                ins_tem_content = common_ins_block_dic[ll_ins_id]
                for ele_ins_tem_content in ins_tem_content:
                    ins_tem_content_start = ele_ins_tem_content[0]
                    ins_tem_content_end = ele_ins_tem_content[1]
                    ins_tem_content_bone_trans = trans_(ins_inserted_dic, [ll_ins_id])[ll_ins_id]
                    ins_tem_content_start_trans = ins_tem_content_bone_trans + ins_tem_content_start
                    ins_tem_content_end_trans = ins_tem_content_bone_trans + ins_tem_content_end
                    common_ins_per_df_list.append([ins_tem_content_start_trans,ins_tem_content_end_trans])
        common_ins_per_df_list_df = pd.DataFrame(common_ins_per_df_list)
    #for  ll in range(len(list_block)):
    #    start_block = list_block[ll][0]
    #    end_block = list_block[ll][1]
    #    blobk_base_number = present_clu_dis2[kk][start_block:end_block+1,:]
    #    standard_ref_code_block = final_ref_code_standard[start_block:end_block+1]
        #gamma_ = [0.99786129,0.00152765,0.00061106,0]
        # theta_ =
        #ll = 9
    #    updated_segment = copy.deepcopy(Ref_dic_block[ll])
    #    segment_reads = Dep_dic[ll]
        #start_block = 1327
        #end_block = 1328
        #blobk_base_number = copy.deepcopy(present_clu_dis2[kk][1327:1329,:])
        #dep_s = sip.select_ref_block(1327,1328,Best_ref_collection[8][2])
        #updated_segment = copy.deepcopy(dep_s)
        #segment_reads = copy.deepcopy(dep)
        read_number_tem = len(segment_reads)
        #read_id_tem = list(segment_reads.keys())
        #global gamma1, gamma2, gamma3, gamma4,theta_,del_beta,del_p,ins_beta,ins_p
        gamma_c =  gamma1[kk]
        gamma_mis = gamma2[kk]
        gamma_del = gamma3[kk]
        gamma_ins = gamma4[kk]
        gamma_del_p = del_p[kk]
        gamma_del_beta = del_beta[kk]
        gamma_ins_p = ins_p[kk]
        gamma_ins_beta = ins_beta[kk]
        pcorrect = theta_[0]
        p_mis_loc = theta_[1]
        pins1_loc = theta_[2]
        pdel1_loc = theta_[3]
        lamuna_ins = theta_[4]
        lamuna_del = theta_[5]
        corretion_eps_from_zero = 0.0000000001
        return_result = []
        # Indel_bound
        # sreening
        state_mapping = []
        original_true_state_ref = []
        start_block_id = int(start_block_id)
        end_block_id = int(end_block_id)
        var_site_collection_abs = []
        var_site_collection_ori = []
        for basss_index in range(start_block_id, end_block_id+1):
            var_site_collection_abs.append(infer_var_site[basss_index])
            var_site_collection_ori.append(infer_var_ori_bone[basss_index])
            if infer_var_ori_state[basss_index]<0:
                state_mapping.append(2)
            else:
                state_mapping.append(1)
            #state_mapping.append(infer_var_ori_state[basss_index])
            if infer_var_ori_state[basss_index] > -1:
                original_true_state_ref.append(-1)
                #state_mapping.append(1)
                
            else:
                #if blobk_base_number[int(basss_index-start_block_id), 4] > 0:
                #    state_mapping.append(2)
                #else:
                #    state_mapping.append(0)
                original_true_state_ref.append(0)
        #state_mapping = [1,1,2,2,2,1,1,2,2]
        # linking
        state_dep_dic = {}
        for basss_index in range(start_block_id, end_block_id+1):
             
            # left-ward direction
            state_now = state_mapping[basss_index-start_block_id]
            present_var_msa = infer_var_site[basss_index]
            if basss_index == start_block_id:
                #left_collection_ins = []
                #left_collection_del = []
                left_collection = []
                #if state_now ==2:
                #    if len(common_del_per)>0:
                #        for dep_del_id_common in range(len(common_del_per)):
                #            if common_del_per[dep_del_id_common][0]<0:
                #                if common_del_per[dep_del_id_common][1] == basss_index:
                #                    id_hang_common = common_del_per[dep_del_id_common][2] 
                #                    common_del_length = common_del_block_df.iloc[id_hang_common,1] - common_del_block_df.iloc[id_hang_common,0] + 1
                #                    left_collection_del = [[1,common_del_length]] # common is 1 
                #else:### consider the common ins
                #    common_ins_per_df_sub = common_ins_per_df[common_ins_per_df.iloc[:,0]==present_var_msa]
                #    ins_bone_original = common_ins_per_df_sub.iloc[0,1]
                #    ins_sub_pos_original = common_ins_per_df_sub.iloc[0,2]
                #    if ins_sub_pos_original>0:
                #        common_ins_length = ins_sub_pos_original
                #        left_collection_ins = [[1,common_ins_length]]
               
            else:
                left_collection = []
                if state_now == 2:
                    shif_l = 1
                    while basss_index-shif_l >= start_block_id:
                        if state_mapping[int(basss_index-shif_l-start_block_id)] == 2:
                            #left_collection.append(int(basss_index-shif_l))
                            break
                        else:
                            if var_site_collection_ori[basss_index-start_block_id] - var_site_collection_ori[int(basss_index-shif_l-start_block_id)]<2:
                                left_collection.append(int(basss_index-shif_l))
                                shif_l = shif_l + 1
                            else:
                                break
    
                else:
                    shif_l = 1
                    while basss_index-shif_l >= start_block_id:
                        if state_mapping[int(basss_index-shif_l-start_block_id)] == 1:
                            if var_site_collection_ori[basss_index-start_block_id] - var_site_collection_ori[int(basss_index-shif_l-start_block_id)]<2:
                                
                                    shif_l = shif_l + 1
                            else:
                                break
                        
                        else:
                            left_collection.append(int(basss_index-shif_l))
                            break
                #left_collection_del = []
                #left_collection_ins = []
                #if state_now == 2:
                #    shif_l = 1
                #    while basss_index-shif_l >= start_block:
                        
                #        if state_mapping[int(basss_index-shif_l-start_block)] == 2:
                            #left_collection.append(int(basss_index-shif_l))
                #            if len(common_del_per)>0:
                #                common_del_per_df_l = common_del_per_df[common_del_per_df.iloc[:,0]<0]
                #                if common_del_per_df_l.iloc[0,1] == basss_index:
                #                    id_hang_common =common_del_per_df_l.iloc[0,2]
                #                    common_del_length = common_del_block_df.iloc[id_hang_common,1] - common_del_block_df.iloc[id_hang_common,0] + 1
                #                    left_collection_del = [[1,common_del_length]] # common is 1
                #            break
                #        else:
                #            common_ins_per_df_sub = common_ins_per_df[common_ins_per_df.iloc[:,0]==present_var_msa]
                #            ins_bone_original = common_ins_per_df_sub.iloc[0,1]
                #            ins_sub_pos_original = common_ins_per_df_sub.iloc[0,2]
                #            if ins_sub_pos_original>0:
                #                common_ins_length = ins_sub_pos_original
                #                left_collection_ins = [[1,common_ins_length]]
                #            left_collection.append(int(basss_index-shif_l))
                #            shif_l = shif_l + 1
    
                #else:
                #    shif_l = 1
                #    while basss_index-shif_l >= start_block:
                #        if state_mapping[int(basss_index-shif_l-start_block)] == 1:
                #            shif_l = shif_l + 1
                #        else:
                #            left_collection.append(int(basss_index-shif_l))
                #            break
    
            # right-ward direction
            if basss_index == end_block_id:
                right_collection = []
            else:
                right_collection = []
                if state_now == 2:
                    shif_l = 1
                    while basss_index+shif_l <= end_block_id:
                        if state_mapping[int(basss_index+shif_l-start_block_id)] == 2:
                            break
                        else:
                            if var_site_collection_ori[basss_index+shif_l-start_block_id] - var_site_collection_ori[basss_index-start_block_id]<2:
                                right_collection.append(int(basss_index+shif_l))
                                shif_l = shif_l + 1
                            else:
                                break
    
                else:
                    shif_l = 1
                    while basss_index+shif_l <= end_block_id:
                        if state_mapping[int(basss_index+shif_l-start_block_id)] == 1:
                            if  var_site_collection_ori[basss_index+shif_l-start_block_id] - var_site_collection_ori[basss_index-start_block_id]<2:
                                shif_l = shif_l + 1
                            else:
                                break
                        else:
                            right_collection.append(int(basss_index+shif_l))
                            break
            final_dep_collec = left_collection+right_collection
            if len(final_dep_collec) > 0:
                state_dep_dic[basss_index] = final_dep_collec
        updated_state_mapping = copy.deepcopy(state_mapping)
        for basss_index in range(start_block_id, end_block_id+1):
            # print(updated_segment)
            #basss_index = start_block_id
            potential_probability = np.zeros(5)
            updated_segment_cp = copy.deepcopy(updated_segment)
            #updated_state_mapping_cp = copy.deepcopy()
            # if
            present_var_msa = var_site_collection_abs[basss_index-start_block_id]
            present_var_ori = var_site_collection_ori[basss_index-start_block_id]
            if infer_var_ori_state[int(basss_index)] > -1:  # ins in the bone
                if basss_index in state_dep_dic:
                    tem_vei = state_dep_dic[basss_index]  # neighbour
                    ll_xun_ = 0
                    case_flag = 1
                    while ll_xun_ < len(tem_vei):
                        if updated_segment_cp[var_site_collection_abs[tem_vei[ll_xun_]-start_block_id]] == 4:
                            case_flag = 0
                            break
                        else:
                            ll_xun_ = ll_xun_ + 1
                    ### check the common part for del
                    if len(common_del_per_df_list_df)>0:
                        left_boundary_del = np.array(common_del_per_df_list_df.iloc[:,1])
                        left_boundary_del_dis = present_var_msa - left_boundary_del
                        if  1 in left_boundary_del_dis:
                            case_flag = 0
                     
                        right_boundary_del = np.array(common_del_per_df_list_df.iloc[:,0])
                        right_boundary_del_dis = right_boundary_del - present_var_msa
                        if 1 in right_boundary_del_dis:
                            case_flag = 0
                    
                    
                else:
                    case_flag = 1
    
                if case_flag == 1:
    
                    for main_base_xun in [0, 1, 2, 3, 5]:
                        if blobk_base_number[int(basss_index-start_block_id),main_base_xun]>0:
                            updated_segment_cp[present_var_msa] = main_base_xun
                            full_full_prob = 0
                            for read_xun_tem in range(read_number_tem):
                                tem_dic = segment_reads[int(read_xun_tem)][2]
                                #begin_id = max(
                                #    segment_reads[int(read_xun_tem)][0], start_block)
                                #end_id = min(
                                #segment_reads[int(read_xun_tem)][1], end_block)
        
                                num_correct = 0
                                num_mis = 0
                                num_del = 0
                                num_ins = 0
                                len_del = {}
                                len_ins = {}
                                start_del_flag = 0
                                end_del_flag = 0
                                #begin_id = 0
                                #end_id = len(var_site_collection_abs)-1
                                #base_xun = int(begin_id)
                                pre_ins_loc = -1
                                for  base_xun  in range(len(var_site_collection_abs)):
                                    xun_base_present_abs  = var_site_collection_abs[int(base_xun)]
                                    if updated_segment_cp[xun_base_present_abs] < 4:
                                        #bone_flag = 1
                                        if xun_base_present_abs in tem_dic:
                                            if tem_dic[xun_base_present_abs] == updated_segment_cp[xun_base_present_abs]:
            
                                                if num_del > 0:
                                                    if num_del not in len_del:
                                                        len_del[num_del] = 1
                                                    else:
                                                        len_del[num_del] = len_del[num_del] + 1
                                                    num_del = 0
            
                                                if num_ins > 0:
                                                    if num_ins not in len_ins:
                                                        len_ins[num_ins] = 1
                                                    else:
                                                        len_ins[num_ins] = len_ins[num_ins] + 1
                                                    num_ins = 0
                                                    pre_ins_loc = -1
                                                num_correct = num_correct + 1
            
                                            elif tem_dic[xun_base_present_abs] < 4 and tem_dic[xun_base_present_abs] != updated_segment_cp[xun_base_present_abs]:
                                                if num_del > 0:
                                                    if num_del not in len_del:
                                                        len_del[num_del] = 1
                                                    else:
                                                        len_del[num_del] = len_del[num_del] + 1
                                                    num_del = 0
                                                if num_ins > 0:
                                                    if num_ins not in len_ins:
                                                        len_ins[num_ins] = 1
                                                    else:
                                                        len_ins[num_ins] = len_ins[num_ins] + 1
                                                    num_ins = 0
                                                    pre_ins_loc = -1
                                                num_mis = num_mis + 1
                                            else:
                                                if num_ins > 0:
                                                    if num_ins not in len_ins:
                                                        len_ins[num_ins] = 1
                                                    else:
                                                        len_ins[num_ins] = len_ins[num_ins] + 1
                                                    num_ins = 0
                                                    pre_ins_loc = -1
                                                    
            
                                                num_del = num_del + 1
                                                #if base_xun == begin_id:
                                                #    start_del_flag = 1
                                                #if base_xun == end_id:
                                                #    end_del_flag = 1
                                    elif updated_segment_cp[xun_base_present_abs] == 4:
                                        if xun_base_present_abs in tem_dic:
                                            if tem_dic[xun_base_present_abs] != 4:
                                                if tem_dic[xun_base_present_abs] < 4:
                                                    if xun_base_present_abs-pre_ins_loc<2:
                                                        num_ins = num_ins + 1
                                                    else:
                                                        if num_ins > 0:
                                                            if num_ins not in len_ins:
                                                                len_ins[num_ins] = 1
                                                            else:
                                                                len_ins[num_ins] = len_ins[num_ins] + 1
                                                            num_ins = 0
                                                        pre_ins_loc = xun_base_present_abs
                                                else:
                                                    print('error10')
                                                    print(basss_index)
                                                    #break
                                    else:
                                        if xun_base_present_abs in tem_dic:
                                            if tem_dic[xun_base_present_abs] == 4:
                                                print('error2')
                                            elif tem_dic[xun_base_present_abs] < 4:
                                                num_ins = num_ins + 1
                                    #base_xun = base_xun + 1
        
                                
                                if num_ins>0:
                                    if num_ins in len_ins:
                                        len_ins[num_ins] = len_ins[num_ins] + 1
                                    else:
                                        len_ins[num_ins] = 1
                                if num_del>0:
                                    if num_del in len_del:
                                        len_del[num_del] = len_del[num_del] + 1
                                    else:
                                        len_del[num_del] = 1
                                tem_probabity_correct = num_correct*np.log(pcorrect+corretion_eps_from_zero)
                                tem_probabity_mis = num_mis*np.log(p_mis_loc+corretion_eps_from_zero)
                                if len(len_ins) > 0:
                                    values_ = len_ins.values()
                                    num_ins_loc = sum(np.array(list(values_)))
                                    tem_probabity_ins = num_ins_loc*np.log(pins1_loc+corretion_eps_from_zero)
                                    for ii in len_ins:
                                        # np.log(stats.expon.pdf(ii,lamuna_ins))
                                        tem_probabity_ins = tem_probabity_ins + \
                                            len_ins[ii] * \
                                            possion_log_probability(ii-1,lamuna_ins)
                                        # print(tem_probabity_ins)
                                        # tem_probabity_ins = tem_probabity_ins + len_ins[ii]*np.log(stats.nbinom.pmf(ii,n_ins,p_ins))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                        # tem_probabity_ins = tem_probabity_ins + len_ins[ii]*np.log(1/(ii*sig_ins*np.sqrt(2*math.pi))*np.exp(-(np.log(ii)-mu_ins)**2/(2*sig_ins**2)))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                else:
                                    tem_probabity_ins = 0
                                    num_ins_loc = 0
        
                                del_whole = 0
                                if len(len_del) > 0:
                                    values_ = len_del.values()
                                    num_del_loc = sum(np.array(list(values_)))
                                    tem_probabity_del = num_del_loc*np.log(pdel1_loc+corretion_eps_from_zero)
                                    for dd in len_del:
                                        # np.log(stats.expon.pdf(ii,lamuna_ins))
                                        tem_probabity_del = tem_probabity_del + \
                                            len_del[dd] * \
                                             possion_log_probability(dd-1,lamuna_del)
                                        # tem_probabity_del = tem_probabity_del + len_del[dd]*np.log(stats.nbinom.pmf(dd,n_del,p_del))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                        # tem_probabity_del = tem_probabity_del + len_del[dd]*np.log(1/(dd*sig_del*np.sqrt(2*math.pi))*np.exp(-(np.log(dd)-mu_del)**2/(2*sig_del**2)))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                        del_whole = del_whole + dd*len_del[dd]
        
                                else:
                                    tem_probabity_del = 0
                                    num_del_loc = 0
        
                                other_gaps = num_correct + num_mis - 1-num_del_loc+start_del_flag + end_del_flag
                                gap_not_ins_pro = other_gaps*np.log(1-pins1_loc+corretion_eps_from_zero)
                                not_del_num = num_correct + num_mis
                                gap_not_del_pro = not_del_num*np.log(1-pdel1_loc+corretion_eps_from_zero)
                                full_prob = gap_not_del_pro+tem_probabity_correct + tem_probabity_mis +\
                                    tem_probabity_ins+tem_probabity_del+gap_not_ins_pro
        
                                full_full_prob = full_full_prob + full_prob
                            
                            ### evo comparision
                            
    
                            ID_tem = [0, 1, 2, 3, 5].index(main_base_xun)
                            ### here we should judge the influence of the variations
                            
                            
                            
                            
                            Tem_ins_evo = {}
                            Tem_del_evo = {}
                            ins_len_evo = 0
                            del_len_evo = 0
                            num_mis_evo = 0
                            num_cor_evo = 0 
                            modify_reads_con_prob = 0
                            for ll_check in range(start_block_id, end_block_id+1):
                                relative_id =  ll_check - start_block_id
                                
                                if original_true_state_ref[relative_id]<0:## meaning this should ins evo fromthe refence
                                    
                                    if del_len_evo>0:
                                        if del_len_evo in Tem_del_evo:
                                            Tem_del_evo[del_len_evo] = Tem_del_evo[del_len_evo] + 1
                                        else:
                                            Tem_del_evo[del_len_evo] = 1
                                        del_len_evo = 0
                                    ### common part for ins
                                    if relative_id<1:
                                        if len(common_ins_per_df_list_df)>0:
                                            left_boundary_ins = np.array(common_ins_per_df_list_df.iloc[:,1])
                                            left_boundary_ins_dis = present_var_msa - left_boundary_ins
                                            if  1 in left_boundary_ins_dis:
                                                common_ins_len_tem_id =  list(left_boundary_ins_dis).index(1)
                                                common_ins_len_ = common_ins_per_df_list_df.iloc[common_ins_len_tem_id,1] - common_ins_per_df_list_df.iloc[common_ins_len_tem_id,0] + 1
                                                ins_len_evo = ins_len_evo
                                                modify_reads_con_prob =  modify_reads_con_prob + common_ins_len_*read_number_tem*(np.log(pcorrect)+np.log(1-p_del_loc)+np.log(1-pins1_loc))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                                           
                                        if len(common_del_per_df_list_df)>0:
                                            left_boundary_del = np.array(common_del_per_df_list_df.iloc[:,1])
                                            left_boundary_del_dis = present_var_msa - left_boundary_del
                                            if  1 in left_boundary_del_dis:
                                                common_del_len_tem_id =  list(left_boundary_del_dis).index(1)
                                                common_del_len_ = common_del_per_df_list_df.iloc[common_del_len_tem_id,1] - common_del_per_df_list_df.iloc[common_del_len_tem_id,0] + 1
                                                del_len_evo = del_len_evo + common_del_len_
                                                modify_reads_con_prob = modify_reads_con_prob + read_number_tem*(np.log(1-pins1_loc))
                                         
                                            
                                    
                                    if updated_state_mapping[relative_id]>0:
                                        ins_len_evo = ins_len_evo + 1
                                    
                                    if len(common_ins_per_df_list_df)>0:
                                        right_boundary_ins = np.array(common_ins_per_df_list_df.iloc[:,0])
                                        right_boundary_ins_dis = present_var_msa - right_boundary_ins
                                        if  1 in right_boundary_ins_dis:
                                            common_ins_len_tem_id =  list(right_boundary_ins_dis).index(1)
                                            common_ins_len_ = common_ins_per_df_list_df.iloc[common_ins_len_tem_id,1] - common_ins_per_df_list_df.iloc[common_ins_len_tem_id,0] + 1    
                                            ins_len_evo = ins_len_evo + common_ins_len_
                                            modify_reads_con_prob = modify_reads_con_prob + common_ins_len_*read_number_tem*(np.log(pcorrect)+np.log(1-p_del_loc)+np.log(1-pins1_loc))                            
                                     
                                    
                                    
                                else:
                                    if ins_len_evo>0:
                                        if ins_len_evo in Tem_ins_evo:
                                            Tem_ins_evo[ins_len_evo] = Tem_ins_evo[ins_len_evo] + 1
                                        else:
                                            Tem_ins_evo[ins_len_evo] =  1
                                        ins_len_evo = 0
                                    ###start position
                                    if relative_id<1:
                                        if len(common_ins_per_df_list_df)>0:
                                            left_boundary_ins = np.array(common_ins_per_df_list_df.iloc[:,1])
                                            left_boundary_ins_dis = present_var_msa - left_boundary_ins
                                            if  1 in left_boundary_ins_dis:
                                                common_ins_len_tem_id =  list(left_boundary_ins_dis).index(1)
                                                common_ins_len_ = common_ins_per_df_list_df.iloc[common_ins_len_tem_id,1] - common_ins_per_df_list_df.iloc[common_ins_len_tem_id,0] + 1
                                                ins_len_evo = ins_len_evo + common_ins_len_
                                                modify_reads_con_prob =  modify_reads_con_prob + common_ins_len_*read_number_tem*(np.log(pcorrect)+np.log(1-p_del_loc)+np.log(1-pins1_loc))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                                           
                                        if len(common_del_per_df_list_df)>0:
                                            left_boundary_del = np.array(common_del_per_df_list_df.iloc[:,1])
                                            left_boundary_del_dis = present_var_msa - left_boundary_del
                                            if  1 in left_boundary_del_dis:
                                                common_del_len_tem_id =  list(left_boundary_del_dis).index(1)
                                                common_del_len_ = common_del_per_df_list_df.iloc[common_del_len_tem_id,1] - common_del_per_df_list_df.iloc[common_del_len_tem_id,0] + 1
                                                del_len_evo = del_len_evo + common_del_len_
                                                modify_reads_con_prob = modify_reads_con_prob + read_number_tem*(np.log(1-pins1_loc))
                                     #### check the case if possible    
                                    
                                    if updated_state_mapping[relative_id]>0:
                                        del_len_evo = del_len_evo + 1
                                        
                                    else:
                                        if ll_check == basss_index:
                                            if main_base_xun == standard_ref_code_block[relative_id]:
                                                num_cor_evo = num_cor_evo + 1
                                            else:
                                                num_mis_evo = num_mis_evo + 1
                                        else:
                                            if updated_segment_cp[var_site_collection_abs[relative_id]] ==  standard_ref_code_block[relative_id]:
                                                num_cor_evo = num_cor_evo + 1
                                            else:
                                                num_mis_evo = num_mis_evo + 1
                                    ### del part
                                    if len(common_del_per_df_list_df)>0:
                                        right_boundary_del = np.array(common_del_per_df_list_df.iloc[:,0])
                                        right_boundary_del_dis = right_boundary_del - present_var_msa
                                        if  1 in right_boundary_del_dis:
                                            common_del_len_tem_id =  list(right_boundary_del_dis).index(1)
                                            common_del_len_ = common_del_per_df_list_df.iloc[common_del_len_tem_id,1] - common_del_per_df_list_df.iloc[common_del_len_tem_id,0] + 1
                                            del_len_evo = del_len_evo + common_del_len_
                            if del_len_evo>0:
                                if del_len_evo in Tem_del_evo:
                                    Tem_del_evo[del_len_evo] = Tem_del_evo[del_len_evo] + 1
                                else:
                                    Tem_del_evo[del_len_evo] = 1
                            if ins_len_evo>0:
                                if ins_len_evo in Tem_ins_evo:
                                    Tem_ins_evo[ins_len_evo] = Tem_ins_evo[ins_len_evo] + 1
                                else:
                                    Tem_ins_evo[ins_len_evo] =  1
                            ### infer the prob
                            tem_evo = 0
                            
                            for del_per in Tem_del_evo:
                                tem_evo = tem_evo + Tem_del_evo[del_per]*(expotential_log_log_probability(del_per-1,del_p[kk],del_beta[kk]) + np.log(gamma_del+corretion_eps_from_zero) )
                            for ins_per in Tem_ins_evo:
                                tem_evo = tem_evo + Tem_ins_evo[ins_per]*(expotential_log_log_probability(ins_per-1,ins_p[kk],ins_beta[kk]) + np.log(gamma_ins+corretion_eps_from_zero) )
                            tem_evo = tem_evo + (num_cor_evo+ num_mis_evo)*np.log((1-gamma_del)) + (num_cor_evo)*np.log(gamma_c) + (num_mis_evo)*np.log(gamma_mis+corretion_eps_from_zero)              
                           
                            
                            potential_probability[int(ID_tem)] = full_full_prob + tem_evo +modify_reads_con_prob
                        else:
                            ID_tem = [0, 1, 2, 3, 5].index(main_base_xun)
                            potential_probability[int(ID_tem)] = float('-inf')
                    #tem_probability = copy.deepcopy(potential_probability/tem_pri)
                    #tem_probability = tem_probability - max(tem_probability)
                    #tem_probability_main = np.exp(tem_probability)
                    #tem_probability_main = tem_probability_main/sum(tem_probability_main)
                    #post_sampling = np.random.multinomial(n= 1,pvals = tem_probability_main)
                    #max_ID = list(post_sampling).index(1)
                    #    n=1, pvals=[6/8, 1/8, 1/16, 1/32, 1/32])
                    #l_i = list(l_i_vec).index(1)+1
                    if max(potential_probability)>float('-inf'):
                        max_ID = list(potential_probability).index(max(potential_probability))
                        
                        if max_ID >3:
                            updated_state_mapping[basss_index- start_block_id] = -1
                        else:
                            updated_state_mapping[basss_index- start_block_id] = 1
                            
                        if max_ID == 4:
                            updated_segment[var_site_collection_abs[basss_index- start_block_id]] = 5
                    else:
                         updated_segment[var_site_collection_abs[basss_index- start_block_id]] = 5
                         updated_state_mapping[basss_index- start_block_id] = -1
                else:
                    updated_segment[var_site_collection_abs[basss_index- start_block_id]] = 5
                    updated_state_mapping[basss_index- start_block_id] = -1
            else:  # del condition
                if basss_index in state_dep_dic:
                    tem_vei = state_dep_dic[basss_index]
                    ll_xun_ = 0
                    case_flag = 1
                    while ll_xun_ < len(tem_vei):
                        if updated_segment_cp[var_site_collection_abs[tem_vei[ll_xun_]-start_block_id]] < 4:## this should be the insertion case
                            case_flag = 0
                            break
                        else:
                            ll_xun_ = ll_xun_ + 1
                    ### check the common part for del
                    if len(common_ins_per_df_list_df)>0:
                        left_boundary_ins = np.array(common_ins_per_df_list_df.iloc[:,1])
                        left_boundary_ins_dis = present_var_msa - left_boundary_ins
                        if  1 in left_boundary_ins_dis:
                            case_flag = 0
                     
                        right_boundary_ins = np.array(common_ins_per_df_list_df.iloc[:,0])
                        right_boundary_ins_dis = right_boundary_ins - present_var_msa
                        if 1 in right_boundary_ins_dis:
                            case_flag = 0
                        ## 0 flag means we have less freedom to decide it
                else:
                    case_flag = 1
    
                if case_flag:
                    for main_base_xun in [0, 1, 2, 3, 4]:
                        if blobk_base_number[int(basss_index-start_block_id),main_base_xun]>0:
                            updated_segment_cp[present_var_msa] = main_base_xun
                            full_full_prob = 0
                            
                            for read_xun_tem in range(read_number_tem):
                                tem_dic = segment_reads[int(read_xun_tem)][2]
                                #begin_id = max(
                                #    segment_reads[int(read_xun_tem)][0], start_block)
                                #end_id = min(
                                #segment_reads[int(read_xun_tem)][1], end_block)
        
                                num_correct = 0
                                num_mis = 0
                                num_del = 0
                                num_ins = 0
                                len_del = {}
                                len_ins = {}
                                start_del_flag = 0
                                end_del_flag = 0
                                #begin_id = 0
                                #end_id = len(var_site_collection_abs)-1
                                #base_xun = int(begin_id)
                                pre_ins_loc = -1
                                #sorted_keys_read = sorted(list(tem_dic.keys()))
                                '''
                                if len(var_site_collection_abs)>1:
                                    ### judge the initial
                                    present_var_msa_xun = var_site_collection_abs[0]
                                    if original_true_state_ref[0]<0:## meaning this should ins evo from the refence
                                        
                                        if del_len>0:
                                            if del_len in Tem_del_evo:
                                                Tem_del_evo[del_len] = Tem_del_evo[del_len] + 1
                                            else:
                                                Tem_del_evo[del_len] = 1
                                            del_len = 0
                                        
                                        
                                        ### common part for ins
                                        #if relative_id<1:
                                        if len(common_ins_per_df_list_df)>0:
                                            left_boundary_ins = np.array(common_ins_per_df_list_df.iloc[:,1])
                                            left_boundary_ins_dis = present_var_msa_xun - left_boundary_ins
                                            if  1 in left_boundary_ins_dis:
                                                common_ins_len_tem_id =  list(left_boundary_ins_dis).index(1)
                                                common_ins_len_ = common_ins_per_df_list_df.iloc[common_ins_len_tem_id,1] - common_ins_per_df_list_df.iloc[common_ins_len_tem_id,0] + 1
                                                ins_len_evo = ins_len_evo + common_ins_len_
                                                #modify_reads_con_prob =  modify_reads_con_prob + common_ins_len_*read_number_tem*(np.log(pcorrect)+np.log(1-p_del_loc)+np.log(1-pins1_loc))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                                           
                                        if len(common_del_per_df_list_df)>0:
                                            left_boundary_del = np.array(common_del_per_df_list_df.iloc[:,1])
                                            left_boundary_del_dis = present_var_msa_xun - left_boundary_del
                                            if  1 in left_boundary_del_dis:
                                                common_del_len_tem_id =  list(left_boundary_del_dis).index(1)
                                                common_del_len_ = common_del_per_df_list_df.iloc[common_del_len_tem_id,1] - common_del_per_df_list_df.iloc[common_del_len_tem_id,0] + 1
                                                del_len = del_len + common_del_len_
                                                #modify_reads_con_prob = modify_reads_con_prob + read_number_tem*(np.log(1-pins1_loc))
                                         
                                        
                                        if updated_state_mapping[relative_id]>0:
                                            ins_len_evo = ins_len_evo + 1
                                        
                                        if len(common_ins_per_df_list_df)>0:
                                            right_boundary_ins = np.array(common_ins_per_df_list_df.iloc[:,0])
                                            right_boundary_ins_dis = present_var_msa_xun - right_boundary_ins
                                            if  1 in right_boundary_ins_dis:
                                                common_ins_len_tem_id =  list(right_boundary_ins_dis).index(1)
                                                common_ins_len_ = common_ins_per_df_list_df.iloc[common_ins_len_tem_id,1] - common_ins_per_df_list_df.iloc[common_ins_len_tem_id,0] + 1    
                                                ins_len_evo = ins_len_evo + common_ins_len_
                                                #modify_reads_con_prob = modify_reads_con_prob + common_ins_len_*read_number_tem*(np.log(pcorrect)+np.log(1-p_del_loc)+np.log(1-pins1_loc))                            
                                                
                                             
                                        
                                    else:
                                        if ins_len_evo>0:
                                            if ins_len_evo in Tem_ins_evo:
                                                Tem_ins_evo[ins_len_evo] = Tem_ins_evo[ins_len_evo] + 1
                                            else:
                                                Tem_ins_evo[ins_len_evo] =  1
                                            ins_len_evo = 0
                                        ###start position
                                        #if relative_id<1:
                                        if len(common_ins_per_df_list_df)>0:
                                            left_boundary_ins = np.array(common_ins_per_df_list_df.iloc[:,1])
                                            left_boundary_ins_dis = present_var_msa_xun - left_boundary_ins
                                            if  1 in left_boundary_ins_dis:
                                                common_ins_len_tem_id =  list(left_boundary_ins_dis).index(1)
                                                common_ins_len_ = common_ins_per_df_list_df.iloc[common_ins_len_tem_id,1] - common_ins_per_df_list_df.iloc[common_ins_len_tem_id,0] + 1
                                                ins_len_evo = ins_len_evo + common_ins_len_
                                                #modify_reads_con_prob =  modify_reads_con_prob + common_ins_len_*read_number_tem*(np.log(pcorrect)+np.log(1-p_del_loc)+np.log(1-pins1_loc))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                                           
                                        if len(common_del_per_df_list_df)>0:
                                            left_boundary_del = np.array(common_del_per_df_list_df.iloc[:,1])
                                            left_boundary_del_dis = present_var_msa_xun - left_boundary_del
                                            if  1 in left_boundary_del_dis:
                                                common_del_len_tem_id =  list(left_boundary_del_dis).index(1)
                                                common_del_len_ = common_del_per_df_list_df.iloc[common_del_len_tem_id,1] - common_del_per_df_list_df.iloc[common_del_len_tem_id,0] + 1
                                                del_len = del_len + common_del_len_
                                                #modify_reads_con_prob = modify_reads_con_prob + read_number_tem*(np.log(1-pins1_loc))
                                         #### check the case if possible 
                                        if updated_state_mapping[0]>0:
                                            del_len_evo = del_len_evo + 1
                                        else:
                                            #relative_id = 
                                            if present_var_msa_xun == infer_var_site[basss_index]:
                                                if main_base_xun == standard_ref_code_block[0]:
                                                    num_cor_evo = num_cor_evo + 1
                                                else:
                                                    num_mis_evo = num_mis_evo + 1
                                            else:
                                                if  updated_segment_cp[present_var_msa_xun] ==  standard_ref_code_block[0]:
                                                    num_cor_evo = num_cor_evo + 1
                                                else:
                                                    num_mis_evo = num_mis_evo + 1
                                        ### del part
                                        if len(common_del_per_df_list_df)>0:
                                            right_boundary_del = np.array(common_del_per_df_list_df.iloc[:,0])
                                            right_boundary_del_dis = right_boundary_del - present_var_msa_xun
                                            if  1 in right_boundary_del_dis:
                                                common_del_len_tem_id =  list(right_boundary_del_dis).index(1)
                                                common_del_len_ = common_del_per_df_list_df.iloc[common_del_len_tem_id,1] - common_del_per_df_list_df.iloc[common_del_len_tem_id,0] + 1
                                                del_len_evo = del_len_evo + common_del_len_
                                                #modify_reads_con_prob = modify_reads_con_prob + read_number_tem*(np.log(1-pins1_loc))
                                       
                                    
                                    
                                    ### >2
                                    for base_xun_id in range(1,len(var_site_collection_abs)):
                                        present_var_msa_xun = var_site_collection_abs[base_xun_id]
                                        if original_true_state_ref[base_xun_id]<0:## meaning this should ins evo from the refence
                                            
                                            if del_len_evo>0:
                                                if del_len_evo in Tem_del_evo:
                                                    Tem_del_evo[del_len_evo] = Tem_del_evo[del_len_evo] + 1
                                                else:
                                                    Tem_del_evo[del_len_evo] = 1
                                                del_len_evo = 0
                                            
                                           
                                            if updated_state_mapping[base_xun_id]>0:
                                                ins_len_evo = ins_len_evo + 1
                                            
                                            if len(common_ins_per_df_list_df)>0:
                                                right_boundary_ins = np.array(common_ins_per_df_list_df.iloc[:,0])
                                                right_boundary_ins_dis = present_var_msa_xun - right_boundary_ins
                                                if  1 in right_boundary_ins_dis:
                                                    common_ins_len_tem_id =  list(right_boundary_ins_dis).index(1)
                                                    common_ins_len_ = common_ins_per_df_list_df.iloc[common_ins_len_tem_id,1] - common_ins_per_df_list_df.iloc[common_ins_len_tem_id,0] + 1    
                                                    ins_len_evo = ins_len_evo + common_ins_len_
                                                    #modify_reads_con_prob = modify_reads_con_prob + common_ins_len_*read_number_tem*(np.log(pcorrect)+np.log(1-p_del_loc)+np.log(1-pins1_loc))                            
                                                    
                                                 
                                            
                                        else:
                                            if ins_len_evo>0:
                                                if ins_len_evo in Tem_ins_evo:
                                                    Tem_ins_evo[ins_len_evo] = Tem_ins_evo[ins_len_evo] + 1
                                                else:
                                                    Tem_ins_evo[ins_len_evo] =  1
                                                ins_len_evo = 0
                                           
                                             #### check the case if possible 
                                            if updated_state_mapping[base_xun_id]>0:
                                                del_len_evo = del_len_evo + 1
                                            else:
                                                #relative_id = 
                                                if present_var_msa_xun == infer_var_site[basss_index]:
                                                    if main_base_xun == standard_ref_code_block[base_xun_id]:
                                                        num_cor_evo = num_cor_evo + 1
                                                    else:
                                                        num_mis_evo = num_mis_evo + 1
                                                else:
                                                    if  updated_segment_cp[present_var_msa_xun] ==  standard_ref_code_block[base_xun_id]:
                                                        num_cor_evo = num_cor_evo + 1
                                                    else:
                                                        num_mis_evo = num_mis_evo + 1
                                            ### del part
                                            if len(common_del_per_df_list_df)>0:
                                                right_boundary_del = np.array(common_del_per_df_list_df.iloc[:,0])
                                                right_boundary_del_dis = right_boundary_del - present_var_msa_xun
                                                if  1 in right_boundary_del_dis:
                                                    common_del_len_tem_id =  list(right_boundary_del_dis).index(1)
                                                    common_del_len_ = common_del_per_df_list_df.iloc[common_del_len_tem_id,1] - common_del_per_df_list_df.iloc[common_del_len_tem_id,0] + 1
                                                    del_len_evo = del_len_evo + common_del_len_
                                                    #modify_reads_con_prob = modify_reads_con_prob + read_number_tem*(np.log(1-pins1_loc))
                                           
                                        
                                    
                                '''
                                for  base_xun  in range(len(var_site_collection_abs)):
                                    xun_base_present_abs  = var_site_collection_abs[base_xun]
                                    if updated_segment_cp[xun_base_present_abs] < 4:
                                        #bone_flag = 1
                                        if xun_base_present_abs in tem_dic:
                                            if tem_dic[xun_base_present_abs] == updated_segment_cp[xun_base_present_abs]:
            
                                                if num_del > 0:
                                                    if num_del not in len_del:
                                                        len_del[num_del] = 1
                                                    else:
                                                        len_del[num_del] = len_del[num_del] + 1
                                                    num_del = 0
            
                                                if num_ins > 0:
                                                    if num_ins not in len_ins:
                                                        len_ins[num_ins] = 1
                                                    else:
                                                        len_ins[num_ins] = len_ins[num_ins] + 1
                                                    num_ins = 0
                                                    pre_ins_loc = -1
                                                num_correct = num_correct + 1
            
                                            elif tem_dic[xun_base_present_abs] < 4 and tem_dic[xun_base_present_abs] != updated_segment_cp[xun_base_present_abs]:
                                                if num_del > 0:
                                                    if num_del not in len_del:
                                                        len_del[num_del] = 1
                                                    else:
                                                        len_del[num_del] = len_del[num_del] + 1
                                                    num_del = 0
                                                if num_ins > 0:
                                                    if num_ins not in len_ins:
                                                        len_ins[num_ins] = 1
                                                    else:
                                                        len_ins[num_ins] = len_ins[num_ins] + 1
                                                    num_ins = 0
                                                    pre_ins_loc = -1
                                                num_mis = num_mis + 1
                                            else:
                                                if num_ins > 0:
                                                    if num_ins not in len_ins:
                                                        len_ins[num_ins] = 1
                                                    else:
                                                        len_ins[num_ins] = len_ins[num_ins] + 1
                                                    num_ins = 0
                                                    pre_ins_loc = -1
                                                    
            
                                                num_del = num_del + 1
                                                #if base_xun == begin_id:
                                                #    start_del_flag = 1
                                                #if base_xun == end_id:
                                                #    end_del_flag = 1
                                    elif updated_segment_cp[xun_base_present_abs] == 4:
                                        if xun_base_present_abs in tem_dic:
                                            if tem_dic[xun_base_present_abs] != 4:
                                                if tem_dic[xun_base_present_abs] < 4:
                                                    if xun_base_present_abs-pre_ins_loc<2:
                                                        num_ins = num_ins + 1
                                                    else:
                                                        if num_ins > 0:
                                                            if num_ins not in len_ins:
                                                                len_ins[num_ins] = 1
                                                            else:
                                                                len_ins[num_ins] = len_ins[num_ins] + 1
                                                            num_ins = 0
                                                        pre_ins_loc = xun_base_present_abs
                                                else:
                                                    print('error11')
                                                    print(basss_index)
                                                    #break
                                    else:
                                        if xun_base_present_abs in tem_dic:
                                            if tem_dic[xun_base_present_abs] == 4:
                                                print('error2')
                                            elif tem_dic[xun_base_present_abs] < 4:
                                                num_ins = num_ins + 1
                                    #base_xun = base_xun + 1
        
                                
                                if num_ins>0:
                                    if num_ins in len_ins:
                                        len_ins[num_ins] = len_ins[num_ins] + 1
                                    else:
                                        len_ins[num_ins] = 1
                                if num_del>0:
                                    if num_del in len_del:
                                        len_del[num_del] = len_del[num_del] + 1
                                    else:
                                        len_del[num_del] = 1
                                tem_probabity_correct = num_correct*np.log(pcorrect+corretion_eps_from_zero)
                                tem_probabity_mis = num_mis*np.log(p_mis_loc+corretion_eps_from_zero)
                                if len(len_ins) > 0:
                                    values_ = len_ins.values()
                                    num_ins_loc = sum(np.array(list(values_)))
                                    tem_probabity_ins = num_ins_loc*np.log(pins1_loc+corretion_eps_from_zero)
                                    for ii in len_ins:
                                        # np.log(stats.expon.pdf(ii,lamuna_ins))
                                        tem_probabity_ins = tem_probabity_ins + \
                                            len_ins[ii] * \
                                            possion_log_probability(ii-1,lamuna_ins)
                                        # print(tem_probabity_ins)
                                        # tem_probabity_ins = tem_probabity_ins + len_ins[ii]*np.log(stats.nbinom.pmf(ii,n_ins,p_ins))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                        # tem_probabity_ins = tem_probabity_ins + len_ins[ii]*np.log(1/(ii*sig_ins*np.sqrt(2*math.pi))*np.exp(-(np.log(ii)-mu_ins)**2/(2*sig_ins**2)))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                else:
                                    tem_probabity_ins = 0
                                    num_ins_loc = 0
        
                                del_whole = 0
                                if len(len_del) > 0:
                                    values_ = len_del.values()
                                    num_del_loc = sum(np.array(list(values_)))
                                    tem_probabity_del = num_del_loc*np.log(pdel1_loc+corretion_eps_from_zero)
                                    for dd in len_del:
                                        # np.log(stats.expon.pdf(ii,lamuna_ins))
                                        tem_probabity_del = tem_probabity_del + \
                                            len_del[dd] * \
                                             possion_log_probability(dd-1,lamuna_del)
                                        # tem_probabity_del = tem_probabity_del + len_del[dd]*np.log(stats.nbinom.pmf(dd,n_del,p_del))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                        # tem_probabity_del = tem_probabity_del + len_del[dd]*np.log(1/(dd*sig_del*np.sqrt(2*math.pi))*np.exp(-(np.log(dd)-mu_del)**2/(2*sig_del**2)))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                        del_whole = del_whole + dd*len_del[dd]
        
                                else:
                                    tem_probabity_del = 0
                                    num_del_loc = 0
        
                                other_gaps = num_correct + num_mis - 1-num_del_loc+start_del_flag + end_del_flag
                                gap_not_ins_pro = other_gaps*np.log(1-pins1_loc+corretion_eps_from_zero)
                                not_del_num = num_correct + num_mis
                                gap_not_del_pro = not_del_num*np.log(1-pdel1_loc+corretion_eps_from_zero)
                                full_prob = gap_not_del_pro+tem_probabity_correct + tem_probabity_mis +\
                                    tem_probabity_ins+tem_probabity_del+gap_not_ins_pro
        
                                full_full_prob = full_full_prob + full_prob
                            
                            ### evo comparision
                            
    
                            ID_tem = [0, 1, 2, 3, 4].index(main_base_xun)
                            
                            Tem_ins_evo = {}
                            Tem_del_evo = {}
                            ins_len_evo = 0
                            del_len_evo = 0
                            num_mis_evo = 0
                            num_cor_evo = 0
                            modify_reads_con_prob = 0
                            for ll_check in range(start_block_id, end_block_id+1):
                                relative_id =  ll_check - start_block_id
                                if original_true_state_ref[relative_id]<0:## meaning this should ins evo from the refence
                                    
                                    if del_len_evo>0:
                                        if del_len_evo in Tem_del_evo:
                                            Tem_del_evo[del_len_evo] = Tem_del_evo[del_len_evo] + 1
                                        else:
                                            Tem_del_evo[del_len_evo] = 1
                                        del_len_evo = 0
                                    
                                    
                                    ### common part for ins
                                    if relative_id<1:
                                        if len(common_ins_per_df_list_df)>0:
                                            left_boundary_ins = np.array(common_ins_per_df_list_df.iloc[:,1])
                                            left_boundary_ins_dis = present_var_msa - left_boundary_ins
                                            if  1 in left_boundary_ins_dis:
                                                common_ins_len_tem_id =  list(left_boundary_ins_dis).index(1)
                                                common_ins_len_ = common_ins_per_df_list_df.iloc[common_ins_len_tem_id,1] - common_ins_per_df_list_df.iloc[common_ins_len_tem_id,0] + 1
                                                ins_len_evo = ins_len_evo + common_ins_len_
                                                modify_reads_con_prob =  modify_reads_con_prob + common_ins_len_*read_number_tem*(np.log(pcorrect)+np.log(1-p_del_loc)+np.log(1-pins1_loc))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                                           
                                        if len(common_del_per_df_list_df)>0:
                                            left_boundary_del = np.array(common_del_per_df_list_df.iloc[:,1])
                                            left_boundary_del_dis = present_var_msa - left_boundary_del
                                            if  1 in left_boundary_del_dis:
                                                common_del_len_tem_id =  list(left_boundary_del_dis).index(1)
                                                common_del_len_ = common_del_per_df_list_df.iloc[common_del_len_tem_id,1] - common_del_per_df_list_df.iloc[common_del_len_tem_id,0] + 1
                                                del_len_evo = del_len_evo + common_del_len_
                                                modify_reads_con_prob = modify_reads_con_prob + read_number_tem*(np.log(1-pins1_loc))
                                         
                                    
                                    if updated_state_mapping[relative_id]>0:
                                        ins_len_evo = ins_len_evo + 1
                                    
                                    if len(common_ins_per_df_list_df)>0:
                                        right_boundary_ins = np.array(common_ins_per_df_list_df.iloc[:,0])
                                        right_boundary_ins_dis = present_var_msa - right_boundary_ins
                                        if  1 in right_boundary_ins_dis:
                                            common_ins_len_tem_id =  list(right_boundary_ins_dis).index(1)
                                            common_ins_len_ = common_ins_per_df_list_df.iloc[common_ins_len_tem_id,1] - common_ins_per_df_list_df.iloc[common_ins_len_tem_id,0] + 1    
                                            ins_len_evo = ins_len_evo + common_ins_len_
                                            modify_reads_con_prob = modify_reads_con_prob + common_ins_len_*read_number_tem*(np.log(pcorrect)+np.log(1-p_del_loc)+np.log(1-pins1_loc))                            
                                            
                                         
                                    
                                else:
                                    if ins_len_evo>0:
                                        if ins_len_evo in Tem_ins_evo:
                                            Tem_ins_evo[ins_len_evo] = Tem_ins_evo[ins_len_evo] + 1
                                        else:
                                            Tem_ins_evo[ins_len_evo] =  1
                                        ins_len_evo = 0
                                    ###start position
                                    if relative_id<1:
                                        if len(common_ins_per_df_list_df)>0:
                                            left_boundary_ins = np.array(common_ins_per_df_list_df.iloc[:,1])
                                            left_boundary_ins_dis = present_var_msa - left_boundary_ins
                                            if  1 in left_boundary_ins_dis:
                                                common_ins_len_tem_id =  list(left_boundary_ins_dis).index(1)
                                                common_ins_len_ = common_ins_per_df_list_df.iloc[common_ins_len_tem_id,1] - common_ins_per_df_list_df.iloc[common_ins_len_tem_id,0] + 1
                                                ins_len_evo = ins_len_evo + common_ins_len_
                                                modify_reads_con_prob =  modify_reads_con_prob + common_ins_len_*read_number_tem*(np.log(pcorrect)+np.log(1-p_del_loc)+np.log(1-pins1_loc))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                                           
                                        if len(common_del_per_df_list_df)>0:
                                            left_boundary_del = np.array(common_del_per_df_list_df.iloc[:,1])
                                            left_boundary_del_dis = present_var_msa - left_boundary_del
                                            if  1 in left_boundary_del_dis:
                                                common_del_len_tem_id =  list(left_boundary_del_dis).index(1)
                                                common_del_len_ = common_del_per_df_list_df.iloc[common_del_len_tem_id,1] - common_del_per_df_list_df.iloc[common_del_len_tem_id,0] + 1
                                                del_len_evo = del_len_evo + common_del_len_
                                                modify_reads_con_prob = modify_reads_con_prob + read_number_tem*(np.log(1-pins1_loc))
                                     #### check the case if possible 
                                    if updated_state_mapping[relative_id]>0:
                                        del_len_evo = del_len_evo + 1
                                    else:
                                        if ll_check == basss_index:
                                            if main_base_xun == standard_ref_code_block[relative_id]:
                                                num_cor_evo = num_cor_evo + 1
                                            else:
                                                num_mis_evo = num_mis_evo + 1
                                        else:
                                            if  updated_segment_cp[var_site_collection_abs[relative_id]] ==  standard_ref_code_block[relative_id]:
                                                num_cor_evo = num_cor_evo + 1
                                            else:
                                                num_mis_evo = num_mis_evo + 1
                                    ### del part
                                    if len(common_del_per_df_list_df)>0:
                                        right_boundary_del = np.array(common_del_per_df_list_df.iloc[:,0])
                                        right_boundary_del_dis = right_boundary_del - present_var_msa
                                        if  1 in right_boundary_del_dis:
                                            common_del_len_tem_id =  list(right_boundary_del_dis).index(1)
                                            common_del_len_ = common_del_per_df_list_df.iloc[common_del_len_tem_id,1] - common_del_per_df_list_df.iloc[common_del_len_tem_id,0] + 1
                                            del_len_evo = del_len_evo + common_del_len_
                                            modify_reads_con_prob = modify_reads_con_prob + read_number_tem*(np.log(1-pins1_loc))
                            
                            if del_len_evo>0:
                                if del_len_evo in Tem_del_evo:
                                    Tem_del_evo[del_len_evo] = Tem_del_evo[del_len_evo] + 1
                                else:
                                    Tem_del_evo[del_len_evo] = 1
                            
                            if ins_len_evo>0:
                                if ins_len_evo in Tem_ins_evo:
                                    Tem_ins_evo[ins_len_evo] = Tem_ins_evo[ins_len_evo] + 1
                                else:
                                    Tem_ins_evo[ins_len_evo] =  1
                            ### infer the prob
                            tem_evo = 0
                            
                            for del_per in Tem_del_evo:
                                tem_evo = tem_evo + Tem_del_evo[del_per]*(expotential_log_log_probability(del_per-1,del_p[kk],del_beta[kk]) + np.log(gamma_del+corretion_eps_from_zero) )
                            
                            for ins_per in Tem_ins_evo:
                                tem_evo = tem_evo + Tem_ins_evo[ins_per]*(expotential_log_log_probability(ins_per-1,ins_p[kk],ins_beta[kk]) + np.log(gamma_ins+corretion_eps_from_zero) )
                            tem_evo = tem_evo + (num_cor_evo+ num_mis_evo)*np.log((1-gamma_del)) + (num_cor_evo)*np.log(gamma_c) + (num_mis_evo)*np.log(gamma_mis+corretion_eps_from_zero)     
                            
                            potential_probability[int(ID_tem)] = full_full_prob + tem_evo + modify_reads_con_prob
                        else:
                            ID_tem = [0, 1, 2, 3, 4].index(main_base_xun)
                            potential_probability[int(ID_tem)] = float('-inf')
                    #tem_probability = copy.deepcopy(potential_probability/tem_pri)
                    #tem_probability = tem_probability - max(tem_probability)
                    #tem_probability_main = np.exp(tem_probability)
                    #tem_probability_main = tem_probability_main/sum(tem_probability_main)
                    #post_sampling = np.random.multinomial(n= 1,pvals = tem_probability_main)
                    #max_ID = list(post_sampling).index(1)
                    #    n=1, pvals=[6/8, 1/8, 1/16, 1/32, 1/32])
                    #l_i = list(l_i_vec).index(1)+1
                    max_ID = list(potential_probability).index(max(potential_probability))
                    if max(potential_probability)>float('-inf'):
                        if max_ID >3:
                            updated_state_mapping[basss_index- start_block_id] = 2
                        else:
                            updated_state_mapping[basss_index- start_block_id] = 0
                        updated_segment[present_var_msa] = max_ID
                    else:
                        updated_state_mapping[basss_index- start_block_id] = 0
                        relative_id = basss_index- start_block_id
                        updated_segment[present_var_msa] = standard_ref_code_block[relative_id]
                    
    
                else:
                    potential_probability = np.zeros(4)
                    for main_base_xun in [0, 1, 2, 3]:
                        if blobk_base_number[int(basss_index-start_block_id),main_base_xun]>0:
                            updated_segment_cp[present_var_msa] = main_base_xun
                            full_full_prob = 0
                            for read_xun_tem in range(read_number_tem):
                                tem_dic = segment_reads[int(read_xun_tem)][2]
                                #begin_id = max(
                                #    segment_reads[int(read_xun_tem)][0], start_block)
                                #end_id = min(
                                #segment_reads[int(read_xun_tem)][1], end_block)
        
                                num_correct = 0
                                num_mis = 0
                                num_del = 0
                                num_ins = 0
                                len_del = {}
                                len_ins = {}
                                start_del_flag = 0
                                end_del_flag = 0
                                #begin_id = 0
                                #end_id = len(var_site_collection_abs)-1
                                #base_xun = int(begin_id)
                                pre_ins_loc = -1
                                for  base_xun  in range(len(var_site_collection_abs)):
                                    xun_base_present_abs  = var_site_collection_abs[int(base_xun)]
                                    if updated_segment_cp[xun_base_present_abs] < 4:
                                        #bone_flag = 1
                                        if xun_base_present_abs in tem_dic:
                                            if tem_dic[xun_base_present_abs] == updated_segment_cp[xun_base_present_abs]:
            
                                                if num_del > 0:
                                                    if num_del not in len_del:
                                                        len_del[num_del] = 1
                                                    else:
                                                        len_del[num_del] = len_del[num_del] + 1
                                                    num_del = 0
            
                                                if num_ins > 0:
                                                    if num_ins not in len_ins:
                                                        len_ins[num_ins] = 1
                                                    else:
                                                        len_ins[num_ins] = len_ins[num_ins] + 1
                                                    num_ins = 0
                                                    pre_ins_loc = -1
                                                num_correct = num_correct + 1
            
                                            elif tem_dic[xun_base_present_abs] < 4 and tem_dic[xun_base_present_abs] != updated_segment_cp[xun_base_present_abs]:
                                                if num_del > 0:
                                                    if num_del not in len_del:
                                                        len_del[num_del] = 1
                                                    else:
                                                        len_del[num_del] = len_del[num_del] + 1
                                                    num_del = 0
                                                if num_ins > 0:
                                                    if num_ins not in len_ins:
                                                        len_ins[num_ins] = 1
                                                    else:
                                                        len_ins[num_ins] = len_ins[num_ins] + 1
                                                    num_ins = 0
                                                    pre_ins_loc = -1
                                                num_mis = num_mis + 1
                                            else:
                                                if num_ins > 0:
                                                    if num_ins not in len_ins:
                                                        len_ins[num_ins] = 1
                                                    else:
                                                        len_ins[num_ins] = len_ins[num_ins] + 1
                                                    num_ins = 0
                                                    pre_ins_loc = -1
                                                    
            
                                                num_del = num_del + 1
                                                #if base_xun == begin_id:
                                                #    start_del_flag = 1
                                                #if base_xun == end_id:
                                                #    end_del_flag = 1
                                    elif updated_segment_cp[xun_base_present_abs] == 4:
                                        if xun_base_present_abs in tem_dic:
                                            if tem_dic[xun_base_present_abs] != 4:
                                                if tem_dic[xun_base_present_abs] < 4:
                                                    if xun_base_present_abs-pre_ins_loc<2:
                                                        num_ins = num_ins + 1
                                                    else:
                                                        if num_ins > 0:
                                                            if num_ins not in len_ins:
                                                                len_ins[num_ins] = 1
                                                            else:
                                                                len_ins[num_ins] = len_ins[num_ins] + 1
                                                            num_ins = 0
                                                        pre_ins_loc = xun_base_present_abs
                                                else:
                                                    print('error1')
                                                    print(basss_index)
                                                    #break
                                    else:
                                        if xun_base_present_abs in tem_dic:
                                            if tem_dic[xun_base_present_abs] == 4:
                                                print('error2')
                                            elif tem_dic[xun_base_present_abs] < 4:
                                                num_ins = num_ins + 1
                                    #base_xun = base_xun + 1
        
                                
                                if num_ins>0:
                                    if num_ins in len_ins:
                                        len_ins[num_ins] = len_ins[num_ins] + 1
                                    else:
                                        len_ins[num_ins] = 1
                                if num_del>0:
                                    if num_del in len_del:
                                        len_del[num_del] = len_del[num_del] + 1
                                    else:
                                        len_del[num_del] = 1
                                tem_probabity_correct = num_correct*np.log(pcorrect+corretion_eps_from_zero)
                                tem_probabity_mis = num_mis*np.log(p_mis_loc+corretion_eps_from_zero)
                                if len(len_ins) > 0:
                                    values_ = len_ins.values()
                                    num_ins_loc = sum(np.array(list(values_)))
                                    tem_probabity_ins = num_ins_loc*np.log(pins1_loc+corretion_eps_from_zero)
                                    for ii in len_ins:
                                        # np.log(stats.expon.pdf(ii,lamuna_ins))
                                        tem_probabity_ins = tem_probabity_ins + \
                                            len_ins[ii] * \
                                            possion_log_probability(ii-1,lamuna_ins)
                                        # print(tem_probabity_ins)
                                        # tem_probabity_ins = tem_probabity_ins + len_ins[ii]*np.log(stats.nbinom.pmf(ii,n_ins,p_ins))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                        # tem_probabity_ins = tem_probabity_ins + len_ins[ii]*np.log(1/(ii*sig_ins*np.sqrt(2*math.pi))*np.exp(-(np.log(ii)-mu_ins)**2/(2*sig_ins**2)))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                else:
                                    tem_probabity_ins = 0
                                    num_ins_loc = 0
        
                                del_whole = 0
                                if len(len_del) > 0:
                                    values_ = len_del.values()
                                    num_del_loc = sum(np.array(list(values_)))
                                    tem_probabity_del = num_del_loc*np.log(pdel1_loc+corretion_eps_from_zero)
                                    for dd in len_del:
                                        # np.log(stats.expon.pdf(ii,lamuna_ins))
                                        tem_probabity_del = tem_probabity_del + \
                                            len_del[dd] * \
                                             possion_log_probability(dd-1,lamuna_del)
                                        # tem_probabity_del = tem_probabity_del + len_del[dd]*np.log(stats.nbinom.pmf(dd,n_del,p_del))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                        # tem_probabity_del = tem_probabity_del + len_del[dd]*np.log(1/(dd*sig_del*np.sqrt(2*math.pi))*np.exp(-(np.log(dd)-mu_del)**2/(2*sig_del**2)))#np.log(stats.expon.pdf(ii,lamuna_ins))
                                        del_whole = del_whole + dd*len_del[dd]
        
                                else:
                                    tem_probabity_del = 0
                                    num_del_loc = 0
        
                                other_gaps = num_correct + num_mis - 1-num_del_loc+start_del_flag + end_del_flag
                                gap_not_ins_pro = other_gaps*np.log(1-pins1_loc+corretion_eps_from_zero)
                                not_del_num = num_correct + num_mis
                                gap_not_del_pro = not_del_num*np.log(1-pdel1_loc+corretion_eps_from_zero)
                                full_prob = gap_not_del_pro+tem_probabity_correct + tem_probabity_mis +\
                                    tem_probabity_ins+tem_probabity_del+gap_not_ins_pro
        
                                full_full_prob = full_full_prob + full_prob
                            
                            ### evo comparision
                            
    
                            ID_tem = [0, 1, 2, 3].index(main_base_xun)
                     
                            updated_state_mapping[basss_index-start_block_id] = 0
                            Tem_ins_evo = {}
                            
                            Tem_del_evo = {}
                            ins_len_evo = 0
                            del_len_evo = 0
                            num_mis_evo = 0
                            num_cor_evo = 0 
                            modify_reads_con_prob = 0
                            for ll_check in range(start_block_id, end_block_id+1):
                                relative_id =  ll_check - start_block_id
                                if original_true_state_ref[relative_id]<0:## meaning this should ins evo fromthe refence
                                    
                                    if del_len_evo>0:
                                        if del_len_evo in Tem_del_evo:
                                            Tem_del_evo[del_len_evo] = Tem_del_evo[del_len_evo] + 1
                                        else:
                                            Tem_del_evo[del_len_evo] = 1
                                        del_len_evo = 0
                                    ###start position
                                    if relative_id<1:
                                        if len(common_ins_per_df_list_df)>0:
                                            left_boundary_ins = np.array(common_ins_per_df_list_df.iloc[:,1])
                                            left_boundary_ins_dis = present_var_msa - left_boundary_ins
                                            if  1 in left_boundary_ins_dis:
                                                common_ins_len_tem_id =  list(left_boundary_ins_dis).index(1)
                                                common_ins_len_ = common_ins_per_df_list_df.iloc[common_ins_len_tem_id,1] - common_ins_per_df_list_df.iloc[common_ins_len_tem_id,0] + 1
                                                ins_len_evo = ins_len_evo + common_ins_len_
                                                modify_reads_con_prob =  modify_reads_con_prob + common_ins_len_*read_number_tem*(np.log(pcorrect)+np.log(1-p_del_loc)+np.log(1-pins1_loc))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                                           
                                        if len(common_del_per_df_list_df)>0:
                                            left_boundary_del = np.array(common_del_per_df_list_df.iloc[:,1])
                                            left_boundary_del_dis = present_var_msa - left_boundary_del
                                            if  1 in left_boundary_del_dis:
                                                common_del_len_tem_id =  list(left_boundary_del_dis).index(1)
                                                common_del_len_ = common_del_per_df_list_df.iloc[common_del_len_tem_id,1] - common_del_per_df_list_df.iloc[common_del_len_tem_id,0] + 1
                                                del_len_evo = del_len_evo + common_del_len_
                                                modify_reads_con_prob = modify_reads_con_prob + read_number_tem*(np.log(1-pins1_loc))
                                    if updated_state_mapping[relative_id]>0:
                                        ins_len_evo = ins_len_evo + 1
                                    
                                    if len(common_ins_per_df_list_df)>0:
                                        right_boundary_ins = np.array(common_ins_per_df_list_df.iloc[:,0])
                                        right_boundary_ins_dis = present_var_msa - right_boundary_ins
                                        if  1 in right_boundary_ins_dis:
                                            common_ins_len_tem_id =  list(right_boundary_ins_dis).index(1)
                                            common_ins_len_ = common_ins_per_df_list_df.iloc[common_ins_len_tem_id,1] - common_ins_per_df_list_df.iloc[common_ins_len_tem_id,0] + 1    
                                            ins_len_evo = ins_len_evo + common_ins_len_
                                            modify_reads_con_prob = modify_reads_con_prob + common_ins_len_*read_number_tem*(np.log(pcorrect)+np.log(1-p_del_loc)+np.log(1-pins1_loc))                            
                                            
                                    
                                else:
                                    if ins_len_evo>0:
                                        if ins_len_evo in Tem_ins_evo:
                                            Tem_ins_evo[ins_len_evo] = Tem_ins_evo[ins_len_evo] + 1
                                        else:
                                            Tem_ins_evo[ins_len_evo] =  1
                                        ins_len_evo = 0
                                    ###start position
                                    if relative_id<1:
                                        if len(common_ins_per_df_list_df)>0:
                                            left_boundary_ins = np.array(common_ins_per_df_list_df.iloc[:,1])
                                            left_boundary_ins_dis = present_var_msa - left_boundary_ins
                                            if  1 in left_boundary_ins_dis:
                                                common_ins_len_tem_id =  list(left_boundary_ins_dis).index(1)
                                                common_ins_len_ = common_ins_per_df_list_df.iloc[common_ins_len_tem_id,1] - common_ins_per_df_list_df.iloc[common_ins_len_tem_id,0] + 1
                                                ins_len_evo = ins_len_evo + common_ins_len_
                                                modify_reads_con_prob =  modify_reads_con_prob + common_ins_len_*read_number_tem*(np.log(pcorrect)+np.log(1-p_del_loc)+np.log(1-pins1_loc))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                                           
                                        if len(common_del_per_df_list_df)>0:
                                            left_boundary_del = np.array(common_del_per_df_list_df.iloc[:,1])
                                            left_boundary_del_dis = present_var_msa - left_boundary_del
                                            if  1 in left_boundary_del_dis:
                                                common_del_len_tem_id =  list(left_boundary_del_dis).index(1)
                                                common_del_len_ = common_del_per_df_list_df.iloc[common_del_len_tem_id,1] - common_del_per_df_list_df.iloc[common_del_len_tem_id,0] + 1
                                                del_len_evo = del_len_evo + common_del_len_
                                                modify_reads_con_prob = modify_reads_con_prob + read_number_tem*(np.log(1-pins1_loc))
                                    if updated_state_mapping[relative_id]>0:
                                        del_len_evo = del_len_evo + 1
                                    else:
                                        if ll_check == basss_index:
                                            if main_base_xun == standard_ref_code_block[relative_id]:
                                                num_cor_evo = num_cor_evo + 1
                                            else:
                                                num_mis_evo = num_mis_evo + 1
                                        else:
                                            if  updated_segment_cp[var_site_collection_abs[relative_id]] ==  standard_ref_code_block[relative_id]:
                                                num_cor_evo = num_cor_evo + 1
                                            else:
                                                num_mis_evo = num_mis_evo + 1
                                    ### del part
                                    if len(common_del_per_df_list_df)>0:
                                        right_boundary_del = np.array(common_del_per_df_list_df.iloc[:,0])
                                        right_boundary_del_dis = right_boundary_del - present_var_msa
                                        if  1 in right_boundary_del_dis:
                                            common_del_len_tem_id =  list(right_boundary_del_dis).index(1)
                                            common_del_len_ = common_del_per_df_list_df.iloc[common_del_len_tem_id,1] - common_del_per_df_list_df.iloc[common_del_len_tem_id,0] + 1
                                            del_len_evo = del_len_evo + common_del_len_
                            
                            if del_len_evo>0:
                                if del_len_evo in Tem_del_evo:
                                    Tem_del_evo[del_len_evo] = Tem_del_evo[del_len_evo] + 1
                                else:
                                    Tem_del_evo[del_len_evo] = 1
                            if ins_len_evo>0:
                                if ins_len_evo in Tem_ins_evo:
                                    Tem_ins_evo[ins_len_evo] = Tem_ins_evo[ins_len_evo] + 1
                                else:
                                    Tem_ins_evo[ins_len_evo] =  1
                            ### infer the prob
                            tem_evo = 0
                            
                            for del_per in Tem_del_evo:
                                tem_evo = tem_evo + Tem_del_evo[del_per]*(expotential_log_log_probability(del_per-1,del_p[kk],del_beta[kk]) + np.log(gamma_del+corretion_eps_from_zero) )
                            for ins_per in Tem_ins_evo:
                                tem_evo = tem_evo + Tem_ins_evo[ins_per]*(expotential_log_log_probability(ins_per-1,ins_p[kk],ins_beta[kk]) + np.log(gamma_ins+corretion_eps_from_zero) )
                            tem_evo = tem_evo + (num_cor_evo+ num_mis_evo)*np.log((1-gamma_del)) + (num_cor_evo)*np.log(gamma_c) + (num_mis_evo)*np.log(gamma_mis+corretion_eps_from_zero)     
                            
                            potential_probability[int(ID_tem)] = full_full_prob + tem_evo + modify_reads_con_prob
                        
                        else:
                            ID_tem = [0, 1, 2, 3].index(main_base_xun)
                            potential_probability[int(ID_tem)] = float('-inf')
                    #tem_probability = copy.deepcopy(potential_probability/tem_pri)
                    #tem_probability = tem_probability - max(tem_probability)
                    #tem_probability_main = np.exp(tem_probability)
                    #tem_probability_main = tem_probability_main/sum(tem_probability_main)
                    #post_sampling = np.random.multinomial(n= 1,pvals = tem_probability_main)
                    #max_ID = list(post_sampling).index(1)
                    #    n=1, pvals=[6/8, 1/8, 1/16, 1/32, 1/32])
                    #l_i = list(l_i_vec).index(1)+1
                    if max(potential_probability)>float('-inf'):
                        max_ID = list(potential_probability).index(max(potential_probability))
                        updated_state_mapping[basss_index - start_block_id] = 0
                        updated_segment[present_var_msa] = max_ID
                    else:
                        updated_state_mapping[basss_index - start_block_id] = 0
                        relative_id = basss_index- start_block_id
                        updated_segment[present_var_msa] = standard_ref_code_block[relative_id]
                    
    
            return_result.append([infer_var_site[basss_index], updated_segment[present_var_msa]])

    
        return(return_result)

def possion_log_probability(l_max,lamuna):
    #K_max = 500
    #lamuna = 0.5
    full_sum = -lamuna + l_max*np.log(lamuna)
    if l_max>0:
        ll_vec = np.array(range(1,int(l_max)+1))
        full_sum = full_sum - sum(np.log(ll_vec))
    return(full_sum)


def expotential_log_log_probability(l_max,p_,beta_):
    #l_max = 0
    #p_ = 1
    #beta_ = 8
    corretion_eps_from_zero = 0.000000000001
    num_return = -np.log(-np.log(min(p_+corretion_eps_from_zero,1-corretion_eps_from_zero))+corretion_eps_from_zero)+(np.log(beta_+corretion_eps_from_zero)+np.log(1-p_+corretion_eps_from_zero)-beta_*l_max)-\
        np.log(1-(1-max(p_,corretion_eps_from_zero))*np.exp(-beta_*l_max)+corretion_eps_from_zero)
    #shiyan = 0
    #shiyan = (1/(-np.log(p_))*beta_*(1-p_)*np.exp(-beta_*l_max)/(1-(1-p_)*np.exp(-beta_*l_max)))
    return(num_return)
    

def convert_sparse_normal(X_sparse, n):
    X_sparse_copy = X_sparse  # X_sparse is the data frmae
    original_matrix = np.zeros([n, n])
    for ll in range(len(X_sparse_copy)):
        if X_sparse_copy.iloc[ll, 2] > 0:
            row_index = int(X_sparse_copy.iloc[ll, 0])
            col_index = int(X_sparse_copy.iloc[ll, 1])
            original_matrix[row_index, col_index] = X_sparse_copy.iloc[ll, 2]
            original_matrix[col_index, row_index] = X_sparse_copy.iloc[ll, 2]
    return(original_matrix)

def core_selection(data_count_vec,main_number_id,ll):

    

    #data_count_vec = full_matrix[ll,:]

    return_location = -1

    if main_number_id<4:## we only casr the refererence with know nucleoside

        main_number = data_count_vec[main_number_id]

        whole_number = sum(data_count_vec)

        correct_ratio = 0.9

        if whole_number>5:

            other_number = whole_number-main_number

            expt_1 = int(whole_number*correct_ratio)

            expt_2 = whole_number - expt_1

            if whole_number>=50 and other_number>=5:

                #print(ll)

                fobs = [main_number,other_number]

                fexp = [expt_1,expt_2]

                s = scipy.stats.chisquare(fobs,fexp)

                if s[1]<0.05:

                   if main_number<expt_1:

                       return_location = ll

                #if s[1]<0.20:

                #   if main_number<expt_1:

                #       coolection_snp11.append(ll)

            else:

                fobs = [expt_2,other_number]

                f_other = [expt_1,main_number]

                #fobs = []

                s = scipy.stats.barnard_exact([fobs,f_other],alternative='less')

                if s.pvalue<0.05:

                   if main_number<expt_1:

                       return_location = ll

    return(return_location)



def collection_var_per_read(read_xun_id):
    global Best_ref,Best_ref_ins,Best_ref_ins_keys,tem_res_collection,keys_ins_common_original_array,\
        only_bone_note_df_keys,ref_code_standard,common_del_block_sites_collection,\
            pcorrect,p_mis_loc,p_del_loc,pins1_loc,lamuna_ins,lamuna_del,tao,portions
    #tem_sum = 0
    #collection = []
    #for ll in Best_ref[0]:
    #    if Best_ref[0][ll] == Best_ref[1][ll]:
    #        tem_sum = tem_sum + 1
    #    else:
    #        collection.append(ll)
    #read_xun_id = 0
#for read_xun_id in range(max_num_reads_refine):
    #print(read_xun_id)
    corretion_eps_from_zero = 0.00000000001
    #read_xun_id = 0
    tem_dic_bone = tem_res_collection[int(read_xun_id)][6]
    begin_id = tem_res_collection[int(read_xun_id)][1]## this is the original mapping start
    end_id = tem_res_collection[int(read_xun_id)][2]## this is the original mapping end
    tem_ins_dic_per_read = tem_res_collection[read_xun_id][5]
    ### common independent bone, we do not have worriies about the del and ins/same for any consensus
    #sub_common_bone = Independet_site_col_no_alg_array[np.where((Independet_site_col_no_alg_array>=begin_id)&(Independet_site_col_no_alg_array<=end_id))]
    keys_ins_common_original_array_sub = keys_ins_common_original_array[np.where((keys_ins_common_original_array>=begin_id)&(keys_ins_common_original_array<=end_id))]
    other_informative_ins_id_per = list(set(tem_ins_dic_per_read.keys()) - set(keys_ins_common_original_array_sub))
    #common_del_sub = common_del_block_sites_collection_array[np.where((common_del_block_sites_collection_array>=begin_id)&(common_del_block_sites_collection_array<=end_id))]

    only_bone_note_df_keys_array = np.array(only_bone_note_df_keys)
    only_bone_note_df_keys_array_sub = sorted(list(only_bone_note_df_keys_array[np.where((only_bone_note_df_keys_array>=begin_id)&(only_bone_note_df_keys_array<=end_id))]))
    remaining_bone_no_del = list(set(range(begin_id,end_id+1))-set(only_bone_note_df_keys_array_sub)-set(common_del_block_sites_collection))
    tem_correct_per_o = len(remaining_bone_no_del)
    
    
    
    ### ins bone collection
    infortive_ins_collection = [[] for x in range(K)]### we need find the ins part
    for kk in range(K):
        Best_ref_ins_keys_sub = Best_ref_ins_keys[kk][np.where((Best_ref_ins_keys[kk]>=begin_id)&(Best_ref_ins_keys[kk]<end_id))]
        infortive_ins_collection[kk] = sorted(list(set(list(Best_ref_ins_keys_sub) + list(tem_ins_dic_per_read.keys())))+only_bone_note_df_keys_array_sub)
    
    
    
    tem_prob_final = np.zeros(K)
    tem_ins_return_dic = []
    tem_del_return_dic = []
    tem_match_resutn_list = []
    tem_mismatch_resutn_list = []
    for kk in range(K):
        #kk = 1
        tem_correct_per = tem_correct_per_o
        tem_mis_per = 0
        tem_ins_content_con = ''
        tem_ins_content_read = ''
        #tem_ins_len = 0
        #tem_del_len = 0
        len_del_group = {}
        len_ins_group = {}
        previous_del_bone = -1
        previous_del_flag = 0
        ll_per_id = 0
        pre_bone_ = -1
        ### do while format
        if len(infortive_ins_collection[kk])>0:
            tem_new_bone_p_o = infortive_ins_collection[kk][ll_per_id]
            if tem_new_bone_p_o in only_bone_note_df_keys_array_sub:
                tem_infer_base_con = Best_ref[kk][tem_new_bone_p_o]
                tem_infer_base_read = tem_dic_bone[tem_new_bone_p_o]
                if tem_infer_base_con <4:### we have the true base
                    if tem_infer_base_read in Standardbase[:4]:
                        if len(tem_ins_content_con)>0 and len(tem_ins_content_read)>0: # con and read in the bone
                            sub_string_con_string = tem_ins_content_con
                            sub_string_read_string = tem_ins_content_read
                            recollection = []
                            seq_collection = [sub_string_con_string, sub_string_read_string]
                                                
                            for seq_index in range(len(seq_collection)):
                                tem_seq = seqbio.NucleotideSequence(
                                    seq_collection[seq_index])
                                recollection.append(tem_seq)
                            try:
                                alignment, order, guide_tree, distance_matrix = align.align_multiple(
                                    recollection,
                                    matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                    gap_penalty=(-10, -0.5),
                                    terminal_penalty=False
                                )
                                alignment_trace = alignment.trace
                                length_compare, num_col = np.shape(
                                    alignment_trace)
                                result_mapping = []
                                for col_index in range(num_col):
                                    #tem_result_mapping = [recollection[col_index]]
                                    coverted_seq = converting_code_2_seq(
                                        alignment_trace[:, col_index], recollection[col_index])
                                    # tem_result_mapping.append(coverted_seq)
                                    result_mapping.append(
                                        [seq_collection[col_index], coverted_seq.lower()])
                            except (ValueError, ZeroDivisionError):
                                max_len = len(seq_collection[0])
                                for col_index in range(len(seq_collection)):
                                    if len(seq_collection[col_index]) > max_len:
                                        max_len = len(
                                            seq_collection[col_index])
                                result_mapping = []
                                for col_index in range(len(seq_collection)):
                                    if len(seq_collection[col_index]) < max_len:
                                        refined_ins_ = seq_collection[col_index]+'N'*(
                                            max_len-len(seq_collection[col_index]))
            
                                    else:
                                        refined_ins_ = seq_collection[col_index]
                                    result_mapping.append(
                                        [seq_collection[col_index], refined_ins_])
                            
                            tem_tem_ins_content_con = ''
                            tem_tem_ins_content_read = ''
                            #Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5] = Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5]+len(result_mapping[0][1])
                            for ll_fnei in range(len(result_mapping[0][1])):
                                if result_mapping[0][1][ll_fnei] in Standardbase[:4] and result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                    #if result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                    ## ins part for con but del for read
                                    if len(tem_tem_ins_content_con)>0:
                                        
                                        if len(tem_tem_ins_content_con) in len_del_group:
                                            len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                        else:
                                            len_del_group[len(tem_tem_ins_content_con)] = 1
                                        tem_tem_ins_content_con = ''
                                    ## del part
                                    if len(tem_tem_ins_content_read)>0:
                                        if len(tem_tem_ins_content_read) in len_ins_group:
                                            len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                        else:
                                            len_ins_group[len(tem_tem_ins_content_read)] = 1
                                        tem_tem_ins_content_read = ''
                                    ## bone judge
                                     
                                    if result_mapping[0][1][ll_fnei] == result_mapping[1][1][ll_fnei]:
                                        tem_correct_per = tem_correct_per + 1
                                    else:
                                        tem_mis_per = tem_mis_per + 1
                                        
                                elif result_mapping[0][1][ll_fnei] in Standardbase[:4]:
                                    tem_tem_ins_content_con = tem_tem_ins_content_con + result_mapping[0][1][ll_fnei]
                                    
                                elif result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                    tem_tem_ins_content_read = tem_tem_ins_content_read + result_mapping[1][1][ll_fnei]
                            
                            ## bone judge            
                            ## ins part for con but del for read
                            if len(tem_tem_ins_content_con)>0:
                                
                                if len(tem_tem_ins_content_con) in len_del_group:
                                    len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                else:
                                    len_del_group[len(tem_tem_ins_content_con)] = 1
                                tem_tem_ins_content_con = ''
                            ## del part
                            if len(tem_tem_ins_content_read)>0:
                                if len(tem_tem_ins_content_read) in len_ins_group:
                                    len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                else:
                                    len_ins_group[len(tem_tem_ins_content_read)] = 1
                                tem_tem_ins_content_read = ''
                        
                        elif len(tem_ins_content_read)>0:
                            if len(tem_ins_content_read) in len_ins_group:
                                len_ins_group[len(tem_ins_content_read)] = len_ins_group[len(tem_ins_content_read)] + 1
                            else:
                                len_ins_group[len(tem_ins_content_read)] = 1
                            tem_ins_content_read = ''
                            
                        elif len(tem_ins_content_con)>0:
                            
                            if len(tem_ins_content_con) in len_del_group:
                                len_del_group[len(tem_ins_content_con)] = len_del_group[len(tem_ins_content_con)] + 1
                            else:
                                len_del_group[len(tem_ins_content_con)] = 1
                            
                            tem_ins_content_con = ''
                           
                        
                        if Standardbase[tem_infer_base_con] == tem_infer_base_read:
                            tem_correct_per = tem_correct_per + 1
                        else:
                            tem_mis_per = tem_mis_per + 1
                        previous_del_flag = 0
                        previous_del_bone = -1
                        tem_ins_content_read = ''
                        tem_ins_content_con = ''
                    else:### this should be the deletion/ common deletion is just like the tunnel to link two deletion neighboured with common deletion
                        if previous_del_flag: # neighboured del do not have inserton 
                            if tem_new_bone_p_o - previous_del_bone<2:
                                #tem_del_len = tem_del_len + 1
                                tem_ins_content_con = tem_ins_content_con + Standardbase[tem_infer_base_con]
                                #tem_del_len = 1 ## check the tunnel effect of the common deletion
                                if len(common_del_block_df_cp)>0:
                                    sub_df_common_del_left = common_del_block_df_cp[common_del_block_df_cp.iloc[:,3]<tem_new_bone_p_o]
                                    if len(sub_df_common_del_left)>0:
                                        last_line_bind_del = sub_df_common_del_left.iloc[-1,3]
                                        if tem_new_bone_p_o - last_line_bind_del<2:##left join
                                            #tem_del_length = tem_del_length + sub_left_df_del.iloc[-1,3] - sub_left_df_del.iloc[-1,2] + 1
                                            #common_del_block_df_cp.iloc[len(sub_left_df_del)-1,4] = 1
                                            #index_del_o = len(sub_left_df_del)-1
                                            #covered_common_del_id.append(index_del_o)
                                            previous_del_bone = common_del_block_df_cp.iloc[len(sub_df_common_del_left)-1,3]
                                        else:
                                            previous_del_bone = tem_new_bone_p_o
                                    else:
                                        previous_del_bone = tem_new_bone_p_o
                                else:
                                    previous_del_bone = tem_new_bone_p_o
                                previous_del_flag = 1
                                tem_ins_content_read = ''
                                #tem_ins_content_con = Standardbase[tem_infer_base_con]
                                
                            else:
                                ### it may occur not neighboured middle insertion
                                if len(tem_ins_content_con)>0 and len(tem_ins_content_read)>0: # con and read in the bone
                                    sub_string_con_string = tem_ins_content_con
                                    sub_string_read_string = tem_ins_content_read
                                    recollection = []
                                    seq_collection = [sub_string_con_string, sub_string_read_string]
                                                        
                                    for seq_index in range(len(seq_collection)):
                                        tem_seq = seqbio.NucleotideSequence(
                                            seq_collection[seq_index])
                                        recollection.append(tem_seq)
                                    try:
                                        alignment, order, guide_tree, distance_matrix = align.align_multiple(
                                            recollection,
                                            matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                            gap_penalty=(-10, -0.5),
                                            terminal_penalty=False
                                        )
                                        alignment_trace = alignment.trace
                                        length_compare, num_col = np.shape(
                                            alignment_trace)
                                        result_mapping = []
                                        for col_index in range(num_col):
                                            #tem_result_mapping = [recollection[col_index]]
                                            coverted_seq = converting_code_2_seq(
                                                alignment_trace[:, col_index], recollection[col_index])
                                            # tem_result_mapping.append(coverted_seq)
                                            result_mapping.append(
                                                [seq_collection[col_index], coverted_seq.lower()])
                                    except (ValueError, ZeroDivisionError):
                                        max_len = len(seq_collection[0])
                                        for col_index in range(len(seq_collection)):
                                            if len(seq_collection[col_index]) > max_len:
                                                max_len = len(
                                                    seq_collection[col_index])
                                        result_mapping = []
                                        for col_index in range(len(seq_collection)):
                                            if len(seq_collection[col_index]) < max_len:
                                                refined_ins_ = seq_collection[col_index]+'N'*(
                                                    max_len-len(seq_collection[col_index]))
                    
                                            else:
                                                refined_ins_ = seq_collection[col_index]
                                            result_mapping.append(
                                                [seq_collection[col_index], refined_ins_])
                                    
                                    tem_tem_ins_content_con = ''
                                    tem_tem_ins_content_read = ''
                                    #Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5] = Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5]+len(result_mapping[0][1])
                                    for ll_fnei in range(len(result_mapping[0][1])):
                                        if result_mapping[0][1][ll_fnei] in Standardbase[:4] and result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                            #if result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                            ## ins part for con but del for read
                                            if len(tem_tem_ins_content_con)>0:
                                                
                                                if len(tem_tem_ins_content_con) in len_del_group:
                                                    len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                                else:
                                                    len_del_group[len(tem_tem_ins_content_con)] = 1
                                                tem_tem_ins_content_con = ''
                                            ## del part
                                            if len(tem_tem_ins_content_read)>0:
                                                if len(tem_tem_ins_content_read) in len_ins_group:
                                                    len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                                else:
                                                    len_ins_group[len(tem_tem_ins_content_read)] = 1
                                                tem_tem_ins_content_read = ''
                                            ## bone judge
                                             
                                            if result_mapping[0][1][ll_fnei] == result_mapping[1][1][ll_fnei]:
                                                tem_correct_per = tem_correct_per + 1
                                            else:
                                                tem_mis_per = tem_mis_per + 1
                                                
                                        elif result_mapping[0][1][ll_fnei] in Standardbase[:4]:
                                            tem_tem_ins_content_con = tem_tem_ins_content_con + result_mapping[0][1][ll_fnei]
                                            
                                        elif result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                            tem_tem_ins_content_read = tem_tem_ins_content_read + result_mapping[1][1][ll_fnei]
                                    
                                    ## bone judge            
                                    ## ins part for con but del for read
                                    if len(tem_tem_ins_content_con)>0:
                                        
                                        if len(tem_tem_ins_content_con) in len_del_group:
                                            len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                        else:
                                            len_del_group[len(tem_tem_ins_content_con)] = 1
                                        tem_tem_ins_content_con = ''
                                    ## del part
                                    if len(tem_tem_ins_content_read)>0:
                                        if len(tem_tem_ins_content_read) in len_ins_group:
                                            len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                        else:
                                            len_ins_group[len(tem_tem_ins_content_read)] = 1
                                        tem_tem_ins_content_read = ''
                                
                                elif len(tem_ins_content_read)>0:
                                    if len(tem_ins_content_read) in len_ins_group:
                                        len_ins_group[len(tem_ins_content_read)] = len_ins_group[len(tem_ins_content_read)] + 1
                                    else:
                                        len_ins_group[len(tem_ins_content_read)] = 1
                                    tem_ins_content_read = ''
                                    
                                elif len(tem_ins_content_con)>0:
                                    
                                    if len(tem_ins_content_con) in len_del_group:
                                        len_del_group[len(tem_ins_content_con)] = len_del_group[len(tem_ins_content_con)] + 1
                                    else:
                                        len_del_group[len(tem_ins_content_con)] = 1
                                    
                                    tem_ins_content_con = ''
                                #tem_del_len = 1 ## check the tunnel effect of the common deletion
                                if len(common_del_block_df_cp)>0:
                                    sub_df_common_del_left = common_del_block_df_cp[common_del_block_df_cp.iloc[:,3]<tem_new_bone_p_o]
                                    if len(sub_df_common_del_left)>0:
                                        last_line_bind_del = sub_df_common_del_left.iloc[-1,3]
                                        if tem_new_bone_p_o - last_line_bind_del<2:##left join
                                            #tem_del_length = tem_del_length + sub_left_df_del.iloc[-1,3] - sub_left_df_del.iloc[-1,2] + 1
                                            #common_del_block_df_cp.iloc[len(sub_left_df_del)-1,4] = 1
                                            #index_del_o = len(sub_left_df_del)-1
                                            #covered_common_del_id.append(index_del_o)
                                            previous_del_bone = common_del_block_df_cp.iloc[len(sub_df_common_del_left)-1,3]
                                        else:
                                            previous_del_bone = tem_new_bone_p_o
                                    else:
                                        previous_del_bone = tem_new_bone_p_o
                                else:
                                    previous_del_bone = tem_new_bone_p_o
                                previous_del_flag = 1
                                tem_ins_content_read = ''
                                tem_ins_content_con = Standardbase[tem_infer_base_con]
                                
                        else:
                            ### it may occur not neighboured middle insertion
                            if len(tem_ins_content_con)>0 and len(tem_ins_content_read)>0: # con and read in the bone
                                sub_string_con_string = tem_ins_content_con
                                sub_string_read_string = tem_ins_content_read
                                recollection = []
                                seq_collection = [sub_string_con_string, sub_string_read_string]
                                                    
                                for seq_index in range(len(seq_collection)):
                                    tem_seq = seqbio.NucleotideSequence(
                                        seq_collection[seq_index])
                                    recollection.append(tem_seq)
                                try:
                                    alignment, order, guide_tree, distance_matrix = align.align_multiple(
                                        recollection,
                                        matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                        gap_penalty=(-10, -0.5),
                                        terminal_penalty=False
                                    )
                                    alignment_trace = alignment.trace
                                    length_compare, num_col = np.shape(
                                        alignment_trace)
                                    result_mapping = []
                                    for col_index in range(num_col):
                                        #tem_result_mapping = [recollection[col_index]]
                                        coverted_seq = converting_code_2_seq(
                                            alignment_trace[:, col_index], recollection[col_index])
                                        # tem_result_mapping.append(coverted_seq)
                                        result_mapping.append(
                                            [seq_collection[col_index], coverted_seq.lower()])
                                except (ValueError, ZeroDivisionError):
                                    max_len = len(seq_collection[0])
                                    for col_index in range(len(seq_collection)):
                                        if len(seq_collection[col_index]) > max_len:
                                            max_len = len(
                                                seq_collection[col_index])
                                    result_mapping = []
                                    for col_index in range(len(seq_collection)):
                                        if len(seq_collection[col_index]) < max_len:
                                            refined_ins_ = seq_collection[col_index]+'N'*(
                                                max_len-len(seq_collection[col_index]))
                
                                        else:
                                            refined_ins_ = seq_collection[col_index]
                                        result_mapping.append(
                                            [seq_collection[col_index], refined_ins_])
                                
                                tem_tem_ins_content_con = ''
                                tem_tem_ins_content_read = ''
                                #Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5] = Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5]+len(result_mapping[0][1])
                                for ll_fnei in range(len(result_mapping[0][1])):
                                    if result_mapping[0][1][ll_fnei] in Standardbase[:4] and result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                        #if result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                        ## ins part for con but del for read
                                        if len(tem_tem_ins_content_con)>0:
                                            
                                            if len(tem_tem_ins_content_con) in len_del_group:
                                                len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                            else:
                                                len_del_group[len(tem_tem_ins_content_con)] = 1
                                            tem_tem_ins_content_con = ''
                                        ## del part
                                        if len(tem_tem_ins_content_read)>0:
                                            if len(tem_tem_ins_content_read) in len_ins_group:
                                                len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                            else:
                                                len_ins_group[len(tem_tem_ins_content_read)] = 1
                                            tem_tem_ins_content_read = ''
                                        ## bone judge
                                         
                                        if result_mapping[0][1][ll_fnei] == result_mapping[1][1][ll_fnei]:
                                            tem_correct_per = tem_correct_per + 1
                                        else:
                                            tem_mis_per = tem_mis_per + 1
                                            
                                    elif result_mapping[0][1][ll_fnei] in Standardbase[:4]:
                                        tem_tem_ins_content_con = tem_tem_ins_content_con + result_mapping[0][1][ll_fnei]
                                        
                                    elif result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                        tem_tem_ins_content_read = tem_tem_ins_content_read + result_mapping[1][1][ll_fnei]
                                
                                ## bone judge            
                                ## ins part for con but del for read
                                if len(tem_tem_ins_content_con)>0:
                                    
                                    if len(tem_tem_ins_content_con) in len_del_group:
                                        len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                    else:
                                        len_del_group[len(tem_tem_ins_content_con)] = 1
                                    tem_tem_ins_content_con = ''
                                ## del part
                                if len(tem_tem_ins_content_read)>0:
                                    if len(tem_tem_ins_content_read) in len_ins_group:
                                        len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                    else:
                                        len_ins_group[len(tem_tem_ins_content_read)] = 1
                                    tem_tem_ins_content_read = ''
                            
                            elif len(tem_ins_content_read)>0:
                                if len(tem_ins_content_read) in len_ins_group:
                                    len_ins_group[len(tem_ins_content_read)] = len_ins_group[len(tem_ins_content_read)] + 1
                                else:
                                    len_ins_group[len(tem_ins_content_read)] = 1
                                tem_ins_content_read = ''
                                
                            elif len(tem_ins_content_con)>0:
                                
                                if len(tem_ins_content_con) in len_del_group:
                                    len_del_group[len(tem_ins_content_con)] = len_del_group[len(tem_ins_content_con)] + 1
                                else:
                                    len_del_group[len(tem_ins_content_con)] = 1
                                
                                tem_ins_content_con = ''
                               
                            
                            #tem_del_len = 1 ## check the tunnel effect of the common deletion
                            if len(common_del_block_df_cp)>0:
                                sub_df_common_del_left = common_del_block_df_cp[common_del_block_df_cp.iloc[:,3]<tem_new_bone_p_o]
                                if len(sub_df_common_del_left)>0:
                
                                    last_line_bind_del = sub_df_common_del_left.iloc[-1,3]
                                    if tem_new_bone_p_o - last_line_bind_del<2:##left join
                                        #tem_del_length = tem_del_length + sub_left_df_del.iloc[-1,3] - sub_left_df_del.iloc[-1,2] + 1
                                        #common_del_block_df_cp.iloc[len(sub_left_df_del)-1,4] = 1
                                        #index_del_o = len(sub_left_df_del)-1
                                        #covered_common_del_id.append(index_del_o)
                                        previous_del_bone = common_del_block_df_cp.iloc[len(sub_df_common_del_left)-1,3]
                                        
                                    else:
                                        previous_del_bone = tem_new_bone_p_o
                                else:
                                    previous_del_bone = tem_new_bone_p_o
                            else:
                                previous_del_bone = tem_new_bone_p_o
                            previous_del_flag = 1
                            tem_ins_content_read = ''
                            tem_ins_content_con = Standardbase[tem_infer_base_con]
                else:#not occur in the con but read
                    if tem_infer_base_read in Standardbase[:4]:
                        if tem_new_bone_p_o - pre_bone_<2:
                            tem_ins_content_read = tem_ins_content_read + tem_infer_base_read
                        else:
                            ### it may occur not neighboured middle insertion
                            if len(tem_ins_content_con)>0 and len(tem_ins_content_read)>0: # con and read in the bone
                                sub_string_con_string = tem_ins_content_con
                                sub_string_read_string = tem_ins_content_read
                                recollection = []
                                seq_collection = [sub_string_con_string, sub_string_read_string]
                                                    
                                for seq_index in range(len(seq_collection)):
                                    tem_seq = seqbio.NucleotideSequence(
                                        seq_collection[seq_index])
                                    recollection.append(tem_seq)
                                try:
                                    alignment, order, guide_tree, distance_matrix = align.align_multiple(
                                        recollection,
                                        matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                        gap_penalty=(-10, -0.5),
                                        terminal_penalty=False
                                    )
                                    alignment_trace = alignment.trace
                                    length_compare, num_col = np.shape(
                                        alignment_trace)
                                    result_mapping = []
                                    for col_index in range(num_col):
                                        #tem_result_mapping = [recollection[col_index]]
                                        coverted_seq = converting_code_2_seq(
                                            alignment_trace[:, col_index], recollection[col_index])
                                        # tem_result_mapping.append(coverted_seq)
                                        result_mapping.append(
                                            [seq_collection[col_index], coverted_seq.lower()])
                                except (ValueError, ZeroDivisionError):
                                    max_len = len(seq_collection[0])
                                    for col_index in range(len(seq_collection)):
                                        if len(seq_collection[col_index]) > max_len:
                                            max_len = len(
                                                seq_collection[col_index])
                                    result_mapping = []
                                    for col_index in range(len(seq_collection)):
                                        if len(seq_collection[col_index]) < max_len:
                                            refined_ins_ = seq_collection[col_index]+'N'*(
                                                max_len-len(seq_collection[col_index]))
                
                                        else:
                                            refined_ins_ = seq_collection[col_index]
                                        result_mapping.append(
                                            [seq_collection[col_index], refined_ins_])
                                
                                tem_tem_ins_content_con = ''
                                tem_tem_ins_content_read = ''
                                #Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5] = Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5]+len(result_mapping[0][1])
                                for ll_fnei in range(len(result_mapping[0][1])):
                                    if result_mapping[0][1][ll_fnei] in Standardbase[:4] and result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                        #if result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                        ## ins part for con but del for read
                                        if len(tem_tem_ins_content_con)>0:
                                            
                                            if len(tem_tem_ins_content_con) in len_del_group:
                                                len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                            else:
                                                len_del_group[len(tem_tem_ins_content_con)] = 1
                                            tem_tem_ins_content_con = ''
                                        ## del part
                                        if len(tem_tem_ins_content_read)>0:
                                            if len(tem_tem_ins_content_read) in len_ins_group:
                                                len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                            else:
                                                len_ins_group[len(tem_tem_ins_content_read)] = 1
                                            tem_tem_ins_content_read = ''
                                        ## bone judge
                                         
                                        if result_mapping[0][1][ll_fnei] == result_mapping[1][1][ll_fnei]:
                                            tem_correct_per = tem_correct_per + 1
                                        else:
                                            tem_mis_per = tem_mis_per + 1
                                            
                                    elif result_mapping[0][1][ll_fnei] in Standardbase[:4]:
                                        tem_tem_ins_content_con = tem_tem_ins_content_con + result_mapping[0][1][ll_fnei]
                                        
                                    elif result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                        tem_tem_ins_content_read = tem_tem_ins_content_read + result_mapping[1][1][ll_fnei]
                                
                                ## bone judge            
                                ## ins part for con but del for read
                                if len(tem_tem_ins_content_con)>0:
                                    
                                    if len(tem_tem_ins_content_con) in len_del_group:
                                        len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                    else:
                                        len_del_group[len(tem_tem_ins_content_con)] = 1
                                    tem_tem_ins_content_con = ''
                                ## del part
                                if len(tem_tem_ins_content_read)>0:
                                    if len(tem_tem_ins_content_read) in len_ins_group:
                                        len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                    else:
                                        len_ins_group[len(tem_tem_ins_content_read)] = 1
                                    tem_tem_ins_content_read = ''
                            
                            elif len(tem_ins_content_read)>0:
                                if len(tem_ins_content_read) in len_ins_group:
                                    len_ins_group[len(tem_ins_content_read)] = len_ins_group[len(tem_ins_content_read)] + 1
                                else:
                                    len_ins_group[len(tem_ins_content_read)] = 1
                                tem_ins_content_read = ''
                                
                            elif len(tem_ins_content_con)>0:
                                
                                if len(tem_ins_content_con) in len_del_group:
                                    len_del_group[len(tem_ins_content_con)] = len_del_group[len(tem_ins_content_con)] + 1
                                else:
                                    len_del_group[len(tem_ins_content_con)] = 1
                                
                                tem_ins_content_con = ''
                               
                            
                            #previous_del_flag = 0
                            tem_ins_content_read = tem_infer_base_read
                            #tem_ins_content_con = ''
                    
                        previous_del_bone = -1
                        previous_del_flag = 0
                        tem_ins_content_con = ''
                    #else is the two del d; do nothin
                    
                if tem_new_bone_p_o in tem_ins_dic_per_read:
                    tem_ins_content_read = tem_ins_content_read + tem_ins_dic_per_read[tem_new_bone_p_o]
                if tem_new_bone_p_o in Best_ref_ins[kk]:
                    tem_ins_content_con = tem_ins_content_con + Best_ref_ins[kk][tem_new_bone_p_o]
            else:
                previous_del_flag = 0
                previous_del_bone = -1
                if tem_new_bone_p_o - pre_bone_>=2:
                    ### ins seperated
                    if len(tem_ins_content_con)>0 and len(tem_ins_content_read)>0: # con and read in the bone
                        sub_string_con_string = tem_ins_content_con
                        sub_string_read_string = tem_ins_content_read
                        recollection = []
                        seq_collection = [sub_string_con_string, sub_string_read_string]
                                            
                        for seq_index in range(len(seq_collection)):
                            tem_seq = seqbio.NucleotideSequence(
                                seq_collection[seq_index])
                            recollection.append(tem_seq)
                        try:
                            alignment, order, guide_tree, distance_matrix = align.align_multiple(
                                recollection,
                                matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                gap_penalty=(-10, -0.5),
                                terminal_penalty=False
                            )
                            alignment_trace = alignment.trace
                            length_compare, num_col = np.shape(
                                alignment_trace)
                            result_mapping = []
                            for col_index in range(num_col):
                                #tem_result_mapping = [recollection[col_index]]
                                coverted_seq = converting_code_2_seq(
                                    alignment_trace[:, col_index], recollection[col_index])
                                # tem_result_mapping.append(coverted_seq)
                                result_mapping.append(
                                    [seq_collection[col_index], coverted_seq.lower()])
                        except (ValueError, ZeroDivisionError):
                            max_len = len(seq_collection[0])
                            for col_index in range(len(seq_collection)):
                                if len(seq_collection[col_index]) > max_len:
                                    max_len = len(
                                        seq_collection[col_index])
                            result_mapping = []
                            for col_index in range(len(seq_collection)):
                                if len(seq_collection[col_index]) < max_len:
                                    refined_ins_ = seq_collection[col_index]+'N'*(
                                        max_len-len(seq_collection[col_index]))
        
                                else:
                                    refined_ins_ = seq_collection[col_index]
                                result_mapping.append(
                                    [seq_collection[col_index], refined_ins_])
                        
                        tem_tem_ins_content_con = ''
                        tem_tem_ins_content_read = ''
                        #Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5] = Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5]+len(result_mapping[0][1])
                        for ll_fnei in range(len(result_mapping[0][1])):
                            if result_mapping[0][1][ll_fnei] in Standardbase[:4] and result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                #if result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                ## ins part for con but del for read
                                if len(tem_tem_ins_content_con)>0:
                                    
                                    if len(tem_tem_ins_content_con) in len_del_group:
                                        len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                    else:
                                        len_del_group[len(tem_tem_ins_content_con)] = 1
                                    tem_tem_ins_content_con = ''
                                ## del part
                                if len(tem_tem_ins_content_read)>0:
                                    if len(tem_tem_ins_content_read) in len_ins_group:
                                        len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                    else:
                                        len_ins_group[len(tem_tem_ins_content_read)] = 1
                                    tem_tem_ins_content_read = ''
                                ## bone judge
                                 
                                if result_mapping[0][1][ll_fnei] == result_mapping[1][1][ll_fnei]:
                                    tem_correct_per = tem_correct_per + 1
                                else:
                                    tem_mis_per = tem_mis_per + 1
                                    
                            elif result_mapping[0][1][ll_fnei] in Standardbase[:4]:
                                tem_tem_ins_content_con = tem_tem_ins_content_con + result_mapping[0][1][ll_fnei]
                                
                            elif result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                tem_tem_ins_content_read = tem_tem_ins_content_read + result_mapping[1][1][ll_fnei]
                        
                        ## bone judge            
                        ## ins part for con but del for read
                        if len(tem_tem_ins_content_con)>0:
                            
                            if len(tem_tem_ins_content_con) in len_del_group:
                                len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                            else:
                                len_del_group[len(tem_tem_ins_content_con)] = 1
                            tem_tem_ins_content_con = ''
                        ## del part
                        if len(tem_tem_ins_content_read)>0:
                            if len(tem_tem_ins_content_read) in len_ins_group:
                                len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                            else:
                                len_ins_group[len(tem_tem_ins_content_read)] = 1
                            tem_tem_ins_content_read = ''
                    
                    elif len(tem_ins_content_read)>0:
                        if len(tem_ins_content_read) in len_ins_group:
                            len_ins_group[len(tem_ins_content_read)] = len_ins_group[len(tem_ins_content_read)] + 1
                        else:
                            len_ins_group[len(tem_ins_content_read)] = 1
                        tem_ins_content_read = ''
                        
                    elif len(tem_ins_content_con)>0:
                        
                        if len(tem_ins_content_con) in len_del_group:
                            len_del_group[len(tem_ins_content_con)] = len_del_group[len(tem_ins_content_con)] + 1
                        else:
                            len_del_group[len(tem_ins_content_con)] = 1
                        
                        tem_ins_content_con = ''
                    ### previous bone or not neighboured ins
                    
                    tem_ins_content_read = ''
                    tem_ins_content_con = ''
                    ### conclude the previous
                    if tem_new_bone_p_o in tem_ins_dic_per_read:
                        tem_ins_content_read = tem_ins_content_read + tem_ins_dic_per_read[tem_new_bone_p_o]
                    if tem_new_bone_p_o in Best_ref_ins[kk]:
                        tem_ins_content_con = tem_ins_content_con + Best_ref_ins[kk][tem_new_bone_p_o]
                else:
                    if tem_new_bone_p_o in tem_ins_dic_per_read:
                        tem_ins_content_read = tem_ins_content_read + tem_ins_dic_per_read[tem_new_bone_p_o]
                    if tem_new_bone_p_o in Best_ref_ins[kk]:
                        tem_ins_content_con = tem_ins_content_con + Best_ref_ins[kk][tem_new_bone_p_o]
                        if tem_new_bone_p_o not in tem_ins_dic_per_read:
                            previous_del_flag = 1
                            previous_del_bone = tem_new_bone_p_o### we borrow the previous ones
                    
                
                
                    
            ll_per_id = ll_per_id + 1
            pre_bone_ = tem_new_bone_p_o
            while ll_per_id < len(infortive_ins_collection[kk]):
               
                tem_new_bone_p_o = infortive_ins_collection[kk][ll_per_id]
                
                if tem_new_bone_p_o in only_bone_note_df_keys_array_sub:
                    tem_infer_base_con = Best_ref[kk][tem_new_bone_p_o]
                    tem_infer_base_read = tem_dic_bone[tem_new_bone_p_o]
                    if tem_infer_base_con <4:### we have the true base
                        if tem_infer_base_read in Standardbase[:4]:
                            if len(tem_ins_content_con)>0 and len(tem_ins_content_read)>0: # con and read in the bone
                                sub_string_con_string = tem_ins_content_con
                                sub_string_read_string = tem_ins_content_read
                                recollection = []
                                seq_collection = [sub_string_con_string, sub_string_read_string]
                                                    
                                for seq_index in range(len(seq_collection)):
                                    tem_seq = seqbio.NucleotideSequence(
                                        seq_collection[seq_index])
                                    recollection.append(tem_seq)
                                try:
                                    alignment, order, guide_tree, distance_matrix = align.align_multiple(
                                        recollection,
                                        matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                        gap_penalty=(-10, -0.5),
                                        terminal_penalty=False
                                    )
                                    alignment_trace = alignment.trace
                                    length_compare, num_col = np.shape(
                                        alignment_trace)
                                    result_mapping = []
                                    for col_index in range(num_col):
                                        #tem_result_mapping = [recollection[col_index]]
                                        coverted_seq = converting_code_2_seq(
                                            alignment_trace[:, col_index], recollection[col_index])
                                        # tem_result_mapping.append(coverted_seq)
                                        result_mapping.append(
                                            [seq_collection[col_index], coverted_seq.lower()])
                                except (ValueError, ZeroDivisionError):
                                    max_len = len(seq_collection[0])
                                    for col_index in range(len(seq_collection)):
                                        if len(seq_collection[col_index]) > max_len:
                                            max_len = len(
                                                seq_collection[col_index])
                                    result_mapping = []
                                    for col_index in range(len(seq_collection)):
                                        if len(seq_collection[col_index]) < max_len:
                                            refined_ins_ = seq_collection[col_index]+'N'*(
                                                max_len-len(seq_collection[col_index]))
                
                                        else:
                                            refined_ins_ = seq_collection[col_index]
                                        result_mapping.append(
                                            [seq_collection[col_index], refined_ins_])
                                
                                tem_tem_ins_content_con = ''
                                tem_tem_ins_content_read = ''
                                #Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5] = Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5]+len(result_mapping[0][1])
                                for ll_fnei in range(len(result_mapping[0][1])):
                                    if result_mapping[0][1][ll_fnei] in Standardbase[:4] and result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                        #if result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                        ## ins part for con but del for read
                                        if len(tem_tem_ins_content_con)>0:
                                            
                                            if len(tem_tem_ins_content_con) in len_del_group:
                                                len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                            else:
                                                len_del_group[len(tem_tem_ins_content_con)] = 1
                                            tem_tem_ins_content_con = ''
                                        ## del part
                                        if len(tem_tem_ins_content_read)>0:
                                            if len(tem_tem_ins_content_read) in len_ins_group:
                                                len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                            else:
                                                len_ins_group[len(tem_tem_ins_content_read)] = 1
                                            tem_tem_ins_content_read = ''
                                        ## bone judge
                                         
                                        if result_mapping[0][1][ll_fnei] == result_mapping[1][1][ll_fnei]:
                                            tem_correct_per = tem_correct_per + 1
                                        else:
                                            tem_mis_per = tem_mis_per + 1
                                            
                                    elif result_mapping[0][1][ll_fnei] in Standardbase[:4]:
                                        tem_tem_ins_content_con = tem_tem_ins_content_con + result_mapping[0][1][ll_fnei]
                                        
                                    elif result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                        tem_tem_ins_content_read = tem_tem_ins_content_read + result_mapping[1][1][ll_fnei]
                                
                                ## bone judge            
                                ## ins part for con but del for read
                                if len(tem_tem_ins_content_con)>0:
                                    
                                    if len(tem_tem_ins_content_con) in len_del_group:
                                        len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                    else:
                                        len_del_group[len(tem_tem_ins_content_con)] = 1
                                    tem_tem_ins_content_con = ''
                                ## del part
                                if len(tem_tem_ins_content_read)>0:
                                    if len(tem_tem_ins_content_read) in len_ins_group:
                                        len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                    else:
                                        len_ins_group[len(tem_tem_ins_content_read)] = 1
                                    tem_tem_ins_content_read = ''
                            
                            elif len(tem_ins_content_read)>0:
                                if len(tem_ins_content_read) in len_ins_group:
                                    len_ins_group[len(tem_ins_content_read)] = len_ins_group[len(tem_ins_content_read)] + 1
                                else:
                                    len_ins_group[len(tem_ins_content_read)] = 1
                                tem_ins_content_read = ''
                                
                            elif len(tem_ins_content_con)>0:
                                
                                if len(tem_ins_content_con) in len_del_group:
                                    len_del_group[len(tem_ins_content_con)] = len_del_group[len(tem_ins_content_con)] + 1
                                else:
                                    len_del_group[len(tem_ins_content_con)] = 1
                                
                                tem_ins_content_con = ''
                               
                            
                            if Standardbase[tem_infer_base_con] == tem_infer_base_read:
                                tem_correct_per = tem_correct_per + 1
                            else:
                                tem_mis_per = tem_mis_per + 1
                            previous_del_flag = 0
                            previous_del_bone = -1
                            tem_ins_content_read = ''
                            tem_ins_content_con = ''
                        else:### this should be the deletion/ common deletion is just like the tunnel to link two deletion neighboured with common deletion
                            if previous_del_flag: # neighboured del do not have inserton 
                                if tem_new_bone_p_o - previous_del_bone<2:
                                    #tem_del_len = tem_del_len + 1
                                    tem_ins_content_con = tem_ins_content_con + Standardbase[tem_infer_base_con]
                                    #tem_del_len = 1 ## check the tunnel effect of the common deletion
                                    if len(common_del_block_df_cp)>0:
                                        sub_df_common_del_left = common_del_block_df_cp[common_del_block_df_cp.iloc[:,3]<tem_new_bone_p_o]
                                        if len(sub_df_common_del_left)>0:
                                            last_line_bind_del = sub_df_common_del_left.iloc[-1,3]
                                            if tem_new_bone_p_o - last_line_bind_del<2:##left join
                                                #tem_del_length = tem_del_length + sub_left_df_del.iloc[-1,3] - sub_left_df_del.iloc[-1,2] + 1
                                                #common_del_block_df_cp.iloc[len(sub_left_df_del)-1,4] = 1
                                                #index_del_o = len(sub_left_df_del)-1
                                                #covered_common_del_id.append(index_del_o)
                                                previous_del_bone = common_del_block_df_cp.iloc[len(sub_df_common_del_left)-1,3]
                                            else:
                                                previous_del_bone = tem_new_bone_p_o
                                        else:
                                            previous_del_bone = tem_new_bone_p_o
                                    else:
                                        previous_del_bone = tem_new_bone_p_o
                                    previous_del_flag = 1
                                    tem_ins_content_read = ''
                                    #tem_ins_content_con = Standardbase[tem_infer_base_con]
                                    
                                else:
                                    ### it may occur not neighboured middle insertion
                                    if len(tem_ins_content_con)>0 and len(tem_ins_content_read)>0: # con and read in the bone
                                        sub_string_con_string = tem_ins_content_con
                                        sub_string_read_string = tem_ins_content_read
                                        recollection = []
                                        seq_collection = [sub_string_con_string, sub_string_read_string]
                                                            
                                        for seq_index in range(len(seq_collection)):
                                            tem_seq = seqbio.NucleotideSequence(
                                                seq_collection[seq_index])
                                            recollection.append(tem_seq)
                                        try:
                                            alignment, order, guide_tree, distance_matrix = align.align_multiple(
                                                recollection,
                                                matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                                gap_penalty=(-10, -0.5),
                                                terminal_penalty=False
                                            )
                                            alignment_trace = alignment.trace
                                            length_compare, num_col = np.shape(
                                                alignment_trace)
                                            result_mapping = []
                                            for col_index in range(num_col):
                                                #tem_result_mapping = [recollection[col_index]]
                                                coverted_seq = converting_code_2_seq(
                                                    alignment_trace[:, col_index], recollection[col_index])
                                                # tem_result_mapping.append(coverted_seq)
                                                result_mapping.append(
                                                    [seq_collection[col_index], coverted_seq.lower()])
                                        except (ValueError, ZeroDivisionError):
                                            max_len = len(seq_collection[0])
                                            for col_index in range(len(seq_collection)):
                                                if len(seq_collection[col_index]) > max_len:
                                                    max_len = len(
                                                        seq_collection[col_index])
                                            result_mapping = []
                                            for col_index in range(len(seq_collection)):
                                                if len(seq_collection[col_index]) < max_len:
                                                    refined_ins_ = seq_collection[col_index]+'N'*(
                                                        max_len-len(seq_collection[col_index]))
                        
                                                else:
                                                    refined_ins_ = seq_collection[col_index]
                                                result_mapping.append(
                                                    [seq_collection[col_index], refined_ins_])
                                        
                                        tem_tem_ins_content_con = ''
                                        tem_tem_ins_content_read = ''
                                        #Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5] = Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5]+len(result_mapping[0][1])
                                        for ll_fnei in range(len(result_mapping[0][1])):
                                            if result_mapping[0][1][ll_fnei] in Standardbase[:4] and result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                                #if result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                                ## ins part for con but del for read
                                                if len(tem_tem_ins_content_con)>0:
                                                    
                                                    if len(tem_tem_ins_content_con) in len_del_group:
                                                        len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                                    else:
                                                        len_del_group[len(tem_tem_ins_content_con)] = 1
                                                    tem_tem_ins_content_con = ''
                                                ## del part
                                                if len(tem_tem_ins_content_read)>0:
                                                    if len(tem_tem_ins_content_read) in len_ins_group:
                                                        len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                                    else:
                                                        len_ins_group[len(tem_tem_ins_content_read)] = 1
                                                    tem_tem_ins_content_read = ''
                                                ## bone judge
                                                 
                                                if result_mapping[0][1][ll_fnei] == result_mapping[1][1][ll_fnei]:
                                                    tem_correct_per = tem_correct_per + 1
                                                else:
                                                    tem_mis_per = tem_mis_per + 1
                                                    
                                            elif result_mapping[0][1][ll_fnei] in Standardbase[:4]:
                                                tem_tem_ins_content_con = tem_tem_ins_content_con + result_mapping[0][1][ll_fnei]
                                                
                                            elif result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                                tem_tem_ins_content_read = tem_tem_ins_content_read + result_mapping[1][1][ll_fnei]
                                        
                                        ## bone judge            
                                        ## ins part for con but del for read
                                        if len(tem_tem_ins_content_con)>0:
                                            
                                            if len(tem_tem_ins_content_con) in len_del_group:
                                                len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                            else:
                                                len_del_group[len(tem_tem_ins_content_con)] = 1
                                            tem_tem_ins_content_con = ''
                                        ## del part
                                        if len(tem_tem_ins_content_read)>0:
                                            if len(tem_tem_ins_content_read) in len_ins_group:
                                                len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                            else:
                                                len_ins_group[len(tem_tem_ins_content_read)] = 1
                                            tem_tem_ins_content_read = ''
                                    
                                    elif len(tem_ins_content_read)>0:
                                        if len(tem_ins_content_read) in len_ins_group:
                                            len_ins_group[len(tem_ins_content_read)] = len_ins_group[len(tem_ins_content_read)] + 1
                                        else:
                                            len_ins_group[len(tem_ins_content_read)] = 1
                                        tem_ins_content_read = ''
                                        
                                    elif len(tem_ins_content_con)>0:
                                        
                                        if len(tem_ins_content_con) in len_del_group:
                                            len_del_group[len(tem_ins_content_con)] = len_del_group[len(tem_ins_content_con)] + 1
                                        else:
                                            len_del_group[len(tem_ins_content_con)] = 1
                                        
                                        tem_ins_content_con = ''
                                    #tem_del_len = 1 ## check the tunnel effect of the common deletion
                                    if len(common_del_block_df_cp)>0:
                                        sub_df_common_del_left = common_del_block_df_cp[common_del_block_df_cp.iloc[:,3]<tem_new_bone_p_o]
                                        if len(sub_df_common_del_left)>0:
                                            last_line_bind_del = sub_df_common_del_left.iloc[-1,3]
                                            if tem_new_bone_p_o - last_line_bind_del<2:##left join
                                                #tem_del_length = tem_del_length + sub_left_df_del.iloc[-1,3] - sub_left_df_del.iloc[-1,2] + 1
                                                #common_del_block_df_cp.iloc[len(sub_left_df_del)-1,4] = 1
                                                #index_del_o = len(sub_left_df_del)-1
                                                #covered_common_del_id.append(index_del_o)
                                                previous_del_bone = common_del_block_df_cp.iloc[len(sub_df_common_del_left)-1,3]
                                            else:
                                                previous_del_bone = tem_new_bone_p_o
                                        else:
                                            previous_del_bone = tem_new_bone_p_o
                                    else:
                                        previous_del_bone = tem_new_bone_p_o
                                    previous_del_flag = 1
                                    tem_ins_content_read = ''
                                    tem_ins_content_con = Standardbase[tem_infer_base_con]
                                    
                            else:
                                ### it may occur not neighboured middle insertion
                                if len(tem_ins_content_con)>0 and len(tem_ins_content_read)>0: # con and read in the bone
                                    sub_string_con_string = tem_ins_content_con
                                    sub_string_read_string = tem_ins_content_read
                                    recollection = []
                                    seq_collection = [sub_string_con_string, sub_string_read_string]
                                                        
                                    for seq_index in range(len(seq_collection)):
                                        tem_seq = seqbio.NucleotideSequence(
                                            seq_collection[seq_index])
                                        recollection.append(tem_seq)
                                    try:
                                        alignment, order, guide_tree, distance_matrix = align.align_multiple(
                                            recollection,
                                            matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                            gap_penalty=(-10, -0.5),
                                            terminal_penalty=False
                                        )
                                        alignment_trace = alignment.trace
                                        length_compare, num_col = np.shape(
                                            alignment_trace)
                                        result_mapping = []
                                        for col_index in range(num_col):
                                            #tem_result_mapping = [recollection[col_index]]
                                            coverted_seq = converting_code_2_seq(
                                                alignment_trace[:, col_index], recollection[col_index])
                                            # tem_result_mapping.append(coverted_seq)
                                            result_mapping.append(
                                                [seq_collection[col_index], coverted_seq.lower()])
                                    except (ValueError, ZeroDivisionError):
                                        max_len = len(seq_collection[0])
                                        for col_index in range(len(seq_collection)):
                                            if len(seq_collection[col_index]) > max_len:
                                                max_len = len(
                                                    seq_collection[col_index])
                                        result_mapping = []
                                        for col_index in range(len(seq_collection)):
                                            if len(seq_collection[col_index]) < max_len:
                                                refined_ins_ = seq_collection[col_index]+'N'*(
                                                    max_len-len(seq_collection[col_index]))
                    
                                            else:
                                                refined_ins_ = seq_collection[col_index]
                                            result_mapping.append(
                                                [seq_collection[col_index], refined_ins_])
                                    
                                    tem_tem_ins_content_con = ''
                                    tem_tem_ins_content_read = ''
                                    #Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5] = Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5]+len(result_mapping[0][1])
                                    for ll_fnei in range(len(result_mapping[0][1])):
                                        if result_mapping[0][1][ll_fnei] in Standardbase[:4] and result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                            #if result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                            ## ins part for con but del for read
                                            if len(tem_tem_ins_content_con)>0:
                                                
                                                if len(tem_tem_ins_content_con) in len_del_group:
                                                    len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                                else:
                                                    len_del_group[len(tem_tem_ins_content_con)] = 1
                                                tem_tem_ins_content_con = ''
                                            ## del part
                                            if len(tem_tem_ins_content_read)>0:
                                                if len(tem_tem_ins_content_read) in len_ins_group:
                                                    len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                                else:
                                                    len_ins_group[len(tem_tem_ins_content_read)] = 1
                                                tem_tem_ins_content_read = ''
                                            ## bone judge
                                             
                                            if result_mapping[0][1][ll_fnei] == result_mapping[1][1][ll_fnei]:
                                                tem_correct_per = tem_correct_per + 1
                                            else:
                                                tem_mis_per = tem_mis_per + 1
                                                
                                        elif result_mapping[0][1][ll_fnei] in Standardbase[:4]:
                                            tem_tem_ins_content_con = tem_tem_ins_content_con + result_mapping[0][1][ll_fnei]
                                            
                                        elif result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                            tem_tem_ins_content_read = tem_tem_ins_content_read + result_mapping[1][1][ll_fnei]
                                    
                                    ## bone judge            
                                    ## ins part for con but del for read
                                    if len(tem_tem_ins_content_con)>0:
                                        
                                        if len(tem_tem_ins_content_con) in len_del_group:
                                            len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                        else:
                                            len_del_group[len(tem_tem_ins_content_con)] = 1
                                        tem_tem_ins_content_con = ''
                                    ## del part
                                    if len(tem_tem_ins_content_read)>0:
                                        if len(tem_tem_ins_content_read) in len_ins_group:
                                            len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                        else:
                                            len_ins_group[len(tem_tem_ins_content_read)] = 1
                                        tem_tem_ins_content_read = ''
                                
                                elif len(tem_ins_content_read)>0:
                                    if len(tem_ins_content_read) in len_ins_group:
                                        len_ins_group[len(tem_ins_content_read)] = len_ins_group[len(tem_ins_content_read)] + 1
                                    else:
                                        len_ins_group[len(tem_ins_content_read)] = 1
                                    tem_ins_content_read = ''
                                    
                                elif len(tem_ins_content_con)>0:
                                    
                                    if len(tem_ins_content_con) in len_del_group:
                                        len_del_group[len(tem_ins_content_con)] = len_del_group[len(tem_ins_content_con)] + 1
                                    else:
                                        len_del_group[len(tem_ins_content_con)] = 1
                                    
                                    tem_ins_content_con = ''
                                   
                                
                                #tem_del_len = 1 ## check the tunnel effect of the common deletion
                                if len(common_del_block_df_cp)>0:
                                    sub_df_common_del_left = common_del_block_df_cp[common_del_block_df_cp.iloc[:,3]<tem_new_bone_p_o]
                                    if len(sub_df_common_del_left)>0:
                                        
                                        last_line_bind_del = sub_df_common_del_left.iloc[-1,3]
                                        if tem_new_bone_p_o - last_line_bind_del<2:##left join
                                            #tem_del_length = tem_del_length + sub_left_df_del.iloc[-1,3] - sub_left_df_del.iloc[-1,2] + 1
                                            #common_del_block_df_cp.iloc[len(sub_left_df_del)-1,4] = 1
                                            #index_del_o = len(sub_left_df_del)-1
                                            #covered_common_del_id.append(index_del_o)
                                            previous_del_bone = common_del_block_df_cp.iloc[len(sub_df_common_del_left)-1,3]
                                        
                                        else:
                                            previous_del_bone = tem_new_bone_p_o
                                    else:
                                        previous_del_bone = tem_new_bone_p_o
                                else:
                                    previous_del_bone = tem_new_bone_p_o
                                previous_del_flag = 1
                                tem_ins_content_read = ''
                                tem_ins_content_con = Standardbase[tem_infer_base_con]
                    else:#not occur in the con but read
                        if tem_infer_base_read in Standardbase[:4]:
                            if tem_new_bone_p_o - pre_bone_<2:
                                tem_ins_content_read = tem_ins_content_read + tem_infer_base_read
                            else:
                                ### it may occur not neighboured middle insertion
                                if len(tem_ins_content_con)>0 and len(tem_ins_content_read)>0: # con and read in the bone
                                    sub_string_con_string = tem_ins_content_con
                                    sub_string_read_string = tem_ins_content_read
                                    recollection = []
                                    seq_collection = [sub_string_con_string, sub_string_read_string]
                                                        
                                    for seq_index in range(len(seq_collection)):
                                        tem_seq = seqbio.NucleotideSequence(
                                            seq_collection[seq_index])
                                        recollection.append(tem_seq)
                                    try:
                                        alignment, order, guide_tree, distance_matrix = align.align_multiple(
                                            recollection,
                                            matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                            gap_penalty=(-10, -0.5),
                                            terminal_penalty=False
                                        )
                                        alignment_trace = alignment.trace
                                        length_compare, num_col = np.shape(
                                            alignment_trace)
                                        result_mapping = []
                                        for col_index in range(num_col):
                                            #tem_result_mapping = [recollection[col_index]]
                                            coverted_seq = converting_code_2_seq(
                                                alignment_trace[:, col_index], recollection[col_index])
                                            # tem_result_mapping.append(coverted_seq)
                                            result_mapping.append(
                                                [seq_collection[col_index], coverted_seq.lower()])
                                    except (ValueError, ZeroDivisionError):
                                        max_len = len(seq_collection[0])
                                        for col_index in range(len(seq_collection)):
                                            if len(seq_collection[col_index]) > max_len:
                                                max_len = len(
                                                    seq_collection[col_index])
                                        result_mapping = []
                                        for col_index in range(len(seq_collection)):
                                            if len(seq_collection[col_index]) < max_len:
                                                refined_ins_ = seq_collection[col_index]+'N'*(
                                                    max_len-len(seq_collection[col_index]))
                    
                                            else:
                                                refined_ins_ = seq_collection[col_index]
                                            result_mapping.append(
                                                [seq_collection[col_index], refined_ins_])
                                    
                                    tem_tem_ins_content_con = ''
                                    tem_tem_ins_content_read = ''
                                    #Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5] = Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5]+len(result_mapping[0][1])
                                    for ll_fnei in range(len(result_mapping[0][1])):
                                        if result_mapping[0][1][ll_fnei] in Standardbase[:4] and result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                            #if result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                            ## ins part for con but del for read
                                            if len(tem_tem_ins_content_con)>0:
                                                
                                                if len(tem_tem_ins_content_con) in len_del_group:
                                                    len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                                else:
                                                    len_del_group[len(tem_tem_ins_content_con)] = 1
                                                tem_tem_ins_content_con = ''
                                            ## del part
                                            if len(tem_tem_ins_content_read)>0:
                                                if len(tem_tem_ins_content_read) in len_ins_group:
                                                    len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                                else:
                                                    len_ins_group[len(tem_tem_ins_content_read)] = 1
                                                tem_tem_ins_content_read = ''
                                            ## bone judge
                                             
                                            if result_mapping[0][1][ll_fnei] == result_mapping[1][1][ll_fnei]:
                                                tem_correct_per = tem_correct_per + 1
                                            else:
                                                tem_mis_per = tem_mis_per + 1
                                                
                                        elif result_mapping[0][1][ll_fnei] in Standardbase[:4]:
                                            tem_tem_ins_content_con = tem_tem_ins_content_con + result_mapping[0][1][ll_fnei]
                                            
                                        elif result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                            tem_tem_ins_content_read = tem_tem_ins_content_read + result_mapping[1][1][ll_fnei]
                                    
                                    ## bone judge            
                                    ## ins part for con but del for read
                                    if len(tem_tem_ins_content_con)>0:
                                        
                                        if len(tem_tem_ins_content_con) in len_del_group:
                                            len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                        else:
                                            len_del_group[len(tem_tem_ins_content_con)] = 1
                                        tem_tem_ins_content_con = ''
                                    ## del part
                                    if len(tem_tem_ins_content_read)>0:
                                        if len(tem_tem_ins_content_read) in len_ins_group:
                                            len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                        else:
                                            len_ins_group[len(tem_tem_ins_content_read)] = 1
                                        tem_tem_ins_content_read = ''
                                
                                elif len(tem_ins_content_read)>0:
                                    if len(tem_ins_content_read) in len_ins_group:
                                        len_ins_group[len(tem_ins_content_read)] = len_ins_group[len(tem_ins_content_read)] + 1
                                    else:
                                        len_ins_group[len(tem_ins_content_read)] = 1
                                    tem_ins_content_read = ''
                                    
                                elif len(tem_ins_content_con)>0:
                                    
                                    if len(tem_ins_content_con) in len_del_group:
                                        len_del_group[len(tem_ins_content_con)] = len_del_group[len(tem_ins_content_con)] + 1
                                    else:
                                        len_del_group[len(tem_ins_content_con)] = 1
                                    
                                    tem_ins_content_con = ''
                                   
                                
                                #previous_del_flag = 0
                                tem_ins_content_read = tem_infer_base_read
                                #tem_ins_content_con = ''
                        
                            previous_del_bone = -1
                            previous_del_flag = 0
                            tem_ins_content_con = ''
                        #else is the two del d; do nothin
                        
                    if tem_new_bone_p_o in tem_ins_dic_per_read:
                        tem_ins_content_read = tem_ins_content_read + tem_ins_dic_per_read[tem_new_bone_p_o]
                    if tem_new_bone_p_o in Best_ref_ins[kk]:
                        tem_ins_content_con = tem_ins_content_con + Best_ref_ins[kk][tem_new_bone_p_o]
                else:
                    previous_del_flag = 0
                    previous_del_bone = -1
                    if tem_new_bone_p_o - pre_bone_>=2:
                        ### ins seperated
                        if len(tem_ins_content_con)>0 and len(tem_ins_content_read)>0: # con and read in the bone
                            sub_string_con_string = tem_ins_content_con
                            sub_string_read_string = tem_ins_content_read
                            recollection = []
                            seq_collection = [sub_string_con_string, sub_string_read_string]
                                                
                            for seq_index in range(len(seq_collection)):
                                tem_seq = seqbio.NucleotideSequence(
                                    seq_collection[seq_index])
                                recollection.append(tem_seq)
                            try:
                                alignment, order, guide_tree, distance_matrix = align.align_multiple(
                                    recollection,
                                    matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                    gap_penalty=(-10, -0.5),
                                    terminal_penalty=False
                                )
                                alignment_trace = alignment.trace
                                length_compare, num_col = np.shape(
                                    alignment_trace)
                                result_mapping = []
                                for col_index in range(num_col):
                                    #tem_result_mapping = [recollection[col_index]]
                                    coverted_seq = converting_code_2_seq(
                                        alignment_trace[:, col_index], recollection[col_index])
                                    # tem_result_mapping.append(coverted_seq)
                                    result_mapping.append(
                                        [seq_collection[col_index], coverted_seq.lower()])
                            except (ValueError, ZeroDivisionError):
                                max_len = len(seq_collection[0])
                                for col_index in range(len(seq_collection)):
                                    if len(seq_collection[col_index]) > max_len:
                                        max_len = len(
                                            seq_collection[col_index])
                                result_mapping = []
                                for col_index in range(len(seq_collection)):
                                    if len(seq_collection[col_index]) < max_len:
                                        refined_ins_ = seq_collection[col_index]+'N'*(
                                            max_len-len(seq_collection[col_index]))
            
                                    else:
                                        refined_ins_ = seq_collection[col_index]
                                    result_mapping.append(
                                        [seq_collection[col_index], refined_ins_])
                            
                            tem_tem_ins_content_con = ''
                            tem_tem_ins_content_read = ''
                            #Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5] = Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5]+len(result_mapping[0][1])
                            for ll_fnei in range(len(result_mapping[0][1])):
                                if result_mapping[0][1][ll_fnei] in Standardbase[:4] and result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                    #if result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                    ## ins part for con but del for read
                                    if len(tem_tem_ins_content_con)>0:
                                        
                                        if len(tem_tem_ins_content_con) in len_del_group:
                                            len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                        else:
                                            len_del_group[len(tem_tem_ins_content_con)] = 1
                                        tem_tem_ins_content_con = ''
                                    ## del part
                                    if len(tem_tem_ins_content_read)>0:
                                        if len(tem_tem_ins_content_read) in len_ins_group:
                                            len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                        else:
                                            len_ins_group[len(tem_tem_ins_content_read)] = 1
                                        tem_tem_ins_content_read = ''
                                    ## bone judge
                                     
                                    if result_mapping[0][1][ll_fnei] == result_mapping[1][1][ll_fnei]:
                                        tem_correct_per = tem_correct_per + 1
                                    else:
                                        tem_mis_per = tem_mis_per + 1
                                        
                                elif result_mapping[0][1][ll_fnei] in Standardbase[:4]:
                                    tem_tem_ins_content_con = tem_tem_ins_content_con + result_mapping[0][1][ll_fnei]
                                    
                                elif result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                                    tem_tem_ins_content_read = tem_tem_ins_content_read + result_mapping[1][1][ll_fnei]
                            
                            ## bone judge            
                            ## ins part for con but del for read
                            if len(tem_tem_ins_content_con)>0:
                                
                                if len(tem_tem_ins_content_con) in len_del_group:
                                    len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                                else:
                                    len_del_group[len(tem_tem_ins_content_con)] = 1
                                tem_tem_ins_content_con = ''
                            ## del part
                            if len(tem_tem_ins_content_read)>0:
                                if len(tem_tem_ins_content_read) in len_ins_group:
                                    len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                                else:
                                    len_ins_group[len(tem_tem_ins_content_read)] = 1
                                tem_tem_ins_content_read = ''
                        
                        elif len(tem_ins_content_read)>0:
                            if len(tem_ins_content_read) in len_ins_group:
                                len_ins_group[len(tem_ins_content_read)] = len_ins_group[len(tem_ins_content_read)] + 1
                            else:
                                len_ins_group[len(tem_ins_content_read)] = 1
                            tem_ins_content_read = ''
                            
                        elif len(tem_ins_content_con)>0:
                            
                            if len(tem_ins_content_con) in len_del_group:
                                len_del_group[len(tem_ins_content_con)] = len_del_group[len(tem_ins_content_con)] + 1
                            else:
                                len_del_group[len(tem_ins_content_con)] = 1
                            
                            tem_ins_content_con = ''
                        ### previous bone or not neighboured ins
                        
                        tem_ins_content_read = ''
                        tem_ins_content_con = ''
                        ### conclude the previous
                        if tem_new_bone_p_o in tem_ins_dic_per_read:
                            tem_ins_content_read = tem_ins_content_read + tem_ins_dic_per_read[tem_new_bone_p_o]
                        if tem_new_bone_p_o in Best_ref_ins[kk]:
                            tem_ins_content_con = tem_ins_content_con + Best_ref_ins[kk][tem_new_bone_p_o]
                    else:
                        if tem_new_bone_p_o in tem_ins_dic_per_read:
                            tem_ins_content_read = tem_ins_content_read + tem_ins_dic_per_read[tem_new_bone_p_o]
                        if tem_new_bone_p_o in Best_ref_ins[kk]:
                            tem_ins_content_con = tem_ins_content_con + Best_ref_ins[kk][tem_new_bone_p_o]
                            if tem_new_bone_p_o not in tem_ins_dic_per_read:
                                previous_del_flag = 1
                                previous_del_bone = tem_new_bone_p_o### we borrow the previous ones
                        
                    
                    
                pre_bone_ = tem_new_bone_p_o
                ll_per_id = ll_per_id + 1
            
            
            ### tail checking
            if len(tem_ins_content_con)>0 and len(tem_ins_content_read)>0: # con and read in the bone
                sub_string_con_string = tem_ins_content_con
                sub_string_read_string = tem_ins_content_read
                recollection = []
                seq_collection = [sub_string_con_string, sub_string_read_string]
                                    
                for seq_index in range(len(seq_collection)):
                    tem_seq = seqbio.NucleotideSequence(
                        seq_collection[seq_index])
                    recollection.append(tem_seq)
                try:
                    alignment, order, guide_tree, distance_matrix = align.align_multiple(
                        recollection,
                        matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                        gap_penalty=(-10, -0.5),
                        terminal_penalty=False
                    )
                    alignment_trace = alignment.trace
                    length_compare, num_col = np.shape(
                        alignment_trace)
                    result_mapping = []
                    for col_index in range(num_col):
                        #tem_result_mapping = [recollection[col_index]]
                        coverted_seq = converting_code_2_seq(
                            alignment_trace[:, col_index], recollection[col_index])
                        # tem_result_mapping.append(coverted_seq)
                        result_mapping.append(
                            [seq_collection[col_index], coverted_seq.lower()])
                except (ValueError, ZeroDivisionError):
                    max_len = len(seq_collection[0])
                    for col_index in range(len(seq_collection)):
                        if len(seq_collection[col_index]) > max_len:
                            max_len = len(
                                seq_collection[col_index])
                    result_mapping = []
                    for col_index in range(len(seq_collection)):
                        if len(seq_collection[col_index]) < max_len:
                            refined_ins_ = seq_collection[col_index]+'N'*(
                                max_len-len(seq_collection[col_index]))
    
                        else:
                            refined_ins_ = seq_collection[col_index]
                        result_mapping.append(
                            [seq_collection[col_index], refined_ins_])
                
                tem_tem_ins_content_con = ''
                tem_tem_ins_content_read = ''
                #Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5] = Table_true_variations.iloc[int(col_o3[ll])+num_rep*int(number_k),5]+len(result_mapping[0][1])
                for ll_fnei in range(len(result_mapping[0][1])):
                    if result_mapping[0][1][ll_fnei] in Standardbase[:4] and result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                        #if result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                        ## ins part for con but del for read
                        if len(tem_tem_ins_content_con)>0:
                            
                            if len(tem_tem_ins_content_con) in len_del_group:
                                len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                            else:
                                len_del_group[len(tem_tem_ins_content_con)] = 1
                            tem_tem_ins_content_con = ''
                        ## del part
                        if len(tem_tem_ins_content_read)>0:
                            if len(tem_tem_ins_content_read) in len_ins_group:
                                len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                            else:
                                len_ins_group[len(tem_tem_ins_content_read)] = 1
                            tem_tem_ins_content_read = ''
                        ## bone judge
                         
                        if result_mapping[0][1][ll_fnei] == result_mapping[1][1][ll_fnei]:
                            tem_correct_per = tem_correct_per + 1
                        else:
                            tem_mis_per = tem_mis_per + 1
                            
                    elif result_mapping[0][1][ll_fnei] in Standardbase[:4]:
                        tem_tem_ins_content_con = tem_tem_ins_content_con + result_mapping[0][1][ll_fnei]
                        
                    elif result_mapping[1][1][ll_fnei] in Standardbase[:4]:
                        tem_tem_ins_content_read = tem_tem_ins_content_read + result_mapping[1][1][ll_fnei]
                
                ## bone judge            
                ## ins part for con but del for read
                if len(tem_tem_ins_content_con)>0:
                    
                    if len(tem_tem_ins_content_con) in len_del_group:
                        len_del_group[len(tem_tem_ins_content_con)] = len_del_group[len(tem_tem_ins_content_con)] + 1
                    else:
                        len_del_group[len(tem_tem_ins_content_con)] = 1
                    tem_tem_ins_content_con = ''
                ## del part
                if len(tem_tem_ins_content_read)>0:
                    if len(tem_tem_ins_content_read) in len_ins_group:
                        len_ins_group[len(tem_tem_ins_content_read)] = len_ins_group[len(tem_tem_ins_content_read)] + 1
                    else:
                        len_ins_group[len(tem_tem_ins_content_read)] = 1
                    tem_tem_ins_content_read = ''
            
            elif len(tem_ins_content_read)>0:
                if len(tem_ins_content_read) in len_ins_group:
                    len_ins_group[len(tem_ins_content_read)] = len_ins_group[len(tem_ins_content_read)] + 1
                else:
                    len_ins_group[len(tem_ins_content_read)] = 1
                tem_ins_content_read = ''
                
            elif len(tem_ins_content_con)>0:
                
                if len(tem_ins_content_con) in len_del_group:
                    len_del_group[len(tem_ins_content_con)] = len_del_group[len(tem_ins_content_con)] + 1
                else:
                    len_del_group[len(tem_ins_content_con)] = 1
                
                tem_ins_content_con = ''
            
               
        
        ### previous bone or not neighboured ins
        
        tem_probabity_correct = tem_correct_per*np.log(pcorrect+corretion_eps_from_zero)
        tem_probabity_mis = tem_mis_per*np.log(p_mis_loc+corretion_eps_from_zero)
        if len(len_ins_group) > 0:
            values_ = len_ins_group.values()
            num_ins_loc = sum(np.array(list(values_)))
            tem_probabity_ins = num_ins_loc*np.log(pins1_loc+corretion_eps_from_zero)
            for ii in len_ins_group:
                # np.log(stats.expon.pdf(ii,lamuna_ins))
                tem_probabity_ins = tem_probabity_ins + \
                    len_ins_group[ii] * \
                    possion_log_probability(ii-1,lamuna_ins)
                # tem_probabity_ins = tem_probabity_ins + len_ins[ii]*np.log(stats.nbinom.pmf(ii,n_ins,p_ins))#np.log(stats.expon.pdf(ii,lamuna_ins))
                # tem_probabity_ins = tem_probabity_ins + len_ins[ii]*np.log(1/(ii*sig_ins*np.sqrt(2*math.pi))*np.exp(-(np.log(ii)-mu_ins)**2/(2*sig_ins**2)))#np.log(stats.expon.pdf(ii,lamuna_ins))
        else:
            tem_probabity_ins = 0
            num_ins_loc = 0

        del_whole = 0
        if len(len_del_group) > 0:
            values_ = len_del_group.values()
            num_del_loc = sum(np.array(list(values_)))
            tem_probabity_del = num_del_loc*np.log(pdel1_loc+corretion_eps_from_zero)
            for dd in len_del_group:
                # np.log(stats.expon.pdf(ii,lamuna_ins))
                tem_probabity_del = tem_probabity_del + \
                    len_del_group[dd] * \
                    possion_log_probability(dd-1,lamuna_del)
                # tem_probabity_del = tem_probabity_del + len_del[dd]*np.log(stats.nbinom.pmf(dd,n_del,p_del))#np.log(stats.expon.pdf(ii,lamuna_ins))
                # tem_probabity_del = tem_probabity_del + len_del[dd]*np.log(1/(dd*sig_del*np.sqrt(2*math.pi))*np.exp(-(np.log(dd)-mu_del)**2/(2*sig_del**2)))#np.log(stats.expon.pdf(ii,lamuna_ins))
                del_whole = del_whole + dd*len_del_group[dd]

        else:
            tem_probabity_del = 0
            num_del_loc = 0
        
        tem_match_resutn_list.append(tem_correct_per)
        tem_mismatch_resutn_list.append(tem_mis_per)
        
        tem_ins_return_dic.append(copy.deepcopy(len_ins_group))
        tem_del_return_dic.append(copy.deepcopy(len_del_group))
        
        other_gaps = tem_correct_per + tem_mis_per - 1-num_del_loc
        gap_not_ins_pro = other_gaps*np.log(1-pins1_loc+corretion_eps_from_zero)
        not_del_num = tem_correct_per + tem_mis_per
        gap_not_del_pro = not_del_num*np.log(1-pdel1_loc+corretion_eps_from_zero)
        full_prob = gap_not_del_pro+tem_probabity_correct + tem_probabity_mis +\
            tem_probabity_ins+tem_probabity_del+gap_not_ins_pro
        
        tem_prob_final[kk] = (full_prob+np.log(portions[kk]+corretion_eps_from_zero))/tao
    
    tem_prob_final = tem_prob_final - max(tem_prob_final)
    tem_prob_final = np.exp(tem_prob_final)
    tem_prob_final = tem_prob_final/sum(tem_prob_final)
    ### sampling
    l_s_vec = np.random.multinomial(
        n=1, pvals=tem_prob_final)
    group_id = list(l_s_vec).index(1)
    
    max_len_return = max([1,len(tem_ins_return_dic[group_id]),len(tem_del_return_dic[group_id])])
    final_return_result_per_read = [[-1,0,0,-1,-1,-1,-1]+[0]*K+[-1] for x in range(max_len_return)]
    final_return_result_per_read[0][1] = tem_match_resutn_list[group_id]
    final_return_result_per_read[0][2] = tem_mismatch_resutn_list[group_id]
    final_return_result_per_read[0][0] = read_xun_id
    final_return_result_per_read[0][7:(7+K)] = list(tem_prob_final)
    final_return_result_per_read[0][7+K] = group_id
    
    ### saving the del sequencing error
    tem_del_return_dic_sel_keys = list(tem_del_return_dic[group_id].keys())
    for seq_id_del in range(len(tem_del_return_dic_sel_keys)):
        final_return_result_per_read[seq_id_del][3:5] = [tem_del_return_dic_sel_keys[seq_id_del],tem_del_return_dic[group_id][tem_del_return_dic_sel_keys[seq_id_del]]]
    
    ### saving the ins sequencing error
    tem_ins_return_dic_sel_keys = list(tem_ins_return_dic[group_id].keys())
    for seq_id_ins in range(len(tem_ins_return_dic_sel_keys)):
        final_return_result_per_read[seq_id_ins][5:7] = [tem_ins_return_dic_sel_keys[seq_id_ins],tem_ins_return_dic[group_id][tem_ins_return_dic_sel_keys[seq_id_ins]]]
    
    
    
    return(final_return_result_per_read)
    

def comparision_between_con_and_ref(id_block):
    #global Result_infer_ , seg_collection_ini,full_infer_var_site_o_df,only_bone_note_df_keys,\
    #    common_ins_block_dic,common_ins_block_dic_pure_num ,common_del_block_df_cp,ref_code_standard,other_not_bone_df
    
#for id_block in range(len(seg_collection_ini)):
    #id_block = 8839
    tem_c = 0
    tem_mis = 0
    tem_del_length = 0
    tem_num_string = ''
    tem_num_string_full = ''
    #tem_del_collection = []
    #tem_ins_collection = []
    start_block_tem = seg_collection_ini[id_block][0]
    end_block_tem = seg_collection_ini[id_block][1]
    Best_ref_ins_list = []
    Best_ref_ins_list_full = []
    Best_ref_list = []
    Del_length = []
    
    covered_common_del_id = []
    ll_per_id = start_block_tem
    last_del_pos = -1
    while ll_per_id <= end_block_tem:
        ll_per = Result_infer_[ll_per_id]
        tem_new_bone_p = ll_per[0]
        tem_new_bone_p_sub_df =  full_infer_var_site_o_df[full_infer_var_site_o_df.iloc[:,0]==tem_new_bone_p]
        tem_new_bone_p_o = tem_new_bone_p_sub_df.iloc[0,1]
        tem_infer_bone_base = ll_per[1]
        tem_infer_state = tem_new_bone_p_sub_df.iloc[0,2]
        if tem_infer_state<0:
            Best_ref_list.append([tem_new_bone_p_o,tem_infer_bone_base])
            ll_per_id = ll_per_id + 1
            
            if tem_infer_bone_base==4: #deletion
                if last_del_pos>0 and tem_new_bone_p_o - last_del_pos>1:
                   
                    if tem_del_length>0:
                        
                        Del_length.append(tem_del_length)
                    
                    tem_del_length = 1
                    last_del_pos = tem_new_bone_p_o
                    ### we need find the neighboured del in the common part the part is left/right join in with the del
                    if len(common_del_block_df_cp)>0:
                        sub_left_df_del = common_del_block_df_cp[common_del_block_df_cp.iloc[:,3]<tem_new_bone_p_o]
                        if len(sub_left_df_del)>0:
                            last_line_bind_del = sub_left_df_del.iloc[-1,3]
                            if tem_new_bone_p_o - last_line_bind_del<2 and sub_left_df_del.iloc[-1,4]<1:##left join
                                tem_del_length = tem_del_length + sub_left_df_del.iloc[-1,3] - sub_left_df_del.iloc[-1,2] + 1
                                #common_del_block_df_cp.iloc[len(sub_left_df_del)-1,4] = 1
                                index_del_o = len(sub_left_df_del)-1
                                covered_common_del_id.append(index_del_o)
                                
                        sub_right_df_del = common_del_block_df_cp[common_del_block_df_cp.iloc[:,2]>tem_new_bone_p_o]
                        if len(sub_right_df_del)>0:
                            first_line_bind_del = sub_right_df_del.iloc[0,2]
                            if first_line_bind_del-tem_new_bone_p_o<2:
                                tem_del_length = tem_del_length + sub_right_df_del.iloc[0,3] - sub_right_df_del.iloc[0,2] + 1
                                #common_del_block_df_cp.iloc[len(sub_right_df_del)-1,4] = 1
                                last_del_pos =  sub_right_df_del.iloc[0,3]
                                index_del_o = list(common_del_block_df_cp.iloc[:,0]).index(sub_right_df_del.iloc[0,0]) 
                                covered_common_del_id.append(index_del_o)
                            
                    ### end of the possble common del
                 
                    
                else:### continus del
                    if last_del_pos>0:## here second judge can not achieve
                        last_del_pos = tem_new_bone_p
                        tem_del_length = tem_del_length + 1
                        if len(common_del_block_df_cp)>0:
                            sub_right_df_del = common_del_block_df_cp[common_del_block_df_cp.iloc[:,2]>tem_new_bone_p_o]
                            if len(sub_right_df_del)>0:
                                first_line_bind_del = sub_right_df_del.iloc[0,2]
                                if first_line_bind_del-tem_new_bone_p_o<2:
                                    tem_del_length = tem_del_length + sub_right_df_del.iloc[0,3] - sub_right_df_del.iloc[0,2] + 1
                                    #common_del_block_df_cp.iloc[len(sub_right_df_del)-1,4] = 1
                                    last_del_pos =  sub_right_df_del.iloc[0,3]
                                    index_del_o = list(common_del_block_df_cp.iloc[:,0]).index(sub_right_df_del.iloc[0,0]) 
                                    covered_common_del_id.append(index_del_o)
                    else:### first
                        tem_del_length = 1
                        ### we need find the neighboured del in the common part the part is left/right join in with the del
                        if len(common_del_block_df_cp)>0:
                            sub_left_df_del = common_del_block_df_cp[common_del_block_df_cp.iloc[:,3]<tem_new_bone_p_o]
                            if len(sub_left_df_del)>0:
                                last_line_bind_del = sub_left_df_del.iloc[-1,3]
                                if tem_new_bone_p_o - last_line_bind_del<2 and sub_left_df_del.iloc[-1,4]<1:##left join
                                    tem_del_length = tem_del_length + sub_left_df_del.iloc[-1,3] - sub_left_df_del.iloc[-1,2] + 1
                                    #common_del_block_df_cp.iloc[len(sub_left_df_del)-1,4] = 1
                                    index_del_o = len(sub_left_df_del)-1
                                    covered_common_del_id.append(index_del_o)
                                    
                            sub_right_df_del = common_del_block_df_cp[common_del_block_df_cp.iloc[:,2]>tem_new_bone_p_o]
                            if len(sub_right_df_del)>0:
                                first_line_bind_del = sub_right_df_del.iloc[0,2]
                                if first_line_bind_del-tem_new_bone_p_o<2:
                                    tem_del_length = tem_del_length + sub_right_df_del.iloc[0,3] - sub_right_df_del.iloc[0,2] + 1
                                    #common_del_block_df_cp.iloc[len(sub_right_df_del)-1,4] = 1
                                    last_del_pos =  sub_right_df_del.iloc[0,3]
                                    index_del_o = list(common_del_block_df_cp.iloc[:,0]).index(sub_right_df_del.iloc[0,0]) 
                                    covered_common_del_id.append(index_del_o)
                            ### end of the possble common del
                    
                
                ### we need find the neighboured del in the common part the part is left/right join in with the del
                if len(common_del_block_df_cp)>0:
                    sub_left_df_del = common_del_block_df_cp[common_del_block_df_cp.iloc[:,3]<tem_new_bone_p_o]
                    if len(sub_left_df_del)>0:
                       last_line_bind_del = sub_left_df_del.iloc[-1,3]
                       if tem_new_bone_p_o - last_line_bind_del<2 and sub_left_df_del.iloc[-1,4]<1:##left join
                           tem_del_length = tem_del_length + sub_left_df_del.iloc[-1,3] - sub_left_df_del.iloc[-1,2] + 1
                           #common_del_block_df_cp.iloc[len(sub_left_df_del)-1,4] = 1
                           index_del_o = len(sub_left_df_del)-1
                           covered_common_del_id.append(index_del_o)
                           
                    sub_right_df_del = common_del_block_df_cp[common_del_block_df_cp.iloc[:,2]>tem_new_bone_p_o]
                    if len(sub_right_df_del)>0:
                       first_line_bind_del = sub_right_df_del.iloc[0,2]
                       if first_line_bind_del-tem_new_bone_p_o<2:
                           tem_del_length = tem_del_length + sub_right_df_del.iloc[0,3] - sub_right_df_del.iloc[0,2] + 1
                           #common_del_block_df_cp.iloc[len(sub_right_df_del)-1,4] = 1
                           last_del_pos =  sub_right_df_del.iloc[0,3]
                           index_del_o = list(common_del_block_df_cp.iloc[:,0]).index(sub_right_df_del.iloc[0,0]) 
                           covered_common_del_id.append(index_del_o)
                ### end of the possble common del
            
            ### get the informative bone inference 
            else:
                #### judge whether they are equal
                if tem_del_length>0:
                    Del_length.append(tem_del_length)
                    tem_del_length = 0
                if tem_infer_bone_base==ref_code_standard[tem_new_bone_p_o]:
                    tem_c = tem_c + 1
                else:
                    tem_mis = tem_mis + 1
                    
                
        #covered_bone_ins
        if tem_new_bone_p_o in ins_inserted_dic:
            bone_k_num = ins_inserted_dic[tem_new_bone_p_o]
            if tem_new_bone_p_o in common_ins_block_dic:
                ins_segments = common_ins_block_dic[tem_new_bone_p_o]
                if len(ins_segments)>1:
                    nei_xun_common_s =  ins_segments[0][0]
                    nei_xun_common_e =  ins_segments[0][1]
                    ### ge the binding position of the ins
                    tem_num_string_com = ins_segments[0][2]
                    tem_num_string_com_num = common_ins_block_dic_pure_num[tem_new_bone_p_o][0][2]
                    #if bone_k_num -(len(tem_num_string_com))>=1:
                    #### the ins string may be divied into neighboured seg as ins_common_ins parts
                    tem_num_string = ''
                    tem_num_string_full = ''
                    for ll_pre_ins in range(nei_xun_common_s):
                        tem_num_string_full = tem_num_string_full + str(Result_infer_[ll_per_id][1])
                        if Result_infer_[ll_per_id][1] <4:
                            tem_num_string = tem_num_string + Standardbase[Result_infer_[ll_per_id][1]]
                            
                            
                            
                        ll_per_id = ll_per_id + 1
                        
                    tem_num_string = tem_num_string + tem_num_string_com
                    tem_num_string_full = tem_num_string_full +tem_num_string_com_num
                    
                    #for ll_aft_ins in range(nei_xun_common_e+1,bone_k_num):
                    #    tem_num_string_full = tem_num_string_full + str(Result_infer_[ll_per_id][1])
                    #    if Result_infer_[ll_per_id][1] <4:
                    #        tem_num_string = tem_num_string + Standardbase[Result_infer_[ll_per_id][1]]
                    #    ll_per_id = ll_per_id + 1
                    pre_end = nei_xun_common_e    
                    for sub_ins_segment_id in range(1,len(ins_segments)):
                        nei_xun_common_s =  ins_segments[sub_ins_segment_id][0]
                        nei_xun_common_e =  ins_segments[sub_ins_segment_id][1]
                        ### mid part between 2 segments
                        for mid_xun_ins in range(pre_end+1,nei_xun_common_s):
                            tem_num_string_full = tem_num_string_full + str(Result_infer_[ll_per_id][1])
                            if Result_infer_[ll_per_id][1] <4:
                                tem_num_string = tem_num_string + Standardbase[Result_infer_[ll_per_id][1]]
                                 
                            ll_per_id = ll_per_id + 1  
                                
                        
                        ### end of mid part
                        
                        ### ge the binding position of the ins
                        tem_num_string_com = ins_segments[sub_ins_segment_id][2]
                        tem_num_string_com_num = common_ins_block_dic_pure_num[tem_new_bone_p_o][sub_ins_segment_id][2]
                        
                        tem_num_string = tem_num_string + tem_num_string_com
                        tem_num_string_full = tem_num_string_full +tem_num_string_com_num
                        
                        pre_end = nei_xun_common_e 
                    
                    for ll_aft_ins in range(nei_xun_common_e+1,bone_k_num):
                        tem_num_string_full = tem_num_string_full + str(Result_infer_[ll_per_id][1])
                        if Result_infer_[ll_per_id][1] <4:
                            tem_num_string = tem_num_string + Standardbase[Result_infer_[ll_per_id][1]]
                        ll_per_id = ll_per_id + 1    
                            
                            
                    if len(tem_num_string)>0:
                        
                        Best_ref_ins_list.append([tem_new_bone_p_o,tem_num_string,len(tem_num_string)])
                        tem_num_string = ''
                            
                            
                        
                    Best_ref_ins_list_full.append([tem_new_bone_p_o,tem_num_string_full,len(tem_num_string_full)])
                    tem_num_string_full = ''
                
                else:
                    nei_xun_common_s =  ins_segments[0][0]
                    nei_xun_common_e =  ins_segments[0][1]
                    ### ge the binding position of the ins
                    tem_num_string_com = ins_segments[0][2]
                    tem_num_string_com_num = common_ins_block_dic_pure_num[tem_new_bone_p_o][0][2]
                    if bone_k_num -(len(tem_num_string_com))>=1:
                        tem_num_string = ''
                        tem_num_string_full = ''
                        for ll_pre_ins in range(nei_xun_common_s):
                            tem_num_string_full = tem_num_string_full + str(Result_infer_[ll_per_id][1])
                            if Result_infer_[ll_per_id][1] <4:
                                tem_num_string = tem_num_string + Standardbase[Result_infer_[ll_per_id][1]]
                                
                                
                                
                            ll_per_id = ll_per_id + 1
                        # mid
                        tem_num_string = tem_num_string + tem_num_string_com
                        tem_num_string_full = tem_num_string_full +tem_num_string_com_num
                        
                        # mid
                        
                        for ll_aft_ins in range(nei_xun_common_e+1,bone_k_num):
                            tem_num_string_full = tem_num_string_full + str(Result_infer_[ll_per_id][1])
                            if Result_infer_[ll_per_id][1] <4:
                                tem_num_string = tem_num_string + Standardbase[Result_infer_[ll_per_id][1]]
                            ll_per_id = ll_per_id + 1  
                        
                    else: ### this means the exact match the common ins
                        tem_num_string = tem_num_string_com
                        tem_num_string_full = tem_num_string_com_num
                        
                        Best_ref_ins_list.append([tem_new_bone_p_o,tem_num_string,len(tem_num_string)])
                        tem_num_string = ''
                    
                    Best_ref_ins_list_full.append([tem_new_bone_p_o,tem_num_string_full,len(tem_num_string_full)])
                    tem_num_string_full = ''
                    
                        
            else:
                tem_num_string = ''
                tem_num_string_full = ''
                for ll_pre_ins in range(bone_k_num):
                    tem_num_string_full = tem_num_string_full + str(Result_infer_[ll_per_id][1])
                    if Result_infer_[ll_per_id][1] <4:
                        tem_num_string = tem_num_string + Standardbase[Result_infer_[ll_per_id][1]]
                    ll_per_id = ll_per_id + 1
                
                if len(tem_num_string)>0:
                    
                    Best_ref_ins_list.append([tem_new_bone_p_o,tem_num_string,len(tem_num_string)])
                    tem_num_string = ''
                Best_ref_ins_list_full.append([tem_new_bone_p_o,tem_num_string_full,len(tem_num_string_full)])
                tem_num_string_full = ''
            #tem_touched_ins.append(tem_new_bone_p_o)#record the inhomo sites insertion
        
        
    if tem_del_length>0:
        Del_length.append(tem_del_length)
    if len(tem_num_string)>0:
        Best_ref_ins_list.append([tem_new_bone_p_o,tem_num_string,len(tem_num_string)])
    if len(tem_num_string_full)>0:
        Best_ref_ins_list_full.append([tem_new_bone_p_o,tem_num_string_full,len(tem_num_string_full)])
        #print(ll_per_id)
    ### end of the last deletion
    max_len_id = max([len(Best_ref_list),len(Best_ref_ins_list),len(Del_length),len(covered_common_del_id),len(Best_ref_ins_list_full)])
    if max_len_id>0:
        full_return_result = [[id_block,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1] for x in range(max_len_id)]
        full_return_result[0][1] = tem_c
        full_return_result[0][2] = tem_mis
        if len(Best_ref_ins_list)>0:
            for ins_per_id in range(len(Best_ref_ins_list)):
                full_return_result[ins_per_id][3:6] = Best_ref_ins_list[ins_per_id]
        
        if len(Best_ref_ins_list_full)>0:
            for ins_per_id_full in range(len(Best_ref_ins_list_full)):
                full_return_result[ins_per_id_full][6:9] = Best_ref_ins_list_full[ins_per_id_full]
        
        
        if len(Del_length)>0:
            for del_per_id in range(len(Del_length)):
                full_return_result[del_per_id][9] = Del_length[del_per_id]
                
        if len(Best_ref_list)>0:
            for bone_id_check in range(len(Best_ref_list)):
                full_return_result[bone_id_check][10:12] = Best_ref_list[bone_id_check]
        
        if len(covered_common_del_id)>0:
            for flag_common_ins_id in range(len(covered_common_del_id)):
                full_return_result[flag_common_ins_id][12] = covered_common_del_id[flag_common_ins_id]
    else:
        full_return_result = [[id_block,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]]
    return(full_return_result)
    



if __name__ == '__main__':
    time_save_list = []
    replication_times = 1
    print(sys.argv)
    t1=time.time()
    
    
    
    
    proposed_k = int(sys.argv[1])
    num_cores = int(sys.argv[2])
    main_url_input_bam = sys.argv[3]
    main_url_input_ref = sys.argv[4]
   
    main_url_save =  sys.argv[5]
    start_region = int(sys.argv[6])
    end_region = int(sys.argv[6])
    if not os.path.exists(main_url_save):
        os.makedirs(main_url_save)
    
    
    
    sum_prob_vec = np.zeros(1)
    
    corretion_eps_from_zero = 0.0000000001
    noninformative_constant = 0.0000000001
    portions_var_fail = 0.0001
    path_exon_bam=main_url_input_bam
    inbam = pysam.AlignmentFile(path_exon_bam, "rb" )
   
    
    Result_table = pd.DataFrame(columns=['Num', 'seq', 'start_mapping', 'end_mapping',
                                 'len_seq', 'Cigar_information'])
    Result_table['Cigar_information'] = Result_table['Cigar_information'].astype(
        'object')
    
    #Final_result_interested = pd.DataFrame(np.zeros([7]))
    
    
    collection_of_mutilple_al = []
    tem_res_collection = []
    #label_table = pd.read_table(main_url_input_label,sep=',')
    #collection_reads = list(label_table.iloc[:,1])
    
    ### it can be done in a parrllel way
    
    #for s_name_num in range(max_num_reads):
    
    s_name_num = 0
    collection_SV_var_read_del = []
    collection_SV_var_read_ins = []
    collection_SV_var_read_clip = []
    for read_xun in inbam.fetch():
        #print(read_xun)
        ## flag_infor
        #read_flag = read_xun.flag
        #if not read_xun.is_secondary:
        if not read_xun.is_secondary and not read_xun.is_unmapped and not read_xun.is_supplementary:
            read_name = read_xun.query_name
            read_start = read_xun.reference_start - start_region
            read_seq = read_xun.query_sequence
            read_cigar = read_xun.cigar
            if len(read_cigar)>0:
                bone_base_count = 0
                pre_clip_len = 0
                after_clip_len = 0
                seq_base_count = 0
                for cigar_ in read_cigar:
                    if cigar_[0] in [0,2,7,8]:
                        if cigar_[0] in [2]:
                            if  cigar_[1]>20:
                                collection_SV_var_read_del.append([s_name_num,read_start+bone_base_count,read_start+bone_base_count+cigar_[1]])
                        else:
                            seq_base_count = seq_base_count + cigar_[1]
                        bone_base_count = bone_base_count + cigar_[1]
                    elif cigar_[0] in [1] and cigar_[1]>20:
                        collection_SV_var_read_ins.append([s_name_num,read_start+bone_base_count-1,cigar_[1],read_seq[seq_base_count:seq_base_count+cigar_[1]]])
                        seq_base_count = seq_base_count + cigar_[1]
                    elif cigar_[0] in [4]:
                        seq_base_count = seq_base_count + cigar_[1]
                read_end = read_start+bone_base_count-1            
                if read_cigar[0][0] in [4]:
                    pre_clip_len = int(read_cigar[0][1])
                    if pre_clip_len>20:
                        collection_SV_var_read_clip.append([s_name_num,read_start-1,pre_clip_len])
                if read_cigar[-1][0] in [4]:
                    after_clip_len = int(read_cigar[-1][1])
                    if after_clip_len>20:
                        collection_SV_var_read_clip.append([s_name_num,read_end,after_clip_len])
                tem_res = bam_2_msa_ins(read_start,read_seq,read_cigar)
                tem_res_collection.append([read_name,read_start,tem_res[0],read_seq,read_cigar,tem_res[1],tem_res[2]])
                Result_table.loc[s_name_num,:] = [s_name_num,read_seq,read_start,read_start+bone_base_count-1,bone_base_count-1,0]
                Result_table.at[s_name_num,'Cigar_information'] = read_cigar
                s_name_num = s_name_num + 1
                print(s_name_num)
            #else:
            #    tem_res_collection.append([])
            #    collection_of_unmapped_index.append(s_name_num)
            #    collection_of_unmapped.append('reads'+str(s_name_num))
            #    Result_table.iloc[s_name_num,:] = [s_name_num,read_seq,-1,-1,-1,-1,0]
        #else:
        #    collection_of_mutilple_al.append([read_xun.query_name,read_xun.reference_start])
        #### multiple ins
    
    ### here we want to test the range of the overlapping region
    other_informative_sites = []
    for del_content in collection_SV_var_read_del:
        other_informative_sites = other_informative_sites + list(range(del_content[1],del_content[2]+1))
    t2=time.time()
    time_save_list.append(['import data and transformation',t2-t1])
    
            
    whole_insertion = {}
    pos_read_index = {}
    for i in range(len(tem_res_collection)):
        if len(tem_res_collection[i])>0:
            tem_insert_dic = tem_res_collection[i][5]
            len_insertion = len(tem_insert_dic)
            if len_insertion>0:
                for keys_ in tem_insert_dic:
                    if keys_ in whole_insertion:
                        whole_insertion[keys_].append(tem_insert_dic[keys_])
                    else:
                        whole_insertion[keys_] = [tem_insert_dic[keys_]]   
                    if keys_ in pos_read_index:
                        pos_read_index[keys_].append(tem_res_collection[i][0])
                    else:
                        pos_read_index[keys_] = [tem_res_collection[i][0]]
    
    #### select remaining multiple sequce alignment
    dic_multiple_ali_ins = {}
    ins_inserted_dic = {}
    ins_inserted_dic_df = []
    tem_sum_trans_ = 0
    whole_insertion_keys_list = sorted(list(whole_insertion.keys()))
    pool = mp.Pool(num_cores)
    dic_multiple_ali_ins_list = pool.map(multiple_seq_alignment2_1,whole_insertion_keys_list)
    pool.close()
    pool.join()
    
    
    for key_id in range(len(whole_insertion_keys_list)):
        dic_multiple_ali_ins[whole_insertion_keys_list[key_id]] = dic_multiple_ali_ins_list[key_id]
        max_len_ =  len(dic_multiple_ali_ins_list[key_id][0][1])
        for sssss_index in range(len(dic_multiple_ali_ins_list[key_id])):
            if len(dic_multiple_ali_ins_list[key_id][sssss_index][1])>max_len_:
                max_len_ = len(dic_multiple_ali_ins_list[key_id][1])
        ins_inserted_dic[whole_insertion_keys_list[key_id]] = max_len_
        tem_sum_trans_ = tem_sum_trans_ + max_len_ 
        ins_inserted_dic_df.append([whole_insertion_keys_list[key_id],max_len_,tem_sum_trans_+whole_insertion_keys_list[key_id]])
    
    ### modified to here we need to filter the unmapped reads from the Result_table
    ins_inserted_dic_df = pd.DataFrame(ins_inserted_dic_df) 
    Result_table_cp = copy.deepcopy(Result_table)
    #Result_table_cp = Result_table_cp.drop(collection_of_unmapped_index)
    Result_table_cp['updated_index'] = range(len(Result_table_cp))
    
    
    ### reference_part
    #post_dir_reference_fasta = str(id_num)+'_0_reads_reference.fasta'
    fasta_file_path = main_url_input_ref #main_url_save + post_dir_reference_fasta
   
    fasta_file = pysam.FastaFile(fasta_file_path)
    seqname = fasta_file.references[0]
    BestRefSeq = (fasta_file.fetch(seqname,start_region,end_region+1))
    n_c = len(BestRefSeq)#stabndard bone
    ref_code_standard = string_2_code(BestRefSeq)
    
    ### refine the reference
    refine_best_ref  = refine_the_fasta(ins_inserted_dic,BestRefSeq)
    
    final_ref_standard = copy.deepcopy(refine_best_ref)
    # finalized_ref_coll.append(final_ref)
    final_ref_code_standard = string_2_code(final_ref_standard)
    n_c_refine = len(final_ref_code_standard)
    
    pool = mp.Pool(num_cores)
    out_put_df_list = pool.map(adding_ins_to_bone_per,list(range(len(Result_table_cp))))
    pool.close()
    pool.join()
    ### finish the change of the sequence
    pool = mp.Pool(num_cores)
    final_full_result = pool.map(coding_mismatch_all_insertion_error2_per,list(range(len(Result_table_cp))))
    pool.close()
    pool.join()
    t3 = time.time()
    time_save_list.append(['MSA time for insertion', t3-t2])
    
    max_num_reads_refine = len(out_put_df_list)
    
    full_matrix = np.zeros([n_c_refine,10],dtype=int)
    for ele_id in range(len(final_full_result)):
        start_tem = final_full_result[ele_id][0]
        end_tem = final_full_result[ele_id][1]
        tem_matrix_ = final_full_result[ele_id][2]
        full_matrix[start_tem:end_tem+1,1:7] = full_matrix[start_tem:end_tem+1,1:7] + tem_matrix_[:,1:]

    

    full_matrix[:,0] = list(range(n_c_refine))
    #full_matrix_full =  np.zeros([len(loc_bam_list),7])
    #full_matrix_full[:,:4] =full_matrix[:,:4]
    full_matrix[:,7] = np.sum(full_matrix[:,1:7],axis=1)
    full_matrix[:,8] = np.max(full_matrix[:,1:7],axis=1)
    full_matrix[:,9] = full_matrix[:,7] - full_matrix[:,8]
    full_matrix_sub = full_matrix[full_matrix[:,9]>0]
    
    ### normal_constant_part
    infer_var_site = sorted(list(full_matrix_sub[:,0]))
    K = proposed_k
    saving_para = np.zeros(K*8+6)
    if len(infer_var_site)>0:## we have the sites for inference 
        
        
        
        constant_loc_coll = sorted(list(set(list(range(n_c_refine)))-set(infer_var_site)))
        num_filter_noinfor_base = len(constant_loc_coll)
        
        full_num_base_c_seq = 0
        Backgroud_evo_c = 0
        Backgroud_evo_m = 0
        Backgroud_evo_cor_suspicious = []
        Backgroud_evo_mis_suspicious = {}
        Backgroud_evo_ins_suspicious = {}
        Backgroud_evo_del_suspicious = []
        Common_infer_site = {}
        ### 0 base / ins and del only record the pos of the ins and del
        min_num_per = 2
        min_num_full = min_num_per*proposed_k
        full_matrix_sub_df = pd.DataFrame(full_matrix_sub)
        for ll_sum_base in range(num_filter_noinfor_base):
            present_loc_tem = constant_loc_coll[ll_sum_base]
            full_num_base_c_seq = full_num_base_c_seq + full_matrix[present_loc_tem,7]
            standard_base = final_ref_code_standard[present_loc_tem]
            if full_matrix[present_loc_tem,7] - full_matrix[present_loc_tem,1+standard_base]<1:
                ### bone match
                if standard_base<5:
                    Backgroud_evo_c = Backgroud_evo_c + 1
                    Backgroud_evo_cor_suspicious.append(present_loc_tem)
                Common_infer_site[present_loc_tem] = standard_base
                #full_num_base_c_seq = full_num_base_c_seq + full_matrix[ll_sum_base,7]
            else:
                if full_matrix[present_loc_tem,8]>=min_num_full:
                    infer_base = list(full_matrix[present_loc_tem,1:7]).index(full_matrix[present_loc_tem,7])
                    if standard_base <4:### this is the bone
                        if infer_base>3:
                            Backgroud_evo_del_suspicious.append(present_loc_tem)
                            Common_infer_site[present_loc_tem] = 4
                        else:
                            Backgroud_evo_m = Backgroud_evo_m + 1
                            Backgroud_evo_mis_suspicious[present_loc_tem] = infer_base
                            Common_infer_site[present_loc_tem] = infer_base 
                            #full_num_base_c_seq = full_num_base_c_seq + full_matrix[ll_sum_base,7]
                    else:### this is the gap between the bones/ins
                        Backgroud_evo_ins_suspicious[present_loc_tem] = infer_base
                        Common_infer_site[present_loc_tem] = infer_base
                        #full_num_base_c_seq = full_num_base_c_seq + full_matrix[ll_sum_base,7]
                else:
                    infer_var_site.append(present_loc_tem)
                    full_matrix_sub_df.loc[len(full_matrix_sub_df),:] = full_matrix[present_loc_tem,:]
            #full_num_base_c_seq = full_num_base_c_seq + full_matrix[ll_sum_base,7]
        Independet_site_col_no_alg = Backgroud_evo_cor_suspicious + list(Backgroud_evo_mis_suspicious.keys())
        Indel_site_col_no_alg = Backgroud_evo_del_suspicious + list(Backgroud_evo_ins_suspicious.keys())
        full_non_informative_site_collection = sorted(Independet_site_col_no_alg + Indel_site_col_no_alg)
        Independet_site_col_no_alg_array = np.array(Independet_site_col_no_alg)
        ### we only filter the mis common and normal commoon
        
        infer_var_site = sorted(list(infer_var_site))
        num_remain_infor_base = len(infer_var_site)
        ## contant ref
        contant_ref_code_var = np.zeros(num_remain_infor_base,dtype=int)
        for infer_var_id  in range(num_remain_infor_base):
            contant_ref_code_var[infer_var_id] = final_ref_code_standard[int(infer_var_site[infer_var_id])]
        ##
        
        full_matrix_sub_df = full_matrix_sub_df.sort_values(by=[0])
        full_matrix_sub = np.array(full_matrix_sub_df)
        #if len(infer_var_site)>1000:
        pool = mp.Pool(num_cores)
        infer_var_site_o = pool.map(reverse_trans2_per,infer_var_site)
        pool.close()
        pool.join()
        #else:
        #    infer_var_site_o = reverse_trans2(ins_inserted_dic_df, infer_var_site,0)
        
        full_infer_var_site_o_df = pd.DataFrame(infer_var_site_o)
        
        only_bone_note_df = full_infer_var_site_o_df[full_infer_var_site_o_df.iloc[:,2]==-1]
        only_bone_note_df_keys = sorted(list(only_bone_note_df.iloc[:,1]))
        other_not_bone_df = full_infer_var_site_o_df[full_infer_var_site_o_df.iloc[:,3]>-1]
        other_not_bone_df_keys = sorted(list(set(only_bone_note_df.iloc[:,0])))
        full_id_bone_list = sorted(only_bone_note_df_keys+other_not_bone_df_keys)
        infer_var_ori_bone = list(full_infer_var_site_o_df.iloc[:,1])
        infer_var_ori_state = list(full_infer_var_site_o_df.iloc[:,2])
        
        
        
        ### common_del_block
        common_del_block = []
        common_del_block_sites_collection = []
        
        if len(Backgroud_evo_del_suspicious)>0:
            present_start_del = Backgroud_evo_del_suspicious[0]
            present_start_del_o = reverse_trans2_per(present_start_del)[1]
            for d_del in range(1,len(Backgroud_evo_del_suspicious)):
                if Backgroud_evo_del_suspicious[d_del]-Backgroud_evo_del_suspicious[d_del-1]>1:
                    
                    common_del_block.append([present_start_del,Backgroud_evo_del_suspicious[d_del-1],\
                                             present_start_del_o,present_start_del_o+Backgroud_evo_del_suspicious[d_del-1]-present_start_del])
                    common_del_block_sites_collection = common_del_block_sites_collection + list(range(common_del_block[-1][2],common_del_block[-1][3]+1))
                    present_start_del = Backgroud_evo_del_suspicious[d_del]
                    present_start_del_o = reverse_trans2_per(present_start_del)[1]
            common_del_block.append([present_start_del,Backgroud_evo_del_suspicious[-1],\
                                     present_start_del_o,present_start_del_o+Backgroud_evo_del_suspicious[-1]-present_start_del])
            common_del_block_sites_collection= common_del_block_sites_collection + list(range(common_del_block[-1][2],common_del_block[-1][3]+1))
        common_del_block_df = pd.DataFrame(common_del_block)
        common_del_block_sites_collection = sorted(common_del_block_sites_collection)
        common_del_block_sites_collection_array = np.array(common_del_block_sites_collection)
        common_del_block_df[4] = 0 ### denote whether it is used
        
        ### common_ins_block
        common_ins_block = []
        common_ins_block_dic = {}
        common_ins_block_dic_pure_num = {}
        keys_ins_common = []
       
        if len(Backgroud_evo_ins_suspicious)>0:
           keys_ins_common = sorted(list(Backgroud_evo_ins_suspicious.keys()))##abs location
           present_start_ins = keys_ins_common[0]
           bones_dic_tem_ins_s = reverse_trans(ins_inserted_dic, [present_start_ins])
           bones_dic_tem_bone = bones_dic_tem_ins_s[present_start_ins][0]
           sub_ins_start = bones_dic_tem_ins_s[present_start_ins][1]
           for i_ins in range(1,len(keys_ins_common)):
               
               if keys_ins_common[i_ins]-keys_ins_common[i_ins-1]>1:
                   bones_dic_tem_ins_e = reverse_trans(ins_inserted_dic, [keys_ins_common[i_ins-1]])
                   sub_ins_end =  bones_dic_tem_ins_e[keys_ins_common[i_ins-1]][1]
                   tem_string_ins = ''
                   tem_string_ins_num = ''
                   for ll_ins_per in range(present_start_ins,keys_ins_common[i_ins-1]+1):
                       tem_string_ins = tem_string_ins + Standardbase[Backgroud_evo_ins_suspicious[ll_ins_per]]
                       tem_string_ins_num = tem_string_ins_num + str(Backgroud_evo_ins_suspicious[ll_ins_per])
                   if bones_dic_tem_bone in  common_ins_block_dic:
                       common_ins_block_dic[bones_dic_tem_bone].append( [sub_ins_start,sub_ins_end,tem_string_ins])
                       common_ins_block_dic_pure_num[bones_dic_tem_bone].append( [sub_ins_start,sub_ins_end,tem_string_ins_num])
                   else:
                       common_ins_block_dic[bones_dic_tem_bone] = [[sub_ins_start,sub_ins_end,tem_string_ins]]
                       common_ins_block_dic_pure_num[bones_dic_tem_bone] = [ [sub_ins_start,sub_ins_end,tem_string_ins_num]]
                   tem_string_ins = ''
                   tem_string_ins_num = ''
                   
                   common_ins_block.append([present_start_ins,keys_ins_common[i_ins-1]])
                   
                   ### update the ins pos
                   present_start_ins = keys_ins_common[i_ins]
                   bones_dic_tem_ins_s = reverse_trans(ins_inserted_dic, [present_start_ins])
                   bones_dic_tem_bone = bones_dic_tem_ins_s[present_start_ins][0]
                   sub_ins_start = bones_dic_tem_ins_s[present_start_ins][1]
               
           if (len(common_ins_block)>0 and keys_ins_common[-1] !=   common_ins_block[-1][1]) or (len(common_ins_block)<1): 
              present_end_ins = keys_ins_common[-1]
              bones_dic_tem_ins_s = reverse_trans(ins_inserted_dic, [present_end_ins])
              bones_dic_tem_bone = bones_dic_tem_ins_s[present_end_ins][0]
              sub_ins_end = bones_dic_tem_ins_s[present_end_ins][1]
              tem_string_ins = ''
              tem_string_ins_num = ''
              for ll_ins_per in range(present_start_ins,keys_ins_common[-1]+1):
                   tem_string_ins = tem_string_ins + Standardbase[Backgroud_evo_ins_suspicious[ll_ins_per]]
                   tem_string_ins_num = tem_string_ins_num + str(Backgroud_evo_ins_suspicious[ll_ins_per])
              if bones_dic_tem_bone in  common_ins_block_dic:
                   common_ins_block_dic[bones_dic_tem_bone].append( [sub_ins_start,sub_ins_end,tem_string_ins])
                   common_ins_block_dic_pure_num[bones_dic_tem_bone].append( [sub_ins_start,sub_ins_end,tem_string_ins_num])
              else:
                   common_ins_block_dic[bones_dic_tem_bone] = [[sub_ins_start,sub_ins_end,tem_string_ins]]
                   common_ins_block_dic_pure_num[bones_dic_tem_bone] = [ [sub_ins_start,sub_ins_end,tem_string_ins_num]]
                   
              common_ins_block.append([present_start_ins,keys_ins_common[-1]])
        keys_ins_common_original = list(common_ins_block_dic.keys())
        keys_ins_common_original_array = np.array(keys_ins_common_original)
        common_ins_block_df = pd.DataFrame(common_ins_block)
       
        
       
        
        
        ### sgementation of the block
        seg_collection_ini = []
        
        previous_del_flag = 0
        last_bone = reverse_trans(ins_inserted_dic, [infer_var_site[0]])[infer_var_site[0]][0]#reverse_trans2_per(infer_var_site[0])[1]
        tem_start_block_id = 0
        if full_matrix_sub[0,5]>0:
            previous_del_flag = 1
            #tem_start_block_id = 0
        #else:
        #    if last_bone not in other_not_bone_df_keys:
        #        seg_collection_ini.append([0,0])## block margin
        
        ll_seg_id =  1
        while  ll_seg_id <len(infer_var_site):
            present_bone =  reverse_trans(ins_inserted_dic, [infer_var_site[ll_seg_id]])[infer_var_site[ll_seg_id]][0]
            present_bone_state = reverse_trans(ins_inserted_dic, [infer_var_site[ll_seg_id]])[infer_var_site[ll_seg_id]][1]
            if previous_del_flag:
                if present_bone_state<0:# this means we can have the true bone
                    if len(common_del_block_df)>0:
                        tem_del_left_boundary = np.array(common_del_block_df.iloc[:,2])
                        tem_del_left_boundary_dis = tem_del_left_boundary - last_bone
                        
                        
                        if 1 in tem_del_left_boundary_dis:## neighboured common_del
                            id_common_del = list(tem_del_left_boundary_dis).index(1)
                            end_possible_del = common_del_block_df.iloc[id_common_del,3]
                            if full_matrix_sub[ll_seg_id,5]<1:
                                tem_end_block_id = ll_seg_id - 1
                                seg_collection_ini.append([tem_start_block_id,tem_end_block_id,2])
                                previous_del_flag = 0
                                tem_start_block_id = ll_seg_id
                            else:
                                if present_bone - end_possible_del<2:
                                    previous_del_flag = 1
                                    #present_bone = 
                                    #last_bone = copy.deepcopy(present_bone)
                                else:
                                    tem_end_block_id = ll_seg_id - 1
                                    seg_collection_ini.append([tem_start_block_id,tem_end_block_id,2])
                                    tem_start_block_id = ll_seg_id
                                    
                                    previous_del_flag = 1
                        else:
                            
                            tem_end_block_id = ll_seg_id - 1
                            seg_collection_ini.append([tem_start_block_id,tem_end_block_id,2])
                            tem_start_block_id = ll_seg_id
                            if full_matrix_sub[ll_seg_id,5]>0:
                                previous_del_flag = 0
                                #tem_end_block_id = ll_seg_id
                                #seg_collection_ini.append([tem_start_block_id,tem_end_block_id])
                               
                            else:
                                previous_del_flag = 1
                    else:
                        tem_end_block_id = ll_seg_id - 1
                        seg_collection_ini.append([tem_start_block_id,tem_end_block_id,2])
                        tem_start_block_id = ll_seg_id
                        if full_matrix_sub[ll_seg_id,5]>0:
                            previous_del_flag = 0
                            #tem_end_block_id = ll_seg_id
                            #seg_collection_ini.append([tem_start_block_id,tem_end_block_id])
                           
                        else:
                            previous_del_flag = 1
                    
                else:
                    tem_end_block_id = ll_seg_id - 1
                    seg_collection_ini.append([tem_start_block_id,tem_end_block_id,2])
                    tem_start_block_id = ll_seg_id
                    previous_del_flag = 0
            
            else:
                if present_bone - last_bone>0:
                    tem_end_block_id = ll_seg_id - 1
                    previous_bone_state =  reverse_trans2_per(infer_var_site[ll_seg_id-1])[2]
                    if previous_bone_state<0:
                        seg_collection_ini.append([tem_start_block_id,tem_end_block_id,0])
                    else:
                        seg_collection_ini.append([tem_start_block_id,tem_end_block_id,1])
                    tem_start_block_id = ll_seg_id
                    if full_matrix_sub[ll_seg_id,5]>0:
                        previous_del_flag = 1
                        
                    else:
                        previous_del_flag = 0
                else:
                    
                    if full_matrix_sub[ll_seg_id,5]>0:
                        previous_del_flag = 1
                    else:
                        previous_del_flag = 0
                    
            last_bone = copy.deepcopy(present_bone)
            ll_seg_id = ll_seg_id + 1
            print(ll_seg_id)
                    
        if tem_end_block_id < tem_start_block_id:
            tem_end_block_id = len(infer_var_site) - 1
            present_bone_state = reverse_trans2_per(infer_var_site[ll_seg_id-1])[2]
            if previous_del_flag:
                seg_collection_ini.append([tem_start_block_id,tem_end_block_id,2])
            else:
                if present_bone_state<0:
                   seg_collection_ini.append([tem_start_block_id,tem_end_block_id,0]) 
                else:
                    seg_collection_ini.append([tem_start_block_id,tem_end_block_id,1]) 
           
        seg_collection_ini_df = pd.DataFrame(seg_collection_ini)
        seg_collection_ini_df_no_del_bone = seg_collection_ini_df[seg_collection_ini_df.iloc[:,-1]==0]
        seg_collection_ini_df_no_del_bone_id = sorted(list(seg_collection_ini_df_no_del_bone.iloc[:,0]))
        
        seg_collection_ini_df_indel_part = seg_collection_ini_df[seg_collection_ini_df.iloc[:,-1]>0]
        seg_collection_ini_df_indel_part_cp = copy.deepcopy(seg_collection_ini_df_indel_part)
        seg_collection_ini_df_indel_part_cp[3] = seg_collection_ini_df_indel_part_cp.iloc[:,1] - seg_collection_ini_df_indel_part_cp.iloc[:,0]
        seg_collection_ini_df_indel_part_cp_sub = seg_collection_ini_df_indel_part_cp[seg_collection_ini_df_indel_part_cp.iloc[:,3]<1]
        seg_collection_ini_df_indel_part_cp_sub_link = seg_collection_ini_df_indel_part_cp[seg_collection_ini_df_indel_part_cp.iloc[:,3]>0]
        
        dep_spe_infor_del_col_main = []
        dep_spe_infor_ins_col_main = []
        for ll_dep_check_block_id in range(len(seg_collection_ini_df_indel_part_cp_sub_link)):
            ## left common_del check
            tem_block_start = seg_collection_ini_df_indel_part_cp_sub_link.iloc[ll_dep_check_block_id,0]
            tem_block_end = seg_collection_ini_df_indel_part_cp_sub_link.iloc[ll_dep_check_block_id,1]
            dep_spe_infor_del = []
            dep_spe_infor_ins = []
            for ll_dep_check_block_id_nei in range(tem_block_start,tem_block_end+1):
                present_var_tem =  infer_var_site[int(ll_dep_check_block_id_nei)]
                if len(common_del_block_df) > 0:
                    common_del_block_df_right_abs = np.array(common_del_block_df.iloc[:,1])
                    tem_left_dis_arr = present_var_tem - common_del_block_df_right_abs
                    dep_flag = 0
                    #dep_spe_infor_del = []
                    #dep_spe_infor_ins = []
                    if 1 in tem_left_dis_arr:
                        #collection_true_dependent_common_indel.append(ll_indep_check_id)
                        #dep_flag = 1
                        dep_spe_infor_del.append([-1,ll_dep_check_block_id_nei,list(tem_left_dis_arr).index(1)])
                    
                    ## right common_del check
                    common_del_block_df_left_abs = np.array(common_del_block_df.iloc[:,0])
                    tem_right_dis_arr = common_del_block_df_left_abs-present_var_tem
                    
                    if 1 in tem_right_dis_arr:
                        #collection_true_dependent_common_indel.append(ll_indep_check_id)
                        #dep_flag = 1
                        dep_spe_infor_del.append([1,ll_dep_check_block_id_nei,list(tem_right_dis_arr).index(1)])
                    
                ## common_ins_prob
                original_flag_ins_content = reverse_trans2_per(present_var_tem)
                original_flag_ins_pos = original_flag_ins_content[1]
                original_flag_ins_pos_sub_pos = original_flag_ins_content[2]
                if original_flag_ins_pos in common_ins_block_dic:
                    if original_flag_ins_pos not in dep_spe_infor_ins:
                        dep_spe_infor_ins.append(original_flag_ins_pos)
                        
               
            dep_spe_infor_del_col_main.append(dep_spe_infor_del)
            dep_spe_infor_ins_col_main.append(dep_spe_infor_ins)
            
        collection_true_independent_indel = [] ## this means they do not reply on the neighbour just single indel not with common ins and del
        #collection_true_dependent_indel = []
        collection_true_dependent_common_indel = []
        dep_spe_infor_ins_col_mk = []
        dep_spe_infor_del_col_mk = []
        seg_collection_ini_df_indel_part_cp_sub_link_full = copy.deepcopy(seg_collection_ini_df_indel_part_cp_sub_link)
        seg_collection_ini_df_indel_part_cp_sub_link_full = seg_collection_ini_df_indel_part_cp_sub_link_full.reset_index(drop=True)
        for ll_indep_check_id in range(len(seg_collection_ini_df_indel_part_cp_sub)):
            ## left common_del check
            present_var_tem =  infer_var_site[int(seg_collection_ini_df_indel_part_cp_sub.iloc[ll_indep_check_id,1])]
            dep_flag = 0
            dep_spe_infor_del = []
            dep_spe_infor_ins = []
            if len(common_del_block_df)>0:
                common_del_block_df_right_abs = np.array(common_del_block_df.iloc[:,1])
                tem_left_dis_arr = present_var_tem - common_del_block_df_right_abs
               
                if 1 in tem_left_dis_arr:
                    #collection_true_dependent_common_indel.append(ll_indep_check_id)
                    dep_flag = 1
                    dep_spe_infor_del.append([-1,seg_collection_ini_df_indel_part_cp_sub.iloc[ll_indep_check_id,0],list(tem_left_dis_arr).index(1)])
                
                ## right common_del check
                common_del_block_df_left_abs = np.array(common_del_block_df.iloc[:,0])
                tem_right_dis_arr = common_del_block_df_left_abs-present_var_tem
                
                if 1 in tem_right_dis_arr:
                    #collection_true_dependent_common_indel.append(ll_indep_check_id)
                    dep_flag = 1
                    dep_spe_infor_del.append([1,seg_collection_ini_df_indel_part_cp_sub.iloc[ll_indep_check_id,0],list(tem_right_dis_arr).index(1)])
                
            ## common_ins_prob
            original_flag_ins_content = reverse_trans2_per(present_var_tem)
            original_flag_ins_pos = original_flag_ins_content[1]
            original_flag_ins_pos_sub_pos = original_flag_ins_content[2]
            if original_flag_ins_pos in ins_inserted_dic:
                full_ins_num = ins_inserted_dic[original_flag_ins_pos]
                if full_ins_num>1:
                    dep_flag = 1
                    #collection_true_dependent_common_indel.append(ll_indep_check_id)
                    if original_flag_ins_pos in common_ins_block_dic and original_flag_ins_pos not in dep_spe_infor_ins:
                        dep_spe_infor_ins.append(original_flag_ins_pos)
                    
                    
            if not dep_flag:
                collection_true_independent_indel.append(seg_collection_ini_df_indel_part_cp_sub.iloc[ll_indep_check_id,0])
            else:
                seg_collection_ini_df_indel_part_cp_sub_link_full.loc[len(seg_collection_ini_df_indel_part_cp_sub_link_full),:] = seg_collection_ini_df_indel_part_cp_sub.iloc[ll_indep_check_id,:]
                dep_spe_infor_del_col_mk.append(dep_spe_infor_del)
                dep_spe_infor_ins_col_mk.append(dep_spe_infor_ins)
        seg_collection_ini_df_dep_id = sorted(collection_true_independent_indel + seg_collection_ini_df_no_del_bone_id)
    
                
        dep_spe_infor_del_col_full = dep_spe_infor_del_col_main + dep_spe_infor_del_col_mk
        dep_spe_infor_ins_col_full = dep_spe_infor_ins_col_main + dep_spe_infor_ins_col_mk
        #seg_collection_ini_df_ins_bone = seg_collection_ini_df[seg_collection_ini_df.iloc[:,-1]==1]
        t4 = time.time()
        time_save_list.append(['Segmentation of dependent and independent', t4-t3])
        ### 
        #collectionerr = []
        #for ll in range(len(Result_infer_)):
        #    if Result_infer_[ll][0] != infer_var_site[ll]:
        #        collectionerr.append(ll)
                
        
        
        
        read_full_bases_filter = []
        other_remaining_num_matrix = np.zeros([len(infer_var_site),7],dtype=int)
        for ele_id in range(len(final_full_result)):
            #ele_id = 2184
            tem_start = final_full_result[ele_id][0]
            tem_end = final_full_result[ele_id][1]
            left_usefull_col_tem = sorted(list(set(list(range(tem_start,tem_end+1))) - set(full_non_informative_site_collection)))
            tem_original_matrix = final_full_result[ele_id][2]
            if len(left_usefull_col_tem)>0:
                tem_saving_matrix_filter = np.zeros([len(left_usefull_col_tem),7],dtype=int)
                for left_id in range(len(left_usefull_col_tem)):
                    relative_pos_tem_original = left_usefull_col_tem[left_id] - tem_start
                    tem_saving_matrix_filter[left_id,:] = tem_original_matrix[relative_pos_tem_original,:7]
                read_full_bases_filter.append(tem_saving_matrix_filter)
            else:
                read_full_bases_filter.append([])
        
        
        #full_num_base_c_evo = len(constant_loc_coll)
        
        pool = mp.Pool(num_cores)
        infer_result_snv = pool.starmap(core_selection,[(full_matrix_sub[re_id,1:7],final_ref_code_standard[int(full_matrix_sub[re_id,0])],\
    
                                           int(full_matrix_sub[re_id,0])) for re_id in range(len(full_matrix_sub))])
    
        pool.close()
        pool.join()
    
        num_id_vec_re = np.array(infer_result_snv)
        num_id_vec_re_sub = list(num_id_vec_re[num_id_vec_re>-1])
        #num_id_vec_re_sub_df = pd.DataFrame(num_id_vec_re_sub)
        
        
        # simulation of making up the missing value
        combination_list2 = pd.DataFrame(
            it.combinations(list(range(max_num_reads_refine)), 2))
        gamma_r = 0.9
        theta_r = 0.9
        labels_initial0 = np.zeros(max_num_reads_refine)
        #labels_initial0 =np.array( Result_table.iloc[:,-3])
        labels_initial0_matrix = change_label_matrix(labels_initial0, 1)
        #proba_list = labels_initial0_matrix
        
        # labels_initial0
        # analysis the error
        #reads_flag = pd.DataFrame(np.zeros([max_num_reads_refine,2]))
        #reads_flag.iloc[:,0] = range(max_num_reads_refine)
        #reads_flag.iloc[:,1] = Result_table.iloc[:,-3]
        
        standard_ref_code = final_ref_code_standard
        
        
        
        coolection_snp1 = copy.deepcopy(num_id_vec_re_sub)
        t5 = time.time()
        time_save_list.append(['Collection of the informative sites', t5-t4])
        transed_var_sites = trans_(ins_inserted_dic,other_informative_sites)
        other_informative_sites_trans = list(transed_var_sites.values())
        for ins_content in collection_SV_var_read_ins:
            pre_ins_bone = trans_(ins_inserted_dic,[ins_content[1]])[ins_content[1]]
            other_informative_sites_trans = other_informative_sites_trans + list(range(pre_ins_bone+1,pre_ins_bone+ins_inserted_dic[ins_content[1]]))
        print(coolection_snp1)
        print(other_informative_sites_trans)
        coolection_snp1 = list(set((coolection_snp1 + other_informative_sites_trans)))
        if  len(coolection_snp1)>0:        
            #alter_one2 = trans_(ins_inserted_dic, coolection_snp1)
            mutation_list = copy.deepcopy(coolection_snp1)
        
                
            with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
         
             
                 Jaccard_matrix_df3=pool.starmap(Jaccard_pair9,[(combination_list2.iloc[index,0],combination_list2.iloc[index,1])\
                                            for index in range(len(combination_list2))])
                 
            
        else:
            print('warning, no sig snp')
            with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                
                Jaccard_matrix_df3=pool.starmap(Jaccard_pair8,[(combination_list2.iloc[index,0],combination_list2.iloc[index,1]) for index in range(len(combination_list2))])
                
                
            #final_inferred=pool.starmap(sampling2,[(pp_random[index,0],pp_random[index,1],combination_list2.iloc[int(index/num_chains),0],combination_list2.iloc[int(index/num_chains),1],out_put_df,collection_of_sites,final_ref_code) for index in range((len(combination_list2)*num_chains))])
            pool.close()
            pool.join()
            
        #Jaccard_matri
        Jaccard_matrix_df3 = pd.DataFrame(Jaccard_matrix_df3)
        #full_T_eq_coll = np.zeros([max_num_reads_refine, max_num_reads_refine],dtype=int)
                          
        #full_T_not_eq_coll = [np.zeros([max_num_reads_refine,max_num_reads_refine]) for xx in range(num_k-1)]
        
        #full_inferred_res = []
        #S_matrix = np.zeros([max_num_reads_refine, num_k-1])
        #same_count_col = []
        #not_same_count_col = []
        #full_label = []
        
        #for K in range(1, num_k+1):
        #K = 3
        # same graph perform different community detection
        # Here we need the carefull splitting for the compoments
        
        
        
        Label_collection = [np.zeros(max_num_reads_refine,dtype=int) for x in range(4)]
        
        X_Jaccard_full = convert_sparse_normal(Jaccard_matrix_df3, max_num_reads_refine)
        graph_nx_o = nx.from_numpy_array(X_Jaccard_full)
        connected_components = list(nx.connected_components(graph_nx_o))
        random_index_read_index = []
        random_method_indicator = np.zeros([4,len(connected_components)],dtype=int)
        for connected_components_sub_id in range(len(connected_components)):
            connected_components_sub_tem = sorted(list(connected_components[connected_components_sub_id]))
            if len(connected_components_sub_tem)<=K:
                random_index_read_index = random_index_read_index + connected_components_sub_tem
            else:
                failing_split_flag = np.zeros(4,dtype = int)
                tem_num_sub_compoments = len(connected_components_sub_tem)
                X_Jaccard = X_Jaccard_full[connected_components_sub_tem]
                X_Jaccard = X_Jaccard[:,connected_components_sub_tem]
                X_Jaccard_dis = np.ones([tem_num_sub_compoments, tem_num_sub_compoments]) - X_Jaccard-np.diag(np.ones(tem_num_sub_compoments))
                X_Jaccard_dis[X_Jaccard_dis < 0] = 0
                graph_nx_o_sub = nx.from_numpy_array(X_Jaccard)
               
                # 1. spectural clustering
                try:
                    clustering = SpectralClustering(
                        n_clusters=proposed_k, assign_labels='discretize', affinity='precomputed_nearest_neighbors').fit(X_Jaccard_dis)
                    # K=3
                    
                    labels_initial1 = clustering.labels_
                    Label_collection[0][connected_components_sub_tem] = labels_initial1
                except:
                    tem_prob_array = np.random.multinomial(1, [1/proposed_k]*proposed_k, size=tem_num_sub_compoments)
                    for ll_tem_ in range(tem_num_sub_compoments):
                        spe_clu = list(tem_prob_array[ll_tem_,:]).index(1)
                        Label_collection[0][int(connected_components_sub_tem[ll_tem_])] = spe_clu
                    random_method_indicator[0,connected_components_sub_id] = 1
                
                        
        
               
                # 2. greedy_modularity_communities
                try:
                    communities2 = greedy_modularity_communities(
                        graph_nx_o_sub, weight='weight', best_n=proposed_k)
                    labels_initial2 = np.zeros(tem_num_sub_compoments)
                    for kk in range(1, K):
                        tem_mem_list = list(communities2[kk])
                        for list_index in tem_mem_list:
                            labels_initial2[int(list_index)] = kk
                    Label_collection[1][connected_components_sub_tem] = labels_initial2
                except:
                    tem_prob_array = np.random.multinomial(1, [1/proposed_k]*proposed_k, size=tem_num_sub_compoments)
                    for ll_tem_ in range(tem_num_sub_compoments):
                        spe_clu = list(tem_prob_array[ll_tem_,:]).index(1)
                        Label_collection[1][int(connected_components_sub_tem[ll_tem_])] = spe_clu
                    random_method_indicator[1,connected_components_sub_id] = 1
                
                
                # 3 Birth
                try:
                    labels_initial8 = np.zeros(tem_num_sub_compoments)
                    pca = PCA()
                    transform_8 = pca.fit_transform(X_Jaccard)  
                    model8 =Birch(
                    n_clusters=  proposed_k).fit(transform_8)
                    labels_initial8 =model8.labels_#usaged bu
                    Label_collection[2][connected_components_sub_tem] = labels_initial8
                except:
                   tem_prob_array = np.random.multinomial(1, [1/proposed_k]*proposed_k, size=tem_num_sub_compoments)
                   for ll_tem_ in range(tem_num_sub_compoments):
                       spe_clu = list(tem_prob_array[ll_tem_,:]).index(1)
                       Label_collection[2][int(connected_components_sub_tem[ll_tem_])] = spe_clu
                   random_method_indicator[2,connected_components_sub_id] = 1
                
                # hierarchical_clustering
                try:
                    labels_initial9 = np.zeros(tem_num_sub_compoments)
                    communities9 = hierarchical_clustering(X_Jaccard, n=proposed_k)
                   
                    for kk in range(1, len(communities9)):
                        tem_mem_list = list(communities9[kk])
                        for list_index in tem_mem_list:
                            labels_initial9[int(list_index)] = kk
                    Label_collection[3][connected_components_sub_tem] = labels_initial9
                except:
                    tem_prob_array = np.random.multinomial(1, [1/proposed_k]*proposed_k, size=tem_num_sub_compoments)
                    for ll_tem_ in range(tem_num_sub_compoments):
                        spe_clu = list(tem_prob_array[ll_tem_,:]).index(1)
                        Label_collection[3][int(connected_components_sub_tem[ll_tem_])] = spe_clu
                    random_method_indicator[3,connected_components_sub_id] = 1
        
        tem_prob_array2 = np.random.multinomial(1, [1/proposed_k]*proposed_k, size=len(random_index_read_index))
        for rand_index_read_ in random_index_read_index:
            
            for ll_tem_ in range(len(random_index_read_index)):
                spe_clu = list(tem_prob_array2[ll_tem_,:]).index(1)
                for kk_ in range(proposed_k):
                    Label_collection[kk_][int(random_index_read_index[ll_tem_])] = spe_clu
        
        #### randomness checking
        no_use_clu = []
        for i_me in range(len(Label_collection)):
            if sum(random_method_indicator[i_me,:])==len(connected_components):
                no_use_clu.append(i_me)
        if len(no_use_clu) == len(Label_collection):
            print('full random initial')
            rabdom_ini_flag = 1
        else:
            rabdom_ini_flag = 0
        
        t6 = time.time()
        time_save_list.append(['Initialization of the group', t6-t5])
        
        ### perfect
        #Label_collection.append(copy.fdeepcopy(Result_table.iloc[:, -1]))
        
        Collection_label = []
        sum_of_prob = []
        method_index_c = []
        Best_ref_collection = {}
        Best_ref_ins_collection = {}
        parameters_saving = {}
        
        for method_index in range(len(Label_collection)):
            print('method_index', method_index)
            if not rabdom_ini_flag:
                #if len(set(Label_collection[method_index])) == K and method_index not in no_use_clu:
                    #break
                    #method_index = 6
                    t_pre = time.time()
                    time_save_list.append(['Initialization of the group', t_pre-t6])
                    
                    reads_infer_clu = pd.DataFrame(np.zeros([max_num_reads_refine, 2],dtype=int))
                    reads_infer_clu.iloc[:, 0] = range(max_num_reads_refine)
                    reads_infer_clu.iloc[:, 1] = copy.deepcopy(Label_collection[method_index])
                    
                    last_probability = np.zeros([max_num_reads_refine, K])
                    for read_index in range(max_num_reads_refine):
                        last_probability[read_index, int(
                            Label_collection[method_index][read_index])] = 1
                    last_portions = np.ones(K)/K
                    present_clu_dis2 = [
                        np.zeros([num_remain_infor_base, 6],dtype = int) for x in range(K)]
                    #present_clu_dis2_pro = [np.zeros([n_c_refine,5]) for x in range(K)]
                    Best_ref = [{} for x in range(K)]
                    Best_ref_ins = [{} for x in range(K)]
                    Best_ref_ins_keys = [[] for x in range(K)]
                    Best_ref_ins_full = [{} for x in range(K)]
                    last_informative_sites_code = [ contant_ref_code_var for x in range(K)]
                    informative_sites_code = [ contant_ref_code_var for x in range(K)]
                    
                    portions = np.zeros(K)
                    Basic_ref = ref_code_standard#final_ref_code_standard
                    #con_ref_prob = np.zeros([max_num_reads_refine,K])
                    con_ref_prob = np.zeros(K)
                    con_ref_prob_full = np.zeros([K, 6])
                    gamma1 = 0.9*np.ones(K)  # correct
                    gamma2 = 0.05*np.ones(K)  # mis
                    gamma3 = 0.05*np.ones(K) # del
                    gamma4 = 0.05*np.ones(K) # ins
                    del_beta = 0.05*np.ones(K)  # del
                    del_p = 0.05*np.ones(K)
                    ins_beta= 0.05*np.ones(K)  # del
                    ins_p = 0.05*np.ones(K)
                    
                    p_mis_loc = np.random.rand(1)/10
                    pins1_loc = np.random.rand(1)/10
                    lamuna_ins = int(np.random.rand(1)*6)+1
                    
                    tao = 1
                    #tao2 = 1
                    tao_step = 0.97
                    #tao_step2 = 0.99
                    pdel1_loc = np.random.rand(1)/10
                    lamuna_del = int(np.random.rand(1)*6)+1
                    
            
                    pcorrect = np.random.rand(1)/10+0.9
                    theta_ = [pcorrect, p_mis_loc, pins1_loc,
                              pdel1_loc, lamuna_ins, lamuna_del]
                    evo_ins_len_dic = []### this is to save the dic of the length of the ins length and del length
                    evo_del_len_dic = []
                    for kk in range(K):
                        ele_list = list(reads_infer_clu[reads_infer_clu.iloc[:, 1] == kk].iloc[:, 0])
                        #ele_list = list(
                        #    reads_infer_cluf[reads_infer_cluf.iloc[:, 1] == kk].iloc[:, 0])
                        portions[kk] = len(ele_list)
                        #tem_eq_num = 0
                        #tem_neq_num = 0
                        for ele in ele_list:
                            #tem_ele_start = final_full_result[ele][0]
                            #tem_ele_end = final_full_result[ele][1]
                            tem_matrix_refine = read_full_bases_filter[ele]
                            if len(tem_matrix_refine)>0:
                                tem_start_refine = infer_var_site.index(tem_matrix_refine[0,0])
                                tem_end_refine = tem_start_refine + len(tem_matrix_refine)-1
                                present_clu_dis2[kk][tem_start_refine:tem_end_refine+1,:] = present_clu_dis2[kk][tem_start_refine:tem_end_refine+1,:] + \
                                    tem_matrix_refine[:,1:]
                        #max_ID = list(present_clu_dis2[kk][int(basss_index),:]).index(max(present_clu_dis2[kk][int(basss_index),:]))
                        tem_del_len_dic = {}
                        tem_ins_len_dic = {}
                        common_del_block_df_cp = copy.deepcopy(common_del_block_df)
                        keys_ins_common_original_cp = copy.deepcopy(keys_ins_common_original)
                        tem_touched_ins = []
                        if len(ele_list) > 0:
                            # 1.seperate the room for seperate zoom for block indel
            
                            remaining_independent = copy.deepcopy(infer_var_site)
                            with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
            
                            # Result_infer_ind=pool.starmap(signle_sites_infer_max,[(remaining_independent[index],present_clu_dis2[kk][int(remaining_independent[index]),:],standard_ref_code[int(remaining_independent[index])],\
                            # gamma[kk],theta_) for index in range(len(remaining_independent))])
            
                                Result_infer_ind = pool.starmap(pre_max_, [(index, present_clu_dis2[kk][(index), :]) for index in range(len(remaining_independent))])
                                #Result_infer_ind =  Result_infer_ind_iter.get()
                            #final_inferred=pool.starmap(sampling2,[(pp_random[index,0],pp_random[index,1],combination_list2.iloc[int(index/num_chains),0],combination_list2.iloc[int(index/num_chains),1],out_put_df,collection_of_sites,final_ref_code) for index in range((len(combination_list2)*num_chains))])
                            #pool.close()
                            #pool.join()
            
                            # making up full result
                            Result_infer_ = copy.deepcopy(Result_infer_ind)
                            Result_infer_df = pd.DataFrame(Result_infer_)
                            last_informative_sites_code[kk] = np.array(Result_infer_df.iloc[:,1])
                            
                            with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                        
                             # Result_infer_ind=pool.starmap(signle_sites_infer_max,[(remaining_independent[index],present_clu_dis2[kk][int(remaining_independent[index]),:],standard_ref_code[int(remaining_independent[index])],\
                             # gamma[kk],theta_) for index in range(len(remaining_independent))])
                        
                                 final_compare_ref_result = pool.map(comparision_between_con_and_ref, list(range(len(seg_collection_ini))))
                                 #Result_infer_ind =  Result_infer_ind_iter.get()
                             #final_inferred=pool.starmap(sampling2,[(pp_random[index,0],pp_random[index,1],combination_list2.iloc[int(index/num_chains),0],combination_list2.iloc[int(index/num_chains),1],out_put_df,collection_of_sites,final_ref_code) for index in range((len(combination_list2)*num_chains))])
                             #pool.close()
                             #pool.join()
                            full_final_compare_ref_result = []
                             
                            for ll_final_id in range(len(final_compare_ref_result)):
                                 full_final_compare_ref_result = full_final_compare_ref_result + final_compare_ref_result[ll_final_id]
                             
                            full_final_compare_ref_result_df = pd.DataFrame(full_final_compare_ref_result)
                            
                            
                            con_ref_prob_full[kk,0] = sum(full_final_compare_ref_result_df.iloc[:,1]) + Backgroud_evo_c
                            con_ref_prob_full[kk,1] = sum(full_final_compare_ref_result_df.iloc[:,2]) + Backgroud_evo_m
                            #sub_sub_result = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,-1]==1]
                            
                            ### ins part
                            full_final_compare_ref_result_df_ins = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,3]>-1]
                            for ll_ins_id in range(len(full_final_compare_ref_result_df_ins)):
                                tem_ins_loc = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,3]
                                tem_ins_content = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,4]
                                tem_ins_length = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,5]
                                Best_ref_ins[kk][tem_ins_loc] = tem_ins_content
                                if tem_ins_length in tem_ins_len_dic:
                                    tem_ins_len_dic[tem_ins_length] = tem_ins_len_dic[tem_ins_length] + 1
                                else:
                                    tem_ins_len_dic[tem_ins_length] = 1
                            
                            collection_id_ins_covered = list(Best_ref_ins[kk].keys())
                            remaining_ins_id = list(set(list(common_ins_block_dic.keys()))-set(collection_id_ins_covered))
                            
                            for ll_ins_id_mk in remaining_ins_id:### common ins part
                                full_ins_tem_list = common_ins_block_dic[ll_ins_id_mk]
                                #for sub_tem_ins_ele in full_ins_tem_list:
                                tem_ins_content = full_ins_tem_list[0][2]
                                tem_ins_length = len(tem_ins_content)
                                Best_ref_ins[kk][ll_ins_id_mk] = tem_ins_content
                                if tem_ins_length in tem_ins_len_dic:
                                    tem_ins_len_dic[tem_ins_length] = tem_ins_len_dic[tem_ins_length] + 1
                                else:
                                    tem_ins_len_dic[tem_ins_length] = 1
                            ### full ins part
                            full_final_compare_ref_result_df_ins_full = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,6]>-1]
                            for ll_ins_id in range(len(full_final_compare_ref_result_df_ins_full)):
                                tem_ins_loc = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,6]
                                tem_ins_content = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,7]
                                tem_ins_length = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,8]
                                Best_ref_ins_full[kk][tem_ins_loc] = tem_ins_content
                            
                            ##Del
                            full_final_compare_ref_result_df_del = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,9]>-1]
                            for ll_del_id in range(len(full_final_compare_ref_result_df_del)):
                                #tem_del_loc = full_final_compare_ref_result_df_del.iloc[ll_del_id,3]
                                #tem_del_content = full_final_compare_ref_result_df_del.iloc[ll_del_id,4]
                                tem_del_length = full_final_compare_ref_result_df_del.iloc[ll_del_id,9]
                                #Best_ref_ins[kk][tem_ins_loc] = tem_ins_content
                                if tem_del_length in tem_del_len_dic:
                                    tem_del_len_dic[tem_del_length] = tem_del_len_dic[tem_del_length] + 1
                                else:
                                    tem_del_len_dic[tem_del_length] = 1
                            
                            collection_id_del_covered = list(set(full_final_compare_ref_result_df.iloc[:,12]))
                            collection_id_del_covered.remove(-1)
                            remaining_del_id = list(set(range(len(common_del_block_df_cp)))-set(collection_id_del_covered))
                            
                            for ll_remain_del_id in remaining_del_id:
                                tem_del_length = common_del_block_df_cp.iloc[ll_remain_del_id,1] - common_del_block_df_cp.iloc[ll_remain_del_id,0] + 1
                                if tem_del_length>0:
                                    if tem_del_length not in tem_del_len_dic:
                                        tem_del_len_dic[tem_del_length] = 1
                                    else:
                                        tem_del_len_dic[tem_del_length] = tem_del_len_dic[tem_del_length] + 1
        
                            con_ref_prob_full[kk,5] =  con_ref_prob_full[kk,0] + con_ref_prob_full[kk,1]
                            for del_key in tem_del_len_dic:
                                con_ref_prob_full[kk,2] = con_ref_prob_full[kk,2] + tem_del_len_dic[del_key]
                               
                            for ins_key in tem_ins_len_dic:
                                con_ref_prob_full[kk,3] = con_ref_prob_full[kk,3] + tem_ins_len_dic[ins_key]
        
                            con_ref_prob_full[kk,4] = con_ref_prob_full[kk,5] - con_ref_prob_full[kk,2]-len(tem_del_len_dic) -1
                            
                            ### update of the bone
                            full_final_compare_ref_result_df_bone = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,10]>-1]
                            for bone_id_df in range(len(full_final_compare_ref_result_df_bone)):
                                bone_num_tem = full_final_compare_ref_result_df_bone.iloc[bone_id_df,10]
                                bone_base_infer = full_final_compare_ref_result_df_bone.iloc[bone_id_df,11]
                                Best_ref[kk][bone_num_tem] = bone_base_infer
                            for ll_ins in ins_inserted_dic:
                                if ll_ins not in Best_ref_ins_full[kk]:
                                    Best_ref_ins_full[kk][ll_ins] = '5'*ins_inserted_dic[ll_ins]
                             
                            
                        else:
                            for basss_index in only_bone_note_df_keys:
                                Best_ref[kk][int(basss_index)] = copy.deepcopy(
                                    Basic_ref[int(basss_index)])
                            con_ref_prob_full[kk, :] = [n_c, 0, 0, 0, n_c-1,n_c]
                            for ll_ins in ins_inserted_dic:
                                Best_ref_ins_full[kk][ll_ins] = '5'*ins_inserted_dic[ll_ins]
                            
                        
                        
                        evo_ins_len_dic.append(copy.deepcopy(tem_ins_len_dic))
                        evo_del_len_dic.append(copy.deepcopy(tem_del_len_dic))
                        
                        Best_ref_ins_keys[kk] = np.array(list(Best_ref_ins[kk].keys()))
                        
                    
            
                    portions = portions/sum(portions)
            
                    conditional_probalility_ = np.zeros([max_num_reads_refine, K])
                    conditional_probalility_o = np.zeros([max_num_reads_refine, K])
                    update_fre = [np.zeros([n_c_refine, 6]) for x in range(K)]
                    whole_ins_dic = {}
                    whole_del_dic = {}
            
                    with  mp.Pool(num_cores) as pool:
                
                     
                         final_compare_con_result = pool.map(collection_var_per_read, list(range(max_num_reads_refine)))
                    
                    ### making up the content
                    final_compare_con_result_full = []
                    for ll_per_read_id in range(max_num_reads_refine):
                        final_compare_con_result_full = final_compare_con_result_full + final_compare_con_result[ll_per_read_id]
                    
                    final_compare_con_result_full_df = pd.DataFrame(final_compare_con_result_full)
                    final_compare_con_result_df_sub = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,-1]>-1]
                    
                    reads_infer_clu.iloc[:,1] = final_compare_con_result_df_sub.iloc[:,-1]
                    ### initial to get the seperation of the Ref/Ins iid then for base; most important is the reads label 
                    ### and the informative parameters
                    for kkk in range(K):
                        
                        tem_gamma_ = np.array([gamma1[kkk], gamma2[kkk], gamma3[kkk], gamma4[kkk], 1-gamma4[kkk]])
                        tem_gamma_ = tem_gamma_ + corretion_eps_from_zero 
                        tem_gamma_[:2] = tem_gamma_[:2]/sum(tem_gamma_[:2])
                        tem_gamma_[2] = min(1,tem_gamma_[2])
                        tem_gamma_[3:] = tem_gamma_[3:]/sum(tem_gamma_[3:])
                        con_ref_prob[kkk] = con_ref_prob[kkk] + con_ref_prob_full[kkk, 0] * np.log(tem_gamma_[0])### correct prob
                        con_ref_prob[kkk] = con_ref_prob[kkk] + con_ref_prob_full[kkk, 1] * np.log(tem_gamma_[1])### mis prob
                        con_ref_prob[kkk] = con_ref_prob[kkk] + (con_ref_prob_full[kkk, 5]) * np.log(1-tem_gamma_[2])### not del prob
                        
                        con_ref_prob[kkk] = con_ref_prob[kkk] + (con_ref_prob_full[kkk, 2]) * np.log(tem_gamma_[2])### not del prob
                        
                        for del_dic_key in evo_del_len_dic[kkk]:
                            con_ref_prob[kkk] = con_ref_prob[kkk] + evo_del_len_dic[kkk][del_dic_key]*expotential_log_log_probability(del_dic_key,del_p[kkk],del_beta[kkk])
                        
                        ### ins part 
                        con_ref_prob[kkk] = con_ref_prob[kkk] + con_ref_prob_full[kkk, 3]*np.log(tem_gamma_[3])
                        for ins_dic_key in evo_ins_len_dic[kkk]:
                            con_ref_prob[kkk] = con_ref_prob[kkk] + evo_ins_len_dic[kkk][ins_dic_key]*expotential_log_log_probability(ins_dic_key,ins_p[kkk],ins_beta[kkk])
                        con_ref_prob[kkk] = con_ref_prob[kkk] +  con_ref_prob_full[kkk, 4]*np.log(tem_gamma_[4])
                        
                        
                        
                    # liklihood_.append(((conditional_probalility_o*conditional_probalility_).sum()+sum(con_ref_prob)))
            
                    creatierstop = 0.0001
                    #correct_num= np.zeros(K)
                    up_gamma1 = np.ones(K)*0.9
                    up_gamma2 = np.ones(K)*0.05
                    up_gamma3 = np.ones(K)*0.05
                    up_del_beta = np.zeros(K)
                    up_del_p = np.zeros(K)
                    up_gamma4 = np.ones(K)*0.0001
                    up_ins_beta = np.zeros(K)
                    up_ins_p = np.zeros(K)
                    #up_gamma5=  np.ones(K)*0.05
                    present_clu_dis2 = [
                        np.zeros([num_remain_infor_base, 6]) for x in range(K)]
                    con_ref_prob_full = np.zeros([K, 6])
                    evo_ins_len_dic = {}### this is to save the dic of the length of the ins length and del length
                    evo_del_len_dic = {}
                    #Best_ref = [{} for x in range(K)]
                    #Best_ref_ins = [{} for x in range(K)]
                    #Best_ref_ins_keys = [[] for x in range(K)]
                    #Best_ref_ins_full = [{} for x in range(K)]
                    for kk in range(K):
                        ele_list = list(reads_infer_clu[reads_infer_clu.iloc[:, 1] == kk].iloc[:, 0])
                        #ele_list = list(
                        #    reads_infer_cluf[reads_infer_cluf.iloc[:, 1] == kk].iloc[:, 0])
                        portions[kk] = len(ele_list)
                        #tem_eq_num = 0
                        #tem_neq_num = 0
                        tem_del_len_dic = {}
                        tem_ins_len_dic = {}
                        if len(ele_list)>0:
                            for ele in ele_list:
                                #tem_ele_start = final_full_result[ele][0]
                                #tem_ele_end = final_full_result[ele][1]
                                tem_matrix_refine = read_full_bases_filter[ele]
                                if len(tem_matrix_refine)>0:
                                    tem_start_refine = infer_var_site.index(tem_matrix_refine[0,0])
                                    tem_end_refine = tem_start_refine + len(tem_matrix_refine)-1
                                    present_clu_dis2[kk][tem_start_refine:tem_end_refine+1,:] = present_clu_dis2[kk][tem_start_refine:tem_end_refine+1,:] + \
                                        tem_matrix_refine[:,1:]
                            #max_ID = list(present_clu_dis2[kk][int(basss_index),:]).index(max(present_clu_dis2[kk][int(basss_index),:]))
                            
                            common_del_block_df_cp = copy.deepcopy(common_del_block_df)
                            keys_ins_common_original_cp = copy.deepcopy(keys_ins_common_original)
                            tem_touched_ins = []
                            
                            ### using the common part of the non-indel as the independent part
                            remaining_independent = copy.deepcopy(seg_collection_ini_df_dep_id)
                            with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
            
                                Result_infer_ind = pool.starmap(signle_sites_infer_max4, [(remaining_independent[index], present_clu_dis2[kk][int(remaining_independent[index]),:], standard_ref_code[infer_var_site[int(remaining_independent[index])]],kk
                                                                                       ) for index in range(len(remaining_independent))])
                            remaining_dependent = copy.deepcopy(seg_collection_ini_df_indel_part_cp_sub_link_full)
                            related_commom_del = copy.deepcopy(dep_spe_infor_del_col_full)
                            related_commom_ins = copy.deepcopy(dep_spe_infor_ins_col_full)
                            ## next we check the dependent part / we should spilt find the read overlapping with the region and common indel
                            with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                                Dep_dic = pool.starmap(select_read_blocks2, [(
                                remaining_dependent.iloc[index,0], remaining_dependent.iloc[index,1]) for index in range(len(remaining_dependent))])
                                #Dep_dic = Dep_dic_iter.get()
                            with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                                Ref_dic_block = pool.starmap(select_ref_block2, [(
                                 remaining_dependent.iloc[index,0], remaining_dependent.iloc[index,1],kk) for index in range(len(remaining_dependent))])
                            with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                                Result_infer_dep = pool.starmap(block_sites_infer_ECM3, [(int(remaining_dependent.iloc[index,0]), int(remaining_dependent.iloc[index,1]), present_clu_dis2[kk][int(remaining_dependent.iloc[index,0]):int(remaining_dependent.iloc[index,1])+1, :],
                                                                                      final_ref_code_standard[int(remaining_dependent.iloc[index,0]):int(
                                                                                          remaining_dependent.iloc[index,1])+1], copy.deepcopy(Ref_dic_block[index]), Dep_dic[index],related_commom_del[index],related_commom_ins[index],kk) for index in range(len(remaining_dependent))])
                                
                            Best_ref[kk] = {}
                            Best_ref_ins[kk] = {} 
                            Best_ref_ins_keys[kk] = []
                            Best_ref_ins_full[kk] = {}
                            
                            
                            # making up full result
                            Result_infer_ = copy.deepcopy(Result_infer_ind)
                            #dep_Result_infer_ = []
                            for ll_xun_ in Result_infer_dep:
                                Result_infer_ = Result_infer_ + ll_xun_### this is to make the reserved parts outsides from signs
                                #dep_Result_infer_ = dep_Result_infer_ + ll_xun_
                            Result_infer_ = sorted(Result_infer_)
                            Result_infer_df = pd.DataFrame(Result_infer_)
                            informative_sites_code[kk] = np.array(Result_infer_df.iloc[:,1])
                            
                            #dep_Result_infer_ = sorted(dep_Result_infer_)
                            
                            
                            #for ll in range(len(seg_collection_ini)):
                            #    if seg_collection_ini[ll][1]==6439:
                            #        print(seg_collection_ini[ll])
                            #        print(ll)
                            
                            #for ll2 in range(len(seg_collection_ini_df_indel_part_cp_sub_link_full)):
                            #    if seg_collection_ini_df_indel_part_cp_sub_link_full.iloc[ll2,1]==6439:
                            #        print(ll2)
                            #        print(seg_collection_ini_df_indel_part_cp_sub_link_full.iloc[ll2,:])
                            
                            with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                        
                             # Result_infer_ind=pool.starmap(signle_sites_infer_max,[(remaining_independent[index],present_clu_dis2[kk][int(remaining_independent[index]),:],standard_ref_code[int(remaining_independent[index])],\
                             # gamma[kk],theta_) for index in range(len(remaining_independent))])
                        
                                 final_compare_ref_result = pool.map(comparision_between_con_and_ref, list(range(len(seg_collection_ini))))
                                 #Result_infer_ind =  Result_infer_ind_iter.get()
                             #final_inferred=pool.starmap(sampling2,[(pp_random[index,0],pp_random[index,1],combination_list2.iloc[int(index/num_chains),0],combination_list2.iloc[int(index/num_chains),1],out_put_df,collection_of_sites,final_ref_code) for index in range((len(combination_list2)*num_chains))])
                             #pool.close()
                             #pool.join()
                            full_final_compare_ref_result = []
                             
                            for ll_final_id in range(len(final_compare_ref_result)):
                                 full_final_compare_ref_result = full_final_compare_ref_result + final_compare_ref_result[ll_final_id]
                             
                            full_final_compare_ref_result_df = pd.DataFrame(full_final_compare_ref_result)
                            
                            con_ref_prob_full[kk,0] = sum(full_final_compare_ref_result_df.iloc[:,1]) + Backgroud_evo_c
                            con_ref_prob_full[kk,1] = sum(full_final_compare_ref_result_df.iloc[:,2]) + Backgroud_evo_m
                            #sub_sub_result = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,-1]==1]
                            
                            ### ins part
                            full_final_compare_ref_result_df_ins = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,3]>-1]
                            for ll_ins_id in range(len(full_final_compare_ref_result_df_ins)):
                                tem_ins_loc = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,3]
                                tem_ins_content = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,4]
                                tem_ins_length = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,5]
                                Best_ref_ins[kk][tem_ins_loc] = tem_ins_content
                                if tem_ins_length in tem_ins_len_dic:
                                    tem_ins_len_dic[tem_ins_length] = tem_ins_len_dic[tem_ins_length] + 1
                                else:
                                    tem_ins_len_dic[tem_ins_length] = 1
                            
                            collection_id_ins_covered = list(Best_ref_ins[kk].keys())
                            remaining_ins_id = list(set(list(common_ins_block_dic.keys()))-set(collection_id_ins_covered))
                            
                            for ll_ins_id_mk in remaining_ins_id:
                                tem_ins_content = common_ins_block_dic[ll_ins_id_mk][0][2]
                                tem_ins_length = len(tem_ins_content)
                                Best_ref_ins[kk][ll_ins_id_mk] = tem_ins_content
                                if tem_ins_length in tem_ins_len_dic:
                                    tem_ins_len_dic[tem_ins_length] = tem_ins_len_dic[tem_ins_length] + 1
                                else:
                                    tem_ins_len_dic[tem_ins_length] = 1
                            
                            ### full ins part
                            full_final_compare_ref_result_df_ins_full = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,6]>-1]
                            for ll_ins_id in range(len(full_final_compare_ref_result_df_ins_full)):
                                tem_ins_loc = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,6]
                                tem_ins_content = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,7]
                                tem_ins_length = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,8]
                                Best_ref_ins_full[kk][tem_ins_loc] = tem_ins_content
                            
                            ##Del
                            full_final_compare_ref_result_df_del = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,9]>-1]
                            for ll_del_id in range(len(full_final_compare_ref_result_df_del)):
                                #tem_del_loc = full_final_compare_ref_result_df_del.iloc[ll_del_id,3]
                                #tem_del_content = full_final_compare_ref_result_df_del.iloc[ll_del_id,4]
                                tem_del_length = full_final_compare_ref_result_df_del.iloc[ll_del_id,9]
                                #Best_ref_ins[kk][tem_ins_loc] = tem_ins_content
                                if tem_del_length in tem_del_len_dic:
                                    tem_del_len_dic[tem_del_length] = tem_del_len_dic[tem_del_length] + 1
                                else:
                                    tem_del_len_dic[tem_del_length] = 1
                            
                            collection_id_del_covered = list(set(full_final_compare_ref_result_df.iloc[:,12]))
                            collection_id_del_covered.remove(-1)
                            remaining_del_id = list(set(range(len(common_del_block_df_cp)))-set(collection_id_del_covered))
                            
                            for ll_remain_del_id in remaining_del_id:
                                tem_del_length = common_del_block_df_cp.iloc[ll_remain_del_id,1] - common_del_block_df_cp.iloc[ll_remain_del_id,0] + 1
                                if tem_del_length>0:
                                    if tem_del_length not in tem_del_len_dic:
                                        tem_del_len_dic[tem_del_length] = 1
                                    else:
                                        tem_del_len_dic[tem_del_length] = tem_del_len_dic[tem_del_length] + 1
                            
                            
            
                            con_ref_prob_full[kk,5] =  con_ref_prob_full[kk,0] + con_ref_prob_full[kk,1]
                            for del_key in tem_del_len_dic:
                                con_ref_prob_full[kk,2] = con_ref_prob_full[kk,2] + tem_del_len_dic[del_key]
                                con_ref_prob_full[kk,5] = con_ref_prob_full[kk,5] + tem_del_len_dic[del_key]
                               
                            for ins_key in tem_ins_len_dic:
                                con_ref_prob_full[kk,3] = con_ref_prob_full[kk,3] + tem_ins_len_dic[ins_key]
            
                            con_ref_prob_full[kk,4] = con_ref_prob_full[kk,5] - con_ref_prob_full[kk,2]-len(tem_del_len_dic) -1
                            
                            ### update of the bone
                            full_final_compare_ref_result_df_bone = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,10]>-1]
                            for bone_id_df in range(len(full_final_compare_ref_result_df_bone)):
                                bone_num_tem = full_final_compare_ref_result_df_bone.iloc[bone_id_df,10]
                                bone_base_infer = full_final_compare_ref_result_df_bone.iloc[bone_id_df,11]
                                Best_ref[kk][bone_num_tem] = bone_base_infer
                            
                            tem_eq_num = con_ref_prob_full[kk, 0]
                            tem_neq_num = con_ref_prob_full[kk, 1]
                            tem_del_loc_num = con_ref_prob_full[kk, 2]
                            tem_del_loc_num_full = con_ref_prob_full[kk,5] 
                            tem_del_loc_prob = tem_del_loc_num/tem_del_loc_num_full
                            tem_ins_num = con_ref_prob_full[kk, 3]
                            tem_no_ins_num = con_ref_prob_full[kk, 4]
                            
                            up_gamma4[kk] = (tem_ins_num) / \
                                (tem_ins_num+tem_no_ins_num)
                            up_gamma1[kk] = tem_eq_num / \
                                (tem_eq_num+tem_neq_num)
                            up_gamma3[kk] = tem_del_loc_prob
                            up_gamma2[kk] = 1-up_gamma1[kk]
                        
                            if len(tem_del_len_dic)>0:
                                tem_dep_c = 0
                                for dd_k in  tem_del_len_dic:
                                    up_del_beta[kk] = up_del_beta[kk] + (dd_k/(1-(1-del_p[kk])*np.exp(-del_beta[kk]*dd_k)))*tem_del_len_dic[dd_k]
                                    tem_dep_c = tem_dep_c + 1/(1-(1-del_p[kk])*np.exp(-del_beta[kk]*dd_k))*tem_del_len_dic[dd_k]
                                import_data_p = (tem_dep_c,tem_del_loc_num)
                                p_00 = 0
                                up_del_p[kk] = fsolve(extimation_p, p_00,args=import_data_p)
                               
                                if up_del_p[kk]>1:
                                    up_del_p[kk] = 1
                                if up_del_p[kk]<0:
                                    up_del_p[kk] = 0
                                #up_del_p[kk] = -tem_del_loc_num*(1-del_p[kk])/(up_del_p[kk] * np.log(del_p[kk]))
                                
                                
                                up_del_beta[kk] = tem_del_loc_num /up_del_beta[kk]
                            else:
                                up_del_beta[kk] = 0.5
                                up_del_p[kk] = 0.5
                            
                            if len(tem_ins_len_dic)>0:
                                tem_dep_c = 0
                                for ii_k in  tem_ins_len_dic:
                                    up_ins_beta[kk] = up_ins_beta[kk] + (ii_k/(1-(1-ins_p[kk])*np.exp(-ins_beta[kk]*ii_k)))*tem_ins_len_dic[ii_k]
                                    #up_ins_p[kk] = up_ins_p[kk] + 1/(1-(1-ins_p[kk])*np.exp(-ins_beta[kk]*ii_k))*tem_ins_len_dic[ii_k]
                                    tem_dep_c = tem_dep_c + 1/(1-(1-ins_p[kk])*np.exp(-ins_beta[kk]*ii_k))*tem_ins_len_dic[ii_k]
                                import_data_p = (tem_dep_c,tem_ins_num)
                                p_00 = 0
                                up_ins_p[kk] = fsolve(extimation_p, p_00,args=import_data_p)
                                if up_ins_p[kk]<0:
                                    up_ins_p[kk] = 0
                                if up_ins_p[kk]>1:
                                    up_ins_p[kk] = 1
                               
                                #up_ins_p[kk] = -tem_ins_num*(1-ins_p[kk])/(up_ins_p[kk] * np.log(ins_p[kk]))
                                up_ins_beta[kk] = tem_ins_num /up_ins_beta[kk]
                            else:
                                up_ins_beta[kk] = 0.5
                                up_ins_p[kk] = 0.5
                            
                            #up_ins_p[kk] = -tem_ins_num*(1-up_ins_p[kk])/(up_ins_p[kk] * np.log(ins_p[kk]))
                            #up_ins_beta[kk] = tem_ins_num /up_ins_beta[kk]
                            #up_beta_a = tem_eq_num + 48
                            #up_beta_b = tem_neq_num + 4
                            #up_gamma[kk] =  np.random.beta(a = up_beta_a,b = up_beta_b,size=1)
                            #up_gamma[kk] =  tem_eq_num/(tem_neq_num+tem_eq_num)
                            if up_gamma1[kk] < 0.9:
                                up_gamma1[kk] = 0.9
                                up_gamma2[kk] = 1-up_gamma1[kk]
                            for ll_ins in ins_inserted_dic:
                                if ll_ins not in Best_ref_ins_full[kk]:
                                    Best_ref_ins_full[kk][ll_ins] = '5'*ins_inserted_dic[ll_ins]
                             
                               
                        else:
                            Best_ref[kk] = {}
                            Best_ref_ins[kk] = {} 
                            Best_ref_ins_keys[kk] = []
                            Best_ref_ins_full[kk] = {}
                            for basss_index in only_bone_note_df_keys:
                                Best_ref[kk][int(basss_index)] = copy.deepcopy(
                                    Basic_ref[int(basss_index)])
                            con_ref_prob_full[kk, :] = [n_c, 0, 0, 0, n_c-1,n_c]
                            up_gamma1[kk] = np.random.beta(a=48, b=4, size=1)
                            up_gamma2[kk] = 1-up_gamma1[kk]
                            up_gamma3[kk] = np.random.rand(1)*0.01
                            up_gamma4[kk] = np.random.rand(1)*0.01
                            
                            up_del_beta[kk] = 0.5
                            up_del_p[kk] = 0.5
                            up_ins_beta[kk] = 0.5
                            up_ins_p[kk] = 0.5
                            for ll_ins in ins_inserted_dic:
                                Best_ref_ins_full[kk][ll_ins] = '5'*ins_inserted_dic[ll_ins]
                            
                            
                        
                        evo_ins_len_dic[kk] = copy.deepcopy(tem_ins_len_dic)
                        evo_del_len_dic[kk] = copy.deepcopy(tem_del_len_dic)
                            
                        Best_ref_ins_keys[kk] = np.array(list(Best_ref_ins[kk].keys()))
                        
                    portions = portions + 0.001
                    portions = portions/sum(portions)
                    check_dis2 = abs(last_portions-portions).max()
            
                    # parameters
                    conditional_probalility_ = np.zeros([max_num_reads_refine, K])
                    conditional_probalility_o = np.zeros([max_num_reads_refine, K])
                    update_fre = [np.zeros([n_c_refine, 6]) for x in range(K)]
                    wei_del = 0
                    wei_ins = 0
                    whole_ins_dic = {}
                    whole_del_dic = {}
            
                    with  mp.Pool(num_cores) as pool:
                
                     
                         final_compare_con_result = pool.map(collection_var_per_read, list(range(max_num_reads_refine)))
                    
                    ### making up the content
                    final_compare_con_result_full = []
                    for ll_per_read_id in range(max_num_reads_refine):
                        final_compare_con_result_full = final_compare_con_result_full + final_compare_con_result[ll_per_read_id]
                    
                    final_compare_con_result_full_df = pd.DataFrame(final_compare_con_result_full)
                    final_compare_con_result_df_sub = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,-1]>-1]
                    
                    reads_infer_clu.iloc[:,1] = final_compare_con_result_df_sub.iloc[:,-1]
                    
                    ### correct_num
                    final_compare_con_result_df_sub_correct = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,1]>-1]
                    wei_cor = sum(final_compare_con_result_df_sub_correct.iloc[:,1])
                    
                    ### mis num
                    final_compare_con_result_df_sub_mis = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,2]>-1]
                    wei_mis = sum(final_compare_con_result_df_sub_mis.iloc[:,2])
                    
                    ### del
                    final_compare_con_result_df_sub_del = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,3]>-1]
                    for del_hang_check in range(len(final_compare_con_result_df_sub_del)):
                        tem_len_key_del = final_compare_con_result_df_sub_del.iloc[del_hang_check,3]
                        tem_len_num_del = final_compare_con_result_df_sub_del.iloc[del_hang_check,4]
                        if  tem_len_key_del in whole_del_dic:
                            whole_del_dic[tem_len_key_del] = whole_del_dic[tem_len_key_del] + tem_len_num_del
                        else:
                            whole_del_dic[tem_len_key_del] = tem_len_num_del
                        wei_del = wei_del + tem_len_num_del
                    
                    ### ins
                    final_compare_con_result_df_sub_ins = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,5]>-1]
                    for ins_hang_check in range(len(final_compare_con_result_df_sub_ins)):
                        tem_len_key_ins = final_compare_con_result_df_sub_ins.iloc[ins_hang_check,5]
                        tem_len_num_ins = final_compare_con_result_df_sub_ins.iloc[ins_hang_check,6]
                        if  tem_len_key_ins in whole_ins_dic:
                            whole_ins_dic[tem_len_key_ins] = whole_ins_dic[tem_len_key_ins] + tem_len_num_ins
                        else:
                            whole_ins_dic[tem_len_key_ins] = tem_len_num_ins
                        wei_ins = wei_ins + tem_len_num_ins
                    
                    ### prob
                    final_compare_con_result_df_sub_prob = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,0]>-1]
                    for ll_read_per in range(len(final_compare_con_result_df_sub_prob)):
                        conditional_probalility_[ll_read_per,:] = final_compare_con_result_df_sub_prob.iloc[ll_read_per,7:(7+K)]
                    
            
                    
                    cum_del = 0
                    for kk_del in whole_del_dic:
                        cum_del = cum_del+whole_del_dic[kk_del]*kk_del
            
                    cum_ins = 0
                    for kk_ins in whole_ins_dic:
                        cum_ins = cum_ins+whole_ins_dic[kk_ins] * kk_ins
                    up_lamuna_del = cum_del/(wei_del+corretion_eps_from_zero) -1 
                    up_lamuna_ins = cum_ins/(wei_ins+corretion_eps_from_zero) -1
                    if up_lamuna_del<0:
                        up_lamuna_del = corretion_eps_from_zero
                    if up_lamuna_ins<0:
                        up_lamuna_ins = corretion_eps_from_zero
                    
                    up_pcorrect = wei_cor/(wei_cor+wei_mis)
                    up_p_mis_loc = wei_mis/(wei_cor+wei_mis)
                    up_pins1_loc = wei_ins/(wei_cor+wei_mis)
                    up_pdel1_loc = wei_del/(wei_cor+wei_mis+cum_del)
                    
                    Differ_matrix = np.zeros([K,K],dtype=int)
                    for i_last in range(K):
                        for j_last in range(K):
                            differ_num_tem_array = abs(last_informative_sites_code[i_last]-informative_sites_code[j_last])
                            differ_num_ = len(differ_num_tem_array[differ_num_tem_array<1])
                            Differ_matrix[i_last,j_last] = differ_num_
                    
                    row_o3_mec_itera, col_o3_mec_itera = linear_sum_assignment(Differ_matrix)
                    differ_num_ =  Differ_matrix[row_o3_mec_itera,col_o3_mec_itera].sum()
                    if differ_num_/num_remain_infor_base<portions_var_fail:
                        check_dis1 = creatierstop/2
                    else:
                        check_dis1 = 1
                    
                    
                    
                    check_dis1 = abs(last_probability -
                                     conditional_probalility_).max()
                    check_dis3 = max(abs(up_lamuna_del-lamuna_del), abs(lamuna_ins-up_lamuna_ins),
                                     abs(up_pcorrect -
                                         pcorrect), abs(up_p_mis_loc-p_mis_loc),
                                     abs(up_pins1_loc -
                                         pins1_loc), abs(up_pdel1_loc-pdel1_loc),
                                     abs(up_gamma1 -
                                         gamma1).max(), abs(up_gamma2-gamma2).max(),
                                     abs(up_gamma3-gamma3).max(), abs(up_gamma4-gamma4).max(),
                                     abs(up_del_beta-del_beta).max(),abs(up_del_p-del_p).max(),
                                     abs(up_ins_beta-ins_beta).max(),abs(up_ins_p-ins_p).max())
            
                    # check_dis3 = max(abs(up_mu_del-mu_del),abs(up_sig_del-sig_del),\
                    #                 abs(up_mu_ins-mu_ins),abs(up_sig_ins-sig_ins),\
                    #                 abs(up_pcorrect-pcorrect),abs(up_p_mis_loc-p_mis_loc),\
                    #                abs(up_pins1_loc-pins1_loc),abs(up_pdel1_loc-pdel1_loc))
                    check_dis = max(check_dis1, check_dis2, check_dis3)
                    max_iter = 10
                    min_iter = 50
                    
                    loop_ = 0
                    check_dis_flag = 1
                    time_vec = []
                    
                    ## initial end itera
                    t_pre_now = time.time()
                    time_save_list.append(['Time for the first iteration', t_pre_now-t_pre])
                    #print(up_gamma4)
                    
                    while check_dis > creatierstop:
                        time_vec.append(time.time())
                        if loop_ < max_iter:
                            t_pre = copy.deepcopy(t_pre_now)
                            tao = tao * tao_step
                            #tao2 = tao2 * tao_step2
                            last_probability = copy.deepcopy(
                                conditional_probalility_)
                            last_fre = copy.deepcopy(update_fre)
                            conditional_probalility_ = np.zeros(
                                [max_num_reads_refine, K])
                            update_fre = [np.zeros([n_c_refine, 6])
                                          for x in range(K)]
                            last_portions = copy.deepcopy(portions)
                            p_mis_loc = copy.deepcopy(up_p_mis_loc)
                            pins1_loc = copy.deepcopy(up_pins1_loc)
                            lamuna_ins = copy.deepcopy(up_lamuna_ins)
                            #mu_ins = copy.deepcopy(up_mu_ins)
                            #sig_ins = copy.deepcopy(up_sig_ins)
                            gamma1 = copy.deepcopy(up_gamma1)
                            gamma2 = copy.deepcopy(up_gamma2)
                            gamma3 = copy.deepcopy(up_gamma3)
                            gamma4 = copy.deepcopy(up_gamma4)
                            del_beta = copy.deepcopy(up_del_beta)  # del
                            del_p = copy.deepcopy(up_del_p)
                            ins_beta= copy.deepcopy(up_ins_beta)  # del
                            ins_p = copy.deepcopy(up_ins_p)
                            
                            pdel1_loc = copy.deepcopy(up_pdel1_loc)
                            lamuna_del = copy.deepcopy(up_lamuna_del)
                            #mu_del = copy.deepcopy(up_mu_del)
                            #sig_del = copy.deepcopy(up_sig_del)
                            pcorrect = copy.deepcopy(up_pcorrect)
                            theta_ = [pcorrect, p_mis_loc, pins1_loc,
                                      pdel1_loc, lamuna_ins, lamuna_del]
                            con_ref_prob = np.zeros(K)
                            #Best_ref = [{} for x in range(K)]
                            #Best_ref_ins = [{} for x in range(K)]
                            #Best_ref_ins_keys = [[] for x in range(K)]
                            #Best_ref_ins_full = [{} for x in range(K)]
                            whole_ins_dic = {}
                            whole_del_dic = {}
                            wei_del = 0
                            wei_ins = 0
                            last_informative_sites_code = copy.deepcopy(informative_sites_code)
                            informative_sites_code = [contant_ref_code_var for x in range(K)]
                            #covering_sites_msa = [[] for x in range(K)]
                            #covering_sites_bone = [[] for x in range(K)]
                            #covering_sites_ins = [[] for x in range(K)]
                            for kkk in range(K):
                                
                                tem_gamma_ = np.array([gamma1[kkk], gamma2[kkk], gamma3[kkk], gamma4[kkk], 1-gamma4[kkk]])
                                tem_gamma_ = tem_gamma_ + corretion_eps_from_zero 
                                tem_gamma_[:2] = tem_gamma_[:2]/sum(tem_gamma_[:2])
                                tem_gamma_[2] = min(1,tem_gamma_[2])
                                tem_gamma_[3:] = tem_gamma_[3:]/sum(tem_gamma_[3:])
                                con_ref_prob[kkk] = con_ref_prob[kkk] + con_ref_prob_full[kkk, 0] * np.log(tem_gamma_[0])### correct prob
                                con_ref_prob[kkk] = con_ref_prob[kkk] + con_ref_prob_full[kkk, 1] * np.log(tem_gamma_[1])### mis prob
                                con_ref_prob[kkk] = con_ref_prob[kkk] + (con_ref_prob_full[kkk, 5]) * np.log(1-tem_gamma_[2])### not del prob
                                
                                con_ref_prob[kkk] = con_ref_prob[kkk] + (con_ref_prob_full[kkk, 2]) * np.log(tem_gamma_[2])### not del prob
                                
                                for del_dic_key in evo_del_len_dic[kkk]:
                                    con_ref_prob[kkk] = con_ref_prob[kkk] + evo_del_len_dic[kkk][del_dic_key]*expotential_log_log_probability(del_dic_key,del_p[kkk],del_beta[kkk])
                                
                                ### ins part 
                                con_ref_prob[kkk] = con_ref_prob[kkk] + con_ref_prob_full[kkk, 3]*np.log(tem_gamma_[3])
                                for ins_dic_key in evo_ins_len_dic[kkk]:
                                    con_ref_prob[kkk] = con_ref_prob[kkk] + evo_ins_len_dic[kkk][ins_dic_key]*expotential_log_log_probability(ins_dic_key,ins_p[kkk],ins_beta[kkk])
                                con_ref_prob[kkk] = con_ref_prob[kkk] +  con_ref_prob_full[kkk, 4]*np.log(tem_gamma_[4])
                                
                                
                                
                            # liklihood_.append(((conditional_probalility_o*conditional_probalility_).sum()+sum(con_ref_prob)))
            
                            #creatierstop = 0.0001
                            #correct_num= np.zeros(K)
                            up_gamma1 = np.ones(K)*0.9
                            up_gamma2 = np.ones(K)*0.05
                            up_gamma3 = np.ones(K)*0.05
                            up_del_beta = np.zeros(K)
                            up_del_p = np.zeros(K)
                            up_gamma4 = np.ones(K)*0.0001
                            up_ins_beta = np.zeros(K)
                            up_ins_p = np.zeros(K)
                            #up_gamma5=  np.ones(K)*0.05
                            present_clu_dis2 = [
                                np.zeros([num_remain_infor_base, 6]) for x in range(K)]
                            con_ref_prob_full = np.zeros([K, 6])
                            
                            evo_ins_len_dic = {}### this is to save the dic of the length of the ins length and del length
                            evo_del_len_dic = {}
                            for kk in range(K):
                                ele_list = list(reads_infer_clu[reads_infer_clu.iloc[:, 1] == kk].iloc[:, 0])
                                #ele_list = list(
                                #    reads_infer_cluf[reads_infer_cluf.iloc[:, 1] == kk].iloc[:, 0])
                                portions[kk] = len(ele_list)
                                #tem_eq_num = 0
                                #tem_neq_num = 0
                                tem_del_len_dic = {}
                                tem_ins_len_dic = {}
                                if len(ele_list)>0:
                                    for ele in ele_list:
                                        #tem_ele_start = final_full_result[ele][0]
                                        #tem_ele_end = final_full_result[ele][1]
                                        tem_matrix_refine = read_full_bases_filter[ele]
                                        if len(tem_matrix_refine)>0:
                                            tem_start_refine = infer_var_site.index(tem_matrix_refine[0,0])
                                            tem_end_refine = tem_start_refine + len(tem_matrix_refine)-1
                                            present_clu_dis2[kk][tem_start_refine:tem_end_refine+1,:] = present_clu_dis2[kk][tem_start_refine:tem_end_refine+1,:] + \
                                                tem_matrix_refine[:,1:]
                                    #max_ID = list(present_clu_dis2[kk][int(basss_index),:]).index(max(present_clu_dis2[kk][int(basss_index),:]))
                                    
                                    common_del_block_df_cp = copy.deepcopy(common_del_block_df)
                                    keys_ins_common_original_cp = copy.deepcopy(keys_ins_common_original)
                                    tem_touched_ins = []
                                    
                                    ### using the common part of the non-indel as the independent part
                                    remaining_independent = copy.deepcopy(seg_collection_ini_df_dep_id)
                                    with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                    
                                        Result_infer_ind = pool.starmap(signle_sites_infer_max4, [(remaining_independent[index], present_clu_dis2[kk][int(remaining_independent[index]),:], standard_ref_code[infer_var_site[int(remaining_independent[index])]],kk
                                                                                               ) for index in range(len(remaining_independent))])
                                    remaining_dependent = copy.deepcopy(seg_collection_ini_df_indel_part_cp_sub_link_full)
                                    related_commom_del = copy.deepcopy(dep_spe_infor_del_col_full)
                                    related_commom_ins = copy.deepcopy(dep_spe_infor_ins_col_full)
                                    ## next we check the dependent part / we should spilt find the read overlapping with the region and common indel
                                    with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                                        Dep_dic = pool.starmap(select_read_blocks2, [(
                                        remaining_dependent.iloc[index,0], remaining_dependent.iloc[index,1]) for index in range(len(remaining_dependent))])
                                        #Dep_dic = Dep_dic_iter.get()
                                    with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                                        Ref_dic_block = pool.starmap(select_ref_block2, [(
                                         remaining_dependent.iloc[index,0], remaining_dependent.iloc[index,1],kk) for index in range(len(remaining_dependent))])
                                    with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                                        Result_infer_dep = pool.starmap(block_sites_infer_ECM3, [(int(remaining_dependent.iloc[index,0]), int(remaining_dependent.iloc[index,1]), present_clu_dis2[kk][int(remaining_dependent.iloc[index,0]):int(remaining_dependent.iloc[index,1])+1, :],
                                                                                              final_ref_code_standard[int(remaining_dependent.iloc[index,0]):int(
                                                                                                  remaining_dependent.iloc[index,1])+1], copy.deepcopy(Ref_dic_block[index]), Dep_dic[index],related_commom_del[index],related_commom_ins[index],kk) for index in range(len(remaining_dependent))])
                                        
                                    Best_ref[kk] = {}
                                    Best_ref_ins[kk] = {} 
                                    Best_ref_ins_keys[kk] = []
                                    Best_ref_ins_full[kk] = {}
                                    
                                    # making up full result
                                    Result_infer_ = copy.deepcopy(Result_infer_ind)
                                    #dep_Result_infer_ = []
                                    for ll_xun_ in Result_infer_dep:
                                        Result_infer_ = Result_infer_ + ll_xun_### this is to make the reserved parts outsides from signs
                                        #dep_Result_infer_ = dep_Result_infer_ + ll_xun_
                                    Result_infer_ = sorted(Result_infer_)
                                    Result_infer_df = pd.DataFrame(Result_infer_)
                                    informative_sites_code[kk] = np.array(Result_infer_df.iloc[:,1])
                                    
                                    
                                    #for ll in range(len(seg_collection_ini)):
                                    #    if seg_collection_ini[ll][1]==6439:
                                    #        print(seg_collection_ini[ll])
                                    #        print(ll)
                                    
                                    #for ll2 in range(len(seg_collection_ini_df_indel_part_cp_sub_link_full)):
                                    #    if seg_collection_ini_df_indel_part_cp_sub_link_full.iloc[ll2,1]==6439:
                                    #        print(ll2)
                                    #        print(seg_collection_ini_df_indel_part_cp_sub_link_full.iloc[ll2,:])
                                    
                                    with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                                
                                     # Result_infer_ind=pool.starmap(signle_sites_infer_max,[(remaining_independent[index],present_clu_dis2[kk][int(remaining_independent[index]),:],standard_ref_code[int(remaining_independent[index])],\
                                     # gamma[kk],theta_) for index in range(len(remaining_independent))])
                                
                                         final_compare_ref_result = pool.map(comparision_between_con_and_ref, list(range(len(seg_collection_ini))))
                                         #Result_infer_ind =  Result_infer_ind_iter.get()
                                     #final_inferred=pool.starmap(sampling2,[(pp_random[index,0],pp_random[index,1],combination_list2.iloc[int(index/num_chains),0],combination_list2.iloc[int(index/num_chains),1],out_put_df,collection_of_sites,final_ref_code) for index in range((len(combination_list2)*num_chains))])
                                     #pool.close()
                                     #pool.join()
                                    full_final_compare_ref_result = []
                                     
                                    for ll_final_id in range(len(final_compare_ref_result)):
                                         full_final_compare_ref_result = full_final_compare_ref_result + final_compare_ref_result[ll_final_id]
                                     
                                    full_final_compare_ref_result_df = pd.DataFrame(full_final_compare_ref_result)
                                    
                                    con_ref_prob_full[kk,0] = sum(full_final_compare_ref_result_df.iloc[:,1]) + Backgroud_evo_c
                                    con_ref_prob_full[kk,1] = sum(full_final_compare_ref_result_df.iloc[:,2]) + Backgroud_evo_m
                                    #sub_sub_result = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,-1]==1]
                                    
                                    ### ins part
                                    full_final_compare_ref_result_df_ins = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,3]>-1]
                                    for ll_ins_id in range(len(full_final_compare_ref_result_df_ins)):
                                        tem_ins_loc = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,3]
                                        tem_ins_content = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,4]
                                        tem_ins_length = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,5]
                                        Best_ref_ins[kk][tem_ins_loc] = tem_ins_content
                                        if tem_ins_length in tem_ins_len_dic:
                                            tem_ins_len_dic[tem_ins_length] = tem_ins_len_dic[tem_ins_length] + 1
                                        else:
                                            tem_ins_len_dic[tem_ins_length] = 1
                                    
                                    collection_id_ins_covered = list(Best_ref_ins[kk].keys())
                                    remaining_ins_id = list(set(list(common_ins_block_dic.keys()))-set(collection_id_ins_covered))
                                    
                                    for ll_ins_id_mk in remaining_ins_id:
                                        tem_ins_content = common_ins_block_dic[ll_ins_id_mk][0][2]
                                        tem_ins_length = len(tem_ins_content)
                                        Best_ref_ins[kk][ll_ins_id_mk] = tem_ins_content
                                        if tem_ins_length in tem_ins_len_dic:
                                            tem_ins_len_dic[tem_ins_length] = tem_ins_len_dic[tem_ins_length] + 1
                                        else:
                                            tem_ins_len_dic[tem_ins_length] = 1
                                    
                                    ### full ins part
                                    full_final_compare_ref_result_df_ins_full = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,6]>-1]
                                    for ll_ins_id in range(len(full_final_compare_ref_result_df_ins_full)):
                                        tem_ins_loc = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,6]
                                        tem_ins_content = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,7]
                                        tem_ins_length = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,8]
                                        Best_ref_ins_full[kk][tem_ins_loc] = tem_ins_content
                                    
                                    ##Del
                                    full_final_compare_ref_result_df_del = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,9]>-1]
                                    for ll_del_id in range(len(full_final_compare_ref_result_df_del)):
                                        #tem_del_loc = full_final_compare_ref_result_df_del.iloc[ll_del_id,3]
                                        #tem_del_content = full_final_compare_ref_result_df_del.iloc[ll_del_id,4]
                                        tem_del_length = full_final_compare_ref_result_df_del.iloc[ll_del_id,9]
                                        #Best_ref_ins[kk][tem_ins_loc] = tem_ins_content
                                        if tem_del_length in tem_del_len_dic:
                                            tem_del_len_dic[tem_del_length] = tem_del_len_dic[tem_del_length] + 1
                                        else:
                                            tem_del_len_dic[tem_del_length] = 1
                                    
                                    collection_id_del_covered = list(set(full_final_compare_ref_result_df.iloc[:,12]))
                                    collection_id_del_covered.remove(-1)
                                    remaining_del_id = list(set(range(len(common_del_block_df_cp)))-set(collection_id_del_covered))
                                    
                                    for ll_remain_del_id in remaining_del_id:
                                        tem_del_length = common_del_block_df_cp.iloc[ll_remain_del_id,1] - common_del_block_df_cp.iloc[ll_remain_del_id,0] + 1
                                        if tem_del_length>0:
                                            if tem_del_length not in tem_del_len_dic:
                                                tem_del_len_dic[tem_del_length] = 1
                                            else:
                                                tem_del_len_dic[tem_del_length] = tem_del_len_dic[tem_del_length] + 1
                                    
                                    
                    
                                    con_ref_prob_full[kk,5] =  con_ref_prob_full[kk,0] + con_ref_prob_full[kk,1]
                                    for del_key in tem_del_len_dic:
                                        con_ref_prob_full[kk,2] = con_ref_prob_full[kk,2] + tem_del_len_dic[del_key]
                                        con_ref_prob_full[kk,5] = con_ref_prob_full[kk,5] + tem_del_len_dic[del_key]
                                       
                                    for ins_key in tem_ins_len_dic:
                                        con_ref_prob_full[kk,3] = con_ref_prob_full[kk,3] + tem_ins_len_dic[ins_key]
                    
                                    con_ref_prob_full[kk,4] = con_ref_prob_full[kk,5] - con_ref_prob_full[kk,2]-len(tem_del_len_dic) -1
                                    
                                    ### update of the bone
                                    full_final_compare_ref_result_df_bone = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,10]>-1]
                                    for bone_id_df in range(len(full_final_compare_ref_result_df_bone)):
                                        bone_num_tem = full_final_compare_ref_result_df_bone.iloc[bone_id_df,10]
                                        bone_base_infer = full_final_compare_ref_result_df_bone.iloc[bone_id_df,11]
                                        Best_ref[kk][bone_num_tem] = bone_base_infer
                                    
                                    tem_eq_num = con_ref_prob_full[kk, 0]
                                    tem_neq_num = con_ref_prob_full[kk, 1]
                                    tem_del_loc_num = con_ref_prob_full[kk, 2]
                                    tem_del_loc_num_full = con_ref_prob_full[kk,5] 
                                    tem_del_loc_prob = tem_del_loc_num/tem_del_loc_num_full
                                    tem_ins_num = con_ref_prob_full[kk, 3]
                                    tem_no_ins_num = con_ref_prob_full[kk, 4]
                                    
                                    up_gamma4[kk] = (tem_ins_num) / \
                                        (tem_ins_num+tem_no_ins_num)
                                    up_gamma1[kk] = tem_eq_num / \
                                        (tem_eq_num+tem_neq_num)
                                    up_gamma3[kk] = tem_del_loc_prob
                                    up_gamma2[kk] = 1-up_gamma1[kk]
                                
                                    if len(tem_del_len_dic)>0:
                                        tem_dep_c = 0
                                        for dd_k in  tem_del_len_dic:
                                            up_del_beta[kk] = up_del_beta[kk] + (dd_k/(1-(1-del_p[kk])*np.exp(-del_beta[kk]*dd_k)))*tem_del_len_dic[dd_k]
                                            tem_dep_c = tem_dep_c + 1/(1-(1-del_p[kk])*np.exp(-del_beta[kk]*dd_k))*tem_del_len_dic[dd_k]
                                        import_data_p = (tem_dep_c,tem_del_loc_num)
                                        p_00 = 0
                                        up_del_p[kk] = fsolve(extimation_p, p_00,args=import_data_p)
                                       
                                        if up_del_p[kk]>1:
                                            up_del_p[kk] = 1
                                        if up_del_p[kk]<0:
                                            up_del_p[kk] = 0
                                        #up_del_p[kk] = -tem_del_loc_num*(1-del_p[kk])/(up_del_p[kk] * np.log(del_p[kk]))
                                        
                                        
                                        up_del_beta[kk] = tem_del_loc_num /up_del_beta[kk]
                                    else:
                                        up_del_beta[kk] = 0.5
                                        up_del_p[kk] = 0.5
                                    
                                    if len(tem_ins_len_dic)>0:
                                        tem_dep_c = 0
                                        for ii_k in  tem_ins_len_dic:
                                            up_ins_beta[kk] = up_ins_beta[kk] + (ii_k/(1-(1-ins_p[kk])*np.exp(-ins_beta[kk]*ii_k)))*tem_ins_len_dic[ii_k]
                                            #up_ins_p[kk] = up_ins_p[kk] + 1/(1-(1-ins_p[kk])*np.exp(-ins_beta[kk]*ii_k))*tem_ins_len_dic[ii_k]
                                            tem_dep_c = tem_dep_c + 1/(1-(1-ins_p[kk])*np.exp(-ins_beta[kk]*ii_k))*tem_ins_len_dic[ii_k]
                                        import_data_p = (tem_dep_c,tem_ins_num)
                                        p_00 = 0
                                        up_ins_p[kk] = fsolve(extimation_p, p_00,args=import_data_p)
                                        if up_ins_p[kk]<0:
                                            up_ins_p[kk] = 0
                                        if up_ins_p[kk]>1:
                                            up_ins_p[kk] = 1
                                       
                                        #up_ins_p[kk] = -tem_ins_num*(1-ins_p[kk])/(up_ins_p[kk] * np.log(ins_p[kk]))
                                        up_ins_beta[kk] = tem_ins_num /up_ins_beta[kk]
                                    else:
                                        up_ins_beta[kk] = 0.5
                                        up_ins_p[kk] = 0.5
                                    
                                    #up_ins_p[kk] = -tem_ins_num*(1-up_ins_p[kk])/(up_ins_p[kk] * np.log(ins_p[kk]))
                                    #up_ins_beta[kk] = tem_ins_num /up_ins_beta[kk]
                                    #up_beta_a = tem_eq_num + 48
                                    #up_beta_b = tem_neq_num + 4
                                    #up_gamma[kk] =  np.random.beta(a = up_beta_a,b = up_beta_b,size=1)
                                    #up_gamma[kk] =  tem_eq_num/(tem_neq_num+tem_eq_num)
                                    if up_gamma1[kk] < 0.9:
                                        up_gamma1[kk] = 0.9
                                        up_gamma2[kk] = 1-up_gamma1[kk]
                                    
                                    for ll_ins in ins_inserted_dic:
                                        if ll_ins not in Best_ref_ins_full[kk]:
                                            Best_ref_ins_full[kk][ll_ins] = '5'*ins_inserted_dic[ll_ins]
                                     
                                       
                                else:
                                    Best_ref[kk] = {}
                                    Best_ref_ins[kk] = {} 
                                    Best_ref_ins_keys[kk] = []
                                    Best_ref_ins_full[kk] = {}
                                    
                                    for basss_index in only_bone_note_df_keys:
                                        Best_ref[kk][int(basss_index)] = copy.deepcopy(
                                            Basic_ref[int(basss_index)])
                                    con_ref_prob_full[kk, :] = [n_c, 0, 0, 0, n_c-1,n_c]
                                    up_gamma1[kk] = np.random.beta(a=48, b=4, size=1)
                                    up_gamma2[kk] = 1-up_gamma1[kk]
                                    up_gamma3[kk] = np.random.rand(1)*0.01
                                    up_gamma4[kk] = np.random.rand(1)*0.01
                                    
                                    up_del_beta[kk] = 0.5
                                    up_del_p[kk] = 0.5
                                    up_ins_beta[kk] = 0.5
                                    up_ins_p[kk] = 0.5
                                    for ll_ins in ins_inserted_dic:
                                        Best_ref_ins_full[kk][ll_ins] = '5'*ins_inserted_dic[ll_ins]
                                    
                                    
                                
                                evo_ins_len_dic[kk] = copy.deepcopy(tem_ins_len_dic)
                                evo_del_len_dic[kk] = copy.deepcopy(tem_del_len_dic)
                                    
                                Best_ref_ins_keys[kk] = np.array(list(Best_ref_ins[kk].keys()))
                                
            
                            portions = portions + 0.001
                            portions = portions/sum(portions)
                            check_dis2 = abs(last_portions-portions).max()
                    
                            # parameters
                            conditional_probalility_ = np.zeros([max_num_reads_refine, K])
                            conditional_probalility_o = np.zeros([max_num_reads_refine, K])
                            update_fre = [np.zeros([n_c_refine, 6]) for x in range(K)]
                            wei_del = 0
                            wei_ins = 0
                            whole_ins_dic = {}
                            whole_del_dic = {}
                    
                            with  mp.Pool(num_cores) as pool:
                        
                             
                                 final_compare_con_result = pool.map(collection_var_per_read, list(range(max_num_reads_refine)))
                            
                            ### making up the content
                            final_compare_con_result_full = []
                            for ll_per_read_id in range(max_num_reads_refine):
                                final_compare_con_result_full = final_compare_con_result_full + final_compare_con_result[ll_per_read_id]
                            
                            final_compare_con_result_full_df = pd.DataFrame(final_compare_con_result_full)
                            final_compare_con_result_df_sub = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,-1]>-1]
                            
                            reads_infer_clu.iloc[:,1] = final_compare_con_result_df_sub.iloc[:,-1]
                            
                            ### correct_num
                            final_compare_con_result_df_sub_correct = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,1]>-1]
                            wei_cor = sum(final_compare_con_result_df_sub_correct.iloc[:,1])
                            
                            ### mis num
                            final_compare_con_result_df_sub_mis = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,2]>-1]
                            wei_mis = sum(final_compare_con_result_df_sub_mis.iloc[:,2])
                            
                            ### del
                            final_compare_con_result_df_sub_del = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,3]>-1]
                            for del_hang_check in range(len(final_compare_con_result_df_sub_del)):
                                tem_len_key_del = final_compare_con_result_df_sub_del.iloc[del_hang_check,3]
                                tem_len_num_del = final_compare_con_result_df_sub_del.iloc[del_hang_check,4]
                                if  tem_len_key_del in whole_del_dic:
                                    whole_del_dic[tem_len_key_del] = whole_del_dic[tem_len_key_del] + tem_len_num_del
                                else:
                                    whole_del_dic[tem_len_key_del] = tem_len_num_del
                                wei_del = wei_del + tem_len_num_del
                            
                            ### ins
                            final_compare_con_result_df_sub_ins = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,5]>-1]
                            for ins_hang_check in range(len(final_compare_con_result_df_sub_ins)):
                                tem_len_key_ins = final_compare_con_result_df_sub_ins.iloc[ins_hang_check,5]
                                tem_len_num_ins = final_compare_con_result_df_sub_ins.iloc[ins_hang_check,6]
                                if  tem_len_key_ins in whole_ins_dic:
                                    whole_ins_dic[tem_len_key_ins] = whole_ins_dic[tem_len_key_ins] + tem_len_num_ins
                                else:
                                    whole_ins_dic[tem_len_key_ins] = tem_len_num_ins
                                wei_ins = wei_ins + tem_len_num_ins
                            
                            ### prob
                            final_compare_con_result_df_sub_prob = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,0]>-1]
                            for ll_read_per in range(len(final_compare_con_result_df_sub_prob)):
                                conditional_probalility_[ll_read_per,:] = final_compare_con_result_df_sub_prob.iloc[ll_read_per,7:(7+K)]
                            
                    
                            
                            cum_del = 0
                            for kk_del in whole_del_dic:
                                cum_del = cum_del+whole_del_dic[kk_del]*kk_del
                    
                            cum_ins = 0
                            for kk_ins in whole_ins_dic:
                                cum_ins = cum_ins+whole_ins_dic[kk_ins] * kk_ins
                            up_lamuna_del = cum_del/(wei_del+corretion_eps_from_zero) -1 
                            
                            up_lamuna_ins = cum_ins/(wei_ins+corretion_eps_from_zero) -1
                            if up_lamuna_del<0:
                                up_lamuna_del = corretion_eps_from_zero
                            if up_lamuna_ins<0:
                                up_lamuna_ins = corretion_eps_from_zero
                            
                            up_pcorrect = wei_cor/(wei_cor+wei_mis)
                            up_p_mis_loc = wei_mis/(wei_cor+wei_mis)
                            up_pins1_loc = wei_ins/(wei_cor+wei_mis)
                            up_pdel1_loc = wei_del/(wei_cor+wei_mis+cum_del)
                            
                            ### estimation of the consensus difference between iterations
                            #check_dis1 = 1
                            Differ_matrix = np.zeros([K,K],dtype=int)
                            for i_last in range(K):
                                for j_last in range(K):
                                    differ_num_tem_array = abs(last_informative_sites_code[i_last]-informative_sites_code[j_last])
                                    differ_num_ = len(differ_num_tem_array[differ_num_tem_array>0])
                                    Differ_matrix[i_last,j_last] = differ_num_
                            
                            row_o3_mec_itera, col_o3_mec_itera = linear_sum_assignment(Differ_matrix)
                            differ_num_ =  Differ_matrix[row_o3_mec_itera,col_o3_mec_itera].sum()
                            differ_num_ratio = differ_num_/(K*num_remain_infor_base)
                            if differ_num_ratio<portions_var_fail:
                                check_dis1 = creatierstop/2
                            else:
                                check_dis1 = 1
                            print("difference ratio of con in iterations", differ_num_ratio)    
                            
                            #check_dis1 = abs(last_probability -
                            #                 conditional_probalility_).max()
                            check_dis3 = max(abs(up_lamuna_del-lamuna_del), abs(lamuna_ins-up_lamuna_ins),
                                             abs(up_pcorrect -
                                                 pcorrect), abs(up_p_mis_loc-p_mis_loc),
                                             abs(up_pins1_loc -
                                                 pins1_loc), abs(up_pdel1_loc-pdel1_loc),
                                             abs(up_gamma1 -
                                                 gamma1).max(), abs(up_gamma2-gamma2).max(),
                                             abs(up_gamma3-gamma3).max(), abs(up_gamma4-gamma4).max(),
                                             abs(up_del_beta-del_beta).max(),abs(up_del_p-del_p).max(),
                                             abs(up_ins_beta-ins_beta).max(),abs(up_ins_p-ins_p).max())
                    
                            # check_dis3 = max(abs(up_mu_del-mu_del),abs(up_sig_del-sig_del),\
                            #                 abs(up_mu_ins-mu_ins),abs(up_sig_ins-sig_ins),\
                            #                 abs(up_pcorrect-pcorrect),abs(up_p_mis_loc-p_mis_loc),\
                            #                abs(up_pins1_loc-pins1_loc),abs(up_pdel1_loc-pdel1_loc))
                            check_dis = max(check_dis1, check_dis2, check_dis3)
                            print(check_dis)
                            t_pre_now = time.time()
                            time_save_list.append(['Time for the '+str(loop_)+' iteration', t_pre_now-t_pre])
                            
                            loop_ = loop_ + 1
                            # if  and loop_>min_iter:
                            #    check_dis_flag = 0
                        else:
                            break
            
                    tem_cluster_index = np.argmax(conditional_probalility_, axis=1)
                    Collection_label.append(tem_cluster_index)
                    #conditional_probalility_o_up = copy.deepcopy(conditional_probalility_o)
                    # for kk_xun in range(K):
                    #    conditional_probalility_o_up[:,kk_xun] = conditional_probalility_o_up[:,kk_xun] +
                    sum_of_prob.append(
                        ((conditional_probalility_o*conditional_probalility_).sum()+sum(con_ref_prob)))
                    method_index_c.append(method_index)
                    Best_ref_collection[method_index] = copy.deepcopy(Best_ref)
                    Best_ref_ins_collection[method_index] = copy.deepcopy(Best_ref_ins)
                    parameters_saving[method_index] = [gamma1,gamma2,gamma3,gamma4,del_beta,del_p,\
                                                       ins_beta, ins_p, pcorrect,p_mis_loc,pins1_loc,lamuna_ins,pdel1_loc,lamuna_del]
                    
                # else:
                #    Collection_label.append([])
            else:
                #if len(set(Label_collection[method_index])) == K and method_index not in no_use_clu:
                    #break
                    #method_index = 6
                    t_pre = time.time()
                    time_save_list.append(['Initialization of the group', t_pre-t6])
                    
                    reads_infer_clu = pd.DataFrame(np.zeros([max_num_reads_refine, 2],dtype=int))
                    reads_infer_clu.iloc[:, 0] = range(max_num_reads_refine)
                    reads_infer_clu.iloc[:, 1] = copy.deepcopy(Label_collection[method_index])
                    
                    last_probability = np.zeros([max_num_reads_refine, K])
                    for read_index in range(max_num_reads_refine):
                        last_probability[read_index, int(
                            Label_collection[method_index][read_index])] = 1
                    last_portions = np.ones(K)/K
                    present_clu_dis2 = [
                        np.zeros([num_remain_infor_base, 6],dtype = int) for x in range(K)]
                    #present_clu_dis2_pro = [np.zeros([n_c_refine,5]) for x in range(K)]
                    Best_ref = [{} for x in range(K)]
                    Best_ref_ins = [{} for x in range(K)]
                    Best_ref_ins_keys = [[] for x in range(K)]
                    Best_ref_ins_full = [{} for x in range(K)]
                    last_informative_sites_code = [ contant_ref_code_var for x in range(K)]
                    informative_sites_code = [ contant_ref_code_var for x in range(K)]
                    
                    portions = np.zeros(K)
                    Basic_ref = ref_code_standard#final_ref_code_standard
                    #con_ref_prob = np.zeros([max_num_reads_refine,K])
                    con_ref_prob = np.zeros(K)
                    con_ref_prob_full = np.zeros([K, 6])
                    gamma1 = 0.9*np.ones(K)  # correct
                    gamma2 = 0.05*np.ones(K)  # mis
                    gamma3 = 0.05*np.ones(K) # del
                    gamma4 = 0.05*np.ones(K) # ins
                    del_beta = 0.05*np.ones(K)  # del
                    del_p = 0.05*np.ones(K)
                    ins_beta= 0.05*np.ones(K)  # del
                    ins_p = 0.05*np.ones(K)
                    
                    p_mis_loc = np.random.rand(1)/10
                    pins1_loc = np.random.rand(1)/10
                    lamuna_ins = int(np.random.rand(1)*6)+1
                    
                    tao = 1
                    #tao2 = 1
                    tao_step = 0.97
                    #tao_step2 = 0.99
                    pdel1_loc = np.random.rand(1)/10
                    lamuna_del = int(np.random.rand(1)*6)+1
                    
            
                    pcorrect = np.random.rand(1)/10+0.9
                    theta_ = [pcorrect, p_mis_loc, pins1_loc,
                              pdel1_loc, lamuna_ins, lamuna_del]
                    evo_ins_len_dic = []### this is to save the dic of the length of the ins length and del length
                    evo_del_len_dic = []
                    for kk in range(K):
                        ele_list = list(reads_infer_clu[reads_infer_clu.iloc[:, 1] == kk].iloc[:, 0])
                        #ele_list = list(
                        #    reads_infer_cluf[reads_infer_cluf.iloc[:, 1] == kk].iloc[:, 0])
                        portions[kk] = len(ele_list)
                        #tem_eq_num = 0
                        #tem_neq_num = 0
                        for ele in ele_list:
                            #tem_ele_start = final_full_result[ele][0]
                            #tem_ele_end = final_full_result[ele][1]
                            tem_matrix_refine = read_full_bases_filter[ele]
                            if len(tem_matrix_refine)>0:
                                tem_start_refine = infer_var_site.index(tem_matrix_refine[0,0])
                                tem_end_refine = tem_start_refine + len(tem_matrix_refine)-1
                                present_clu_dis2[kk][tem_start_refine:tem_end_refine+1,:] = present_clu_dis2[kk][tem_start_refine:tem_end_refine+1,:] + \
                                    tem_matrix_refine[:,1:]
                        #max_ID = list(present_clu_dis2[kk][int(basss_index),:]).index(max(present_clu_dis2[kk][int(basss_index),:]))
                        tem_del_len_dic = {}
                        tem_ins_len_dic = {}
                        common_del_block_df_cp = copy.deepcopy(common_del_block_df)
                        keys_ins_common_original_cp = copy.deepcopy(keys_ins_common_original)
                        tem_touched_ins = []
                        if len(ele_list) > 0:
                            # 1.seperate the room for seperate zoom for block indel
            
                            remaining_independent = copy.deepcopy(infer_var_site)
                            with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
            
                            # Result_infer_ind=pool.starmap(signle_sites_infer_max,[(remaining_independent[index],present_clu_dis2[kk][int(remaining_independent[index]),:],standard_ref_code[int(remaining_independent[index])],\
                            # gamma[kk],theta_) for index in range(len(remaining_independent))])
            
                                Result_infer_ind = pool.starmap(pre_max_, [(index, present_clu_dis2[kk][(index), :]) for index in range(len(remaining_independent))])
                                #Result_infer_ind =  Result_infer_ind_iter.get()
                            #final_inferred=pool.starmap(sampling2,[(pp_random[index,0],pp_random[index,1],combination_list2.iloc[int(index/num_chains),0],combination_list2.iloc[int(index/num_chains),1],out_put_df,collection_of_sites,final_ref_code) for index in range((len(combination_list2)*num_chains))])
                            #pool.close()
                            #pool.join()
            
                            # making up full result
                            Result_infer_ = copy.deepcopy(Result_infer_ind)
                            Result_infer_df = pd.DataFrame(Result_infer_)
                            last_informative_sites_code[kk] = np.array(Result_infer_df.iloc[:,1])
                            
                            with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                        
                             # Result_infer_ind=pool.starmap(signle_sites_infer_max,[(remaining_independent[index],present_clu_dis2[kk][int(remaining_independent[index]),:],standard_ref_code[int(remaining_independent[index])],\
                             # gamma[kk],theta_) for index in range(len(remaining_independent))])
                        
                                 final_compare_ref_result = pool.map(comparision_between_con_and_ref, list(range(len(seg_collection_ini))))
                                 #Result_infer_ind =  Result_infer_ind_iter.get()
                             #final_inferred=pool.starmap(sampling2,[(pp_random[index,0],pp_random[index,1],combination_list2.iloc[int(index/num_chains),0],combination_list2.iloc[int(index/num_chains),1],out_put_df,collection_of_sites,final_ref_code) for index in range((len(combination_list2)*num_chains))])
                             #pool.close()
                             #pool.join()
                            full_final_compare_ref_result = []
                             
                            for ll_final_id in range(len(final_compare_ref_result)):
                                 full_final_compare_ref_result = full_final_compare_ref_result + final_compare_ref_result[ll_final_id]
                             
                            full_final_compare_ref_result_df = pd.DataFrame(full_final_compare_ref_result)
                            
                            
                            con_ref_prob_full[kk,0] = sum(full_final_compare_ref_result_df.iloc[:,1]) + Backgroud_evo_c
                            con_ref_prob_full[kk,1] = sum(full_final_compare_ref_result_df.iloc[:,2]) + Backgroud_evo_m
                            #sub_sub_result = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,-1]==1]
                            
                            ### ins part
                            full_final_compare_ref_result_df_ins = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,3]>-1]
                            for ll_ins_id in range(len(full_final_compare_ref_result_df_ins)):
                                tem_ins_loc = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,3]
                                tem_ins_content = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,4]
                                tem_ins_length = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,5]
                                Best_ref_ins[kk][tem_ins_loc] = tem_ins_content
                                if tem_ins_length in tem_ins_len_dic:
                                    tem_ins_len_dic[tem_ins_length] = tem_ins_len_dic[tem_ins_length] + 1
                                else:
                                    tem_ins_len_dic[tem_ins_length] = 1
                            
                            collection_id_ins_covered = list(Best_ref_ins[kk].keys())
                            remaining_ins_id = list(set(list(common_ins_block_dic.keys()))-set(collection_id_ins_covered))
                            
                            for ll_ins_id_mk in remaining_ins_id:
                                tem_ins_content = common_ins_block_dic[ll_ins_id_mk][0][2]
                                tem_ins_length = len(tem_ins_content)
                                Best_ref_ins[kk][ll_ins_id_mk] = tem_ins_content
                                if tem_ins_length in tem_ins_len_dic:
                                    tem_ins_len_dic[tem_ins_length] = tem_ins_len_dic[tem_ins_length] + 1
                                else:
                                    tem_ins_len_dic[tem_ins_length] = 1
                            ### full ins part
                            full_final_compare_ref_result_df_ins_full = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,6]>-1]
                            for ll_ins_id in range(len(full_final_compare_ref_result_df_ins_full)):
                                tem_ins_loc = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,6]
                                tem_ins_content = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,7]
                                tem_ins_length = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,8]
                                Best_ref_ins_full[kk][tem_ins_loc] = tem_ins_content
                            
                            ##Del
                            full_final_compare_ref_result_df_del = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,9]>-1]
                            for ll_del_id in range(len(full_final_compare_ref_result_df_del)):
                                #tem_del_loc = full_final_compare_ref_result_df_del.iloc[ll_del_id,3]
                                #tem_del_content = full_final_compare_ref_result_df_del.iloc[ll_del_id,4]
                                tem_del_length = full_final_compare_ref_result_df_del.iloc[ll_del_id,9]
                                #Best_ref_ins[kk][tem_ins_loc] = tem_ins_content
                                if tem_del_length in tem_del_len_dic:
                                    tem_del_len_dic[tem_del_length] = tem_del_len_dic[tem_del_length] + 1
                                else:
                                    tem_del_len_dic[tem_del_length] = 1
                            
                            collection_id_del_covered = list(set(full_final_compare_ref_result_df.iloc[:,12]))
                            collection_id_del_covered.remove(-1)
                            remaining_del_id = list(set(range(len(common_del_block_df_cp)))-set(collection_id_del_covered))
                            
                            for ll_remain_del_id in remaining_del_id:
                                tem_del_length = common_del_block_df_cp.iloc[ll_remain_del_id,1] - common_del_block_df_cp.iloc[ll_remain_del_id,0] + 1
                                if tem_del_length>0:
                                    if tem_del_length not in tem_del_len_dic:
                                        tem_del_len_dic[tem_del_length] = 1
                                    else:
                                        tem_del_len_dic[tem_del_length] = tem_del_len_dic[tem_del_length] + 1
        
                            con_ref_prob_full[kk,5] =  con_ref_prob_full[kk,0] + con_ref_prob_full[kk,1]
                            for del_key in tem_del_len_dic:
                                con_ref_prob_full[kk,2] = con_ref_prob_full[kk,2] + tem_del_len_dic[del_key]
                               
                            for ins_key in tem_ins_len_dic:
                                con_ref_prob_full[kk,3] = con_ref_prob_full[kk,3] + tem_ins_len_dic[ins_key]
        
                            con_ref_prob_full[kk,4] = con_ref_prob_full[kk,5] - con_ref_prob_full[kk,2]-len(tem_del_len_dic) -1
                            
                            ### update of the bone
                            full_final_compare_ref_result_df_bone = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,10]>-1]
                            for bone_id_df in range(len(full_final_compare_ref_result_df_bone)):
                                bone_num_tem = full_final_compare_ref_result_df_bone.iloc[bone_id_df,10]
                                bone_base_infer = full_final_compare_ref_result_df_bone.iloc[bone_id_df,11]
                                Best_ref[kk][bone_num_tem] = bone_base_infer
                            
                            for ll_ins in ins_inserted_dic:
                                if ll_ins not in Best_ref_ins_full[kk]:
                                    Best_ref_ins_full[kk][ll_ins] = '5'*ins_inserted_dic[ll_ins]
                             
                            
                        else:
                            for basss_index in only_bone_note_df_keys:
                                Best_ref[kk][int(basss_index)] = copy.deepcopy(
                                    Basic_ref[int(basss_index)])
                            con_ref_prob_full[kk, :] = [n_c, 0, 0, 0, n_c-1,n_c]
                            for ll_ins in ins_inserted_dic:
                                Best_ref_ins_full[kk][ll_ins] = '5'*ins_inserted_dic[ll_ins]
                            
                        
                        
                        evo_ins_len_dic.append(copy.deepcopy(tem_ins_len_dic))
                        evo_del_len_dic.append(copy.deepcopy(tem_del_len_dic))
                        
                        Best_ref_ins_keys[kk] = np.array(list(Best_ref_ins[kk].keys()))
                        
                    
            
                    portions = portions/sum(portions)
            
                    conditional_probalility_ = np.zeros([max_num_reads_refine, K])
                    conditional_probalility_o = np.zeros([max_num_reads_refine, K])
                    update_fre = [np.zeros([n_c_refine, 6]) for x in range(K)]
                    whole_ins_dic = {}
                    whole_del_dic = {}
            
                    with  mp.Pool(num_cores) as pool:
                
                     
                         final_compare_con_result = pool.map(collection_var_per_read, list(range(max_num_reads_refine)))
                    
                    ### making up the content
                    final_compare_con_result_full = []
                    for ll_per_read_id in range(max_num_reads_refine):
                        final_compare_con_result_full = final_compare_con_result_full + final_compare_con_result[ll_per_read_id]
                    
                    final_compare_con_result_full_df = pd.DataFrame(final_compare_con_result_full)
                    final_compare_con_result_df_sub = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,-1]>-1]
                    
                    reads_infer_clu.iloc[:,1] = final_compare_con_result_df_sub.iloc[:,-1]
                    ### initial to get the seperation of the Ref/Ins iid then for base; most important is the reads label 
                    ### and the informative parameters
                    for kkk in range(K):
                        
                        tem_gamma_ = np.array([gamma1[kkk], gamma2[kkk], gamma3[kkk], gamma4[kkk], 1-gamma4[kkk]])
                        tem_gamma_ = tem_gamma_ + corretion_eps_from_zero 
                        tem_gamma_[:2] = tem_gamma_[:2]/sum(tem_gamma_[:2])
                        tem_gamma_[2] = min(1,tem_gamma_[2])
                        tem_gamma_[3:] = tem_gamma_[3:]/sum(tem_gamma_[3:])
                        con_ref_prob[kkk] = con_ref_prob[kkk] + con_ref_prob_full[kkk, 0] * np.log(tem_gamma_[0])### correct prob
                        con_ref_prob[kkk] = con_ref_prob[kkk] + con_ref_prob_full[kkk, 1] * np.log(tem_gamma_[1])### mis prob
                        con_ref_prob[kkk] = con_ref_prob[kkk] + (con_ref_prob_full[kkk, 5]) * np.log(1-tem_gamma_[2])### not del prob
                        
                        con_ref_prob[kkk] = con_ref_prob[kkk] + (con_ref_prob_full[kkk, 2]) * np.log(tem_gamma_[2])### not del prob
                        
                        for del_dic_key in evo_del_len_dic[kkk]:
                            con_ref_prob[kkk] = con_ref_prob[kkk] + evo_del_len_dic[kkk][del_dic_key]*expotential_log_log_probability(del_dic_key,del_p[kkk],del_beta[kkk])
                        
                        ### ins part 
                        con_ref_prob[kkk] = con_ref_prob[kkk] + con_ref_prob_full[kkk, 3]*np.log(tem_gamma_[3])
                        for ins_dic_key in evo_ins_len_dic[kkk]:
                            con_ref_prob[kkk] = con_ref_prob[kkk] + evo_ins_len_dic[kkk][ins_dic_key]*expotential_log_log_probability(ins_dic_key,ins_p[kkk],ins_beta[kkk])
                        con_ref_prob[kkk] = con_ref_prob[kkk] +  con_ref_prob_full[kkk, 4]*np.log(tem_gamma_[4])
                        
                        
                        
                    # liklihood_.append(((conditional_probalility_o*conditional_probalility_).sum()+sum(con_ref_prob)))
            
                    creatierstop = 0.0001
                    #correct_num= np.zeros(K)
                    up_gamma1 = np.ones(K)*0.9
                    up_gamma2 = np.ones(K)*0.05
                    up_gamma3 = np.ones(K)*0.05
                    up_del_beta = np.zeros(K)
                    up_del_p = np.zeros(K)
                    up_gamma4 = np.ones(K)*0.0001
                    up_ins_beta = np.zeros(K)
                    up_ins_p = np.zeros(K)
                    #up_gamma5=  np.ones(K)*0.05
                    present_clu_dis2 = [
                        np.zeros([num_remain_infor_base, 6]) for x in range(K)]
                    con_ref_prob_full = np.zeros([K, 6])
                    evo_ins_len_dic = {}### this is to save the dic of the length of the ins length and del length
                    evo_del_len_dic = {}
                    #Best_ref = [{} for x in range(K)]
                    #Best_ref_ins = [{} for x in range(K)]
                    #Best_ref_ins_keys = [[] for x in range(K)]
                    #Best_ref_ins_full = [{} for x in range(K)]
                    for kk in range(K):
                        ele_list = list(reads_infer_clu[reads_infer_clu.iloc[:, 1] == kk].iloc[:, 0])
                        #ele_list = list(
                        #    reads_infer_cluf[reads_infer_cluf.iloc[:, 1] == kk].iloc[:, 0])
                        portions[kk] = len(ele_list)
                        #tem_eq_num = 0
                        #tem_neq_num = 0
                        tem_del_len_dic = {}
                        tem_ins_len_dic = {}
                        if len(ele_list)>0:
                            for ele in ele_list:
                                #tem_ele_start = final_full_result[ele][0]
                                #tem_ele_end = final_full_result[ele][1]
                                tem_matrix_refine = read_full_bases_filter[ele]
                                if len(tem_matrix_refine)>0:
                                    tem_start_refine = infer_var_site.index(tem_matrix_refine[0,0])
                                    tem_end_refine = tem_start_refine + len(tem_matrix_refine)-1
                                    present_clu_dis2[kk][tem_start_refine:tem_end_refine+1,:] = present_clu_dis2[kk][tem_start_refine:tem_end_refine+1,:] + \
                                        tem_matrix_refine[:,1:]
                            #max_ID = list(present_clu_dis2[kk][int(basss_index),:]).index(max(present_clu_dis2[kk][int(basss_index),:]))
                            
                            common_del_block_df_cp = copy.deepcopy(common_del_block_df)
                            keys_ins_common_original_cp = copy.deepcopy(keys_ins_common_original)
                            tem_touched_ins = []
                            
                            ### using the common part of the non-indel as the independent part
                            remaining_independent = copy.deepcopy(seg_collection_ini_df_dep_id)
                            with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
            
                                Result_infer_ind = pool.starmap(signle_sites_infer_max4, [(remaining_independent[index], present_clu_dis2[kk][int(remaining_independent[index]),:], standard_ref_code[infer_var_site[int(remaining_independent[index])]],kk
                                                                                       ) for index in range(len(remaining_independent))])
                            remaining_dependent = copy.deepcopy(seg_collection_ini_df_indel_part_cp_sub_link_full)
                            related_commom_del = copy.deepcopy(dep_spe_infor_del_col_full)
                            related_commom_ins = copy.deepcopy(dep_spe_infor_ins_col_full)
                            ## next we check the dependent part / we should spilt find the read overlapping with the region and common indel
                            with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                                Dep_dic = pool.starmap(select_read_blocks2, [(
                                remaining_dependent.iloc[index,0], remaining_dependent.iloc[index,1]) for index in range(len(remaining_dependent))])
                                #Dep_dic = Dep_dic_iter.get()
                            with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                                Ref_dic_block = pool.starmap(select_ref_block2, [(
                                 remaining_dependent.iloc[index,0], remaining_dependent.iloc[index,1],kk) for index in range(len(remaining_dependent))])
                            with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                                Result_infer_dep = pool.starmap(block_sites_infer_ECM3, [(int(remaining_dependent.iloc[index,0]), int(remaining_dependent.iloc[index,1]), present_clu_dis2[kk][int(remaining_dependent.iloc[index,0]):int(remaining_dependent.iloc[index,1])+1, :],
                                                                                      final_ref_code_standard[int(remaining_dependent.iloc[index,0]):int(
                                                                                          remaining_dependent.iloc[index,1])+1], copy.deepcopy(Ref_dic_block[index]), Dep_dic[index],related_commom_del[index],related_commom_ins[index],kk) for index in range(len(remaining_dependent))])
                                
                            Best_ref[kk] = {}
                            Best_ref_ins[kk] = {} 
                            Best_ref_ins_keys[kk] = []
                            Best_ref_ins_full[kk] = {}
                            
                            
                            # making up full result
                            Result_infer_ = copy.deepcopy(Result_infer_ind)
                            #dep_Result_infer_ = []
                            for ll_xun_ in Result_infer_dep:
                                Result_infer_ = Result_infer_ + ll_xun_### this is to make the reserved parts outsides from signs
                                #dep_Result_infer_ = dep_Result_infer_ + ll_xun_
                            Result_infer_ = sorted(Result_infer_)
                            Result_infer_df = pd.DataFrame(Result_infer_)
                            informative_sites_code[kk] = np.array(Result_infer_df.iloc[:,1])
                            
                            #dep_Result_infer_ = sorted(dep_Result_infer_)
                            
                            
                            #for ll in range(len(seg_collection_ini)):
                            #    if seg_collection_ini[ll][1]==6439:
                            #        print(seg_collection_ini[ll])
                            #        print(ll)
                            
                            #for ll2 in range(len(seg_collection_ini_df_indel_part_cp_sub_link_full)):
                            #    if seg_collection_ini_df_indel_part_cp_sub_link_full.iloc[ll2,1]==6439:
                            #        print(ll2)
                            #        print(seg_collection_ini_df_indel_part_cp_sub_link_full.iloc[ll2,:])
                            
                            with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                        
                             # Result_infer_ind=pool.starmap(signle_sites_infer_max,[(remaining_independent[index],present_clu_dis2[kk][int(remaining_independent[index]),:],standard_ref_code[int(remaining_independent[index])],\
                             # gamma[kk],theta_) for index in range(len(remaining_independent))])
                        
                                 final_compare_ref_result = pool.map(comparision_between_con_and_ref, list(range(len(seg_collection_ini))))
                                 #Result_infer_ind =  Result_infer_ind_iter.get()
                             #final_inferred=pool.starmap(sampling2,[(pp_random[index,0],pp_random[index,1],combination_list2.iloc[int(index/num_chains),0],combination_list2.iloc[int(index/num_chains),1],out_put_df,collection_of_sites,final_ref_code) for index in range((len(combination_list2)*num_chains))])
                             #pool.close()
                             #pool.join()
                            full_final_compare_ref_result = []
                             
                            for ll_final_id in range(len(final_compare_ref_result)):
                                 full_final_compare_ref_result = full_final_compare_ref_result + final_compare_ref_result[ll_final_id]
                             
                            full_final_compare_ref_result_df = pd.DataFrame(full_final_compare_ref_result)
                            
                            con_ref_prob_full[kk,0] = sum(full_final_compare_ref_result_df.iloc[:,1]) + Backgroud_evo_c
                            con_ref_prob_full[kk,1] = sum(full_final_compare_ref_result_df.iloc[:,2]) + Backgroud_evo_m
                            #sub_sub_result = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,-1]==1]
                            
                            ### ins part
                            full_final_compare_ref_result_df_ins = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,3]>-1]
                            for ll_ins_id in range(len(full_final_compare_ref_result_df_ins)):
                                tem_ins_loc = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,3]
                                tem_ins_content = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,4]
                                tem_ins_length = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,5]
                                Best_ref_ins[kk][tem_ins_loc] = tem_ins_content
                                if tem_ins_length in tem_ins_len_dic:
                                    tem_ins_len_dic[tem_ins_length] = tem_ins_len_dic[tem_ins_length] + 1
                                else:
                                    tem_ins_len_dic[tem_ins_length] = 1
                            
                            collection_id_ins_covered = list(Best_ref_ins[kk].keys())
                            remaining_ins_id = list(set(list(common_ins_block_dic.keys()))-set(collection_id_ins_covered))
                            
                            for ll_ins_id_mk in remaining_ins_id:
                                tem_ins_content = common_ins_block_dic[ll_ins_id_mk][0][2]
                                tem_ins_length = len(tem_ins_content)
                                Best_ref_ins[kk][ll_ins_id_mk] = tem_ins_content
                                if tem_ins_length in tem_ins_len_dic:
                                    tem_ins_len_dic[tem_ins_length] = tem_ins_len_dic[tem_ins_length] + 1
                                else:
                                    tem_ins_len_dic[tem_ins_length] = 1
                            
                            ### full ins part
                            full_final_compare_ref_result_df_ins_full = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,6]>-1]
                            for ll_ins_id in range(len(full_final_compare_ref_result_df_ins_full)):
                                tem_ins_loc = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,6]
                                tem_ins_content = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,7]
                                tem_ins_length = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,8]
                                Best_ref_ins_full[kk][tem_ins_loc] = tem_ins_content
                            
                            ##Del
                            full_final_compare_ref_result_df_del = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,9]>-1]
                            for ll_del_id in range(len(full_final_compare_ref_result_df_del)):
                                #tem_del_loc = full_final_compare_ref_result_df_del.iloc[ll_del_id,3]
                                #tem_del_content = full_final_compare_ref_result_df_del.iloc[ll_del_id,4]
                                tem_del_length = full_final_compare_ref_result_df_del.iloc[ll_del_id,9]
                                #Best_ref_ins[kk][tem_ins_loc] = tem_ins_content
                                if tem_del_length in tem_del_len_dic:
                                    tem_del_len_dic[tem_del_length] = tem_del_len_dic[tem_del_length] + 1
                                else:
                                    tem_del_len_dic[tem_del_length] = 1
                            
                            collection_id_del_covered = list(set(full_final_compare_ref_result_df.iloc[:,12]))
                            collection_id_del_covered.remove(-1)
                            remaining_del_id = list(set(range(len(common_del_block_df_cp)))-set(collection_id_del_covered))
                            
                            for ll_remain_del_id in remaining_del_id:
                                tem_del_length = common_del_block_df_cp.iloc[ll_remain_del_id,1] - common_del_block_df_cp.iloc[ll_remain_del_id,0] + 1
                                if tem_del_length>0:
                                    if tem_del_length not in tem_del_len_dic:
                                        tem_del_len_dic[tem_del_length] = 1
                                    else:
                                        tem_del_len_dic[tem_del_length] = tem_del_len_dic[tem_del_length] + 1
                            
                            
            
                            con_ref_prob_full[kk,5] =  con_ref_prob_full[kk,0] + con_ref_prob_full[kk,1]
                            for del_key in tem_del_len_dic:
                                con_ref_prob_full[kk,2] = con_ref_prob_full[kk,2] + tem_del_len_dic[del_key]
                                con_ref_prob_full[kk,5] = con_ref_prob_full[kk,5] + tem_del_len_dic[del_key]
                               
                            for ins_key in tem_ins_len_dic:
                                con_ref_prob_full[kk,3] = con_ref_prob_full[kk,3] + tem_ins_len_dic[ins_key]
            
                            con_ref_prob_full[kk,4] = con_ref_prob_full[kk,5] - con_ref_prob_full[kk,2]-len(tem_del_len_dic) -1
                            
                            ### update of the bone
                            full_final_compare_ref_result_df_bone = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,10]>-1]
                            for bone_id_df in range(len(full_final_compare_ref_result_df_bone)):
                                bone_num_tem = full_final_compare_ref_result_df_bone.iloc[bone_id_df,10]
                                bone_base_infer = full_final_compare_ref_result_df_bone.iloc[bone_id_df,11]
                                Best_ref[kk][bone_num_tem] = bone_base_infer
                            
                            tem_eq_num = con_ref_prob_full[kk, 0]
                            tem_neq_num = con_ref_prob_full[kk, 1]
                            tem_del_loc_num = con_ref_prob_full[kk, 2]
                            tem_del_loc_num_full = con_ref_prob_full[kk,5] 
                            tem_del_loc_prob = tem_del_loc_num/tem_del_loc_num_full
                            tem_ins_num = con_ref_prob_full[kk, 3]
                            tem_no_ins_num = con_ref_prob_full[kk, 4]
                            
                            up_gamma4[kk] = (tem_ins_num) / \
                                (tem_ins_num+tem_no_ins_num)
                            up_gamma1[kk] = tem_eq_num / \
                                (tem_eq_num+tem_neq_num)
                            up_gamma3[kk] = tem_del_loc_prob
                            up_gamma2[kk] = 1-up_gamma1[kk]
                        
                            if len(tem_del_len_dic)>0:
                                tem_dep_c = 0
                                for dd_k in  tem_del_len_dic:
                                    up_del_beta[kk] = up_del_beta[kk] + (dd_k/(1-(1-del_p[kk])*np.exp(-del_beta[kk]*dd_k)))*tem_del_len_dic[dd_k]
                                    tem_dep_c = tem_dep_c + 1/(1-(1-del_p[kk])*np.exp(-del_beta[kk]*dd_k))*tem_del_len_dic[dd_k]
                                import_data_p = (tem_dep_c,tem_del_loc_num)
                                p_00 = 0
                                up_del_p[kk] = fsolve(extimation_p, p_00,args=import_data_p)
                               
                                if up_del_p[kk]>1:
                                    up_del_p[kk] = 1
                                if up_del_p[kk]<0:
                                    up_del_p[kk] = 0
                                #up_del_p[kk] = -tem_del_loc_num*(1-del_p[kk])/(up_del_p[kk] * np.log(del_p[kk]))
                                
                                
                                up_del_beta[kk] = tem_del_loc_num /up_del_beta[kk]
                            else:
                                up_del_beta[kk] = 0.5
                                up_del_p[kk] = 0.5
                            
                            if len(tem_ins_len_dic)>0:
                                tem_dep_c = 0
                                for ii_k in  tem_ins_len_dic:
                                    up_ins_beta[kk] = up_ins_beta[kk] + (ii_k/(1-(1-ins_p[kk])*np.exp(-ins_beta[kk]*ii_k)))*tem_ins_len_dic[ii_k]
                                    #up_ins_p[kk] = up_ins_p[kk] + 1/(1-(1-ins_p[kk])*np.exp(-ins_beta[kk]*ii_k))*tem_ins_len_dic[ii_k]
                                    tem_dep_c = tem_dep_c + 1/(1-(1-ins_p[kk])*np.exp(-ins_beta[kk]*ii_k))*tem_ins_len_dic[ii_k]
                                import_data_p = (tem_dep_c,tem_ins_num)
                                p_00 = 0
                                up_ins_p[kk] = fsolve(extimation_p, p_00,args=import_data_p)
                                if up_ins_p[kk]<0:
                                    up_ins_p[kk] = 0
                                if up_ins_p[kk]>1:
                                    up_ins_p[kk] = 1
                               
                                #up_ins_p[kk] = -tem_ins_num*(1-ins_p[kk])/(up_ins_p[kk] * np.log(ins_p[kk]))
                                up_ins_beta[kk] = tem_ins_num /up_ins_beta[kk]
                            else:
                                up_ins_beta[kk] = 0.5
                                up_ins_p[kk] = 0.5
                            
                            #up_ins_p[kk] = -tem_ins_num*(1-up_ins_p[kk])/(up_ins_p[kk] * np.log(ins_p[kk]))
                            #up_ins_beta[kk] = tem_ins_num /up_ins_beta[kk]
                            #up_beta_a = tem_eq_num + 48
                            #up_beta_b = tem_neq_num + 4
                            #up_gamma[kk] =  np.random.beta(a = up_beta_a,b = up_beta_b,size=1)
                            #up_gamma[kk] =  tem_eq_num/(tem_neq_num+tem_eq_num)
                            if up_gamma1[kk] < 0.9:
                                up_gamma1[kk] = 0.9
                                up_gamma2[kk] = 1-up_gamma1[kk]
                            
                            for ll_ins in ins_inserted_dic:
                                if ll_ins not in Best_ref_ins_full[kk]:
                                    Best_ref_ins_full[kk][ll_ins] = '5'*ins_inserted_dic[ll_ins]
                             
                               
                        else:
                            Best_ref[kk] = {}
                            Best_ref_ins[kk] = {} 
                            Best_ref_ins_keys[kk] = []
                            Best_ref_ins_full[kk] = {}
                            for basss_index in only_bone_note_df_keys:
                                Best_ref[kk][int(basss_index)] = copy.deepcopy(
                                    Basic_ref[int(basss_index)])
                            con_ref_prob_full[kk, :] = [n_c, 0, 0, 0, n_c-1,n_c]
                            up_gamma1[kk] = np.random.beta(a=48, b=4, size=1)
                            up_gamma2[kk] = 1-up_gamma1[kk]
                            up_gamma3[kk] = np.random.rand(1)*0.01
                            up_gamma4[kk] = np.random.rand(1)*0.01
                            
                            up_del_beta[kk] = 0.5
                            up_del_p[kk] = 0.5
                            up_ins_beta[kk] = 0.5
                            up_ins_p[kk] = 0.5
                            for ll_ins in ins_inserted_dic:
                                Best_ref_ins_full[kk][ll_ins] = '5'*ins_inserted_dic[ll_ins]
                            
                        
                        evo_ins_len_dic[kk] = copy.deepcopy(tem_ins_len_dic)
                        evo_del_len_dic[kk] = copy.deepcopy(tem_del_len_dic)
                            
                        Best_ref_ins_keys[kk] = np.array(list(Best_ref_ins[kk].keys()))
                        
                    portions = portions + 0.001
                    portions = portions/sum(portions)
                    check_dis2 = abs(last_portions-portions).max()
            
                    # parameters
                    conditional_probalility_ = np.zeros([max_num_reads_refine, K])
                    conditional_probalility_o = np.zeros([max_num_reads_refine, K])
                    update_fre = [np.zeros([n_c_refine, 6]) for x in range(K)]
                    wei_del = 0
                    wei_ins = 0
                    whole_ins_dic = {}
                    whole_del_dic = {}
            
                    with  mp.Pool(num_cores) as pool:
                
                     
                         final_compare_con_result = pool.map(collection_var_per_read, list(range(max_num_reads_refine)))
                    
                    ### making up the content
                    final_compare_con_result_full = []
                    for ll_per_read_id in range(max_num_reads_refine):
                        final_compare_con_result_full = final_compare_con_result_full + final_compare_con_result[ll_per_read_id]
                    
                    final_compare_con_result_full_df = pd.DataFrame(final_compare_con_result_full)
                    final_compare_con_result_df_sub = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,-1]>-1]
                    
                    reads_infer_clu.iloc[:,1] = final_compare_con_result_df_sub.iloc[:,-1]
                    
                    ### correct_num
                    final_compare_con_result_df_sub_correct = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,1]>-1]
                    wei_cor = sum(final_compare_con_result_df_sub_correct.iloc[:,1])
                    
                    ### mis num
                    final_compare_con_result_df_sub_mis = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,2]>-1]
                    wei_mis = sum(final_compare_con_result_df_sub_mis.iloc[:,2])
                    
                    ### del
                    final_compare_con_result_df_sub_del = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,3]>-1]
                    for del_hang_check in range(len(final_compare_con_result_df_sub_del)):
                        tem_len_key_del = final_compare_con_result_df_sub_del.iloc[del_hang_check,3]
                        tem_len_num_del = final_compare_con_result_df_sub_del.iloc[del_hang_check,4]
                        if  tem_len_key_del in whole_del_dic:
                            whole_del_dic[tem_len_key_del] = whole_del_dic[tem_len_key_del] + tem_len_num_del
                        else:
                            whole_del_dic[tem_len_key_del] = tem_len_num_del
                        wei_del = wei_del + tem_len_num_del
                    
                    ### ins
                    final_compare_con_result_df_sub_ins = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,5]>-1]
                    for ins_hang_check in range(len(final_compare_con_result_df_sub_ins)):
                        tem_len_key_ins = final_compare_con_result_df_sub_ins.iloc[ins_hang_check,5]
                        tem_len_num_ins = final_compare_con_result_df_sub_ins.iloc[ins_hang_check,6]
                        if  tem_len_key_ins in whole_ins_dic:
                            whole_ins_dic[tem_len_key_ins] = whole_ins_dic[tem_len_key_ins] + tem_len_num_ins
                        else:
                            whole_ins_dic[tem_len_key_ins] = tem_len_num_ins
                        wei_ins = wei_ins + tem_len_num_ins
                    
                    ### prob
                    final_compare_con_result_df_sub_prob = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,0]>-1]
                    for ll_read_per in range(len(final_compare_con_result_df_sub_prob)):
                        conditional_probalility_[ll_read_per,:] = final_compare_con_result_df_sub_prob.iloc[ll_read_per,7:(7+K)]
                    
            
                    
                    cum_del = 0
                    for kk_del in whole_del_dic:
                        cum_del = cum_del+whole_del_dic[kk_del]*kk_del
            
                    cum_ins = 0
                    for kk_ins in whole_ins_dic:
                        cum_ins = cum_ins+whole_ins_dic[kk_ins] * kk_ins
                    up_lamuna_del = cum_del/(wei_del+corretion_eps_from_zero) -1 
                    up_lamuna_ins = cum_ins/(wei_ins+corretion_eps_from_zero) -1
                    if up_lamuna_del<0:
                        up_lamuna_del = corretion_eps_from_zero
                    if up_lamuna_ins<0:
                        up_lamuna_ins = corretion_eps_from_zero
                    
                    up_pcorrect = wei_cor/(wei_cor+wei_mis)
                    up_p_mis_loc = wei_mis/(wei_cor+wei_mis)
                    up_pins1_loc = wei_ins/(wei_cor+wei_mis)
                    up_pdel1_loc = wei_del/(wei_cor+wei_mis+cum_del)
                    
                    Differ_matrix = np.zeros([K,K],dtype=int)
                    for i_last in range(K):
                        for j_last in range(K):
                            differ_num_tem_array = abs(last_informative_sites_code[i_last]-informative_sites_code[j_last])
                            differ_num_ = len(differ_num_tem_array[differ_num_tem_array<1])
                            Differ_matrix[i_last,j_last] = differ_num_
                    
                    row_o3_mec_itera, col_o3_mec_itera = linear_sum_assignment(Differ_matrix)
                    differ_num_ =  Differ_matrix[row_o3_mec_itera,col_o3_mec_itera].sum()
                    if differ_num_/num_remain_infor_base<portions_var_fail:
                        check_dis1 = creatierstop/2
                    else:
                        check_dis1 = 1
                    
                    
                    
                    check_dis1 = abs(last_probability -
                                     conditional_probalility_).max()
                    check_dis3 = max(abs(up_lamuna_del-lamuna_del), abs(lamuna_ins-up_lamuna_ins),
                                     abs(up_pcorrect -
                                         pcorrect), abs(up_p_mis_loc-p_mis_loc),
                                     abs(up_pins1_loc -
                                         pins1_loc), abs(up_pdel1_loc-pdel1_loc),
                                     abs(up_gamma1 -
                                         gamma1).max(), abs(up_gamma2-gamma2).max(),
                                     abs(up_gamma3-gamma3).max(), abs(up_gamma4-gamma4).max(),
                                     abs(up_del_beta-del_beta).max(),abs(up_del_p-del_p).max(),
                                     abs(up_ins_beta-ins_beta).max(),abs(up_ins_p-ins_p).max())
            
                    # check_dis3 = max(abs(up_mu_del-mu_del),abs(up_sig_del-sig_del),\
                    #                 abs(up_mu_ins-mu_ins),abs(up_sig_ins-sig_ins),\
                    #                 abs(up_pcorrect-pcorrect),abs(up_p_mis_loc-p_mis_loc),\
                    #                abs(up_pins1_loc-pins1_loc),abs(up_pdel1_loc-pdel1_loc))
                    check_dis = max(check_dis1, check_dis2, check_dis3)
                    max_iter = 10
                    min_iter = 50
                    
                    loop_ = 0
                    check_dis_flag = 1
                    time_vec = []
                    
                    ## initial end itera
                    t_pre_now = time.time()
                    time_save_list.append(['Time for the first iteration', t_pre_now-t_pre])
                    
                    while check_dis > creatierstop:
                        time_vec.append(time.time())
                        if loop_ < max_iter:
                            t_pre = copy.deepcopy(t_pre_now)
                            tao = tao * tao_step
                            #tao2 = tao2 * tao_step2
                            last_probability = copy.deepcopy(
                                conditional_probalility_)
                            last_fre = copy.deepcopy(update_fre)
                            conditional_probalility_ = np.zeros(
                                [max_num_reads_refine, K])
                            update_fre = [np.zeros([n_c_refine, 6])
                                          for x in range(K)]
                            last_portions = copy.deepcopy(portions)
                            p_mis_loc = copy.deepcopy(up_p_mis_loc)
                            pins1_loc = copy.deepcopy(up_pins1_loc)
                            lamuna_ins = copy.deepcopy(up_lamuna_ins)
                            #mu_ins = copy.deepcopy(up_mu_ins)
                            #sig_ins = copy.deepcopy(up_sig_ins)
                            gamma1 = copy.deepcopy(up_gamma1)
                            gamma2 = copy.deepcopy(up_gamma2)
                            gamma3 = copy.deepcopy(up_gamma3)
                            gamma4 = copy.deepcopy(up_gamma4)
                            del_beta = copy.deepcopy(up_del_beta)  # del
                            del_p = copy.deepcopy(up_del_p)
                            ins_beta= copy.deepcopy(up_ins_beta)  # del
                            ins_p = copy.deepcopy(up_ins_p)
                            
                            pdel1_loc = copy.deepcopy(up_pdel1_loc)
                            lamuna_del = copy.deepcopy(up_lamuna_del)
                            #mu_del = copy.deepcopy(up_mu_del)
                            #sig_del = copy.deepcopy(up_sig_del)
                            pcorrect = copy.deepcopy(up_pcorrect)
                            theta_ = [pcorrect, p_mis_loc, pins1_loc,
                                      pdel1_loc, lamuna_ins, lamuna_del]
                            con_ref_prob = np.zeros(K)
                            #Best_ref = [{} for x in range(K)]
                            #Best_ref_ins = [{} for x in range(K)]
                            #Best_ref_ins_keys = [[] for x in range(K)]
                            #Best_ref_ins_full = [{} for x in range(K)]
                            whole_ins_dic = {}
                            whole_del_dic = {}
                            wei_del = 0
                            wei_ins = 0
                            last_informative_sites_code = copy.deepcopy(informative_sites_code)
                            informative_sites_code = [contant_ref_code_var for x in range(K)]
                            #covering_sites_msa = [[] for x in range(K)]
                            #covering_sites_bone = [[] for x in range(K)]
                            #covering_sites_ins = [[] for x in range(K)]
                            for kkk in range(K):
                                
                                tem_gamma_ = np.array([gamma1[kkk], gamma2[kkk], gamma3[kkk], gamma4[kkk], 1-gamma4[kkk]])
                                tem_gamma_ = tem_gamma_ + corretion_eps_from_zero 
                                tem_gamma_[:2] = tem_gamma_[:2]/sum(tem_gamma_[:2])
                                tem_gamma_[2] = min(1,tem_gamma_[2])
                                tem_gamma_[3:] = tem_gamma_[3:]/sum(tem_gamma_[3:])
                                con_ref_prob[kkk] = con_ref_prob[kkk] + con_ref_prob_full[kkk, 0] * np.log(tem_gamma_[0])### correct prob
                                con_ref_prob[kkk] = con_ref_prob[kkk] + con_ref_prob_full[kkk, 1] * np.log(tem_gamma_[1])### mis prob
                                con_ref_prob[kkk] = con_ref_prob[kkk] + (con_ref_prob_full[kkk, 5]) * np.log(1-tem_gamma_[2])### not del prob
                                
                                con_ref_prob[kkk] = con_ref_prob[kkk] + (con_ref_prob_full[kkk, 2]) * np.log(tem_gamma_[2])### not del prob
                                
                                for del_dic_key in evo_del_len_dic[kkk]:
                                    con_ref_prob[kkk] = con_ref_prob[kkk] + evo_del_len_dic[kkk][del_dic_key]*expotential_log_log_probability(del_dic_key,del_p[kkk],del_beta[kkk])
                                
                                ### ins part 
                                con_ref_prob[kkk] = con_ref_prob[kkk] + con_ref_prob_full[kkk, 3]*np.log(tem_gamma_[3])
                                for ins_dic_key in evo_ins_len_dic[kkk]:
                                    con_ref_prob[kkk] = con_ref_prob[kkk] + evo_ins_len_dic[kkk][ins_dic_key]*expotential_log_log_probability(ins_dic_key,ins_p[kkk],ins_beta[kkk])
                                con_ref_prob[kkk] = con_ref_prob[kkk] +  con_ref_prob_full[kkk, 4]*np.log(tem_gamma_[4])
                                
                                
                                
                            # liklihood_.append(((conditional_probalility_o*conditional_probalility_).sum()+sum(con_ref_prob)))
            
                            #creatierstop = 0.0001
                            #correct_num= np.zeros(K)
                            up_gamma1 = np.ones(K)*0.9
                            up_gamma2 = np.ones(K)*0.05
                            up_gamma3 = np.ones(K)*0.05
                            up_del_beta = np.zeros(K)
                            up_del_p = np.zeros(K)
                            up_gamma4 = np.ones(K)*0.0001
                            up_ins_beta = np.zeros(K)
                            up_ins_p = np.zeros(K)
                            #up_gamma5=  np.ones(K)*0.05
                            present_clu_dis2 = [
                                np.zeros([num_remain_infor_base, 6]) for x in range(K)]
                            con_ref_prob_full = np.zeros([K, 6])
                            
                            evo_ins_len_dic = {}### this is to save the dic of the length of the ins length and del length
                            evo_del_len_dic = {}
                            for kk in range(K):
                                ele_list = list(reads_infer_clu[reads_infer_clu.iloc[:, 1] == kk].iloc[:, 0])
                                #ele_list = list(
                                #    reads_infer_cluf[reads_infer_cluf.iloc[:, 1] == kk].iloc[:, 0])
                                portions[kk] = len(ele_list)
                                #tem_eq_num = 0
                                #tem_neq_num = 0
                                tem_del_len_dic = {}
                                tem_ins_len_dic = {}
                                if len(ele_list)>0:
                                    for ele in ele_list:
                                        #tem_ele_start = final_full_result[ele][0]
                                        #tem_ele_end = final_full_result[ele][1]
                                        tem_matrix_refine = read_full_bases_filter[ele]
                                        if len(tem_matrix_refine)>0:
                                            tem_start_refine = infer_var_site.index(tem_matrix_refine[0,0])
                                            tem_end_refine = tem_start_refine + len(tem_matrix_refine)-1
                                            present_clu_dis2[kk][tem_start_refine:tem_end_refine+1,:] = present_clu_dis2[kk][tem_start_refine:tem_end_refine+1,:] + \
                                                tem_matrix_refine[:,1:]
                                    #max_ID = list(present_clu_dis2[kk][int(basss_index),:]).index(max(present_clu_dis2[kk][int(basss_index),:]))
                                    
                                    common_del_block_df_cp = copy.deepcopy(common_del_block_df)
                                    keys_ins_common_original_cp = copy.deepcopy(keys_ins_common_original)
                                    tem_touched_ins = []
                                    
                                    ### using the common part of the non-indel as the independent part
                                    remaining_independent = copy.deepcopy(seg_collection_ini_df_dep_id)
                                    with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                    
                                        Result_infer_ind = pool.starmap(signle_sites_infer_max4, [(remaining_independent[index], present_clu_dis2[kk][int(remaining_independent[index]),:], standard_ref_code[infer_var_site[int(remaining_independent[index])]],kk
                                                                                               ) for index in range(len(remaining_independent))])
                                    remaining_dependent = copy.deepcopy(seg_collection_ini_df_indel_part_cp_sub_link_full)
                                    related_commom_del = copy.deepcopy(dep_spe_infor_del_col_full)
                                    related_commom_ins = copy.deepcopy(dep_spe_infor_ins_col_full)
                                    ## next we check the dependent part / we should spilt find the read overlapping with the region and common indel
                                    with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                                        Dep_dic = pool.starmap(select_read_blocks2, [(
                                        remaining_dependent.iloc[index,0], remaining_dependent.iloc[index,1]) for index in range(len(remaining_dependent))])
                                        #Dep_dic = Dep_dic_iter.get()
                                    with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                                        Ref_dic_block = pool.starmap(select_ref_block2, [(
                                         remaining_dependent.iloc[index,0], remaining_dependent.iloc[index,1],kk) for index in range(len(remaining_dependent))])
                                    with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                                        Result_infer_dep = pool.starmap(block_sites_infer_ECM3, [(int(remaining_dependent.iloc[index,0]), int(remaining_dependent.iloc[index,1]), present_clu_dis2[kk][int(remaining_dependent.iloc[index,0]):int(remaining_dependent.iloc[index,1])+1, :],
                                                                                              final_ref_code_standard[int(remaining_dependent.iloc[index,0]):int(
                                                                                                  remaining_dependent.iloc[index,1])+1], copy.deepcopy(Ref_dic_block[index]), Dep_dic[index],related_commom_del[index],related_commom_ins[index],kk) for index in range(len(remaining_dependent))])
                                        
                                    Best_ref[kk] = {}
                                    Best_ref_ins[kk] = {} 
                                    Best_ref_ins_keys[kk] = []
                                    Best_ref_ins_full[kk] = {}
                                    
                                    # making up full result
                                    Result_infer_ = copy.deepcopy(Result_infer_ind)
                                    #dep_Result_infer_ = []
                                    for ll_xun_ in Result_infer_dep:
                                        Result_infer_ = Result_infer_ + ll_xun_### this is to make the reserved parts outsides from signs
                                        #dep_Result_infer_ = dep_Result_infer_ + ll_xun_
                                    Result_infer_ = sorted(Result_infer_)
                                    Result_infer_df = pd.DataFrame(Result_infer_)
                                    informative_sites_code[kk] = np.array(Result_infer_df.iloc[:,1])
                                    
                                    
                                    #for ll in range(len(seg_collection_ini)):
                                    #    if seg_collection_ini[ll][1]==6439:
                                    #        print(seg_collection_ini[ll])
                                    #        print(ll)
                                    
                                    #for ll2 in range(len(seg_collection_ini_df_indel_part_cp_sub_link_full)):
                                    #    if seg_collection_ini_df_indel_part_cp_sub_link_full.iloc[ll2,1]==6439:
                                    #        print(ll2)
                                    #        print(seg_collection_ini_df_indel_part_cp_sub_link_full.iloc[ll2,:])
                                    
                                    with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:
                                
                                     # Result_infer_ind=pool.starmap(signle_sites_infer_max,[(remaining_independent[index],present_clu_dis2[kk][int(remaining_independent[index]),:],standard_ref_code[int(remaining_independent[index])],\
                                     # gamma[kk],theta_) for index in range(len(remaining_independent))])
                                
                                         final_compare_ref_result = pool.map(comparision_between_con_and_ref, list(range(len(seg_collection_ini))))
                                         #Result_infer_ind =  Result_infer_ind_iter.get()
                                     #final_inferred=pool.starmap(sampling2,[(pp_random[index,0],pp_random[index,1],combination_list2.iloc[int(index/num_chains),0],combination_list2.iloc[int(index/num_chains),1],out_put_df,collection_of_sites,final_ref_code) for index in range((len(combination_list2)*num_chains))])
                                     #pool.close()
                                     #pool.join()
                                    full_final_compare_ref_result = []
                                     
                                    for ll_final_id in range(len(final_compare_ref_result)):
                                         full_final_compare_ref_result = full_final_compare_ref_result + final_compare_ref_result[ll_final_id]
                                     
                                    full_final_compare_ref_result_df = pd.DataFrame(full_final_compare_ref_result)
                                    
                                    con_ref_prob_full[kk,0] = sum(full_final_compare_ref_result_df.iloc[:,1]) + Backgroud_evo_c
                                    con_ref_prob_full[kk,1] = sum(full_final_compare_ref_result_df.iloc[:,2]) + Backgroud_evo_m
                                    #sub_sub_result = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,-1]==1]
                                    
                                    ### ins part
                                    full_final_compare_ref_result_df_ins = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,3]>-1]
                                    for ll_ins_id in range(len(full_final_compare_ref_result_df_ins)):
                                        tem_ins_loc = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,3]
                                        tem_ins_content = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,4]
                                        tem_ins_length = full_final_compare_ref_result_df_ins.iloc[ll_ins_id,5]
                                        Best_ref_ins[kk][tem_ins_loc] = tem_ins_content
                                        if tem_ins_length in tem_ins_len_dic:
                                            tem_ins_len_dic[tem_ins_length] = tem_ins_len_dic[tem_ins_length] + 1
                                        else:
                                            tem_ins_len_dic[tem_ins_length] = 1
                                    
                                    collection_id_ins_covered = list(Best_ref_ins[kk].keys())
                                    remaining_ins_id = list(set(list(common_ins_block_dic.keys()))-set(collection_id_ins_covered))
                                    
                                    for ll_ins_id_mk in remaining_ins_id:
                                        tem_ins_content = common_ins_block_dic[ll_ins_id_mk][0][2]
                                        tem_ins_length = len(tem_ins_content)
                                        Best_ref_ins[kk][ll_ins_id_mk] = tem_ins_content
                                        if tem_ins_length in tem_ins_len_dic:
                                            tem_ins_len_dic[tem_ins_length] = tem_ins_len_dic[tem_ins_length] + 1
                                        else:
                                            tem_ins_len_dic[tem_ins_length] = 1
                                    
                                    ### full ins part
                                    full_final_compare_ref_result_df_ins_full = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,6]>-1]
                                    for ll_ins_id in range(len(full_final_compare_ref_result_df_ins_full)):
                                        tem_ins_loc = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,6]
                                        tem_ins_content = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,7]
                                        tem_ins_length = full_final_compare_ref_result_df_ins_full.iloc[ll_ins_id,8]
                                        Best_ref_ins_full[kk][tem_ins_loc] = tem_ins_content
                                    
                                    ##Del
                                    full_final_compare_ref_result_df_del = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,9]>-1]
                                    for ll_del_id in range(len(full_final_compare_ref_result_df_del)):
                                        #tem_del_loc = full_final_compare_ref_result_df_del.iloc[ll_del_id,3]
                                        #tem_del_content = full_final_compare_ref_result_df_del.iloc[ll_del_id,4]
                                        tem_del_length = full_final_compare_ref_result_df_del.iloc[ll_del_id,9]
                                        #Best_ref_ins[kk][tem_ins_loc] = tem_ins_content
                                        if tem_del_length in tem_del_len_dic:
                                            tem_del_len_dic[tem_del_length] = tem_del_len_dic[tem_del_length] + 1
                                        else:
                                            tem_del_len_dic[tem_del_length] = 1
                                    
                                    collection_id_del_covered = list(set(full_final_compare_ref_result_df.iloc[:,12]))
                                    collection_id_del_covered.remove(-1)
                                    remaining_del_id = list(set(range(len(common_del_block_df_cp)))-set(collection_id_del_covered))
                                    
                                    for ll_remain_del_id in remaining_del_id:
                                        tem_del_length = common_del_block_df_cp.iloc[ll_remain_del_id,1] - common_del_block_df_cp.iloc[ll_remain_del_id,0] + 1
                                        if tem_del_length>0:
                                            if tem_del_length not in tem_del_len_dic:
                                                tem_del_len_dic[tem_del_length] = 1
                                            else:
                                                tem_del_len_dic[tem_del_length] = tem_del_len_dic[tem_del_length] + 1
                                    
                                    
                    
                                    con_ref_prob_full[kk,5] =  con_ref_prob_full[kk,0] + con_ref_prob_full[kk,1]
                                    for del_key in tem_del_len_dic:
                                        con_ref_prob_full[kk,2] = con_ref_prob_full[kk,2] + tem_del_len_dic[del_key]
                                        con_ref_prob_full[kk,5] = con_ref_prob_full[kk,5] + tem_del_len_dic[del_key]
                                       
                                    for ins_key in tem_ins_len_dic:
                                        con_ref_prob_full[kk,3] = con_ref_prob_full[kk,3] + tem_ins_len_dic[ins_key]
                    
                                    con_ref_prob_full[kk,4] = con_ref_prob_full[kk,5] - con_ref_prob_full[kk,2]-len(tem_del_len_dic) -1
                                    
                                    ### update of the bone
                                    full_final_compare_ref_result_df_bone = full_final_compare_ref_result_df[full_final_compare_ref_result_df.iloc[:,10]>-1]
                                    for bone_id_df in range(len(full_final_compare_ref_result_df_bone)):
                                        bone_num_tem = full_final_compare_ref_result_df_bone.iloc[bone_id_df,10]
                                        bone_base_infer = full_final_compare_ref_result_df_bone.iloc[bone_id_df,11]
                                        Best_ref[kk][bone_num_tem] = bone_base_infer
                                    
                                    tem_eq_num = con_ref_prob_full[kk, 0]
                                    tem_neq_num = con_ref_prob_full[kk, 1]
                                    tem_del_loc_num = con_ref_prob_full[kk, 2]
                                    tem_del_loc_num_full = con_ref_prob_full[kk,5] 
                                    tem_del_loc_prob = tem_del_loc_num/tem_del_loc_num_full
                                    tem_ins_num = con_ref_prob_full[kk, 3]
                                    tem_no_ins_num = con_ref_prob_full[kk, 4]
                                    
                                    up_gamma4[kk] = (tem_ins_num) / \
                                        (tem_ins_num+tem_no_ins_num)
                                    up_gamma1[kk] = tem_eq_num / \
                                        (tem_eq_num+tem_neq_num)
                                    up_gamma3[kk] = tem_del_loc_prob
                                    up_gamma2[kk] = 1-up_gamma1[kk]
                                
                                    if len(tem_del_len_dic)>0:
                                        tem_dep_c = 0
                                        for dd_k in  tem_del_len_dic:
                                            up_del_beta[kk] = up_del_beta[kk] + (dd_k/(1-(1-del_p[kk])*np.exp(-del_beta[kk]*dd_k)))*tem_del_len_dic[dd_k]
                                            tem_dep_c = tem_dep_c + 1/(1-(1-del_p[kk])*np.exp(-del_beta[kk]*dd_k))*tem_del_len_dic[dd_k]
                                        import_data_p = (tem_dep_c,tem_del_loc_num)
                                        p_00 = 0
                                        up_del_p[kk] = fsolve(extimation_p, p_00,args=import_data_p)
                                       
                                        if up_del_p[kk]>1:
                                            up_del_p[kk] = 1
                                        if up_del_p[kk]<0:
                                            up_del_p[kk] = 0
                                        #up_del_p[kk] = -tem_del_loc_num*(1-del_p[kk])/(up_del_p[kk] * np.log(del_p[kk]))
                                        
                                        
                                        up_del_beta[kk] = tem_del_loc_num /up_del_beta[kk]
                                    else:
                                        up_del_beta[kk] = 0.5
                                        up_del_p[kk] = 0.5
                                    
                                    if len(tem_ins_len_dic)>0:
                                        tem_dep_c = 0
                                        for ii_k in  tem_ins_len_dic:
                                            up_ins_beta[kk] = up_ins_beta[kk] + (ii_k/(1-(1-ins_p[kk])*np.exp(-ins_beta[kk]*ii_k)))*tem_ins_len_dic[ii_k]
                                            #up_ins_p[kk] = up_ins_p[kk] + 1/(1-(1-ins_p[kk])*np.exp(-ins_beta[kk]*ii_k))*tem_ins_len_dic[ii_k]
                                            tem_dep_c = tem_dep_c + 1/(1-(1-ins_p[kk])*np.exp(-ins_beta[kk]*ii_k))*tem_ins_len_dic[ii_k]
                                        import_data_p = (tem_dep_c,tem_ins_num)
                                        p_00 = 0
                                        up_ins_p[kk] = fsolve(extimation_p, p_00,args=import_data_p)
                                        if up_ins_p[kk]<0:
                                            up_ins_p[kk] = 0
                                        if up_ins_p[kk]>1:
                                            up_ins_p[kk] = 1
                                       
                                        #up_ins_p[kk] = -tem_ins_num*(1-ins_p[kk])/(up_ins_p[kk] * np.log(ins_p[kk]))
                                        up_ins_beta[kk] = tem_ins_num /up_ins_beta[kk]
                                    else:
                                        up_ins_beta[kk] = 0.5
                                        up_ins_p[kk] = 0.5
                                    
                                    #up_ins_p[kk] = -tem_ins_num*(1-up_ins_p[kk])/(up_ins_p[kk] * np.log(ins_p[kk]))
                                    #up_ins_beta[kk] = tem_ins_num /up_ins_beta[kk]
                                    #up_beta_a = tem_eq_num + 48
                                    #up_beta_b = tem_neq_num + 4
                                    #up_gamma[kk] =  np.random.beta(a = up_beta_a,b = up_beta_b,size=1)
                                    #up_gamma[kk] =  tem_eq_num/(tem_neq_num+tem_eq_num)
                                    if up_gamma1[kk] < 0.9:
                                        up_gamma1[kk] = 0.9
                                        up_gamma2[kk] = 1-up_gamma1[kk]
                                    
                                    for ll_ins in ins_inserted_dic:
                                        if ll_ins not in Best_ref_ins_full[kk]:
                                            Best_ref_ins_full[kk][ll_ins] = '5'*ins_inserted_dic[ll_ins]
                                        
                                       
                                else:
                                    Best_ref[kk] = {}
                                    Best_ref_ins[kk] = {} 
                                    Best_ref_ins_keys[kk] = []
                                    Best_ref_ins_full[kk] = {}
                                    
                                    for basss_index in only_bone_note_df_keys:
                                        Best_ref[kk][int(basss_index)] = copy.deepcopy(
                                            Basic_ref[int(basss_index)])
                                    con_ref_prob_full[kk, :] = [n_c, 0, 0, 0, n_c-1,n_c]
                                    up_gamma1[kk] = np.random.beta(a=48, b=4, size=1)
                                    up_gamma2[kk] = 1-up_gamma1[kk]
                                    up_gamma3[kk] = np.random.rand(1)*0.01
                                    up_gamma4[kk] = np.random.rand(1)*0.01
                                    
                                    up_del_beta[kk] = 0.5
                                    up_del_p[kk] = 0.5
                                    up_ins_beta[kk] = 0.5
                                    up_ins_p[kk] = 0.5
                                    for ll_ins in ins_inserted_dic:
                                        Best_ref_ins_full[kk][ll_ins] = '5'*ins_inserted_dic[ll_ins]
                                    
                                    
                                
                                evo_ins_len_dic[kk] = copy.deepcopy(tem_ins_len_dic)
                                evo_del_len_dic[kk] = copy.deepcopy(tem_del_len_dic)
                                    
                                Best_ref_ins_keys[kk] = np.array(list(Best_ref_ins[kk].keys()))
                                
            
                            portions = portions + 0.001
                            portions = portions/sum(portions)
                            check_dis2 = abs(last_portions-portions).max()
                    
                            # parameters
                            conditional_probalility_ = np.zeros([max_num_reads_refine, K])
                            conditional_probalility_o = np.zeros([max_num_reads_refine, K])
                            update_fre = [np.zeros([n_c_refine, 6]) for x in range(K)]
                            wei_del = 0
                            wei_ins = 0
                            whole_ins_dic = {}
                            whole_del_dic = {}
                    
                            with  mp.Pool(num_cores) as pool:
                        
                             
                                 final_compare_con_result = pool.map(collection_var_per_read, list(range(max_num_reads_refine)))
                            
                            ### making up the content
                            final_compare_con_result_full = []
                            for ll_per_read_id in range(max_num_reads_refine):
                                final_compare_con_result_full = final_compare_con_result_full + final_compare_con_result[ll_per_read_id]
                            
                            final_compare_con_result_full_df = pd.DataFrame(final_compare_con_result_full)
                            final_compare_con_result_df_sub = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,-1]>-1]
                            
                            reads_infer_clu.iloc[:,1] = final_compare_con_result_df_sub.iloc[:,-1]
                            
                            ### correct_num
                            final_compare_con_result_df_sub_correct = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,1]>-1]
                            wei_cor = sum(final_compare_con_result_df_sub_correct.iloc[:,1])
                            
                            ### mis num
                            final_compare_con_result_df_sub_mis = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,2]>-1]
                            wei_mis = sum(final_compare_con_result_df_sub_mis.iloc[:,2])
                            
                            ### del
                            final_compare_con_result_df_sub_del = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,3]>-1]
                            for del_hang_check in range(len(final_compare_con_result_df_sub_del)):
                                tem_len_key_del = final_compare_con_result_df_sub_del.iloc[del_hang_check,3]
                                tem_len_num_del = final_compare_con_result_df_sub_del.iloc[del_hang_check,4]
                                if  tem_len_key_del in whole_del_dic:
                                    whole_del_dic[tem_len_key_del] = whole_del_dic[tem_len_key_del] + tem_len_num_del
                                else:
                                    whole_del_dic[tem_len_key_del] = tem_len_num_del
                                wei_del = wei_del + tem_len_num_del
                            
                            ### ins
                            final_compare_con_result_df_sub_ins = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,5]>-1]
                            for ins_hang_check in range(len(final_compare_con_result_df_sub_ins)):
                                tem_len_key_ins = final_compare_con_result_df_sub_ins.iloc[ins_hang_check,5]
                                tem_len_num_ins = final_compare_con_result_df_sub_ins.iloc[ins_hang_check,6]
                                if  tem_len_key_ins in whole_ins_dic:
                                    whole_ins_dic[tem_len_key_ins] = whole_ins_dic[tem_len_key_ins] + tem_len_num_ins
                                else:
                                    whole_ins_dic[tem_len_key_ins] = tem_len_num_ins
                                wei_ins = wei_ins + tem_len_num_ins
                            
                            ### prob
                            final_compare_con_result_df_sub_prob = final_compare_con_result_full_df[final_compare_con_result_full_df.iloc[:,0]>-1]
                            for ll_read_per in range(len(final_compare_con_result_df_sub_prob)):
                                conditional_probalility_[ll_read_per,:] = final_compare_con_result_df_sub_prob.iloc[ll_read_per,7:(7+K)]
                            
                    
                            
                            cum_del = 0
                            for kk_del in whole_del_dic:
                                cum_del = cum_del+whole_del_dic[kk_del]*kk_del
                    
                            cum_ins = 0
                            for kk_ins in whole_ins_dic:
                                cum_ins = cum_ins+whole_ins_dic[kk_ins] * kk_ins
                            up_lamuna_del = cum_del/(wei_del+corretion_eps_from_zero) -1 
                            
                            up_lamuna_ins = cum_ins/(wei_ins+corretion_eps_from_zero) -1
                            if up_lamuna_del<0:
                                up_lamuna_del = corretion_eps_from_zero
                            if up_lamuna_ins<0:
                                up_lamuna_ins = corretion_eps_from_zero
                            
                            up_pcorrect = wei_cor/(wei_cor+wei_mis)
                            up_p_mis_loc = wei_mis/(wei_cor+wei_mis)
                            up_pins1_loc = wei_ins/(wei_cor+wei_mis)
                            up_pdel1_loc = wei_del/(wei_cor+wei_mis+cum_del)
                            
                            ### estimation of the consensus difference between iterations
                            #check_dis1 = 1
                            Differ_matrix = np.zeros([K,K],dtype=int)
                            for i_last in range(K):
                                for j_last in range(K):
                                    differ_num_tem_array = abs(last_informative_sites_code[i_last]-informative_sites_code[j_last])
                                    differ_num_ = len(differ_num_tem_array[differ_num_tem_array>0])
                                    Differ_matrix[i_last,j_last] = differ_num_
                            
                            row_o3_mec_itera, col_o3_mec_itera = linear_sum_assignment(Differ_matrix)
                            differ_num_ =  Differ_matrix[row_o3_mec_itera,col_o3_mec_itera].sum()
                            differ_num_ratio = differ_num_/(K*num_remain_infor_base)
                            if differ_num_ratio<portions_var_fail:
                                check_dis1 = creatierstop/2
                            else:
                                check_dis1 = 1
                            print("difference ratio of con in iterations", differ_num_ratio)    
                            
                            #check_dis1 = abs(last_probability -
                            #                 conditional_probalility_).max()
                            check_dis3 = max(abs(up_lamuna_del-lamuna_del), abs(lamuna_ins-up_lamuna_ins),
                                             abs(up_pcorrect -
                                                 pcorrect), abs(up_p_mis_loc-p_mis_loc),
                                             abs(up_pins1_loc -
                                                 pins1_loc), abs(up_pdel1_loc-pdel1_loc),
                                             abs(up_gamma1 -
                                                 gamma1).max(), abs(up_gamma2-gamma2).max(),
                                             abs(up_gamma3-gamma3).max(), abs(up_gamma4-gamma4).max(),
                                             abs(up_del_beta-del_beta).max(),abs(up_del_p-del_p).max(),
                                             abs(up_ins_beta-ins_beta).max(),abs(up_ins_p-ins_p).max())
                    
                            # check_dis3 = max(abs(up_mu_del-mu_del),abs(up_sig_del-sig_del),\
                            #                 abs(up_mu_ins-mu_ins),abs(up_sig_ins-sig_ins),\
                            #                 abs(up_pcorrect-pcorrect),abs(up_p_mis_loc-p_mis_loc),\
                            #                abs(up_pins1_loc-pins1_loc),abs(up_pdel1_loc-pdel1_loc))
                            check_dis = max(check_dis1, check_dis2, check_dis3)
                            print(check_dis)
                            t_pre_now = time.time()
                            time_save_list.append(['Time for the '+str(loop_)+' iteration', t_pre_now-t_pre])
                            
                            loop_ = loop_ + 1
                            # if  and loop_>min_iter:
                            #    check_dis_flag = 0
                        else:
                            break
            
                    tem_cluster_index = np.argmax(conditional_probalility_, axis=1)
                    Collection_label.append(tem_cluster_index)
                    #conditional_probalility_o_up = copy.deepcopy(conditional_probalility_o)
                    # for kk_xun in range(K):
                    #    conditional_probalility_o_up[:,kk_xun] = conditional_probalility_o_up[:,kk_xun] +
                    sum_of_prob.append(
                        ((conditional_probalility_o*conditional_probalility_).sum()+sum(con_ref_prob)))
                    method_index_c.append(method_index)
                    Best_ref_collection[method_index] = copy.deepcopy(Best_ref)
                    Best_ref_ins_collection[method_index] = copy.deepcopy(Best_ref_ins)
                    parameters_saving[method_index] = [gamma1,gamma2,gamma3,gamma4,del_beta,del_p,\
                                                       ins_beta, ins_p, pcorrect,p_mis_loc,pins1_loc,lamuna_ins,pdel1_loc,lamuna_del]
                # else:
                #    Collection_label.append([])
        ####
        sum_of_prob_max_ = max(sum_of_prob)
        sum_of_prob_max_index = sum_of_prob.index(sum_of_prob_max_)
        method_select_id = method_index_c[int(sum_of_prob_max_index)]
        Selected_Best_ref = Best_ref_collection[method_select_id]
        Selected_Best_ref_ins = Best_ref_ins_collection[method_select_id]
        
        #sum_of_prob_max_index = 0
        #method_select_id = 0
        # finish the estimation
        cluster_index = Collection_label[int(sum_of_prob_max_index)]
        #cluster_index = Label_collection[0]
        #cluster_index = label
        parameters_f = copy.deepcopy(parameters_saving[method_select_id])
        saving_para[:K] = parameters_f[0]
        saving_para[K:2*K] = parameters_f[1]
        saving_para[2*K:3*K] = parameters_f[2]
        saving_para[3*K:4*K] = parameters_f[3]
        saving_para[4*K:5*K] = parameters_f[4]
        saving_para[5*K:6*K] = parameters_f[5]
        saving_para[6*K:7*K] = parameters_f[6]
        saving_para[7*K:8*K] = parameters_f[7]
        saving_para[8*K:] = parameters_f[8:]
        #cluster_index = Result_table.iloc[:,8]
        #cluster_index = labels_initial1
        reads_infer_cluf = pd.DataFrame(np.zeros([max_num_reads_refine, 2]))
        reads_infer_cluf.iloc[:,0] = range(max_num_reads_refine)
        reads_infer_cluf.iloc[:, 1] = copy.deepcopy(cluster_index)
        
        
        para_ = 5+K+len(final_ref_code_standard)*K+K-1
        BICo2 = -2*sum_of_prob_max_ + np.log(max_num_reads_refine)*para_
        AICo2 = -2*sum_of_prob_max_ + 2*para_
        HBICo2 = -2*sum_of_prob_max_ + \
            np.log(max_num_reads_refine)*(5+K+K-1)
        sum_prob_vec = sum_of_prob_max_
        for kk in range(K):
            num_com = len(reads_infer_cluf.iloc[:, 1] == kk)+0.5
            HBICo2 = HBICo2 + \
                np.log(num_com)*len(final_ref_code_standard)
        
        
        infer_var_table =pd.DataFrame(columns=['chrom','pos','type','content'])
        infer_loop = 0
        for s_con in range(proposed_k):
            for keys_bone in  Selected_Best_ref[s_con]:
                if Selected_Best_ref[s_con][keys_bone]!=ref_code_standard[keys_bone]:
                    if Selected_Best_ref[s_con][keys_bone]<4:
                        infer_var_table.loc[infer_loop] = [s_con, keys_bone, 'Mis', Standardbase[Selected_Best_ref[s_con][keys_bone]]]
                        #infer_loop = infer_loop + 1
                    else:
                        infer_var_table.loc[infer_loop] = [s_con, keys_bone, 'Del', '-']
                    infer_loop = infer_loop + 1
            #collection_mis_o = []
            for con_mis_key in Backgroud_evo_mis_suspicious:
                revers_original_bone = reverse_trans2_per(con_mis_key)[1]
                infer_var_table.loc[infer_loop] = [s_con, revers_original_bone, 'Mis', Standardbase[Backgroud_evo_mis_suspicious[con_mis_key]]]
                infer_loop = infer_loop + 1
            for keys_ins in Selected_Best_ref_ins[s_con]:
                    infer_var_table.loc[infer_loop] = [s_con, keys_ins, 'Ins', Selected_Best_ref_ins[s_con][keys_ins]]
                    infer_loop = infer_loop + 1
            for com_del_hang in range(len(common_del_block_df)):
                start_common_del = common_del_block_df.iloc[com_del_hang,2]
                end_common_del = common_del_block_df.iloc[com_del_hang,3]
                for del_base_bone in range(start_common_del,end_common_del+1):
                    infer_var_table.loc[infer_loop] = [s_con, del_base_bone, 'Del', '-']
                    infer_loop = infer_loop + 1
           
            
                
        
       
        # making uop the matrix of output
        # 1. final label
        #id_num = '0'
        #id_num = sys.argv[1]
        Label_matrxi = pd.DataFrame(np.zeros([max_num_reads_refine, 2]))
        #for kk in range(num_k):
        tem_res_collection_df = pd.DataFrame(tem_res_collection)
        Label_matrxi.iloc[:, 0] = copy.deepcopy(tem_res_collection_df.iloc[:, 0])
        Label_matrxi.iloc[:, 1] = cluster_index
        Label_matrxi.to_csv(main_url_save+str(id_num) +
                            '_label_inference.txt', sep='\t', index=0)
        # 2. MEC score/CPR score/AC/Pre/Recal/F1s1/ARI/AIC/BIC/HBIC/4smean
        Other_result = pd.DataFrame(np.zeros([3, 1]))
        #Other_result.iloc[0, :] = ACo2
        #Other_result.iloc[1, :] = Pre1o2
        #Other_result.iloc[2, :] = Recal1o2
        #Other_result.iloc[3, :] = F1s1o2
        #Other_result.iloc[4, :] = ARIo2
        Other_result.iloc[0, :] = AICo2
        Other_result.iloc[1, :] = BICo2
        Other_result.iloc[2, :] = HBICo2
        #Other_result.iloc[10, 1:] = S_means[:, 0]
        Other_result.to_csv(main_url_save+str(id_num)+'_Othertest.txt', sep='\t', index=0)
        # saving the bone and Ins
        #df_tempt_bone.to_csv(main_url_save+id_num+'_Bone.txt', sep='\t', index=0)
        #df_tempt_ins2.to_csv(main_url_save+id_num+'_Ins.txt', sep='\t', index=0)
        #tem_df_kk_ins = pd.DataFrame(columns=['pos', 'ins_len'])
        #loop_l_ins = 0
        #for ins_index in ins_inserted_dic:
        #    tem_df_kk_ins.loc[loop_l_ins] = [ins_index,ins_inserted_dic[ins_index]]
        #    loop_l_ins = loop_l_ins + 1
        #tem_df_kk_ins.to_csv(main_url_save+id_num+'_read_ins.txt', sep='\t', index=0)
        
        #tem_df_kk = pd.DataFrame(columns=['Consensus_id', 'pos', 'ins_content'])
        
        #loop_l = 0
        #for kk in range(len(mutation_rate_list)):
        #    kk_ins_dic = collection_ins_dic[kk]
        #    for xx in kk_ins_dic:
        #        tem_df_kk.loc[loop_l] = [kk, xx, kk_ins_dic[xx]]
        #        loop_l = loop_l + 1
        #tem_df_kk.to_csv(main_url_save+id_num+'_true_ins.txt', sep='\t', index=0)
        ###parameter saving
        time_save_list_df = pd.DataFrame(time_save_list)
        time_save_list_df.to_csv(main_url_save+'running_time.txt', sep='\t', index=0)
    
        #Table_true_variations.to_csv(main_url_save+str(id_num)+'_correct_information_sites.txt', sep='\t', index=0)
        infer_var_table.to_csv(main_url_save+'infer_information_sites.txt', sep='\t', index=0)
        #present_clu_dis2_1 = pd.DataFrame(present_clu_dis2[0])
        ### unmapping
        #collection_of_unmapped_df = pd.DataFrame(collection_of_unmapped)
        #collection_of_unmapped_df.to_csv(main_url_save+str(id_num)+'_unmapped.txt', sep='\t', index=0)

    
    else: #we have the common hap in the present region 
        full_var_site = sorted(list(full_matrix[:,0]))
        remaining_independent = copy.deepcopy(full_var_site)
        Selected_Best_ref = {}
        Selected_Best_ref_ins = {}
        up_del_beta = 0.05  # del
        up_del_p = 0.05
        up_ins_beta= 0.05  # del
        up_ins_p = 0.05
        with  mp.Pool(num_cores,maxtasksperchild=1000) as pool:

        # Result_infer_ind=pool.starmap(signle_sites_infer_max,[(remaining_independent[index],present_clu_dis2[kk][int(remaining_independent[index]),:],standard_ref_code[int(remaining_independent[index])],\
        # gamma[kk],theta_) for index in range(len(remaining_independent))])

            Result_infer_ind = pool.starmap(pre_max_, [(index, full_matrix[index,1:7]) for index in range(len(remaining_independent))])
            #Result_infer_ind =  Result_infer_ind_iter.get()
        #final_inferred=pool.starmap(sampling2,[(pp_random[index,0],pp_random[index,1],combination_list2.iloc[int(index/num_chains),0],combination_list2.iloc[int(index/num_chains),1],out_put_df,collection_of_sites,final_ref_code) for index in range((len(combination_list2)*num_chains))])
        #pool.close()
        #pool.join()

        # making up full result
        Result_infer_ = copy.deepcopy(Result_infer_ind)
        
        original_information_collection = reverse_trans2(ins_inserted_dic_df, full_var_site,1)
        tem_ins_ = ''
        pre_bone = -1
        ins_len_dic = {}
        del_len_dic = {}
        del_len_evo = 0
        correct_evo_num = 0
        mis_evo_num = 0
        infer_full_time = time.time()
        time_save_list.append(['Inferecne of the consensus', infer_full_time-t3])
        for ll_infer_cell_id  in range(len(Result_infer_)):
            ll_original_pos_ = original_information_collection[ll_infer_cell_id][1]
            ll_original_pos_state = original_information_collection[ll_infer_cell_id][2]
            present_base_infer = Result_infer_[ll_infer_cell_id][1]
            if ll_original_pos_state<0:
                if len(tem_ins_)>0:
                    Selected_Best_ref_ins[pre_bone] = tem_ins_
                    tem_ins_ = ''
                Selected_Best_ref[ll_original_pos_] = present_base_infer
                pre_bone = ll_original_pos_
                
                if present_base_infer<4:
                    if del_len_evo >0:
                        if del_len_evo in del_len_dic:
                            del_len_dic[del_len_evo] = del_len_dic[del_len_evo] + 1
                        else:
                            del_len_dic[del_len_evo] = 1
                    
                        del_len_evo = 0  
                    if present_base_infer == ref_code_standard[ll_original_pos_]:
                        correct_evo_num = correct_evo_num + 1
                    else:
                        mis_evo_num = mis_evo_num + 1
                
                    
                else:
                    del_len_evo = del_len_evo + 1
                    
            else:
                if del_len_evo >0:
                    if del_len_evo in del_len_dic:
                        del_len_dic[del_len_evo] = del_len_dic[del_len_evo] + 1
                    else:
                        del_len_dic[del_len_evo] = 1
                
                    del_len_evo = 0  
                tem_ins_ = tem_ins_ + Standardbase[present_base_infer]
        
        if len(tem_ins_)>0:
            Selected_Best_ref_ins[pre_bone] = tem_ins_
            
        if del_len_evo >0:
            if del_len_evo in del_len_dic:
                del_len_dic[del_len_evo] = del_len_dic[del_len_evo] + 1
            else:
                del_len_dic[del_len_evo] = 1
        
        if len(Selected_Best_ref_ins)>0:
            for ins_id in Selected_Best_ref_ins:
                if len(Selected_Best_ref_ins[ins_id]) in ins_len_dic:
                    ins_len_dic[len(Selected_Best_ref_ins[ins_id])] = ins_len_dic[len(Selected_Best_ref_ins[ins_id])] + 1
                else:
                    ins_len_dic[len(Selected_Best_ref_ins[ins_id])] = 1
        
        pcorrect = 1
        p_mis_loc = 0
        pins1_loc = 0
        pdel1_loc = 0
        lamuna_ins = 0
        lamuna_del = 0
        gamma1_per = correct_evo_num/(correct_evo_num+mis_evo_num)
        gamma2_per = mis_evo_num/(correct_evo_num+mis_evo_num)
        
        num_del = 0
        for key_del in del_len_dic:
            num_del = num_del + del_len_dic[key_del]
        gamma3_per = num_del/(num_del+correct_evo_num+mis_evo_num)
        
        num_ins = len(Selected_Best_ref_ins)
        gamma4_per = num_ins/(correct_evo_num+mis_evo_num)
        max_dis_para = 1
        max_createria_dis = 0.00001
        while max_dis_para > max_createria_dis:
            del_p = copy.deepcopy(up_del_p)
            del_beta = copy.deepcopy(up_del_beta)
            ins_p = copy.deepcopy(up_ins_p)
            ins_beta = copy.deepcopy(up_ins_beta)
            up_del_p = 0
            up_del_beta = 0
            up_ins_p = 0
            up_ins_beta = 0
            ## change the code
            if len(del_len_dic)>0:
                tem_dep_c = 0
                tem_del_loc_num = 0
                for dd_k in  del_len_dic:
                    up_del_beta = up_del_beta + (dd_k/(1-(1-del_p)*np.exp(-del_beta*dd_k)))*del_len_dic[dd_k]
                    tem_dep_c = tem_dep_c + 1/(1-(1-del_p)*np.exp(-del_beta*dd_k))*del_len_dic[dd_k]
                    tem_del_loc_num = tem_del_loc_num + del_len_dic[dd_k]
                import_data_p = (tem_dep_c,tem_del_loc_num)
                p_00 = 0
                up_del_p = fsolve(extimation_p, p_00,args=import_data_p)
               
                
                #up_del_p[kk] = -tem_del_loc_num*(1-del_p[kk])/(up_del_p[kk] * np.log(del_p[kk]))
                
                if up_del_p>1:
                   up_del_p = 1
                   
                if up_del_p<0:
                    up_del_p = 0
                  
                
                up_del_beta = tem_del_loc_num /up_del_beta
            else:
                up_del_beta = 0.5
                up_del_p = 0.5
            
            if len(ins_len_dic)>0:
                tem_dep_c = 0
                tem_ins_num = 0
                for ii_k in ins_len_dic:
                    up_ins_beta = up_ins_beta + (ii_k/(1-(1-ins_p)*np.exp(-ins_beta*ii_k)))*ins_len_dic[ii_k]
                    #up_ins_p[kk] = up_ins_p[kk] + 1/(1-(1-ins_p[kk])*np.exp(-ins_beta[kk]*ii_k))*tem_ins_len_dic[ii_k]
                    tem_dep_c = tem_dep_c + 1/(1-(1-ins_p)*np.exp(-ins_beta*ii_k))*ins_len_dic[ii_k]
                    tem_ins_num = tem_ins_num + ins_len_dic[ii_k]
                import_data_p = (tem_dep_c,tem_ins_num)
                p_00 = 0
                up_ins_p = fsolve(extimation_p, p_00,args=import_data_p)
                
                if up_ins_p>1:
                   up_ins_p = 1
                   
                if up_ins_p<0:
                    up_ins_p = 0
                    
                #up_ins_p[kk] = -tem_ins_num*(1-ins_p[kk])/(up_ins_p[kk] * np.log(ins_p[kk]))
                up_ins_beta = tem_ins_num /up_ins_beta
            else:
                up_ins_beta = 0.5
                up_ins_p = 0.5
            
            max_dis_para = max([abs(up_del_p-del_p),abs(up_del_beta-del_beta),\
                       abs(up_ins_p-ins_p),abs(up_ins_beta-ins_beta)])
        
        ### finish the parameter estimation
        parameters_f = [gamma1_per,gamma2_per , gamma3_per , gamma4_per ,\
            del_beta , del_p, ins_beta , ins_p ,pcorrect,p_mis_loc,pins1_loc,lamuna_ins,pdel1_loc,lamuna_del]
        finalized_time = time.time()
        time_save_list.append(['Final estimation parameters', finalized_time-infer_full_time])
        saving_para[:K] = parameters_f[0]
        saving_para[K:2*K] = parameters_f[1]
        saving_para[2*K:3*K] = parameters_f[2]
        saving_para[3*K:4*K] = parameters_f[3]
        saving_para[4*K:5*K] = parameters_f[4]
        saving_para[5*K:6*K] = parameters_f[5]
        saving_para[6*K:7*K] = parameters_f[6]
        saving_para[7*K:8*K] = parameters_f[7]
        saving_para[8*K:] = parameters_f[8:]
        
        infer_var_table =pd.DataFrame(columns=['chrom','pos','type','content'])
        infer_loop = 0
        
        for keys_bone in  Selected_Best_ref:
            if Selected_Best_ref[keys_bone]!=ref_code_standard[keys_bone]:
                if Selected_Best_ref[keys_bone]<4:
                    for s_con in range(proposed_k):
                        infer_var_table.loc[infer_loop] = [s_con, keys_bone, 'Mis', Standardbase[Selected_Best_ref[keys_bone]]]
                        infer_loop = infer_loop + 1
                    #infer_loop = infer_loop + 1
                else:
                    for s_con in range(proposed_k):
                        infer_var_table.loc[infer_loop] = [s_con, keys_bone, 'Del', '-']
                        infer_loop = infer_loop + 1
        for keys_ins in Selected_Best_ref_ins:
            for s_con in range(proposed_k):
                infer_var_table.loc[infer_loop] = [s_con, keys_ins, 'Ins', Selected_Best_ref_ins[keys_ins]]
                infer_loop = infer_loop + 1
        
        
        
        
        #saving_para_df = pd.DataFrame(saving_para)
        #saving_para_df.to_csv(main_url_save+str(id_num)+'_true_para.txt', sep='\t', index=0)
        #Table_true_variations.to_csv(main_url_save+str(id_num)+'_correct_information_sites.txt', sep='\t', index=0)
        infer_var_table.to_csv(main_url_save+'infer_information_sites.txt', sep='\t', index=0)
        
        reads_infer_cluf = pd.DataFrame(np.zeros([max_num_reads_refine, 2]))
        reads_infer_cluf.iloc[:,0] = range(max_num_reads_refine)
        reads_infer_cluf.iloc[:, 1] = np.zeros(max_num_reads_refine,dtype=int)
        time_save_list_df = pd.DataFrame(time_save_list)
        time_save_list_df.to_csv(main_url_save+'running_time.txt', sep='\t', index=0)
        
        
    
        
    
