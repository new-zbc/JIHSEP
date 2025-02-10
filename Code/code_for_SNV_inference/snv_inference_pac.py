#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  1 10:32:18 2023

@author: zhen
"""
### this is to save the code for inference
# finalize version
#from scipy.stats import fisher_exact
import scipy.stats as stats
#from scipy.stats import chisquare
#import sys
import numpy as np
#import re
import pandas as pd
#from pandas.core.frame import DataFrame
#from random import sample, choice
#from Bio import Align
#from Bio.Seq import Seq
#from Bio import SeqIO
#from Bio.SeqRecord import SeqRecord
#import pysam
#from scipy import stats
import biotite.sequence as seqbio
import biotite.sequence.align as align
#import itertools as it
#import math
#from multiprocessing import Pool
#import networkx as nx
import copy
#import seaborn as sns
#import itertools
#import matplotlib.pyplot as plt
#from fitter import Fitter, get_common_distributions, get_distributions
#from distfit import distfit
#from scipy.special import comb, perm
#from random import sample
Standardbase = ['a', 't', 'c', 'g', '-', 'n']
def bam_2_msa_ins(start,seq,cigar):
    
    #start = Result_table.iloc[2,2]
    #seq = Result_table.iloc[2,1]
    #cigar = Result_table.iloc[2,5]
    df_tempt_ins2_xun_dic= {}
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
        final_ref = final_ref + seg_ment_bone_seq[kk_key] + ins_num_dic[kk_key]*'n'
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
                        finalized_ins_dic[kkey] = 'n'
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
                        finalized_ins_dic[kkey] = 'n'*ins_inserted_dic[kkey]

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




def extract_s_e(s, e):
    return([s, e])


def Jaccard_pair8(read_index1, read_index2, final_indicator_table):
    # here we should refine the condition
    #threshold = 0
    #read_index1 = 2
    #read_index2 = 3
    #ref_code = finalized_ref_code_coll
    #coverage_df = reads_coverage
    #theta_r = 0.9
    #proba_list = labels_initial0_matrix
    #reads_infer_clu = reads_infer_clu
    start1 = final_indicator_table.iloc[read_index1, 1]
    end1 = final_indicator_table.iloc[read_index1, 2]
    start2 = final_indicator_table.iloc[read_index2, 1]
    end2 = final_indicator_table.iloc[read_index2, 2]
    # read_flag1=reads_infer_clu.iloc[read_index1,1]
    # read_flag2=reads_infer_clu.iloc[read_index2,1]
    #l_min =threshold*(end1-start1+1)

    l_min = 30
    max_s = max([start1, start2])
    min_e = min([end1, end2])
    return_coe = 0

    if min_e - max_s > l_min:

        #cover_s = max_s
        #cover_e = min_e
        #sub_cover = coverage_df.iloc[int(cover_s):int(cover_e+1),:]
        #sub_whole_cover = sum(sub_cover.iloc[:,1])
        mismatched_col_whole1 = final_indicator_table.iloc[read_index1, 3]
        matched_col_whole1 = final_indicator_table.iloc[read_index1, 4]
        mismatched_col_whole2 = final_indicator_table.iloc[read_index2, 3]
        matched_col_whole2 = final_indicator_table.iloc[read_index2, 4]
        whole1 = np.array(sorted(mismatched_col_whole1+matched_col_whole1))
        whole2 = np.array(sorted(mismatched_col_whole2+matched_col_whole2))
        #code_s = max_s*num_base - 1
        #code_e = (min_e+1)*num_base
        #selected_base1 = whole1[np.where((whole1>code_s)&((whole1<code_e)))]
        selected_base1_dic = {int(xx/6): xx % 6 for xx in whole1}
        #selected_base2 = whole2[np.where((whole2>code_s)&((whole2<code_e)))]
        selected_base2_dic = {int(xx/6): xx % 6 for xx in whole2}
        # for read_flag1 in range(len(proba_list_1)):
        #ref1 = ref_code[int(read_flag1)]
        

        Jaccard_coe = 0

        fenmu_ = 0
        for base_r in range(int(max_s), int(min_e+1)):
            if selected_base2_dic[base_r] != 5:
                fenmu_ = fenmu_ + 1
                if selected_base2_dic[base_r] == selected_base1_dic[base_r]:
                    Jaccard_coe = Jaccard_coe + 1
                    #Jaccard_coe = Jaccard_coe + 1/(cover_e-cover_s+1)

        # return_coe = np.exp((cover_1_c+cover_2_c)/(2*(min_e-max_s+1))*np.log(p_correct)+\
        #    (cover_1_w+cover_2_w)/(2*(min_e-max_s+1))*np.log(p_wrong+adjust_r))*Jaccard_coe
        return_coe = Jaccard_coe/fenmu_
    return ([read_index1, read_index2, return_coe])


def Jaccard_pair9(read_index1, read_index2, final_indicator_table, mutation_list):
    # here we should refine the condition
    #threshold = 0.5
    #read_index1 = 0
    #read_index2 = 2
    #ref_code = finalized_ref_code_coll
    #coverage_df = reads_coverage
    #theta_r = 0.9
    #proba_list = labels_initial0_matrix
    #reads_infer_clu = reads_infer_clu
    start1 = final_indicator_table.iloc[read_index1, 1]
    end1 = final_indicator_table.iloc[read_index1, 2]
    start2 = final_indicator_table.iloc[read_index2, 1]
    end2 = final_indicator_table.iloc[read_index2, 2]
    # read_flag1=reads_infer_clu.iloc[read_index1,1]
    # read_flag2=reads_infer_clu.iloc[read_index2,1]
    #l_min =threshold*min((end1-start1+1),(end2-start2+1))
    #NN= 5
    l_min = 30
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
            mismatched_col_whole1 = final_indicator_table.iloc[read_index1, 3]
            matched_col_whole1 = final_indicator_table.iloc[read_index1, 4]
            mismatched_col_whole2 = final_indicator_table.iloc[read_index2, 3]
            matched_col_whole2 = final_indicator_table.iloc[read_index2, 4]
            whole1 = np.array(sorted(mismatched_col_whole1+matched_col_whole1))
            whole2 = np.array(sorted(mismatched_col_whole2+matched_col_whole2))

            #selected_base1 = whole1[np.where((whole1>code_s)&((whole1<code_e)))]
            selected_base1_dic = {int(xx/6): xx % 6 for xx in whole1}
            #selected_base2 = whole2[np.where((whole2>code_s)&((whole2<code_e)))]
            selected_base2_dic = {int(xx/6): xx % 6 for xx in whole2}
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
            if fenmu_ > 0:
                return_coe = Jaccard_coe/fenmu_
    return ([read_index1, read_index2, return_coe])


def change_label_matrix(label_array, K_):
    #K_ = len(set(label_array))
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
            tem_seq = tem_seq + 'n'
        else:
            tem_seq = tem_seq + seq_o[int(code_unit)]
    return(tem_seq)


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
                [seq_collection[col_index], coverted_seq.lower()])
    except (ValueError, ZeroDivisionError):
        max_len = len(seq_collection[0])
        for col_index in range(len(seq_collection)):
            if len(seq_collection[col_index]) > max_len:
                max_len = len(seq_collection[col_index])
        result_mapping = []
        for col_index in range(len(seq_collection)):
            if len(seq_collection[col_index]) < max_len:
                refined_ins_ = seq_collection[col_index] + \
                    'n'*(max_len-len(seq_collection[col_index]))

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
        print(sites)
        ins_num = 0
        #sites = 6
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

    return(original_sites_var)



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



def signle_sites_infer_max2(site_id, base_number, standard_ref_code_site, gamma_, theta_,tem_pri):
    #kk = 2
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
    #    base_number = present_clu_dis2[0][899,:]
        #gamma_ = [0.98838986,0.01145738,0.00015277,0]
        #theta_ = copy.deepcopy(parameters_f[4:])
        gamma_c =  gamma_[0]
        gamma_mis = gamma_[1]
        gamma_ins = gamma_[3]
        gamma_del = gamma_[2]
        pcorrect = theta_[0]
        p_mis_loc = theta_[1]
        pins1_loc = theta_[2]
        pdel1_loc = theta_[3]
        lamuna_ins = theta_[4]
        lamuna_del = theta_[5]
        potential_probability = np.zeros(5)
        main_ref = standard_ref_code_site
        corretion_eps_from_zero = 0.000001
        #base_number = [0,0,0,2,0,0]
        #main_ref = 2
        if main_ref > 4:  # ins in the bone
            for main_base_xun in range(4):
                same_base = base_number[main_base_xun]
                not_same_base = sum(base_number[:4]) - same_base
                del_num = base_number[5]
                potential_probability[main_base_xun] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero) +\
                    (same_base+not_same_base)*np.log(1-pdel1_loc+corretion_eps_from_zero) + del_num*np.log(pdel1_loc+corretion_eps_from_zero) +\
                    del_num*possion_log_probability(0,lamuna_del) + \
                    np.log(gamma_ins+corretion_eps_from_zero)
            potential_probability[4] = (del_num)*np.log(1-pins1_loc+corretion_eps_from_zero) +\
                sum(base_number[:4])*(np.log(pins1_loc+corretion_eps_from_zero) +
                                      possion_log_probability(0,lamuna_ins))#+np.log(1-gamma_ins+corretion_eps_from_zero)
            #tem_probability = copy.deepcopy(potential_probability/tem_pri)
            #tem_probability = tem_probability - max(tem_probability)
            #tem_probability_main = np.exp(tem_probability)
            #tem_probability_main = tem_probability_main/sum(tem_probability_main)
            #post_sampling = np.random.multinomial(n= 1,pvals = tem_probability_main)
            #max_ID = list(post_sampling).index(1)
            #    n=1, pvals=[6/8, 1/8, 1/16, 1/32, 1/32])
            #l_i = list(l_i_vec).index(1)+1
            max_ID = list(potential_probability).index(max(potential_probability))
            if max_ID == 4:
                max_ID = 5
        else:
            if sum(base_number) > 0:
                for main_base_xun in range(4):  # del condtion: ref base in the 0-4
                    #main_base_xun = 1
                    same_base = base_number[main_base_xun]
                    not_same_base = sum(base_number[:4]) - same_base
                    del_num = base_number[4]
                    if main_base_xun == main_ref:
                        potential_probability[main_base_xun] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero) +\
                            (same_base+not_same_base)*np.log(1-pdel1_loc+corretion_eps_from_zero) + del_num*np.log(pdel1_loc+corretion_eps_from_zero) +\
                            del_num*possion_log_probability(0,lamuna_del)+np.log(gamma_c-corretion_eps_from_zero)
                    else:
                        potential_probability[main_base_xun] = same_base*np.log(pcorrect-corretion_eps_from_zero)+not_same_base*np.log(p_mis_loc+corretion_eps_from_zero) +\
                            (same_base+not_same_base)*np.log(1-pdel1_loc+corretion_eps_from_zero) + del_num*np.log(pdel1_loc+corretion_eps_from_zero) +\
                            del_num*possion_log_probability(0,lamuna_del)+np.log(gamma_mis+corretion_eps_from_zero)
    
                potential_probability[4] = (del_num)*np.log(1-pins1_loc+corretion_eps_from_zero) +\
                    sum(base_number[:4])*(np.log(pins1_loc+corretion_eps_from_zero) +
                                         possion_log_probability(0,lamuna_ins))+np.log(gamma_del+corretion_eps_from_zero)
    
                #tem_probability = copy.deepcopy(potential_probability/tem_pri)
                #tem_probability = tem_probability - max(tem_probability)
                #tem_probability_main = np.exp(tem_probability)
                #tem_probability_main = tem_probability_main/sum(tem_probability_main)
                #post_sampling = np.random.multinomial(n= 1,pvals = tem_probability_main)
                #max_ID = list(post_sampling).index(1)
                #    n=1, pvals=[6/8, 1/8, 1/16, 1/32, 1/32])
                #l_i = list(l_i_vec).index(1)+1
                max_ID = list(potential_probability).index(max(potential_probability))
            else:
                max_ID = main_ref
        #max_ID = list(potential_probability).index(max(potential_probability))
        return([site_id, max_ID])


def select_read_blocks(start_block, end_block, ele_list_tem, reads_info):
    #start_block = list_block_nei[0][0]
    #end_block = list_block[0][1]
    tem_return_res = []
    #reads_info = rselct_result
    #ele_list_tem = ele_list
    
    for read_label_xun in ele_list_tem:
        read_label_xun = int(read_label_xun)
        read_start = reads_info[read_label_xun][0]
        read_end = reads_info[read_label_xun][1]

        if not (start_block > read_end) and not (read_start >= end_block):
            tem_dic = reads_info[read_label_xun][2]
            select_dic = {k: v for k, v in tem_dic.items() if (
                (k >= start_block) and (k <= end_block))}
            tem_return_res.append([read_start, read_end, select_dic])
    return(tem_return_res)


def select_ref_block(start_block, end_block, standard_ref_code_dic):
    #tem_return_res = []

    select_dic = {k: v for k, v in standard_ref_code_dic.items() if (
        (k >= start_block) and (k <= end_block))}
    return(select_dic)


def block_sites_infer_ECM2(start_block, end_block, blobk_base_number, standard_ref_code_block, gamma_, theta_, updated_segment, segment_reads,tem_pri):
    
#for  ll in range(len(list_block)):
    #start_block = list_block[ll][0]
    #end_block = list_block[ll][1]
    #blobk_base_number = present_clu_dis2[kk][start_block:end_block+1,:]
    #standard_ref_code_block = final_ref_code_standard[start_block:end_block+1]
    #gamma_ = [0.99786129,0.00152765,0.00061106,0]
    # theta_ =
    #ll = 9
    #updated_segment = copy.deepcopy(Ref_dic_block[ll])
    #segment_reads = Dep_dic[ll]
    #start_block = 1327
    #end_block = 1328
    #blobk_base_number = copy.deepcopy(present_clu_dis2[kk][1327:1329,:])
    #dep_s = sip.select_ref_block(1327,1328,Best_ref_collection[8][2])
    #updated_segment = copy.deepcopy(dep_s)
    #segment_reads = copy.deepcopy(dep)
    read_number_tem = len(segment_reads)
    #read_id_tem = list(segment_reads.keys())
    gamma_c = gamma_[0]
    gamma_mis = gamma_[1]
    gamma_ins = gamma_[3]
    gamma_del = gamma_[2]
    return_result = []
    pcorrect = theta_[0]
    p_mis_loc = theta_[1]
    pins1_loc = theta_[2]
    pdel1_loc = theta_[3]
    lamuna_ins = theta_[4]
    lamuna_del = theta_[5]
    corretion_eps_from_zero = 0.0000000001

    # Indel_bound
    # sreening
    state_mapping = []
    for basss_index in range(start_block, end_block+1):
        if standard_ref_code_block[int(basss_index-start_block)] > 4:
            state_mapping.append(1)
        else:
            state_mapping.append(2)
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

    for basss_index in range(start_block, end_block+1):
        # print(updated_segment)
        #basss_index = 899
        potential_probability = np.zeros(5)
        updated_segment_cp = copy.deepcopy(updated_segment)
        # if

        if state_mapping[int(basss_index-start_block)] < 2:  # ins in the bone
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
                                    #np.log(stats.poisson.pmf(ii-1, lamuna_ins))
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
                                    #np.log(stats.poisson.pmf(dd-1, lamuna_del))
                                    
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

                    ID_tem = [0, 1, 2, 3, 5].index(main_base_xun)
                    if main_base_xun < 5:
                        full_full_prob = full_full_prob + np.log(gamma_ins+corretion_eps_from_zero)
                    else:
                        full_full_prob = full_full_prob + np.log(1-gamma_ins+corretion_eps_from_zero)

                    potential_probability[int(ID_tem)] = full_full_prob
                #tem_probability = copy.deepcopy(potential_probability/tem_pri)
                #tem_probability = tem_probability - max(tem_probability)
                #tem_probability_main = np.exp(tem_probability)
                #tem_probability_main = tem_probability_main/sum(tem_probability_main)
                #post_sampling = np.random.multinomial(n= 1,pvals = tem_probability_main)
                #max_ID = list(post_sampling).index(1)
                #    n=1, pvals=[6/8, 1/8, 1/16, 1/32, 1/32])
                #l_i = list(l_i_vec).index(1)+1
                max_ID = list(potential_probability).index(max(potential_probability))
                
                if max_ID == 4:
                    updated_segment[basss_index] = 5
            else:
                updated_segment[basss_index] = 5
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
                                    #np.log(stats.poisson.pmf(ii-1, lamuna_ins))
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
                                    #np.log(stats.poisson.pmf(dd-1, lamuna_del))
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

                    ID_tem = [0, 1, 2, 3, 4].index(main_base_xun)
                    if standard_ref_code_block[int(basss_index-start_block)] == main_base_xun:
                        full_full_prob = full_full_prob + np.log(gamma_c+corretion_eps_from_zero)
                    else:
                        if main_base_xun < 4:
                            full_full_prob = full_full_prob + np.log(gamma_mis+corretion_eps_from_zero)
                        else:
                            full_full_prob = full_full_prob + np.log(gamma_del+corretion_eps_from_zero)
                    potential_probability[int(ID_tem)] = full_full_prob
                #tem_probability = copy.deepcopy(potential_probability/tem_pri)
                #tem_probability = tem_probability - max(tem_probability)
                #tem_probability_main = np.exp(tem_probability)
                #tem_probability_main = tem_probability_main/sum(tem_probability_main)
                #post_sampling = np.random.multinomial(n= 1,pvals = tem_probability_main)
                #max_ID = list(post_sampling).index(1)
                #    n=1, pvals=[6/8, 1/8, 1/16, 1/32, 1/32])
                #l_i = list(l_i_vec).index(1)+1
                max_ID = list(potential_probability).index(max(potential_probability))
                

            else:
                for main_base_xun in [0, 1, 2, 3]:
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
                                    #np.log(stats.poisson.pmf(ii-1, lamuna_ins))
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
                                    #np.log(stats.poisson.pmf(dd-1, lamuna_del))
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
                    if standard_ref_code_block[int(basss_index-start_block)] == main_base_xun:
                        full_full_prob = full_full_prob + np.log(gamma_c+corretion_eps_from_zero)
                    else:
                        full_full_prob = full_full_prob + np.log(gamma_mis+corretion_eps_from_zero)
                    potential_probability[int(ID_tem)] = full_full_prob

                #tem_probability = copy.deepcopy(potential_probability/tem_pri)
                #tem_probability = tem_probability - max(tem_probability)
                #tem_probability_main = np.exp(tem_probability)
                #tem_probability_main = tem_probability_main/sum(tem_probability_main)
                #post_sampling = np.random.multinomial(n= 1,pvals = tem_probability_main)
                #max_ID = list(post_sampling).index(1)
                #    n=1, pvals=[6/8, 1/8, 1/16, 1/32, 1/32])
                #l_i = list(l_i_vec).index(1)+1
                max_ID = list(potential_probability).index(max(potential_probability))
                

            updated_segment[basss_index] = max_ID

        return_result.append([basss_index, updated_segment[basss_index]])

    
    return(return_result)

def possion_log_probability(l_max,lamuna):
    #K_max = 500
    #lamuna = 0.5
    full_sum = -lamuna + l_max*np.log(lamuna)
    if l_max>0:
        ll_vec = np.array(range(1,int(l_max)+1))
        full_sum = full_sum - sum(np.log(ll_vec))
    return(full_sum)

