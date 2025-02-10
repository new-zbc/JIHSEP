#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 09:24:31 2021

@author: zhenzhang
"""
# finalize version
import scipy.stats as stats
from scipy.stats import chisquare
import sys
import numpy as np
import os
import pandas as pd
import pysam
import biotite.sequence as seqbio
import biotite.sequence.align as align
import itertools as it
import multiprocessing as mp
import networkx as nx
import copy
import itertools
from sklearn.cluster import SpectralClustering,DBSCAN,OPTICS,Birch
from scipy.optimize import linear_sum_assignment
from communities.algorithms import louvain_method, girvan_newman, hierarchical_clustering
from networkx.algorithms.community import greedy_modularity_communities, naive_greedy_modularity_communities,\
    lukes_partitioning, asyn_lpa_communities, label_propagation_communities, asyn_fluidc
sys.path.append('/path/of/code_for_SNV_inference/') # this should be the path of the sub-folder of the packages of the this main file
import snv_inference_pac as sip
from sklearn.decomposition import PCA
#import pickle


Standardbase = ['a', 't', 'c', 'g', '-', 'n']




if __name__ == '__main__':
    print(sys.argv)
    proposed_k = int(sys.argv[1])
    main_url_input_bam = sys.argv[2]
    main_url_input_ref = sys.argv[3]
    main_url_save = sys.argv[4]

  
    corretion_eps_from_zero = 0.0000000001

    path_exon_bam=main_url_input_bam
    inbam = pysam.AlignmentFile(path_exon_bam, "rb" )
    
    
    Result_table = pd.DataFrame(columns=['Num', 'seq', 'start_mapping', 'end_mapping',
                                 'len_seq', 'Cigar_information'])
    Result_table['Cigar_information'] = Result_table['Cigar_information'].astype(
        'object')
    
    Final_result_interested = pd.DataFrame(np.zeros([7]))
    

    
    tem_res_collection = []
    s_name_num = 0
    for read_xun in inbam.fetch():
            if not read_xun.is_secondary and not read_xun.is_unmapped and not read_xun.is_supplementary:
                read_name = read_xun.query_name
                read_start = read_xun.reference_start
                read_seq = (read_xun.query_sequence).lower()
                read_cigar = read_xun.cigar
                if len(read_cigar)>0:
                    bone_base_count = 0
                    pre_clip_len = 0
                    after_clip_len = 0
                    for cigar_ in read_cigar:
                        if cigar_[0] in [0,2,7,8]:
                            bone_base_count = bone_base_count + cigar_[1]
                    if read_cigar[0][0] in [4]:
                        pre_clip_len = int(read_cigar[0][1])
                    if read_cigar[-1][0] in [4]:
                        after_clip_len = int(read_cigar[-1][1])
                    tem_res = sip.bam_2_msa_ins(read_start,read_seq,read_cigar)
                    tem_res_collection.append([read_name,read_start,tem_res[0],read_seq,read_cigar,tem_res[1]])
                    Result_table.loc[s_name_num,:] = [read_name,read_seq,read_start,read_start+bone_base_count-1,bone_base_count-1,0]
                    Result_table.iat[s_name_num,-1] = read_cigar
                    s_name_num = s_name_num + 1
                

    
            
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
    for keys_ in whole_insertion:
        
        dic_multiple_ali_ins[keys_] = sip.multiple_seq_alignment(whole_insertion[keys_])
        max_len_ =  len(dic_multiple_ali_ins[keys_][0][1])
        for sssss_index in range(len(dic_multiple_ali_ins[keys_])):
            if len(dic_multiple_ali_ins[keys_][sssss_index][1])>max_len_:
                max_len_ = len(dic_multiple_ali_ins[keys_])
        ins_inserted_dic[keys_] = max_len_
    ### modified to here we need to filter the unmapped reads from the Result_table
    Result_table_cp = copy.deepcopy(Result_table)
    Result_table_cp['updated_index'] = range(len(Result_table_cp))
    
    ## construct the dictionary
    maping_dic = {}
    for ll in range(len(Result_table_cp)):
        maping_dic[ll] = Result_table_cp.iloc[ll,0]
    
    ### reference_part
   
    fasta_file_path = main_url_input_ref 
    fasta_file = pysam.FastaFile(fasta_file_path)
    seqname = fasta_file.references[0]
    BestRefSeq = (fasta_file.fetch(seqname)).lower()
    n_c = len(BestRefSeq)#stabndard bone
    
    ### refine the reference
    refine_best_ref  = sip.refine_the_fasta(ins_inserted_dic,BestRefSeq)
    
    final_ref_standard = copy.deepcopy(refine_best_ref)
    final_ref_code_standard = sip.string_2_code(final_ref_standard)
    
    [out_put_df, collection_of_siteso] = sip.adding_ins_to_bone(
        Result_table_cp, ins_inserted_dic, dic_multiple_ali_ins, final_ref_code_standard)
    [refined_count_table, final_indicator_table] = sip.coding_mismatch_all_insertion_error2(
        out_put_df, final_ref_standard)
    
    full_bone = list(range(len(final_ref_code_standard)))
    rselct_result = sip.trans_to_dic(final_indicator_table, full_bone)
    ### The selection of important sites
    
    full_list_consider = list(range(len(full_bone)))
       
    
    max_num_reads_refine = len(out_put_df)
    read_full_bases = [np.zeros([len(final_ref_code_standard), 6])
                       for x in range(max_num_reads_refine)]
    for ele in range(max_num_reads_refine):
        for bases_xun in rselct_result[int(ele)][2]:
            read_full_bases[int(ele)][int(bases_xun), int(rselct_result[int(ele)][2][bases_xun])] = \
                read_full_bases[int(ele)][int(bases_xun), int(
                    rselct_result[int(ele)][2][bases_xun])] + 1
           
    
    # simulation of making up the missing value
    combination_list2 = pd.DataFrame(
        it.combinations(list(range(max_num_reads_refine)), 2))

    labels_initial0 = np.zeros(max_num_reads_refine)
    labels_initial0_matrix = sip.change_label_matrix(labels_initial0, 1)
    
    
    standard_ref_code = final_ref_code_standard
    
        
    correct_ratio = 0.9
    coolection_snp1 = []
    #coolection_snp11 = []
    for ll in range(len(refined_count_table)):
        #ll = 756
        main_ref_index = final_ref_code_standard[ll]#Standardbase.index(real_con1[ll])
        main_number = refined_count_table.iloc[ll,2+int(main_ref_index)]
        whole_number = sum(refined_count_table.iloc[ll,2:8])
        #ll = 6460
        if whole_number>5:
            other_number = whole_number-main_number
            expt_1 = int(whole_number*correct_ratio)
            expt_2 = whole_number - expt_1
            if whole_number>=50 and other_number>=5:
                #print(ll)
                fobs = [main_number,other_number]
                fexp = [expt_1,expt_2]
                s = chisquare(fobs,fexp)
                if s[1]<0.05:
                   if main_number<expt_1:
                       coolection_snp1.append(ll)
               
            else:
                fobs = [expt_2,other_number]
                f_other = [expt_1,main_number]
                s = stats.barnard_exact([fobs,f_other],alternative='less')
                if s.pvalue<0.05:
                   if main_number<expt_1:
                       coolection_snp1.append(ll)
                
   
    
    
    if  len(coolection_snp1)>0:        

        new_keys2 = copy.deepcopy(coolection_snp1)
       
        with  mp.get_context("spawn").Pool(5,maxtasksperchild=1000) as pool:
    
        
            Jaccard_matrix_df3=pool.starmap(sip.Jaccard_pair9,[(combination_list2.iloc[index,0],combination_list2.iloc[index,1],final_indicator_table,\
                                     new_keys2) for index in range(len(combination_list2))])
    
        pool.close()
        pool.join()
    else:
        print('warning, no sig snp')
        with  mp.get_context("spawn").Pool(5,maxtasksperchild=1000) as pool:
    
        
    
            Jaccard_matrix_df3 = pool.starmap(sip.Jaccard_pair8, [(combination_list2.iloc[index, 0], combination_list2.iloc[index, 1], final_indicator_table,
                                                           ) for index in range(len(combination_list2))])
    
        pool.close()
        pool.join()
        
    #Jaccard_matri
    Jaccard_matrix_df3 = pd.DataFrame(Jaccard_matrix_df3)
    full_T_eq_coll = np.zeros([max_num_reads_refine, max_num_reads_refine])
                      
    n_trail_ele = np.ones([max_num_reads_refine, max_num_reads_refine])* 0.0001 
    full_inferred_res = []
    same_count_col = []
    not_same_count_col = []
    full_label = []
    K = proposed_k
    saving_para = np.zeros(K*4+6)
    
    Label_collection = []
    
    X_Jaccard = sip.convert_sparse_normal(Jaccard_matrix_df3, max_num_reads_refine)
    X_Jaccard_col_sum = sum(X_Jaccard)
    collection_informative_line = []
    collection_uninformative_line = []
    for ll in range(len(X_Jaccard_col_sum)):
        if X_Jaccard_col_sum[ll]>0:
            collection_informative_line.append(ll)
        else:
            collection_uninformative_line.append(ll)
    X_Jaccard_dis = np.ones(
        [max_num_reads_refine, max_num_reads_refine]) - X_Jaccard-np.diag(np.ones(max_num_reads_refine))
    X_Jaccard_dis[X_Jaccard_dis < 0] = 0
   
    
    
    # 1. spectural clustering
    clustering = SpectralClustering(
        n_clusters=proposed_k, assign_labels='discretize', affinity='precomputed_nearest_neighbors').fit(X_Jaccard_dis)
    # K=3
    labels_initial1 = clustering.labels_
    if len(set(labels_initial1)) != K:
        labels_initial1 = np.zeros(max_num_reads_refine)
    Label_collection.append(labels_initial1)
    
    graph_nx_o = nx.from_numpy_matrix(X_Jaccard)
    graph_nx_d = nx.from_numpy_matrix(X_Jaccard_dis)
    
    # 2. greedy_modularity_communities
    communities2 = greedy_modularity_communities(
        graph_nx_o, weight='weight', cutoff=K,best_n=K)
    
    labels_initial2 = np.zeros(max_num_reads_refine)
    if len(communities2) == K:
        for kk in range(1, K):
            tem_mem_list = list(communities2[kk])
            for list_index in tem_mem_list:
                labels_initial2[int(list_index)] = kk
    Label_collection.append(labels_initial2)
    
    
    
    
    # 4.asyn_lpa_communities/no number of communities
    labels_initial4 = np.zeros(max_num_reads_refine)
    result4 = asyn_lpa_communities(graph_nx_o, weight='weight')
    communities4 = [frozenset(c) for c in result4]
    if len(communities4) == K:
        for kk in range(1, len(communities4)):
            tem_mem_list = list(communities4[kk])
            for list_index in tem_mem_list:
                labels_initial4[int(list_index)] = kk
    Label_collection.append(labels_initial4)
    
    # 5.label_propagation_communities/no number of communities
    labels_initial5 = np.zeros(max_num_reads_refine)
    communities5 = list(label_propagation_communities(graph_nx_o))
    if len(communities5) == K:
        for kk in range(1, len(communities5)):
            tem_mem_list = list(communities5[kk])
            for list_index in tem_mem_list:
                labels_initial5[int(list_index)] = kk
    Label_collection.append(labels_initial5)
    
    
    # 7.asyn_fluidc
    labels_initial7 = np.zeros(max_num_reads_refine)
    try:
        result7 = asyn_fluidc(graph_nx_o, k=K)
        communities7 = [frozenset(c) for c in result7]
        if len(communities7) == K:
            for kk in range(1, len(communities7)):
                tem_mem_list = list(communities7[kk])
                for list_index in tem_mem_list:
                    labels_initial7[int(list_index)] = kk
            Label_collection.append(labels_initial7)
        else:
            Label_collection.append(labels_initial7)
    except Exception:
        Label_collection.append(labels_initial7)
    
    # 8 Birth
    labels_initial8 = np.zeros(max_num_reads_refine)
    pca = PCA()
    transform_8 = pca.fit_transform(X_Jaccard)  
    model8 =Birch(
    n_clusters=  proposed_k).fit(transform_8)
    labels_initial8 =model8.labels_#usaged bu
   
    Label_collection.append(labels_initial8)
    
    # hierarchical_clustering
    labels_initial9 = np.zeros(max_num_reads_refine)
    communities9 = hierarchical_clustering(X_Jaccard, n=K)
    if len(communities9) == K:
        for kk in range(1, len(communities9)):
            tem_mem_list = list(communities9[kk])
            for list_index in tem_mem_list:
                labels_initial9[int(list_index)] = kk
    Label_collection.append(labels_initial9)
    
    ### perfect
    #Label_collection.append(copy.deepcopy(Result_table.iloc[:, -1]))
    
    Collection_label = []
    sum_of_prob = []
    method_index_c = []
    Best_ref_collection = {}
    parameters_saving = {}
    for method_index in range(1):#range(len(Label_collection)):
        print('method_index', method_index)
        if len(set(Label_collection[method_index])) == K:
            #method_index = 6
            reads_infer_clu = pd.DataFrame(np.zeros([max_num_reads_refine, 2]))
            reads_infer_clu.iloc[:, 0] = range(max_num_reads_refine)
            reads_infer_clu.iloc[:, 1] = copy.deepcopy(Label_collection[method_index])
  
            last_probability = np.zeros([max_num_reads_refine, K])
            for read_index in range(max_num_reads_refine):
                last_probability[read_index, int(
                    Label_collection[method_index][read_index])] = 1
            last_portions = np.ones(K)/K
            present_clu_dis2 = [
                np.zeros([len(full_bone), 6]) for x in range(K)]
            #present_clu_dis2_pro = [np.zeros([len(full_bone),5]) for x in range(K)]
            Best_ref = [{} for x in range(K)]
    
            portions = np.zeros(K)
            Basic_ref = final_ref_code_standard
            #con_ref_prob = np.zeros([max_num_reads_refine,K])
            con_ref_prob = np.zeros(K)
            con_ref_prob_full = np.zeros([K, 5])
            gamma1 = 0.9*np.ones(K)  # correct
            gamma2 = 0.05*np.ones(K)  # mis
            gamma3 = 0.05*np.ones(K)  # del
            gamma4 = 0.05*np.ones(K)  # ins
            
            p_mis_loc = np.random.rand(1)/10
            pins1_loc = np.random.rand(1)/10
            lamuna_ins = int(np.random.rand(1)*6)+1
            
            tao = 1
            tao2 = 1
            tao_step = 0.97
            tao_step2 = 0.99
            pdel1_loc = np.random.rand(1)/10
            lamuna_del = int(np.random.rand(1)*6)+1
            
    
            pcorrect = np.random.rand(1)/10+0.9
            theta_ = [pcorrect, p_mis_loc, pins1_loc,
                      pdel1_loc, lamuna_ins, lamuna_del]
            for kk in range(K):
                ele_list = list(
                    reads_infer_clu[reads_infer_clu.iloc[:, 1] == kk].iloc[:, 0])
                
                portions[kk] = len(ele_list)
                
                for ele in ele_list:
                    present_clu_dis2[kk] = present_clu_dis2[kk] + \
                        read_full_bases[int(ele)]
                if len(ele_list) > 0:
    
                    remaining_independent = copy.deepcopy(full_bone)
                    with  mp.get_context("spawn").Pool(5,maxtasksperchild=1000) as pool:
    
                    
    
                        Result_infer_ind = pool.starmap(sip.signle_sites_infer_max2, [(remaining_independent[index], present_clu_dis2[kk][int(remaining_independent[index]), :], standard_ref_code[int(remaining_independent[index])],
                                                                               [gamma1[kk], gamma2[kk], gamma3[kk], gamma4[kk]], theta_,tao2) for index in range(len(remaining_independent))])
    
                    pool.close()
                    pool.join()
    
                    # making up full result
                    Result_infer_ = copy.deepcopy(Result_infer_ind)
    
                    num_del_evo_block = 0
                    first_del_flag = 0
                    end_del_flag = 0
                    previous_del_flag = 0
                    tem_del_length = 0
                    for ll_per in Result_infer_:
                        Best_ref[kk][ll_per[0]] = copy.deepcopy(ll_per[1])
                    bone_k = 0
    
                    if Basic_ref[int(Result_infer_[0][0])] == Result_infer_[0][1]:
                        if Basic_ref[int(Result_infer_[0][0])] not in [4,5]:
                            con_ref_prob_full[kk,
                                              0] = con_ref_prob_full[kk, 0] + 1
                    else:
                        if Result_infer_[0][1] == 4:
                            con_ref_prob_full[kk,
                                              2] = con_ref_prob_full[kk, 2] + 1
                            #num_del_evo_block = 1
                            previous_del_flag = 1
                            first_del_flag = 1
                            tem_del_length = tem_del_length + 1
    
                        else:
                            con_ref_prob_full[kk,
                                              1] = con_ref_prob_full[kk, 1] + 1
    
                    seq_idex = 1
                    if bone_k in ins_inserted_dic:
                        bone_k_num = ins_inserted_dic[bone_k]
                        tem_num_string = ''
                        for kk_ins in range(bone_k_num):
                            tem_num_string = tem_num_string + \
                                str(Result_infer_[seq_idex+kk_ins][1])
                        tem_num_string_sub = tem_num_string.replace(
                            '5', '')
                        if len(tem_num_string_sub) > 0:
                            con_ref_prob_full[kk, 3] = con_ref_prob_full[kk,
                                                                         3] + 1
    
                        seq_idex = seq_idex + ins_inserted_dic[bone_k]
    
                    bone_k = bone_k + 1
    
                    while bone_k < n_c:
                        # if Basic_ref[int(seq_idex)] == 5:
                        #    print(seq_idex)
    
                        if Basic_ref[int(seq_idex)] == Best_ref[kk][seq_idex]:
                            if Basic_ref[int(Result_infer_[0][0])] not in [4,5]:
                                con_ref_prob_full[kk,
                                                  0] = con_ref_prob_full[kk, 0] + 1
                                if tem_del_length > 0:
                                    num_del_evo_block = num_del_evo_block + 1
                                    previous_del_flag = 0
                                    tem_del_length = 0
                        else:
                            if Best_ref[kk][seq_idex] == 4:
                                con_ref_prob_full[kk,
                                                  2] = con_ref_prob_full[kk, 2] + 1
                                if previous_del_flag:
                                    tem_del_length = tem_del_length + 1
                                else:
                                    previous_del_flag = 1
                                    tem_del_length = 1
    
                                if bone_k == n_c - 1:
                                    end_del_flag = 1
    
                            else:
                                con_ref_prob_full[kk,
                                                  1] = con_ref_prob_full[kk, 1] + 1
                                if tem_del_length > 0:
                                    num_del_evo_block = num_del_evo_block + 1
                                    previous_del_flag = 0
                                    tem_del_length = 0
                        pre_seq_index_bone = copy.deepcopy(seq_idex)
                        seq_idex = seq_idex + 1
                        if bone_k in ins_inserted_dic:
                            bone_k_num = ins_inserted_dic[bone_k]
                            tem_num_string = ''
                            for kk_ins in range(bone_k_num):
                                tem_num_string = tem_num_string + \
                                    str(Best_ref[kk][seq_idex+kk_ins])
                            tem_num_string_sub = tem_num_string.replace(
                                '5', '')
                            if len(tem_num_string_sub) > 0:
                                con_ref_prob_full[kk, 3] = con_ref_prob_full[kk,
                                                                             3] + 1
    
                            seq_idex = seq_idex + ins_inserted_dic[bone_k]
    
                        bone_k = bone_k + 1
                    if tem_del_length > 0:
                        num_del_evo_block = num_del_evo_block + 1
                    con_ref_prob_full[kk, 4] = n_c - 2 * \
                        num_del_evo_block + first_del_flag + end_del_flag
                   
                else:
                    for basss_index in full_bone:
                        Best_ref[kk][int(basss_index)] = copy.deepcopy(
                            Basic_ref[int(basss_index)])
                    con_ref_prob_full[kk, :] = [n_c, 0, 0, 0, n_c-1]
    
                    
    
            portions = portions/sum(portions)
    
            conditional_probalility_ = np.zeros([max_num_reads_refine, K])
            conditional_probalility_o = np.zeros([max_num_reads_refine, K])
            update_fre = [np.zeros([len(full_bone), 6]) for x in range(K)]
            whole_correct = np.zeros([max_num_reads_refine, K])
            whole_mis = np.zeros([max_num_reads_refine, K])
            whole_ins = np.zeros([max_num_reads_refine, K])
            whole_del = np.zeros([max_num_reads_refine, K])
            whole_len_del = {}
            whole_len_ins = {}
            whole_ins_dic = {}
            whole_del_dic = {}
    
           
    
            for read_xun in range(max_num_reads_refine):
                tem_dic = rselct_result[int(read_xun)][2]
                begin_id = final_indicator_table.iloc[read_xun, 1]
                end_id = final_indicator_table.iloc[read_xun, 2]
                for kk in range(K):
                    #kk = 0
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
                        if Best_ref[kk][base_xun] < 4:
                            #bone_flag = 1
                            if tem_dic[base_xun] == Best_ref[kk][base_xun]:
    
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
    
                            elif tem_dic[base_xun] < 4 and tem_dic[base_xun] != Best_ref[kk][base_xun]:
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
                        elif Best_ref[kk][base_xun] == 4:
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
                            tem_probabity_ins = tem_probabity_ins + \
                                len_ins[ii] * \
                                sip.possion_log_probability(ii-1,lamuna_del)
                    else:
                        tem_probabity_ins = 0
                        num_ins_loc = 0
    
                    del_whole = 0
                    if len(len_del) > 0:
                        values_ = len_del.values()
                        num_del_loc = sum(np.array(list(values_)))
                        tem_probabity_del = num_del_loc*np.log(pdel1_loc+corretion_eps_from_zero)
                        for dd in len_del:
                            tem_probabity_del = tem_probabity_del + \
                                len_del[dd] * \
                                sip.possion_log_probability(dd-1,lamuna_del)
                                
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
                    whole_correct[read_xun, kk] = num_correct
                    whole_mis[read_xun, kk] = num_mis
                    whole_ins[read_xun, kk] = num_ins_loc
                    whole_del[read_xun, kk] = num_del_loc
                    whole_len_del[read_xun*K+kk] = len_del
                    whole_len_ins[read_xun*K+kk] = len_ins
                
                    conditional_probalility_[read_xun, kk] = (
                        full_prob + np.log(portions[kk]+corretion_eps_from_zero))/tao
                    conditional_probalility_o[read_xun,
                                              kk] = full_prob + np.log(portions[kk]+corretion_eps_from_zero)
                conditional_probalility_[read_xun, :] = conditional_probalility_[
                    read_xun, :]-max(conditional_probalility_[read_xun, :])
                conditional_probalility_[read_xun, :] = np.exp(
                    conditional_probalility_[read_xun, :])
                conditional_probalility_[read_xun, :] = conditional_probalility_[
                    read_xun, :]/sum(conditional_probalility_[read_xun, :])
                # annealing sampling procedure
                l_s_vec = np.random.multinomial(
                    n=1, pvals=conditional_probalility_[read_xun, :])
                reads_infer_clu.iloc[read_xun, 1] = list(l_s_vec).index(1)
                # Mstep
            for kkk in range(K):
                tem_gamma_ = np.array([gamma1[kkk], gamma2[kkk], gamma3[kkk], gamma4[kkk], 1-gamma4[kkk]])
                tem_gamma_ = tem_gamma_ + corretion_eps_from_zero 
                tem_gamma_[:3] = tem_gamma_[:3]/sum(tem_gamma_[:3])
                tem_gamma_[3:] = tem_gamma_[3:]/sum(tem_gamma_[3:])
                con_ref_prob[kkk] = sum(con_ref_prob_full[kkk, :]*np.log(
                    tem_gamma_))
    
            creatierstop = 0.0001
            #correct_num= np.zeros(K)
            up_gamma1 = np.ones(K)*0.9
            up_gamma2 = np.ones(K)*0.05
            up_gamma3 = np.ones(K)*0.05
            up_gamma4 = np.ones(K)*0.0001
            #up_gamma5=  np.ones(K)*0.05
            present_clu_dis2 = [
                np.zeros([len(full_bone), 6]) for x in range(K)]
            con_ref_prob_full = np.zeros([K, 5])
    
            for kk in range(K):
                ele_list = list(
                    reads_infer_clu[reads_infer_clu.iloc[:, 1] == kk].iloc[:, 0])
                portions[kk] = len(ele_list)
                # if len(ele_list)>0:
                if len(ele_list) > 0:
                    
                    for ele in ele_list:
                        present_clu_dis2[kk] = present_clu_dis2[kk] + \
                            read_full_bases[int(ele)]
                    # 1.seperate the room for seperate zoom for block indel
                    full_index_ = np.array(full_bone)
                    block_index_cluster = full_index_[
                        (present_clu_dis2[kk][:, 4]+present_clu_dis2[kk][:, 5]) > 0]
                    remaining_independent = list(
                        set(full_bone)-set(block_index_cluster))
    
                    if len(block_index_cluster) > 0:
                        #block_idnex = 0
                        start_tem = block_index_cluster[0]
                        list_block = []
                        loop_num = 1
                        while loop_num < len(block_index_cluster):
                            if block_index_cluster[loop_num]-block_index_cluster[loop_num-1] > 1:
                                end_tem = block_index_cluster[loop_num-1]
                                if end_tem - start_tem > 0:
                                    list_block.append([start_tem, end_tem])
                                else:
                                    remaining_independent.append(start_tem)
                                start_tem = block_index_cluster[loop_num]
                                #block_idnex = block_idnex +1
                            loop_num = loop_num + 1
    
                        with  mp.get_context("spawn").Pool(5,maxtasksperchild=1000) as pool:
    
                            Result_infer_ind = pool.starmap(sip.signle_sites_infer_max2, [(remaining_independent[index], present_clu_dis2[kk][int(remaining_independent[index]), :], standard_ref_code[int(remaining_independent[index])],
                                                                                   [gamma1[kk], gamma2[kk], gamma3[kk], gamma4[kk]], theta_,tao) for index in range(len(remaining_independent))])
    
                            Dep_dic = pool.starmap(sip.select_read_blocks, [(
                            list_block[index][0], list_block[index][1], ele_list, rselct_result) for index in range(len(list_block))])
                            Ref_dic_block = pool.starmap(sip.select_ref_block, [(
                            list_block[index][0], list_block[index][1], Best_ref[kk]) for index in range(len(list_block))])
                            Result_infer_dep = pool.starmap(sip.block_sites_infer_ECM2, [(list_block[index][0], list_block[index][1], present_clu_dis2[kk][int(list_block[index][0]):int(list_block[index][1])+1, :],
                                                                                  final_ref_code_standard[int(list_block[index][0]):int(
                                                                                      list_block[index][1])+1], [gamma1[kk], gamma2[kk], gamma3[kk], gamma4[kk]],
                                                                                  theta_, copy.deepcopy(Ref_dic_block[index]), Dep_dic[index],tao) for index in range(len(list_block))])
    
                        pool.close()
                        pool.join()
    
                        
    
                        # making up full result
                        Result_infer_ = copy.deepcopy(Result_infer_ind)
                        for ll_xun_ in Result_infer_dep:
                            Result_infer_ = Result_infer_ + ll_xun_
                        collection_inferred = []
                        for ll_per in Result_infer_:
                            collection_inferred.append(ll_per[0])
                            Best_ref[kk][ll_per[0]] = copy.deepcopy(
                                ll_per[1])
                        num_del_evo_block = 0
                        first_del_flag = 0
                        end_del_flag = 0
                        previous_del_flag = 0
                        tem_del_length = 0
                        for ll_per in Result_infer_:
                            Best_ref[kk][ll_per[0]] = copy.deepcopy(
                                ll_per[1])
                        bone_k = 0
    
                        if Basic_ref[int(Result_infer_[0][0])] == Result_infer_[0][1]:
                            if Basic_ref[int(Result_infer_[0][0])] not in [4,5]:
                                con_ref_prob_full[kk,
                                                  0] = con_ref_prob_full[kk, 0] + 1
                        else:
                            if Result_infer_[0][1] == 4:
                                con_ref_prob_full[kk,
                                                  2] = con_ref_prob_full[kk, 2] + 1
                                #num_del_evo_block = 1
                                previous_del_flag = 1
                                first_del_flag = 1
                                tem_del_length = tem_del_length + 1
    
                            else:
                                con_ref_prob_full[kk,
                                                  1] = con_ref_prob_full[kk, 1] + 1
    
                        seq_idex = 1
                        if bone_k in ins_inserted_dic:
                            bone_k_num = ins_inserted_dic[bone_k]
                            tem_num_string = ''
                            for kk_ins in range(bone_k_num):
                                tem_num_string = tem_num_string + \
                                    str(Result_infer_[seq_idex+kk_ins][1])
                            tem_num_string_sub = tem_num_string.replace(
                                '5', '')
                            if len(tem_num_string_sub) > 0:
                                con_ref_prob_full[kk, 3] = con_ref_prob_full[kk,
                                                                             3] + 1
    
                            seq_idex = seq_idex + ins_inserted_dic[bone_k]
    
                        bone_k = bone_k + 1
    
                        while bone_k < n_c:
                           
    
                            if Basic_ref[int(seq_idex)] == Best_ref[kk][seq_idex]:
                                if Basic_ref[int(Result_infer_[0][0])] not in [4,5]:
                                    con_ref_prob_full[kk,
                                                      0] = con_ref_prob_full[kk, 0] + 1
                                    if tem_del_length > 0:
                                        num_del_evo_block = num_del_evo_block + 1
                                        previous_del_flag = 0
                                        tem_del_length = 0
                            else:
                                if Best_ref[kk][seq_idex] == 4:
                                    con_ref_prob_full[kk,
                                                      2] = con_ref_prob_full[kk, 2] + 1
                                    if previous_del_flag:
                                        tem_del_length = tem_del_length + 1
                                    else:
                                        previous_del_flag = 1
                                        tem_del_length = 1
    
                                    if bone_k == n_c - 1:
                                        end_del_flag = 1
    
                                else:
                                    con_ref_prob_full[kk,
                                                      1] = con_ref_prob_full[kk, 1] + 1
                                    if tem_del_length > 0:
                                        num_del_evo_block = num_del_evo_block + 1
                                        previous_del_flag = 0
                                        tem_del_length = 0
                            pre_seq_index_bone = copy.deepcopy(seq_idex)
                            seq_idex = seq_idex + 1
                            if bone_k in ins_inserted_dic:
                                bone_k_num = ins_inserted_dic[bone_k]
                                tem_num_string = ''
                                for kk_ins in range(bone_k_num):
                                    tem_num_string = tem_num_string + \
                                        str(Best_ref[kk][seq_idex+kk_ins])
                                tem_num_string_sub = tem_num_string.replace(
                                    '5', '')
                                if len(tem_num_string_sub) > 0:
                                    con_ref_prob_full[kk, 3] = con_ref_prob_full[kk, 3] + 1
    
                                seq_idex = seq_idex + \
                                    ins_inserted_dic[bone_k]
    
                            bone_k = bone_k + 1
                        if tem_del_length > 0:
                            num_del_evo_block = num_del_evo_block + 1
                        con_ref_prob_full[kk, 4] = n_c - 2 * \
                            num_del_evo_block + first_del_flag + end_del_flag
                    else:
                        with  mp.get_context("spawn").Pool(5,maxtasksperchild=1000) as pool:
    
                            Result_infer_ind = pool.starmap(sip.signle_sites_infer_max2, [(remaining_independent[index], present_clu_dis2[kk][int(remaining_independent[index]), :], standard_ref_code[int(remaining_independent[index])],
                                                                                   [gamma1[kk], gamma2[kk], gamma3[kk], gamma4[kk]], theta_,tao2) for index in range(len(remaining_independent))])
    
                        pool.close()
                        pool.join()
                        # making up full result
                        Result_infer_ = copy.deepcopy(Result_infer_ind)
    
                        for ll_per in Result_infer_:
                            Best_ref[kk][ll_per[0]] = copy.deepcopy(
                                ll_per[1])
                        num_del_evo_block = 0
                        first_del_flag = 0
                        end_del_flag = 0
                        previous_del_flag = 0
                        tem_del_length = 0
    
                        bone_k = 0
    
                        if Basic_ref[int(Result_infer_[0][0])] == Result_infer_[0][1]:
                            if Basic_ref[int(Result_infer_[0][0])] not in [4,5]:
                                con_ref_prob_full[kk,
                                                  0] = con_ref_prob_full[kk, 0] + 1
                        else:
                            if Result_infer_[0][1] == 4:
                                con_ref_prob_full[kk,
                                                  2] = con_ref_prob_full[kk, 2] + 1
                                #num_del_evo_block = 1
                                previous_del_flag = 1
                                first_del_flag = 1
                                tem_del_length = tem_del_length + 1
    
                            else:
                                con_ref_prob_full[kk,
                                                  1] = con_ref_prob_full[kk, 1] + 1
    
                        seq_idex = 1
                        if bone_k in ins_inserted_dic:
                            bone_k_num = ins_inserted_dic[bone_k]
                            tem_num_string = ''
                            for kk_ins in range(bone_k_num):
                                tem_num_string = tem_num_string + \
                                    str(Result_infer_[seq_idex+kk_ins][1])
                            tem_num_string_sub = tem_num_string.replace(
                                '5', '')
                            if len(tem_num_string_sub) > 0:
                                con_ref_prob_full[kk, 3] = con_ref_prob_full[kk,
                                                                             3] + 1
    
                            seq_idex = seq_idex + ins_inserted_dic[bone_k]
    
                        bone_k = bone_k + 1
    
                        while bone_k < n_c:
                            if Basic_ref[int(seq_idex)] == Result_infer_[seq_idex][1]:
                                if Basic_ref[int(Result_infer_[0][0])] not in [4,5]:
                                    con_ref_prob_full[kk,
                                                      0] = con_ref_prob_full[kk, 0] + 1
                                    if tem_del_length > 0:
                                        num_del_evo_block = num_del_evo_block + 1
                                        previous_del_flag = 0
                                        tem_del_length = 0
                            else:
                                if Result_infer_[seq_idex][1] == 4:
                                    con_ref_prob_full[kk,
                                                      2] = con_ref_prob_full[kk, 2] + 1
                                    if previous_del_flag:
                                        tem_del_length = tem_del_length + 1
                                    else:
                                        previous_del_flag = 1
                                        tem_del_length = 1
    
                                    if bone_k == n_c - 1:
                                        end_del_flag = 1
    
                                else:
                                    con_ref_prob_full[kk,
                                                      1] = con_ref_prob_full[kk, 1] + 1
                                    if tem_del_length > 0:
                                        num_del_evo_block = num_del_evo_block + 1
                                        previous_del_flag = 0
                                        tem_del_length = 0
                            pre_seq_index_bone = copy.deepcopy(seq_idex)
                            seq_idex = seq_idex + 1
                            if bone_k in ins_inserted_dic:
                                bone_k_num = ins_inserted_dic[bone_k]
                                tem_num_string = ''
                                for kk_ins in range(bone_k_num):
                                    tem_num_string = tem_num_string + \
                                        str(Result_infer_[
                                            seq_idex+kk_ins][1])
                                tem_num_string_sub = tem_num_string.replace(
                                    '5', '')
                                if len(tem_num_string_sub) > 0:
                                    con_ref_prob_full[kk, 3] = con_ref_prob_full[kk, 3] + 1
    
                                seq_idex = seq_idex + \
                                    ins_inserted_dic[bone_k]
    
                            bone_k = bone_k + 1
                        if tem_del_length > 0:
                            num_del_evo_block = num_del_evo_block + 1
                        con_ref_prob_full[kk, 4] = n_c - 2 * \
                            num_del_evo_block + first_del_flag + end_del_flag
                    tem_eq_num = con_ref_prob_full[kk, 0]
                    tem_neq_num = con_ref_prob_full[kk, 1]
                    tem_del_num = con_ref_prob_full[kk, 2]
                    tem_ins_num = con_ref_prob_full[kk, 3]
                    tem_no_ins_num = con_ref_prob_full[kk, 4]
                    up_gamma4[kk] = (tem_ins_num) / \
                        (tem_ins_num+tem_no_ins_num)
                    up_gamma1[kk] = tem_eq_num / \
                        (tem_eq_num+tem_neq_num+tem_del_num)
                    up_gamma3[kk] = (
                        tem_del_num)/(tem_eq_num+tem_neq_num+tem_del_num)
                    up_gamma2[kk] = 1-up_gamma1[kk]-up_gamma3[kk]
                    
                    if up_gamma1[kk] < 0.9:
                        up_gamma1[kk] = 0.9
                        up_gamma3[kk] = (
                            tem_del_num)/(tem_eq_num+tem_neq_num+tem_del_num)
                        up_gamma2[kk] = 1-up_gamma1[kk]-up_gamma3[kk]
                else:
    
                    for basss_index in full_bone:
                        Best_ref[kk][int(basss_index)] = copy.deepcopy(
                            Basic_ref[int(basss_index)])
                    con_ref_prob_full[kk, :] = [n_c, 0, 0, 0, n_c-1]
                    up_gamma1[kk] = np.random.beta(a=48, b=4, size=1)
                    up_gamma2[kk] = np.random.rand(1)*0.01
                    up_gamma3[kk] = 1-up_gamma1[kk]-up_gamma2[kk]
                    up_gamma4[kk] = np.random.rand(1)*0.01
    
                    
            check_dis1 = abs(last_probability -
                             conditional_probalility_).max()
    
            portions = portions + 0.001
            portions = portions/sum(portions)
            check_dis2 = abs(last_portions-portions).max()
    
            # parameters
    
            wei_cor = (whole_correct*conditional_probalility_).sum()
            wei_mis = (whole_mis*conditional_probalility_).sum()
            wei_ins = (whole_ins*conditional_probalility_).sum()
            wei_del = (whole_del*conditional_probalility_).sum()
    
            for read_xun in range(max_num_reads_refine):
                if len(whole_len_del[read_xun]) > 0:
                    for kk in range(K):
                        for sub_len in whole_len_del[read_xun*K+kk]:
                            if sub_len not in whole_del_dic:
                                whole_del_dic[sub_len] = whole_len_del[read_xun*K +
                                                                       kk][sub_len]*conditional_probalility_[read_xun, kk]
                            else:
                                whole_del_dic[sub_len] = whole_del_dic[sub_len] + \
                                    whole_len_del[read_xun*K+kk][sub_len] * \
                                    conditional_probalility_[read_xun, kk]
    
                        for sub_len2 in whole_len_ins[read_xun*K+kk]:
                            if sub_len2 not in whole_ins_dic:
                                whole_ins_dic[sub_len2] = whole_len_ins[read_xun*K +
                                                                        kk][sub_len2]*conditional_probalility_[read_xun, kk]
                            else:
                                whole_ins_dic[sub_len2] = whole_ins_dic[sub_len2] + \
                                    whole_len_ins[read_xun*K+kk][sub_len2] * \
                                    conditional_probalility_[read_xun, kk]
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
            check_dis3 = max(abs(up_lamuna_del-lamuna_del), abs(lamuna_ins-up_lamuna_ins),
                             abs(up_pcorrect -
                                 pcorrect), abs(up_p_mis_loc-p_mis_loc),
                             abs(up_pins1_loc -
                                 pins1_loc), abs(up_pdel1_loc-pdel1_loc),
                             abs(up_gamma1 -
                                 gamma1).max(), abs(up_gamma2-gamma2).max(),
                             abs(up_gamma3-gamma3).max(), abs(up_gamma4-gamma4).max())
    
            
            check_dis = max(check_dis1, check_dis2, check_dis3)
            max_iter = 100
            min_iter = 50
            loop_ = 0
            check_dis_flag = 1
            while check_dis > creatierstop:
    
                if loop_ < max_iter:
                    tao = tao * tao_step
                    tao2 = tao2 * tao_step2
                    last_probability = copy.deepcopy(
                        conditional_probalility_)
                    last_fre = copy.deepcopy(update_fre)
                    conditional_probalility_ = np.zeros(
                        [max_num_reads_refine, K])
                    update_fre = [np.zeros([len(full_bone), 6])
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
                    pdel1_loc = copy.deepcopy(up_pdel1_loc)
                    lamuna_del = copy.deepcopy(up_lamuna_del)
                    #mu_del = copy.deepcopy(up_mu_del)
                    #sig_del = copy.deepcopy(up_sig_del)
                    pcorrect = copy.deepcopy(up_pcorrect)
                    theta_ = [pcorrect, p_mis_loc, pins1_loc,
                              pdel1_loc, lamuna_ins, lamuna_del]
                    whole_correct = np.zeros([max_num_reads_refine, K])
                    whole_mis = np.zeros([max_num_reads_refine, K])
                    whole_ins = np.zeros([max_num_reads_refine, K])
                    whole_del = np.zeros([max_num_reads_refine, K])
                    whole_len_del = {}
                    whole_len_ins = {}
                    whole_ins_dic = {}
                    whole_del_dic = {}
                    for read_xun in range(max_num_reads_refine):
                        #read_xun = 0
                        tem_dic = rselct_result[int(read_xun)][2]
                        begin_id = final_indicator_table.iloc[read_xun, 1]
                        end_id = final_indicator_table.iloc[read_xun, 2]
                        for kk in range(K):
                            #kk = 0
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
                                if Best_ref[kk][base_xun] < 4:
                                    #bone_flag = 1
                                    if tem_dic[base_xun] == Best_ref[kk][base_xun]:
    
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
    
                                    elif tem_dic[base_xun] < 4 and tem_dic[base_xun] != Best_ref[kk][base_xun]:
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
                                elif Best_ref[kk][base_xun] == 4:
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
                                    tem_probabity_ins = tem_probabity_ins + \
                                        len_ins[ii] * \
                                        sip.possion_log_probability(ii-1,lamuna_ins)
                                   
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
                                        sip.possion_log_probability(ii-1,lamuna_del)
                                       
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
                            whole_correct[read_xun, kk] = num_correct
                            whole_mis[read_xun, kk] = num_mis
                            whole_ins[read_xun, kk] = num_ins_loc
                            whole_del[read_xun, kk] = num_del_loc
                            whole_len_del[read_xun*K+kk] = len_del
                            whole_len_ins[read_xun*K+kk] = len_ins
                        
                            conditional_probalility_[read_xun, kk] = (
                                full_prob + np.log(portions[kk]+corretion_eps_from_zero))/tao
                            conditional_probalility_o[read_xun,
                                                      kk] = full_prob + np.log(portions[kk]+corretion_eps_from_zero)
                        conditional_probalility_[read_xun, :] = conditional_probalility_[
                            read_xun, :]-max(conditional_probalility_[read_xun, :])
                        conditional_probalility_[read_xun, :] = np.exp(
                            conditional_probalility_[read_xun, :])
                        conditional_probalility_[read_xun, :] = conditional_probalility_[
                            read_xun, :]/sum(conditional_probalility_[read_xun, :])
                        # annealing sampling procedure
                        l_s_vec = np.random.multinomial(
                            n=1, pvals=conditional_probalility_[read_xun, :])
                        reads_infer_clu.iloc[read_xun, 1] = list(l_s_vec).index(1)
                        # Mstep
                    for kkk in range(K):
                        tem_gamma_ = np.array([gamma1[kkk], gamma2[kkk], gamma3[kkk], gamma4[kkk], 1-gamma4[kkk]])
                        tem_gamma_ = tem_gamma_ + corretion_eps_from_zero 
                        tem_gamma_[:3] = tem_gamma_[:3]/sum(tem_gamma_[:3])
                        tem_gamma_[3:] = tem_gamma_[3:]/sum(tem_gamma_[3:])
                        con_ref_prob[kkk] = sum(con_ref_prob_full[kkk, :]*np.log(
                            tem_gamma_))
                   
                    up_gamma1 = np.ones(K)*0.9
                    up_gamma2 = np.ones(K)*0.05
                    up_gamma3 = np.ones(K)*0.05
                    up_gamma4 = np.ones(K)*0.0001
                    #up_gamma5=  np.ones(K)*0.05
                    present_clu_dis2 = [
                        np.zeros([len(full_bone), 6]) for x in range(K)]
                    con_ref_prob_full = np.zeros([K, 5])
                    
    
                    for kk in range(K):
                            ele_list = list(
                                reads_infer_clu[reads_infer_clu.iloc[:, 1] == kk].iloc[:, 0])
                            portions[kk] = len(ele_list)
                            if len(ele_list) > 0:
                                # for basss_index in full_bone:
                                #max_ID = list(present_clu_dis2[kk][int(basss_index),:]).index(max(present_clu_dis2[kk][int(basss_index),:]))
                                for ele in ele_list:
                                    present_clu_dis2[kk] = present_clu_dis2[kk] + \
                                        read_full_bases[int(ele)]
                                # 1.seperate the room for seperate zoom for block indel
                                full_index_ = np.array(full_bone)
                                block_index_cluster = full_index_[
                                    (present_clu_dis2[kk][:, 4]+present_clu_dis2[kk][:, 5]) > 0]
                                remaining_independent = list(
                                    set(full_bone)-set(block_index_cluster))
    
                                if len(block_index_cluster) > 0:
                                    #block_idnex = 0
                                    start_tem = block_index_cluster[0]
                                    list_block = []
                                    loop_num = 1
                                    while loop_num < len(block_index_cluster):
                                        if block_index_cluster[loop_num]-block_index_cluster[loop_num-1] > 1:
                                            end_tem = block_index_cluster[loop_num-1]
                                            if end_tem - start_tem > 0:
                                                list_block.append(
                                                    [start_tem, end_tem])
                                            else:
                                                remaining_independent.append(
                                                    start_tem)
                                            start_tem = block_index_cluster[loop_num]
                                            #block_idnex = block_idnex +1
                                        loop_num = loop_num + 1
    
                                    if block_index_cluster[-1]-start_tem > 0:
                                        end_tem = block_index_cluster[-1]
                                        list_block.append(
                                            [start_tem, end_tem])
                                    else:
                                        remaining_independent.append(
                                            start_tem)
    
                                    with  mp.get_context("spawn").Pool(5,maxtasksperchild=1000) as pool:
                    
                                   
                    
                                        Result_infer_ind = pool.starmap(sip.signle_sites_infer_max2, [(remaining_independent[index], present_clu_dis2[kk][int(remaining_independent[index]), :], standard_ref_code[int(remaining_independent[index])],
                                                                                               [gamma1[kk], gamma2[kk], gamma3[kk], gamma4[kk]], theta_,tao2) for index in range(len(remaining_independent))])
                    
                                    pool.close()
                                    pool.join()
                                    
                                    
    
                                    # select the block
                                    with  mp.get_context("spawn").Pool(5,maxtasksperchild=1000) as pool:
    
                                        Dep_dic = pool.starmap(sip.select_read_blocks, [(
                                        list_block[index][0], list_block[index][1], ele_list, rselct_result) for index in range(len(list_block))])
    
                                    pool.close()
                                    pool.join()
    
                                   
                                    with  mp.get_context("spawn").Pool(5,maxtasksperchild=1000) as pool:
                                        Ref_dic_block = pool.starmap(sip.select_ref_block, [(
                                        list_block[index][0], list_block[index][1], Best_ref[kk]) for index in range(len(list_block))])
                                    
    
                                    pool.close()
                                    pool.join()
    
                                    # perfomr the block selection
                                    with  mp.get_context("spawn").Pool(5,maxtasksperchild=1000) as pool:
    
                                    
                                        Result_infer_dep = pool.starmap(sip.block_sites_infer_ECM2, [(list_block[index][0], list_block[index][1], present_clu_dis2[kk][int(list_block[index][0]):int(list_block[index][1])+1, :],
                                              final_ref_code_standard[int(list_block[index][0]):int(
                                                       list_block[index][1])+1], [gamma1[kk], gamma2[kk], gamma3[kk], gamma4[kk]],
                                                   theta_, copy.deepcopy(Ref_dic_block[index]), Dep_dic[index],tao2) for index in range(len(list_block))])
                                    pool.close()
                                    pool.join()
    
                                    # making up full result
                                    Result_infer_ = copy.deepcopy(
                                        Result_infer_ind)
                                    for ll_xun_ in Result_infer_dep:
                                        Result_infer_ = Result_infer_ + ll_xun_
                                    for ll_per in Result_infer_:
                                        Best_ref[kk][ll_per[0]] = copy.deepcopy(
                                            ll_per[1])
                                    num_del_evo_block = 0
                                    first_del_flag = 0
                                    end_del_flag = 0
                                    previous_del_flag = 0
                                    tem_del_length = 0
                                    for ll_per in Result_infer_:
                                        Best_ref[kk][ll_per[0]] = copy.deepcopy(
                                            ll_per[1])
                                    bone_k = 0
    
                                    if Basic_ref[int(Result_infer_[0][0])] == Result_infer_[0][1]:
                                        if Basic_ref[int(Result_infer_[0][0])] not in [4,5]:
                                            con_ref_prob_full[kk,
                                                              0] = con_ref_prob_full[kk, 0] + 1
                                    else:
                                        if Result_infer_[0][1] == 4:
                                            con_ref_prob_full[kk,
                                                              2] = con_ref_prob_full[kk, 2] + 1
                                            #num_del_evo_block = 1
                                            previous_del_flag = 1
                                            first_del_flag = 1
                                            tem_del_length = tem_del_length + 1
    
                                        else:
                                            con_ref_prob_full[kk,
                                                              1] = con_ref_prob_full[kk, 1] + 1
    
                                    seq_idex = 1
                                    if bone_k in ins_inserted_dic:
                                        bone_k_num = ins_inserted_dic[bone_k]
                                        tem_num_string = ''
                                        for kk_ins in range(bone_k_num):
                                            tem_num_string = tem_num_string + \
                                                str(Result_infer_[
                                                    seq_idex+kk_ins][1])
                                        tem_num_string_sub = tem_num_string.replace(
                                            '5', '')
                                        if len(tem_num_string_sub) > 0:
                                            con_ref_prob_full[kk, 3] = con_ref_prob_full[kk, 3] + 1
    
                                        seq_idex = seq_idex + \
                                            ins_inserted_dic[bone_k]
    
                                    bone_k = bone_k + 1
    
                                    while bone_k < n_c:
                                        # if Basic_ref[int(seq_idex)] == 5:
                                        #    print(seq_idex)
    
                                        if Basic_ref[int(seq_idex)] == Best_ref[kk][seq_idex]:
                                            if Basic_ref[int(Result_infer_[0][0])] not in [4,5]:
                                                con_ref_prob_full[kk,
                                                                  0] = con_ref_prob_full[kk, 0] + 1
                                                if tem_del_length > 0:
                                                    num_del_evo_block = num_del_evo_block + 1
                                                    previous_del_flag = 0
                                                    tem_del_length = 0
                                        else:
                                            if Best_ref[kk][seq_idex] == 4:
                                                con_ref_prob_full[kk,
                                                                  2] = con_ref_prob_full[kk, 2] + 1
                                                if previous_del_flag:
                                                    tem_del_length = tem_del_length + 1
                                                else:
                                                    previous_del_flag = 1
                                                    tem_del_length = 1
    
                                                if bone_k == n_c - 1:
                                                    end_del_flag = 1
    
                                            else:
                                                con_ref_prob_full[kk,
                                                                  1] = con_ref_prob_full[kk, 1] + 1
                                                if tem_del_length > 0:
                                                    num_del_evo_block = num_del_evo_block + 1
                                                    previous_del_flag = 0
                                                    tem_del_length = 0
                                        pre_seq_index_bone = copy.deepcopy(
                                            seq_idex)
                                        seq_idex = seq_idex + 1
                                        if bone_k in ins_inserted_dic:
                                            bone_k_num = ins_inserted_dic[bone_k]
                                            tem_num_string = ''
                                            for kk_ins in range(bone_k_num):
                                                tem_num_string = tem_num_string + \
                                                    str(Best_ref[kk]
                                                        [seq_idex+kk_ins])
                                            tem_num_string_sub = tem_num_string.replace(
                                                '5', '')
                                            if len(tem_num_string_sub) > 0:
                                                con_ref_prob_full[kk, 3] = con_ref_prob_full[kk, 3] + 1
    
                                            seq_idex = seq_idex + \
                                                ins_inserted_dic[bone_k]
    
                                        bone_k = bone_k + 1
                                    if tem_del_length > 0:
                                        num_del_evo_block = num_del_evo_block + 1
                                    con_ref_prob_full[kk, 4] = n_c - 2 * \
                                        num_del_evo_block + first_del_flag + end_del_flag
                                else:
                                    with  mp.get_context("spawn").Pool(5,maxtasksperchild=1000) as pool:
    
                                        Result_infer_ind = pool.starmap(sip.signle_sites_infer_max2, [(remaining_independent[index], present_clu_dis2[kk][int(remaining_independent[index]), :], standard_ref_code[int(remaining_independent[index])],
                                                                                               [gamma1[kk], gamma2[kk], gamma3[kk], gamma4[kk]], theta_,tao2) for index in range(len(remaining_independent))])
    
                                    
                                    pool.close()
                                    pool.join()
                                    
                                    num_del_evo_block = 0
                                    first_del_flag = 0
                                    end_del_flag = 0
                                    previous_del_flag = 0
                                    tem_del_length = 0
                                    for ll_per in Result_infer_:
                                        Best_ref[kk][ll_per[0]] = copy.deepcopy(
                                            ll_per[1])
                                    bone_k = 0
    
                                    if Basic_ref[int(Result_infer_[0][0])] == Result_infer_[0][1]:
                                        if Basic_ref[int(Result_infer_[0][0])] not in [4,5]:
                                            con_ref_prob_full[kk,
                                                              0] = con_ref_prob_full[kk, 0] + 1
                                    else:
                                        if Result_infer_[0][1] == 4:
                                            con_ref_prob_full[kk,
                                                              2] = con_ref_prob_full[kk, 2] + 1
                                            #num_del_evo_block = 1
                                            previous_del_flag = 1
                                            first_del_flag = 1
                                            tem_del_length = tem_del_length + 1
    
                                        else:
                                            con_ref_prob_full[kk,
                                                              1] = con_ref_prob_full[kk, 1] + 1
    
                                    seq_idex = 1
                                    if bone_k in ins_inserted_dic:
                                        bone_k_num = ins_inserted_dic[bone_k]
                                        tem_num_string = ''
                                        for kk_ins in range(bone_k_num):
                                            tem_num_string = tem_num_string + \
                                                str(Result_infer_[
                                                    seq_idex+kk_ins][1])
                                        tem_num_string_sub = tem_num_string.replace(
                                            '5', '')
                                        if len(tem_num_string_sub) > 0:
                                            con_ref_prob_full[kk, 3] = con_ref_prob_full[kk, 3] + 1
    
                                        seq_idex = seq_idex + \
                                            ins_inserted_dic[bone_k]
    
                                    bone_k = bone_k + 1
    
                                    while bone_k < n_c:
                                        # if Basic_ref[int(seq_idex)] == 5:
                                        #    print(seq_idex)
    
                                        if Basic_ref[int(seq_idex)] == Best_ref[kk][seq_idex]:
                                            if Basic_ref[int(Result_infer_[0][0])] not in [4,5]:
                                                con_ref_prob_full[kk,
                                                                  0] = con_ref_prob_full[kk, 0] + 1
                                                if tem_del_length > 0:
                                                    num_del_evo_block = num_del_evo_block + 1
                                                    previous_del_flag = 0
                                                    tem_del_length = 0
                                        else:
                                            if Best_ref[kk][seq_idex] == 4:
                                                con_ref_prob_full[kk,
                                                                  2] = con_ref_prob_full[kk, 2] + 1
                                                if previous_del_flag:
                                                    tem_del_length = tem_del_length + 1
                                                else:
                                                    previous_del_flag = 1
                                                    tem_del_length = 1
    
                                                if bone_k == n_c - 1:
                                                    end_del_flag = 1
    
                                            else:
                                                con_ref_prob_full[kk,
                                                                  1] = con_ref_prob_full[kk, 1] + 1
                                                if tem_del_length > 0:
                                                    num_del_evo_block = num_del_evo_block + 1
                                                    previous_del_flag = 0
                                                    tem_del_length = 0
                                        pre_seq_index_bone = copy.deepcopy(
                                            seq_idex)
                                        seq_idex = seq_idex + 1
                                        if bone_k in ins_inserted_dic:
                                            bone_k_num = ins_inserted_dic[bone_k]
                                            tem_num_string = ''
                                            for kk_ins in range(bone_k_num):
                                                tem_num_string = tem_num_string + \
                                                    str(Best_ref[kk]
                                                        [seq_idex+kk_ins])
                                            tem_num_string_sub = tem_num_string.replace(
                                                '5', '')
                                            if len(tem_num_string_sub) > 0:
                                                con_ref_prob_full[kk, 3] = con_ref_prob_full[kk, 3] +1
    
                                            seq_idex = seq_idex + \
                                                ins_inserted_dic[bone_k]
    
                                        bone_k = bone_k + 1
                                    if tem_del_length > 0:
                                        num_del_evo_block = num_del_evo_block + 1
                                    con_ref_prob_full[kk, 4] = n_c - 2 * \
                                        num_del_evo_block + first_del_flag + end_del_flag-1
                                    # block_sites_infer_ECM(start_block,end_block,blobk_base_number,standard_ref_code_block,gamma_,theta_,updated_segment,segment_reads)
                                tem_eq_num = con_ref_prob_full[kk, 0]
                                tem_neq_num = con_ref_prob_full[kk, 1]
                                tem_del_num = con_ref_prob_full[kk, 2]
                                tem_ins_num = con_ref_prob_full[kk, 3]
                                tem_no_ins_num = con_ref_prob_full[kk, 4]
                                up_gamma4[kk] = (
                                    tem_ins_num)/(tem_ins_num+tem_no_ins_num)
                                up_gamma1[kk] = tem_eq_num / \
                                    (tem_eq_num+tem_neq_num+tem_del_num)
                                up_gamma3[kk] = (
                                    tem_del_num)/(tem_eq_num+tem_neq_num+tem_del_num)
                                up_gamma2[kk] = 1 - \
                                    up_gamma1[kk]-up_gamma3[kk]
                                #up_beta_a = tem_eq_num + 48
                                #up_beta_b = tem_neq_num + 4
                                #up_gamma[kk] =  np.random.beta(a = up_beta_a,b = up_beta_b,size=1)
                                #up_gamma[kk] =  tem_eq_num/(tem_neq_num+tem_eq_num)
                                if up_gamma1[kk] < 0.9:
                                    up_gamma1[kk] = 0.9
                                    up_gamma3[kk] = (
                                        tem_del_num)/(tem_eq_num+tem_neq_num+tem_del_num)
                                    up_gamma2[kk] = 1 - \
                                        up_gamma1[kk]-up_gamma3[kk]
    
                            else:
    
                                for basss_index in full_bone:
                                    Best_ref[kk][int(basss_index)] = copy.deepcopy(
                                        Basic_ref[int(basss_index)])
                                con_ref_prob_full[kk, :] = [
                                    n_c, 0, 0, 0, n_c-1]
                                up_gamma1[kk] = np.random.beta(
                                    a=48, b=4, size=1)
                                up_gamma2[kk] = np.random.rand(1)*0.01
                                up_gamma3[kk] = 1 - \
                                    up_gamma1[kk]-up_gamma2[kk]
                                up_gamma4[kk] = np.random.rand(1)*0.01
                                
    
                    check_dis1 = abs(last_probability -
                                     conditional_probalility_).max()
    
                    #portions = sum(conditional_probalility_)
                    portions = portions + 0.001
                    portions = portions/sum(portions)
    
                    check_dis2 = abs(last_portions-portions).max()
                    # parameters
    
                    wei_cor = (whole_correct *
                               conditional_probalility_).sum()
                    wei_mis = (whole_mis*conditional_probalility_).sum()
                    wei_ins = (whole_ins*conditional_probalility_).sum()
                    wei_del = (whole_del*conditional_probalility_).sum()
    
                    for read_xun in range(max_num_reads_refine):
                        if len(whole_len_del[read_xun]) > 0:
                            for kk in range(K):
                                for sub_len in whole_len_del[read_xun*K+kk]:
                                    if sub_len not in whole_del_dic:
                                        whole_del_dic[sub_len] = whole_len_del[read_xun*K +
                                                                               kk][sub_len]*conditional_probalility_[read_xun, kk]
                                    else:
                                        whole_del_dic[sub_len] = whole_del_dic[sub_len] + \
                                            whole_len_del[read_xun*K+kk][sub_len] * \
                                            conditional_probalility_[
                                                read_xun, kk]
    
                                for sub_len2 in whole_len_ins[read_xun*K+kk]:
                                    if sub_len2 not in whole_ins_dic:
                                        whole_ins_dic[sub_len2] = whole_len_ins[read_xun*K +
                                                                                kk][sub_len2]*conditional_probalility_[read_xun, kk]
                                    else:
                                        whole_ins_dic[sub_len2] = whole_ins_dic[sub_len2] + \
                                            whole_len_ins[read_xun*K+kk][sub_len2] * \
                                            conditional_probalility_[
                                                read_xun, kk]
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
                    check_dis3 = max(abs(up_lamuna_del-lamuna_del), abs(lamuna_ins-up_lamuna_ins),
                                     abs(up_pcorrect -
                                         pcorrect), abs(up_p_mis_loc-p_mis_loc),
                                     abs(up_pins1_loc -
                                         pins1_loc), abs(up_pdel1_loc-pdel1_loc),
                                     abs(up_gamma1 -
                                         gamma1).max(), abs(up_gamma2-gamma2).max(),
                                     abs(up_gamma3-gamma3).max(), abs(up_gamma4-gamma4).max())
                    
                    check_dis = max(check_dis1, check_dis2, check_dis3)
                    print(check_dis)
                    loop_ = loop_ + 1
                   
                else:
                    break
    
            tem_cluster_index = np.argmax(conditional_probalility_, axis=1)
            Collection_label.append(tem_cluster_index)
            
            sum_of_prob.append(
                ((conditional_probalility_o*conditional_probalility_).sum()+sum(con_ref_prob)))
            method_index_c.append(method_index)
            Best_ref_collection[method_index] = copy.deepcopy(Best_ref)
            parameters_saving[method_index] = [gamma1,gamma2,gamma3,gamma4,pcorrect,p_mis_loc,pins1_loc,lamuna_ins,pdel1_loc,lamuna_del]
        
    ####
    sum_of_prob_max = max(sum_of_prob)
    sum_of_prob_max_index = sum_of_prob.index(sum_of_prob_max)
    # finish the estimation
    method_index_select = method_index_c[int(sum_of_prob_max_index)]
    cluster_index = Collection_label[int(sum_of_prob_max_index)]
    
    making_up_con = ['' for x in range(int(proposed_k))]
    for s_con in range(proposed_k):
        for lll_loc in range(len(Best_ref_collection[method_index_select][s_con])):
            making_up_con[s_con] = making_up_con[s_con] + Standardbase[Best_ref_collection[method_index_select][s_con][lll_loc]]
    
    ## split the best ref inferred into bone and ins part
    Stand_bone_infer = [{} for x in range(proposed_k)]
    Stand_ins_infer = [{} for x in range(proposed_k)]
    bone_screen_index = 0
    ll_xun = 0
    while ll_xun <len(full_bone):
        
        for kk in range(proposed_k):
            Stand_bone_infer[kk][bone_screen_index] = copy.deepcopy(Standardbase[int(Best_ref_collection[method_index_select][kk][ll_xun])])
        
        if bone_screen_index in ins_inserted_dic:
            potential_ins_num = int(ins_inserted_dic[bone_screen_index])
            start_ins_tem = ll_xun + 1
            end_ins_tem = ll_xun + potential_ins_num
            for kk in range(proposed_k):
                tem_str = ''
                for ins_index in range(start_ins_tem,end_ins_tem+1):
                    tem_str = tem_str + Standardbase[int(Best_ref_collection[method_index_select][kk][ins_index])]
                Stand_ins_infer[kk][bone_screen_index] = tem_str
            
            ll_xun = ll_xun + potential_ins_num
        ll_xun = ll_xun + 1
        bone_screen_index = bone_screen_index + 1   
    ##### MEC score Traning######
    num_wrong_mec_full = 0
    num_full = 0
    for read_index_xun_mec in range(len(final_indicator_table)):
        tem_dic_mec = rselct_result[int(read_index_xun_mec)][2]
        begin_id_mec = final_indicator_table.iloc[read_index_xun_mec, 1]
        end_id_mec = final_indicator_table.iloc[read_index_xun_mec, 2]
        #tem_values = np.zeros(K)
        kk_read_id = int(reads_infer_clu.iloc[read_index_xun_mec, 1])

        # for kk in range(K):
        #kk = 0
        #num_correct_mec = 0
        #tem_mec = np.zeros(proposed_k)
        num_compare = np.zeros(proposed_k)
        #for kk_read_id in range(proposed_k):
            
        num_wrong_mec = 0

        base_xun = int(begin_id_mec)
        compare_base_number = 0
        while base_xun <= end_id_mec:
            if Best_ref_collection[int(method_index_select)][kk_read_id][base_xun] < 4:
                #bone_flag = 1
                if tem_dic_mec[base_xun] != Best_ref_collection[int(method_index_select)][kk_read_id][base_xun]:

                    num_wrong_mec = num_wrong_mec + 1
                compare_base_number = compare_base_number + 1
            elif Best_ref_collection[int(method_index_select)][kk_read_id][base_xun] == 4:
                if tem_dic_mec[base_xun] != 4:
                    if tem_dic_mec[base_xun] < 4:
                        num_wrong_mec = num_wrong_mec + 1
                        compare_base_number = compare_base_number + 1
                    else:
                        print('error1')
            else:
                if tem_dic_mec[base_xun] == 4:
                    print('error2')
                elif tem_dic_mec[base_xun] < 4:
                    num_wrong_mec = num_wrong_mec + 1
                    compare_base_number = compare_base_number + 1
            base_xun = base_xun + 1
        #tem_mec[kk_read_id] = num_wrong_mec
        num_wrong_mec_full = num_wrong_mec_full+num_wrong_mec
        num_full = num_full+compare_base_number
    print( num_wrong_mec_full/num_full)
   
    infer_var_table =pd.DataFrame(columns=['chrom','pos','type','content'])
    infer_loop = 0
    for kk_infer in range(proposed_k):
        bone_check_loopl = 0
        tem_seq_bone = ''
        loop_now = 0
        tem_infer_ins_dic = {}
        while bone_check_loopl <= n_c-1:
            tem_seq_bone = tem_seq_bone + making_up_con[kk_infer][int(loop_now)]
            if making_up_con[kk_infer][int(loop_now)] != BestRefSeq[bone_check_loopl]:
                if making_up_con[kk_infer][int(loop_now)]=='-':
                    
                    infer_var_table.loc[infer_loop] = [kk_infer,bone_check_loopl,'Del','-']
                else:
                    infer_var_table.loc[infer_loop] = [kk_infer,bone_check_loopl,'Mis',making_up_con[kk_infer][int(loop_now)]]
                infer_loop = infer_loop + 1
           
            if bone_check_loopl in ins_inserted_dic:
                
                tem_ins_string = making_up_con[kk_infer][int(loop_now+1):int(1+loop_now+ins_inserted_dic[bone_check_loopl])]
                tem_ins_string_refined = tem_ins_string.replace('n','')
                if len(tem_ins_string_refined)>0:
                    tem_infer_ins_dic[bone_check_loopl] = tem_ins_string_refined
                    infer_var_table.loc[infer_loop] = [kk_infer,bone_check_loopl,'Ins',tem_ins_string_refined]
                    infer_loop = infer_loop + 1
                loop_now = loop_now + ins_inserted_dic[bone_check_loopl]
            
                    
            
            loop_now = loop_now +1
            bone_check_loopl = bone_check_loopl + 1
       
                
                   

    ### saving the test part
    infer_var_table.to_csv(main_url_save+'inference_'+str(proposed_k)+'.txt',index=0)
