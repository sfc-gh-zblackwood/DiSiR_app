#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 14:56:42 2021

@author: E0463430
"""
import numpy as np

from math import ceil
from math import log
from math import log2

import pandas as pd
import streamlit as st
from scipy.io import mmread
from statsmodels.stats.multitest import multipletests

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as pch
from matplotlib import rcParams
mpl.use('Agg')
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams.update({'font.size': 15})


@st.cache
def load_data(input_mtx, uploaded_metadata):
    '''this function loads all the data'''
    scRNA_data = mmread(input_mtx)
    scRNA_array = scRNA_data.toarray()
    scRNA_array = np.transpose(scRNA_array)
    if '.json' in uploaded_metadata.name:
        meta_in = pd.read_json(uploaded_metadata)
        metadata = {}
        for metakey in sorted(meta_in):
            metadata[metakey] = meta_in[metakey]['label_list']
        metadata = pd.DataFrame(metadata)
        metadata.index.name = 'cell_id'
    else:
        metadata = pd.read_csv(uploaded_metadata, index_col=0)
    return scRNA_array, metadata


def filter_out(scRNA_array,
               metadata,
               select_track,
               exclude_items,
               subset):
    'filter out the data based on the empty tracks (hard-coded) and user-defined category'''
    cell_type_labels = metadata[select_track]
    cell_type_keep = [i[0] for i in enumerate(cell_type_labels) if i[1] not in exclude_items]
    
    if subset and subset != "None":
        subset_track = list(subset.keys())[0]
        subset_labels = metadata[subset_track]
        subset_cat  = subset[subset_track]
        cell_type_keep = [i for i, j in enumerate(subset_labels)
                          if i in cell_type_keep and j in subset_cat]
        
    cell_type_labels = [i[1] for i in enumerate(cell_type_labels) if i[0] in cell_type_keep]
    scRNA_array = scRNA_array[ : , cell_type_keep]
    return scRNA_array, cell_type_labels


def calculate_celltype_average_expressions(expressions_scRNA,
                                           gene_names,
                                           cell_type_labels):
    '''creates relevant data frames'''
    # Read meta data (cell type labels)
    unique_cell_type_labels = np.unique(cell_type_labels)
    noClasses = len(unique_cell_type_labels)
    cell_type_labels_df = pd.DataFrame(cell_type_labels)
    gene_number = len(gene_names)
    # Calculate cell type average expression matrix
    average_expression = np.zeros([gene_number,noClasses])
    number_expression = np.zeros([gene_number,noClasses])
    totalcell = np.zeros([gene_number,noClasses])
    cell_type_numbers = [0]*noClasses
    for k in range(noClasses):
        cell_type_index = cell_type_labels_df.isin([unique_cell_type_labels[k]])
        seriesObj = cell_type_index.any(axis = 1)
        columnNames = list(seriesObj[seriesObj == True].index) 
        cell_type_numbers[k] = len(columnNames)
        for i in range(gene_number):
            selected_cells = np.concatenate([expressions_scRNA[i,columnNames], [0]*50])
            average_expression[i,k] = np.mean(expressions_scRNA[i,columnNames])
            number_expression[i,k] = len(selected_cells[selected_cells > 0])
            totalcell[i,k]  = len(expressions_scRNA[i,columnNames])
    
    average_expression_df = pd.DataFrame(average_expression, index = gene_names, columns = unique_cell_type_labels)
    number_expression_df = pd.DataFrame(number_expression, index = gene_names, columns = unique_cell_type_labels)
    totalcell_df = pd.DataFrame(totalcell, index = gene_names, columns = unique_cell_type_labels)   
    cell_type_specific_expressions_df = pd.DataFrame(average_expression/np.max(average_expression), index = gene_names, columns = unique_cell_type_labels)
    
    return average_expression_df, number_expression_df, totalcell_df, cell_type_specific_expressions_df, unique_cell_type_labels, cell_type_numbers


def calculate_LR_interactions_matrix_prefiltered(average_expression_df,
                            gene_names,
                            unique_cell_type_labels,
                            totalcell_df,
                            cell_type_specific_expressions_df,
                            number_expression_df,
                            threshold_fractions,
                            threshold_expressions,
                            threshold_number):
    '''interactions matrix'''
    gene_number = len(gene_names)
    noClasses = len(unique_cell_type_labels)
    
    interactions = np.zeros(gene_number ** 2 * noClasses ** 2) 
    interactions_class_names = ['a'] * (gene_number ** 2 * noClasses ** 2)
    interactions_gene_names = ['a'] * (gene_number ** 2 * noClasses ** 2)
    counter = 0
    for i_g in gene_names:
        for j_g in gene_names:
            for i_c in unique_cell_type_labels:
                for j_c in unique_cell_type_labels:
                    interactions_class_names[counter] = i_c + ' '+ '|' + ' ' + j_c
                    interactions_gene_names[counter] = i_g + ' '+ '|' + ' '  + j_g
                    if totalcell_df.loc[i_g,i_c] >= threshold_fractions and cell_type_specific_expressions_df.loc[i_g,i_c] >= threshold_expressions and number_expression_df.loc[i_g,i_c] >= threshold_number: 
                        if totalcell_df.loc[j_g,j_c] >= threshold_fractions and cell_type_specific_expressions_df.loc[j_g,j_c] >= threshold_expressions and number_expression_df.loc[j_g,j_c] >= threshold_number:
                            interactions[counter] = average_expression_df.loc[i_g,i_c] * average_expression_df.loc[j_g,j_c]
                    counter = counter + 1
    return interactions, interactions_gene_names, interactions_class_names


def calculate_LR_interactions_pvalues(interactions,
                                   interactions_gene_names,
                                   interactions_class_names,
                                   expressions_scRNA,
                                   cell_type_numbers,
                                   gene_names,
                                   unique_cell_type_labels,
                                   iteration):
    '''perform permutations and get p-values'''
    gene_number = len(gene_names)
    noClasses = len(unique_cell_type_labels)
    # Shuffle cell labels between all cells for "iteration" times 
    interactions_shuffle = np.zeros([gene_number ** 2 * noClasses ** 2, iteration])  
    average_expression_shuffle = np.zeros([gene_number,noClasses])   
    for no_iter in range(iteration):
        for k in range(noClasses):
            columnNames = np.random.choice(range(expressions_scRNA.shape[1]), cell_type_numbers[k]).tolist()
            for i in range(gene_number):
                average_expression_shuffle[i,k] = np.mean(expressions_scRNA[i,columnNames])
        counter = 0
        for i_g in range(gene_number):
            for j_g in range(gene_number):
                for i_c in range(noClasses):
                    for j_c in range(noClasses):
                        # interactions_shuffle[counter,no_iter] = np.mean([average_expression_shuffle[i_g,i_c], average_expression_shuffle[j_g,j_c]])
                        interactions_shuffle[counter,no_iter] = average_expression_shuffle[i_g,i_c] * average_expression_shuffle[j_g,j_c]
                        counter = counter + 1
    p_values_matrix = 1000* np.ones([gene_number ** 2, noClasses ** 2])
    expression_celltype_matrix = np.zeros([gene_number ** 2, noClasses ** 2])
    p_values =  1000*np.ones([gene_number ** 2 * noClasses ** 2])
    counter = 0
    counter_1 = 0
    for i_g in range(gene_number):
        for j_g in range(gene_number):
            counter_2 = 0 
            for i_c in range(noClasses):
                for j_c in range(noClasses):
                    if interactions[counter] > 0:
                        nhigher = interactions_shuffle[counter,:][interactions_shuffle[counter,:] > interactions[counter]]
                        p = len(nhigher) / len(interactions_shuffle[counter,:])
                    else:
                        p = 1
                    p_values_matrix[counter_1,counter_2] = p
                    expression_celltype_matrix[counter_1,counter_2] = interactions[counter]
                    p_values[counter] = p
                    counter = counter + 1
                    counter_2 = counter_2 + 1
            counter_1 = counter_1 + 1
    output_table_results = pd.concat([pd.DataFrame(interactions_gene_names), 
                         pd.DataFrame(interactions_class_names),
                         pd.DataFrame(p_values), 
                         pd.DataFrame(interactions)],
                         axis = 1)    
    p_values_gene_names = ['a'] * (gene_number ** 2)
    counter = 0
    for i_g in range(gene_number):
        for j_g in range(gene_number):
            p_values_gene_names[counter] = gene_names[i_g] + ' '+ '|' + ' '  + gene_names[j_g]
            counter = counter + 1
    p_values_celltype_names = ['a'] * (noClasses ** 2)
    counter = 0
    for i_c in range(noClasses):
        for j_c in range(noClasses):
            p_values_celltype_names[counter] = unique_cell_type_labels[i_c] + ' '+ '|' + ' '  + unique_cell_type_labels[j_c]
            counter = counter + 1
    return pd.DataFrame(p_values_matrix), pd.DataFrame(p_values_gene_names), pd.DataFrame(p_values_celltype_names), pd.DataFrame(expression_celltype_matrix), output_table_results    
 
      
def build_interactions_graph_subunit(p_values_matrix,
                             p_values_gene_names,
                             unique_cell_type_labels,
                             expression_celltype_matrix,
                             selected_interactions,
                             threshold):
    '''builds intercation graph'''
    noClasses = len(unique_cell_type_labels)
    matched_gene_names = p_values_gene_names.iloc[:,0].isin(selected_interactions)
    index_matched_gene_names = list(matched_gene_names[matched_gene_names == True].index)
    selected_matched_gene_names = [i for i in list(p_values_gene_names[0]) if i in selected_interactions]
    p_values_matrix = p_values_matrix.iloc[index_matched_gene_names, :]
    p_values_matrix = p_values_matrix.values
    p_values_matrix = np.array([multipletests(i, method='fdr_bh')[1] for i in p_values_matrix])
    expression_celltype_matrix = expression_celltype_matrix.iloc[index_matched_gene_names, :]
    expression_celltype_matrix = expression_celltype_matrix.values
    node_a_list = []
    node_b_list = []
    link_list = []
    link_weight_list = []
    link_pwt_list = []
    count = 0
    count2 = 0
    for i in range(noClasses):
        for j in range(noClasses):
            all_pvals = p_values_matrix[:,count]
            significant_pvalues = all_pvals <= threshold
            interaction_list = expression_celltype_matrix[significant_pvalues == True, count]
            interaction_gene_names_list = [i for i,j in zip(selected_matched_gene_names, list(significant_pvalues)) if j == True]
            # Check if both sub units are presented:
            both_subunit_presented = []
            both_subunit_presented_exp = []
            both_subunit_presented_pvals = []
            #print(subUnit_interactions_data, 'subUnit_interactions_data')
            if len(set(selected_interactions) & set(interaction_gene_names_list)) == len(selected_interactions):
                both_subunit_presented.append('+'.join(selected_interactions))
                both_subunit_presented_exp.append(np.mean([i for i,j in zip(interaction_list, 
                                                                                interaction_gene_names_list) if j in selected_interactions]))
                both_subunit_presented_pvals.append(-log(np.mean(all_pvals)+1e-87))
            
            number_of_interactions = len(both_subunit_presented)
            for k in range(number_of_interactions):
                node_a_list.append(i + 1)
                node_b_list.append(noClasses + j + 1)
                link_weight_list.append(both_subunit_presented_exp[k])
                link_list.append(both_subunit_presented[k])
                link_pwt_list.append(both_subunit_presented_pvals[k])
                count2 = count2 + 1
            count = count + 1
    
    graph_data_links = pd.concat([pd.DataFrame(node_a_list), 
                         pd.DataFrame(node_b_list),
                         pd.DataFrame(link_weight_list),
                         pd.DataFrame(link_pwt_list),
                         pd.DataFrame(link_list)],
                         axis = 1)
    
    node_data_id = list(range(1,2*noClasses+1))
    node_data_type = [2]*len(unique_cell_type_labels) + [1]*len(unique_cell_type_labels) 
    node_data_name =  unique_cell_type_labels.tolist() + unique_cell_type_labels.tolist() 
    x = list(range(1,noClasses * 20 + 1, 20)) + list(range(1,noClasses * 20 + 1, 20))
    y = node_data_type
    graph_data_nodes = pd.concat([pd.DataFrame(node_data_id),
                         pd.DataFrame(node_data_name), 
                         pd.DataFrame(node_data_type),
                         pd.DataFrame(y),
                         pd.DataFrame(x)],
                         axis = 1)
    return graph_data_nodes, graph_data_links


def calculate_celltype_interactions_heatmap(graph_data_links,
                                            unique_cell_type_labels,
                                            average_expression_df,
                                            selected_interactions):
    '''build heatmap'''
    heatmap_dict = {}
    noClasses = len(unique_cell_type_labels)
    receptors = [i.split(' | ')[1] for i in selected_interactions]
    ligands = [i.split(' | ')[0] for i in selected_interactions]
    heatmap_dict = {}
    heatmap_dict['Heatmap (significant only)'] = np.zeros([noClasses, noClasses])
    if not graph_data_links.empty:
        for i in range(1,noClasses+1):
            for j in range(noClasses+1, 2*noClasses+1):
                index_location_x = np.where(graph_data_links.iloc[:,0] == i)
                index_location_y = np.where(graph_data_links.iloc[:,1] == j)
                index_intersect = np.intersect1d(index_location_x ,index_location_y)
                if len(index_intersect) != 0:
                    heatmap_dict['Heatmap (significant only)'][i - 1,j - (noClasses+1)] = graph_data_links.iloc[index_intersect, 2]      
    heatmap_dict['Heatmap (all_interactions)'] = np.zeros([noClasses, noClasses])
    for i, label_i in enumerate(unique_cell_type_labels):
        for j, label_j in enumerate(unique_cell_type_labels):
            expression_sum = np.mean([average_expression_df.loc[k, label_i]*average_expression_df.loc[m, label_j] for k,m in zip(ligands,receptors)])
            heatmap_dict['Heatmap (all_interactions)'][i,j] = expression_sum
    heatmap_dict['Heatmap (significant only)'] = pd.DataFrame(heatmap_dict['Heatmap (significant only)'], index = unique_cell_type_labels, columns = unique_cell_type_labels)
    heatmap_dict['Heatmap (all_interactions)'] = pd.DataFrame(heatmap_dict['Heatmap (all_interactions)'], index = unique_cell_type_labels, columns = unique_cell_type_labels)
    return heatmap_dict


def plot_heatmaps(heatmaps, cellstates):
    '''plots heatmaps'''
    for ht in heatmaps:
        heatmatrix = heatmaps[ht]
        st.text(ht)
        maxval = max(heatmatrix.max())
        fig = plt.figure(figsize=(7.5, 5.5))
        ax = plt.axes()
        pcm = ax.imshow(heatmatrix, cmap = 'YlGnBu', vmin=0, vmax=maxval)
        ax.set_xticks(range(len(cellstates)))
        ax.set_yticks(range(len(cellstates)))
        ax.set_xticklabels(cellstates, rotation=90)
        ax.set_yticklabels(cellstates)
        ax.set_xlabel('Receptor cell state')
        ax.set_ylabel('Ligand cell state')
        fig.colorbar(pcm, ax=ax, label='LxR expression')
        for val in np.arange(0.5, len(cellstates)):
            ax.axhline(val, color='0.8', lw=1)
            ax.axvline(val, color='0.8', lw=1)
        st.pyplot(fig)
        plt.close('all')


def plot_cell_dist(average_expression_df, 
                   number_expression_df,
                   totalcell_df):
    '''bubble plots'''
    smaller_frac_expression_df = number_expression_df/(totalcell_df+50)
    actual_frac_expression_df = number_expression_df/totalcell_df
    Ny = len(average_expression_df)
    Nx = len(average_expression_df.columns)
    fig = plt.figure(figsize=(Nx*1.2, Ny), dpi=1500)
    ax = plt.axes()
    ylabels = list(average_expression_df.index)
    xlabels = list(average_expression_df.columns)
    max_val = ceil(max(average_expression_df.max()))
    norm = mpl.colors.Normalize( vmin=0, vmax=max_val )
    mapper = cm.ScalarMappable(norm=norm, cmap='Reds')
    circles = []
    for i in range(Nx):
        for j in range(Ny):
            rad = log2(smaller_frac_expression_df.iloc[j, i]*100+1)/13.28
            col = mapper.to_rgba(average_expression_df.iloc[j, i])
            patch = pch.Circle((i+1, j+0.5), rad, fc=col, ec='black')
            text = f'{int(number_expression_df.iloc[j, i])} ({round(actual_frac_expression_df.iloc[j, i]*100, 1)}%)'
            ax.text(x=i+0.55, y=j+0.9, s=text, size=7.5)
            ax.add_patch(patch)
            circles.append(patch)
    ax.set(xticks=np.arange(1, Nx+2), yticks=np.arange(0.5, Ny+0.5), yticklabels=ylabels)
    ax.set_ylim(0, Ny+0.25)
    ax.set_xticklabels(xlabels+[''], rotation=80)
    plt.colorbar(mapper, ax=ax, label='log(expression)', pad=0.025)
    st.pyplot(fig)
    plt.close('all')
    
####################################################################################################
  
