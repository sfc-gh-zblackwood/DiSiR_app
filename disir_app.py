import streamlit as st

import pandas as pd
import numpy as np

import sys

import os
#import subprocess

append = '/'.join(os.path.abspath(__file__).split('/')[:-1]) + '/'
sys.path.append(append)
import disir_app_functions as ds
np.random.seed(0)

st.title('DiSiR')

col1, col2 = st.columns(2)

uploaded_mtx = st.file_uploader("Provide path to matrix.mtx", key="mtx")

uploaded_genes = col1.file_uploader("Upload genes.txt")
uploaded_metadata = col1.file_uploader("Upload metadata (csv or json)")
select_track = col1.text_input('Metadata column with cell annotatation', value='CellStates')

threshold_fractions = col2.number_input('min fraction of cells in a cell type that expresses the gene', min_value = 0.0, max_value=1.0, value=0.0, step=0.1)
threshold_expressions = col2.number_input('threshold on scaled (max-normalized) average expression of each ligand or receptor within a cell type', min_value = 0.0, value=0.0, step=0.1)
threshold_number = col2.number_input('min number of cells in a cell type that expresses the gene', min_value = 0, value=0)
iteration = col2.number_input('number of permutations', min_value=1, max_value=1000000,value=500, step=100)

EXCLUDE_ITEMS = ['Unclassified']
THRESHOLD = 0.05 # p-value, essentially

################### RUN APP HERE ####################

if uploaded_mtx and uploaded_genes and uploaded_metadata::
    gene_names_all = list(pd.read_csv(uploaded_genes, header = None, index_col = None)[0])
    ligand_names = st.text_input('Please choose the ligands of interest in order, separated by commas').replace(' ', '').split(',')
    receptor_names = st.text_input('Please choose the receptor of interest in order, separated by commas').replace(' ', '').split(',')
    
    if 'metadata' not in st.session_state and 'scRNA' not in st.session_state:
        data_load_state = st.text('If you are loading your matrix for the first time, it might take a while...')
        scRNA_array, metadata = ds.load_data(uploaded_mtx, uploaded_metadata)
        st.session_state['metadata'] = metadata
        st.session_state['scRNA'] = scRNA_array
        data_load_state.text('Done!')
    else:
        metadata = st.session_state['metadata']
        scRNA_array = st.session_state['scRNA']
    
    categories = sorted(metadata.columns)
    select_cat = st.selectbox('Choose subsetting category if applicable', categories)
    
    choices = sorted(set(metadata[select_cat]))
    select_choice = st.multiselect('Choose item from the category', choices)
    
    subset = {}
    if select_choice:
        subset = {select_cat : select_choice}

    launch = st.button("Let's go!")
    
    # this is where it really starts to run
    if ((ligand_names and receptor_names)
        and (len(set(ligand_names) & set(gene_names_all))==len(set(ligand_names))) 
        and (len(set(receptor_names) & set(gene_names_all))==len(receptor_names))
        and (len(ligand_names)==len(receptor_names)) 
        and launch):
        
        work_text = st.text('Working....')
      
        
        scRNA_array, cell_type_labels = ds.filter_out(scRNA_array,
                                                   metadata,
                                                   select_track,
                                                   EXCLUDE_ITEMS,
                                                   subset)
        
        gene_names = sorted(set(ligand_names) | set(receptor_names))
        subunit_info = [f'{i} | {j}' for i,j in zip(ligand_names, receptor_names)]
        gene_index = [gene_names_all.index(i) for i in gene_names]
        scRNA_array = scRNA_array[gene_index, : ]

        
        average_expression_df, number_expression_df, totalcell_df, cell_type_specific_expressions_df, unique_cell_type_labels, cell_type_numbers = ds.calculate_celltype_average_expressions(scRNA_array,
                                            gene_names,
                                            cell_type_labels)


        interactions_prefiltered, interactions_gene_names, interactions_class_names = ds.calculate_LR_interactions_matrix_prefiltered(average_expression_df,
                                    gene_names,
                                    unique_cell_type_labels,
                                    totalcell_df,
                                    cell_type_specific_expressions_df,
                                    number_expression_df,
                                    threshold_fractions,
                                    threshold_expressions,
                                    threshold_number)


        p_values_matrix, p_values_gene_names, p_values_celltype_names, expression_celltype_matrix, output_table_results = ds.calculate_LR_interactions_pvalues(interactions_prefiltered,
                                           interactions_gene_names,
                                           interactions_class_names,
                                           scRNA_array,
                                           cell_type_numbers,
                                           gene_names,
                                           unique_cell_type_labels,
                                           iteration)

        
        graph_data_nodes, graph_data_links = ds.build_interactions_graph_subunit(p_values_matrix,
                                     p_values_gene_names,
                                     unique_cell_type_labels,
                                     expression_celltype_matrix,
                                     subunit_info,
                                     THRESHOLD)
        
        
        heatmap_dict = ds.calculate_celltype_interactions_heatmap(graph_data_links,
                                                    unique_cell_type_labels,
                                                    average_expression_df,
                                                    subunit_info)
        
        ds.plot_heatmaps(heatmap_dict, unique_cell_type_labels)


        ds.plot_cell_dist(average_expression_df,
                          number_expression_df,
                          totalcell_df)
        
        work_text.text('Done working!')
        del st.session_state['metadata']
        del st.session_state['scRNA']
        
    elif len( (set(ligand_names) | set(receptor_names)) - set(gene_names_all) ) >= 1:
        not_in_dataset = (set(ligand_names) | set(receptor_names)) - set(gene_names_all)
        if not_in_dataset and not_in_dataset != set(['']):
            col1.text(f'{not_in_dataset} not present in this dataset!')        
                
    elif ligand_names and receptor_names and launch:
        col1.text('Please make sure that the number of ligands and receptors is the same')
        col1.text(f'The ligands are {ligand_names}')
        col1.text(f'The corresponding receptors are {receptor_names}')
 
elif not uploaded_mtx or not uploaded_genes or not uploaded_metadata:
    data_load_state = st.text('Please load the matrix, the genes, and the metadata (csv or json)')


######################################################
