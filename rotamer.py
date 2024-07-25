#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 16 11:38:18 2022

@author: avinashmandaiya
"""
import os
import numpy as np    
import sys
from Bio.PDB import *
import _pickle as pickle
from numpy import random
import matplotlib.pyplot as plt
import math
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.DSSP import DSSP
import warnings
from sklearn.decomposition import PCA
import timeit

warnings.resetwarnings()
warnings.filterwarnings(action="ignore", message="unclosed", category=ResourceWarning,append=True)
warnings.filterwarnings(action="error",append=True)

#%%
def residue_type(AA):
    AA_list = "RNDQEHKCMGPSTWYILVFA"
    result = AA_list.find(AA)
    
    return result

def trans_origin(ini_data,no_rows):
    
    data_p = np.zeros((no_rows,3),dtype = float)
    
    data_p = ini_data[0:no_rows][:]

    avg = data_p.sum(axis=0)/no_rows
    avg_matrix = np.tile(avg, (no_rows,1)) 
    data_p = data_p - avg_matrix
    
    ini_data[0:no_rows][:] = data_p
    
    return ini_data



def motif_analyser(ini_data,num,edge_thres):
    edge_matrix = np.zeros((num,num),int)
    
    no_rows = ini_data[0].shape[0]
    no_cols = ini_data[0].shape[1]
    # temp = np.zeros((no_rows,no_cols),float)
    dis_matrix = np.zeros((num,num),float)
    size = int((num*num-num)/2)
    dis_array = np.zeros(size)
    
    counter = 0
    for i in range(num):
        for j in range(i+1,num):
            # temp = ini_data[i]-ini_data[j]
            dis_matrix[i][j] = motif_dis(ini_data[i],ini_data[j])
            
            dis_array[counter] = dis_matrix[i][j]
            if dis_array[counter] < edge_thres:
                edge_matrix[i][j] = 1
                edge_matrix[j][i] = 1
            counter += 1
            
    fig, ax = plt.subplots(figsize =(10, 7))
    ax.hist(dis_array,bins = 100)
    plt.show()     

    return edge_matrix       

    
def pca_analyser(ini_data,num,state):
    no_rows = ini_data[0].shape[0]
    no_cols = ini_data[0].shape[1]
    
    input_data = np.zeros((num,no_rows*no_cols),float)
    
    for i in range(num):
        for j in range(no_rows):
            for k in range(no_cols):
                input_data[i][j*no_cols+k] = ini_data[i][j][k]

    pca = PCA(n_components=2)
    PCA_data = pca.fit_transform(input_data)  
    
    fig, ax = plt.subplots()
    for i in range(num):
        if state[i] == 2:
            ax.scatter(PCA_data[i,0], PCA_data[i,1],s=2, c ='red')
            
    fig1, ax1 = plt.subplots()
    ax1.scatter(PCA_data[:,0], PCA_data[:,1],s=2, c ='blue')
    
    plt.show()
    
    
def only_pca_analyser(ini_data,num):
    no_rows = ini_data[0].shape[0]
    no_cols = ini_data[0].shape[1]
    
    input_data = np.zeros((num,no_rows*no_cols),float)
    
    for i in range(num):
        for j in range(no_rows):
            for k in range(no_cols):
                input_data[i][j*no_cols+k] = ini_data[i][j][k]

    pca = PCA(n_components=2)
    PCA_data = pca.fit_transform(input_data)  
    
    return PCA_data
            

    
    
def dominating_set(edge_matrix):
    size = edge_matrix.shape[0]
    state = np.zeros(size,int)
    nonds_edges = np.zeros(size,int)
    
    unds = size
    
    while unds!=0:        
  
        unds = size
        for i in range(size):
            
            if state[i] != 2:
                num = 0
                for j in range(size):
                    if i==j:
                        continue
                    if edge_matrix[i][j] == 1 and state[j] == 0:
                        num += 1 
                
                nonds_edges[i] = num
                
                if num == 0 and state[i] == 0:
                    state[i] = 2
                    unds -= 1
                
            else:
                nonds_edges[i] = 0 
                unds -= 1
                
        if unds == 0:
            break
        
        max_nonds_edges = -1
        for i in range(size):
            if state[i] != 2 and nonds_edges[i] > max_nonds_edges:
                max_nonds_edges = nonds_edges[i]
                max_edge_node = i
                
        state[max_edge_node] = 2
        
        for i in range(size):
            if max_edge_node == i:
                continue
            
            elif edge_matrix[max_edge_node,i] == 1 and state[i] == 0:
                state[i] = 1
        
        unds = 0
        for i in range(size):
            if state[i] == 0:
                unds += 1
                
    return state

def motif_transform_mean(ini_data,count):
    
    no_motifs = count
    num_points = ini_data.shape[1]
    mean_block = np.zeros((num_points,3))
    
    for i in range(no_motifs):
        mean_block += (ini_data[i,:,0:3]).astype(float)
        
    mean_block /= no_motifs
    
    mean_block = trans_origin(mean_block,num_points)
    
    for i in range(no_motifs):
        real_data = (ini_data[i,:,0:3]).astype(float)
        real_data = trans_origin(real_data,num_points)
        Q = np.transpose(real_data)
        
        if mean_block.shape  == real_data.shape:
            if mean_block.shape[1] == 3:
                matrix_M = Q@mean_block
            else:
                print("vectors are not 3-dimensional")
                
        U, S, WT = np.linalg.svd(matrix_M,full_matrices= True)
        
        det = np.linalg.det(np.dot(U,np.transpose(WT)))
        S_prime = np.zeros((3, 3), int)
        np.fill_diagonal(S_prime, 1)
        
        if det > 0:
            S_prime[2,2] = 1
        else:
            S_prime[2,2] = -1
        
        final_data = np.transpose(Q)@(U@S_prime@WT)
        ini_data[i,:,0:3] = final_data[:,:]
        
    return ini_data



def motif_dis(m1,m2):
    
    no_rows = m1.shape[0]
    no_cols = m2.shape[1]
    
    motif1 = trans_origin(m1,no_rows)
    motif2 = trans_origin(m2,no_rows)
        
    Q = np.transpose(motif2)
    
    if motif1.shape  == motif2.shape:
        if motif1.shape[1] == 3:
            matrix_M = Q@motif1
        else:
            print("vectors are not 3-dimensional")
            
    U, S, WT = np.linalg.svd(matrix_M,full_matrices= True)
    
    det = np.linalg.det(np.dot(U,np.transpose(WT)))
    S_prime = np.zeros((3, 3), int)
    np.fill_diagonal(S_prime, 1)
    
    if det > 0:
        S_prime[2,2] = 1
    else:
        S_prime[2,2] = -1
    
    proj_motif = np.transpose(Q)@(U@S_prime@WT)
    temp = motif1-proj_motif
    dis = np.linalg.norm(temp)
             
    return dis/np.sqrt(no_rows*no_cols)


def motif_dis_b(m1,m2):
    
    no_rows = m1.shape[0]
    no_cols = m2.shape[1]
    
    mp1 = m1[0:3,:]
    mp2 = m2[0:3,:]
    
    avg = mp1.sum(axis=0)/no_rows
    avg_matrix = np.tile(avg, (no_rows,1)) 
    m1 = m1 - avg_matrix    
    
    avg = mp2.sum(axis=0)/no_rows
    avg_matrix = np.tile(avg, (no_rows,1)) 
    m2 = m2 - avg_matrix  
    
    motif1 = trans_origin(mp1,3)
    motif2 = trans_origin(mp2,3)
        
    Q = np.transpose(motif2)
    
    if motif1.shape  == motif2.shape:
        if motif1.shape[1] == 3:
            matrix_M = Q@motif1
        else:
            print("vectors are not 3-dimensional")
            
    U, S, WT = np.linalg.svd(matrix_M,full_matrices= True)
    
    det = np.linalg.det(np.dot(U,np.transpose(WT)))
    S_prime = np.zeros((3, 3), int)
    np.fill_diagonal(S_prime, 1)
    
    if det > 0:
        S_prime[2,2] = 1
    else:
        S_prime[2,2] = -1
    
    proj_motif = m2@(U@S_prime@WT)
    temp = m1-proj_motif
    dis = np.linalg.norm(temp)
             
    return np.sqrt(dis/(no_rows*no_cols))
#%%

max_length = 14
max_motifs = 10000 
lenAA = np.zeros(20,dtype = int)
    
lenAA[0] = 11
lenAA[1] = 8  
lenAA[2] = 8  
lenAA[3] = 9  
lenAA[4] = 9  
lenAA[5] = 10 
lenAA[6] = 9  
lenAA[7] = 6  
lenAA[8] = 8  
lenAA[9] = 4  
lenAA[10] = 7  
lenAA[11] = 6  
lenAA[12] = 7  
lenAA[13] = 14 
lenAA[14] = 12 
lenAA[15] = 8  
lenAA[16] = 8  
lenAA[17] = 7  
lenAA[18] = 11 
lenAA[19] = 5 

AA_list = "RNDQEHKCMGPSTWYILVFA"

rot_motifs = np.zeros((20,max_motifs,max_length,3),dtype = float)
count_motifs = np.zeros(20,dtype = int)

for i in range(20):
    lenAA[i] -= 1

#%%
full_path = os.path.realpath(__file__)
os.chdir(os.path.dirname(full_path))
file1 = open('list_proteins.txt', 'r')

Line = file1.readline()
Line = Line[:-1]
file1.close()
pdb_ids = Line.split(',')

os.chdir('high-res_proteins')

np.seterr(all='raise')

used_proteins = 0
unused_proteins = 0
count_same = 0

start = timeit.default_timer()
# for i in range(0,len(pdb_ids)):
for i in range(1,2):    
    pdb_id = pdb_ids[i]
    filename = pdb_id+".pdb"     

    if i%100 == 0:
        print(i)
        end = timeit.default_timer()
        print(end-start)         
        print("\n\n\n\n\n\n")    

    try:
        warnings.filterwarnings(action="error", message="Chain")
        mmcif_dict = MMCIF2Dict.MMCIF2Dict(filename)
        parser = PDBParser()
        structure = parser.get_structure(pdb_id, filename)  
        used_proteins += 1
    except:
        unused_proteins += 1
        continue
    
    model = structure[0]
    dssp = DSSP(model, filename, dssp='mkdssp')    # for the dssp module
    
    id_chain = list(model.child_dict.keys())[0]
    main_chain = model[id_chain]
    
    
    all_atoms = Selection.unfold_entities(structure, "A")
    all_residues = Selection.unfold_entities(main_chain, "R")

    for residue1 in all_residues:
 
        if (residue1.get_parent()).get_id() != id_chain:
            continue

        if residue1.get_id()[0] != ' ':
            continue    

        try:
            a_key1 = (id_chain,residue1.get_id())
            a1 = residue_type(dssp[a_key1][1])
        except:
            continue

        if count_motifs[a1] >= 10000:
            continue
            
        num_atoms = 0 
        
        for atom in residue1:
            if atom.get_id()[0] not in "CNOS" or atom.get_id() == 'OXT' or atom.get_id() == 'O':
                continue
                   
            rot_motifs[a1][count_motifs[a1]][num_atoms][:] = atom.get_coord()
            num_atoms = num_atoms+1
            
        if num_atoms == lenAA[a1]:
            same = 0
            # for motif_num in range(count_motifs[a1]):
            #     new_dis = motif_dis(rot_motifs[a1][count_motifs[a1]][0:lenAA[a1]],rot_motifs[a1][motif_num][0:lenAA[a1]])
            #     if new_dis < 0.1:
            #         same = 1
            #         break
                
            if same == 0: 
                trans_origin(rot_motifs[a1][count_motifs[a1]],lenAA[a1])
                count_motifs[a1] += 1
            
#%%
"""PCA analysis"""
os.chdir(os.path.dirname(full_path))
figs, axs = plt.subplots(4,5)
# axs = figs.gca()
# fig, ax = plt.subplots()
# figs = plt.gcf()
# figs.set_size_inches(18.5, 14.5)


for a1 in range(20):
    if count_motifs[a1] == 0:
        print(a1)
        continue
    
    temp_motifs = np.zeros((count_motifs[a1],lenAA[a1],3),dtype = float)
    temp_motifs = rot_motifs[a1,0:count_motifs[a1],0:lenAA[a1],:]

    # motif_transform_mean(temp_motifs,count_motifs[a1])

    # edges = motif_analyser(temp_motifs,count_motifs[a1],0.05)
    # ds_state = dominating_set(edges)  
    
    # pca_analyser(temp_motifs,count_motifs[a1],ds_state) 

    # pca_data = only_pca_analyser(temp_motifs,count_motifs[a1]) 
    
    # fig, ax = plt.subplots()
    # axs[a1//5,a1%5].scatter(pca_data[:,0], pca_data[:,1],s=0.2, c ='blue')
    # axs[a1//5,a1%5].title.set_text('Residue: {}'.format(AA_list[a1]))
    # ax.scatter(pca_data[:,0], pca_data[:,1],s=0.5, c ='blue')
    # plt.title('PCA: AA = {}'.format(AA_list[a1]))
    
# plt.savefig('plot.png', dpi=300, bbox_inches='tight')    
# figs.savefig('rot_PCA.eps')
# figs.savefig( dpi=100)     
# plt.show()

#%%
"""writing backbone-solvent library"""
os.chdir(os.path.dirname(full_path))
main_protein = "5DMA"
AA_list = "RNDQEHKCMGPSTWYILVFA"

filename = main_protein+"_rotV1.txt"

with open(filename, 'w') as fh:
    for i in range(20):
                        
            if count_motifs[i] == 0:
                continue
            
            for z in range(count_motifs[i]):

                fh.write('{}\n'.format(AA_list[i]))
            
                for q in range(lenAA[i]):
                    fh.write('{:.6f} \t {:.6f} \t {:.6f}\n'.format(rot_motifs[i,z,q,0], rot_motifs[i,z,q,1], rot_motifs[i,z,q,2]))                   
                fh.write('\n')             