#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 07:40:30 2021

@author: avinashmandaiya
"""
import numpy as np
import timeit
import os
import pandas as pd
import json
from string import ascii_uppercase
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import _pickle as pickle
from sklearn.decomposition import PCA

#%%
full_path = os.path.realpath(__file__)
os.chdir(os.path.dirname(full_path))
os.chdir('../test_libs')

#%%
def residue_type(AA):
    HPAA = "ACILMFV"
    result = HPAA.find(AA)
    
    if result == -1:
        return 0
    else:
        return result+1
    
def ss_type3(s):
    
    if s == 'H' or s == 'G' or s == 'I':
        return 0
    elif s == 'E' or s == 'B':
        return 1
    else:
        return 2
        
#%%

def retrieve_files(file):
    
    filename = file+".txt"
    f = open(filename, "r")
    
    if file == "MlibBbone_PCA":
        num_points = 8
    elif file == "MlibSolA_PCA" or file == "MlibSolA_PCA":
        num_points = 5
        
    max_rows = 5000
    no_col = 3*num_points
    main_matrix = np.zeros((max_rows,no_col)) 
    ss_type = np.zeros(max_rows,dtype = str)  
    resdiff = np.zeros(max_rows,dtype = int)  
    
    ini = 0
    i = 0 
    
    for x in f:
        l = len(x)
        
        if x == "\n":
            i = i+1
            ini = 0
            continue
        
        elif x != "\n" and l < 10:
            continue
            # if num_points == 8:
            #     rt_types[i][0] = residue_type(x[0])
            #     rt_types[i][1] = residue_type(x[1])
            # else:
            #     rt_types[i][0] = residue_type(x[0])
                
        else:
            numbers = [float(x) for x in x.split()[0:3]] 
            if ini == 6:
                res1 = int(x.split()[3])
                ss_type[i] = x.split()[5]
                # ss_type[i] = ss_type3(x.split()[5])
            if ini == 12:
                resdiff[i] = abs(int(x.split()[3]) - res1)
            main_matrix[i,ini:ini+3] = numbers[0:3]
            ini = ini+3            
            
    
    main_matrix = np.delete(main_matrix,slice(i, max_rows),axis = 0) 
    ss_type = np.delete(ss_type,slice(i, max_rows),axis = 0)  
    resdiff = np.delete(resdiff,slice(i, max_rows),axis = 0) 

    return main_matrix,ss_type,resdiff   
#%%
file = "MlibBbone_PCA"
main_matrix_bb, ss_type, resdiff= retrieve_files(file)

pca = PCA(n_components=2)
PCA_data = pca.fit_transform(main_matrix_bb)

#%%
"""writing the PCA data in text files to plot in mathematica"""
    
def PCA_export(filename, PCA_data):
    textfile = filename + ".txt"
    no_motifs = PCA_data.shape[0] 
    # no_comp = PCA_data.shape[1]
    
    with open(textfile, 'w') as fh:
        
        for k in range(no_motifs):
            fh.write('{:.6f} \t {:.6f} \t {:.6f}\n'.format(PCA_data[k,0],PCA_data[k,1],PCA_data[k,2]))  

#%%
os.chdir(os.path.dirname(full_path))

num_motifs = PCA_data.shape[0]

fig, ax1 = plt.subplots()

color_map = {2: 'cyan', 3: 'violet', 4: 'blue', 5: 'purple', 6: 'red'}
colors_map1 = {'H':'violet','G':'cyan','I':'purple','A':'red','P':'green','B':'orange','S':'black','T':'lime','-':'blue'}
# colors_map2 = {0:'red',1:'blue',2:'blue',3:'blue',4:'blue',5:'blue',6:'blue',7:'blue'}
# ax.scatter(PCA_data[:,0], PCA_data[:,1],c=colors_map2.get(ss_type[i]),,s=2)

for i in range(num_motifs):
    if resdiff[i] > 5:
        resdiff[i] = 6
    ax1.scatter(PCA_data[i,0], PCA_data[i,1],c= color_map.get(resdiff[i]), s=2)

    
# plt.savefig('BB.png', format='png') 
plt.legend()
plt.show()

#%%
fig, ax1 = plt.subplots()
ax1.scatter(PCA_data[np.where(resdiff == 2),0], PCA_data[np.where(resdiff == 2),1], c= "cyan", s=2, label = "Turn")
#%%
fig, ax = plt.subplots()
ax.hist(a, bins = [0, 25, 50, 75, 100])
#%%
# PCA_data_lhc = pca.transform(main_matrix_bb[0:count_local])

# num_motifs = PCA_data_lhc.shape[0]

# fig, ax = plt.subplots()

# colors_map1 = {0:'red',1:'green',2:'blue',3:'yellow',4:'violet',5:'cyan',6:'lime',7:'orange'}
# colors_map2 = {0:'red',1:'blue',2:'blue',3:'blue',4:'blue',5:'blue',6:'blue',7:'blue'}
# # ax.scatter(PCA_data[:,0], PCA_data[:,1])
# for i in range(num_motifs):
#     ax.scatter(PCA_data_lhc[i,0], PCA_data_lhc[i,1], c=colors_map0.get(rt_types_lhc[i][0]),s=2)
#     ax.set_ylim(-7, 7)
#     ax.set_xlim(-6.5,9)
    
# plt.savefig('l_pca.eps', format='eps') 
    
# plt.show()

#%%
# PCA_data_lhc = pca.transform(main_matrix_bb[count_local:count_local+count_nonlocal])

# num_motifs = PCA_data_lhc.shape[0]

# fig, ax = plt.subplots()

# colors_map1 = {0:'red',1:'green',2:'blue',3:'yellow',4:'violet',5:'cyan',6:'lime',7:'orange'}
# colors_map2 = {0:'red',1:'blue',2:'blue',3:'blue',4:'blue',5:'blue',6:'blue',7:'blue'}
# # ax.scatter(PCA_data[:,0], PCA_data[:,1])
# for i in range(num_motifs):
#     ax.scatter(PCA_data_lhc[i,0], PCA_data_lhc[i,1], c=colors_map2.get(rt_types_lhc[i][0]),s=2)
#     ax.set_ylim(-7, 7)
#     ax.set_xlim(-6.5,9)
    
# plt.savefig('nl_pca.eps', format='eps') 
# plt.show()
#%%
# """plotting in 3D"""
# from mpl_toolkits import mplot3d


# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.scatter3D(PCA_data[:,0], PCA_data[:,1], PCA_data[:,2])
# ax = fig.add_subplot(1, 1, 1)
# ax.scatter(PCA_data[:,0], PCA_data[:,1])
# plt.title('PCA of {} library'.format(file))
# pp = PdfPages('no_bond.pdf')
# pp.savefig()
# pp.close()




