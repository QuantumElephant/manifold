# -*- coding: utf-8 -*-
"""
Created on Sat Aug 12 04:14:24 2017

@author: SKY
"""

from sklearn.decomposition import PCA
import numpy as np
 
 
 
n=int(3)               #power of number of samples
m=int(3)               #number of molecules

 
 
 
coordi_internal=np.loadtxt('coordi_3_internal.txt') 
co_internal=coordi_internal[0:10**n]
U, s, V = np.linalg.svd(coordi_internal, full_matrices=True) 
 
pca = PCA(n_components=3)
X_r = pca.fit(coordi_internal).transform(coordi_internal)



coordi_internal=np.loadtxt('coordi_4_internal.txt')
co_internal=coordi_internal[0:10**n]
U, s, V = np.linalg.svd(coordi_internal)


