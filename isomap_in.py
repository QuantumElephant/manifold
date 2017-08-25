# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 05:41:57 2017

@author: SKY
"""
import numpy as np
from scipy.linalg import orthogonal_procrustes
from scipy.spatial import distance
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
from sklearn import manifold, datasets

def procrustes(data1, data2):
    
    mtx1 = np.array(data1, dtype=np.double, copy=True)
    mtx2 = np.array(data2, dtype=np.double, copy=True)

    if mtx1.ndim != 2 or mtx2.ndim != 2:
        raise ValueError("Input matrices must be two-dimensional")
    if mtx1.shape != mtx2.shape:
        raise ValueError("Input matrices must be of same shape")
    if mtx1.size == 0:
        raise ValueError("Input matrices must be >0 rows and >0 cols")

    # translate all the data to the origin
    mtx1 -= np.mean(mtx1, 0)
    mtx2 -= np.mean(mtx2, 0)

    norm1 = np.linalg.norm(mtx1)
    norm2 = np.linalg.norm(mtx2)

    if norm1 == 0 or norm2 == 0:
        raise ValueError("Input matrices must contain >1 unique points")

    # change scaling of data (in rows) such that trace(mtx*mtx') = 1
    #mtx1 /= norm1
    #mtx2 /= norm2

    # transform mtx2 to minimize disparity
    R, s = orthogonal_procrustes(mtx1, mtx2)    #R=The matrix solution of the orthogonal Procrustes problem.
    mtx2 = np.dot(mtx2, R.T)                #s=Sum of the singular values of ``dot(A.T, B)``
    # measure the dissimilarity between the two datasets
    distance_mol = math.sqrt(np.sum(np.square(mtx1 - mtx2)))

    return mtx1, mtx2, distance_mol

internal_coordi=np.loadtxt('coordi_3_internal.txt')
coordi=np.loadtxt('coordi_3_10^5.txt')

#iso = manifold.Isomap(n_neighbors=4, n_components=3)
#internal_coordi=internal_coordi[1:2000,:]
#iso.fit(internal_coordi)
#manifold = iso.transform(internal_coordi)



iso = manifold.Isomap(n_neighbors=4, n_components=3)

coordi=coordi[0:1000*3].reshape(1000,9)
iso.fit(coordi)
manifold = iso.transform(coordi)


