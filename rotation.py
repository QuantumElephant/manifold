# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 10:44:03 2017

@author: SKY
"""

import numpy as np
from scipy.linalg import orthogonal_procrustes
from scipy.spatial import distance
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
from sklearn import manifold, datasets

def distance_mol(molecule1,molecule2):
    a=np.square(molecule1-molecule2)
    dis=math.sqrt(np.sum(a))
    return dis
    

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

#    if norm1 == 0 or norm2 == 0:
#        raise ValueError("Input matrices must contain >1 unique points")

    # change scaling of data (in rows) such that trace(mtx*mtx') = 1
    #mtx1 /= norm1
    #mtx2 /= norm2

    # transform mtx2 to minimize disparity
    R, s = orthogonal_procrustes(mtx1, mtx2)    #R=The matrix solution of the orthogonal Procrustes problem.
    mtx2 = np.dot(mtx2, R.T)                #s=Sum of the singular values of ``dot(A.T, B)``
    # measure the dissimilarity between the two datasets
    distance_mol = math.sqrt(np.sum(np.square(mtx1 - mtx2)))

    return mtx1, mtx2, distance_mol

# Left with 2 dimensions


n=int(2)               #power of number of samples
m=int(3)               #number of molecules

coordi_internal=np.loadtxt('coordi_3_internal.txt')
coordi=np.loadtxt('coordi_3_10^5.txt')
out=np.empty((0,m))

coordi_internal=coordi[0:10**n]





geodesic_dis_matrix=np.empty([10**n,10**n])

for i in range(10**n):
    for j in range(10**n):
        if i<j:
#            print(i)
            array0=coordi[i*m:(i+1)*m]
            array1=coordi[j*m:(j+1)*m]
#            print(array0)
#            print(array1)
            mtx1, mtx2, dis = procrustes(array0, array1)
#            print(distance_mol)
            geodesic_dis_matrix[i][j]= dis
            
np.fill_diagonal(geodesic_dis_matrix, 0)
geodesic_dis_matrix = geodesic_dis_matrix + geodesic_dis_matrix.T
print(geodesic_dis_matrix)  
 

iso = manifold.Isomap(n_neighbors=7, n_components=3)
dis_matrix_mani=np.empty([10**n,10**n])                                      
co_internal=coordi_internal[0:10**n]
iso.fit(co_internal)
mani = iso.transform(co_internal)
dis_matrix_mani=iso.dist_matrix_
minus=geodesic_dis_matrix-dis_matrix_mani

#                                      for 3 molecules



#for i in range(10**n):
#    for j in range(10**n):
#        if i<j:
##            print(i)
#            array0=manifold[i:(i+1)]
#            array1=manifold[j:(j+1)]
#            dis=distance_mol(array0,array1)
##            print(distance_mol)
#            dis_matrix_mani[i][j]= dis
#            
#np.fill_diagonal(dis_matrix_mani, 0)
#dis_matrix_mani = dis_matrix_mani + dis_matrix_mani.T
#print(dis_matrix_mani)       
           
                                      
minus=geodesic_dis_matrix-dis_matrix_mani    
    
              #mtx2 is rotated

#m=4 
#coordi=np.loadtxt('coordi_4_10^5.txt')
#array0=coordi[0:m]
#fig = plt.figure()  
#ax = fig.add_subplot(111, projection='3d')  
#X = array0[:,0] 
#Y = array0[:,1]  
#Z = array0[:,2]  
#ax.scatter(X, Y, Z)
#fig.savefig('G:/Anaconda/code/output/scatter3D_4_%standard.png', dpi=600)
#for i in range(10**n):                                        #for 3 molecules
#    array=coordi[(i+1)*m:(i+2)*m]
#    
#    mtx1, mtx2, disparity = procrustes(array0, array)          #mtx2 is rotated
#    
#    print(mtx2)
#    fig = plt.figure()  
#    ax = fig.add_subplot(111, projection='3d')  
#    X = mtx2[:,0] 
#    Y = mtx2[:,1]  
#    Z = mtx2[:,2]  
#    ax.scatter(X, Y, Z)
#    fig.savefig('G:/Anaconda/code/output/scatter3D_4_%i.png' % i, dpi=600)
#
# 
#m=5
#coordi=np.loadtxt('coordi_5_10^5.txt')
#array0=coordi[0:m]
#fig = plt.figure()  
#ax = fig.add_subplot(111, projection='3d')  
#X = array0[:,0] 
#Y = array0[:,1]  
#Z = array0[:,2]  
#ax.scatter(X, Y, Z)
#fig.savefig('G:/Anaconda/code/output/scatter3D_5_%standard.png', dpi=600)
#for i in range(10**n):                                        #for 3 molecules
#    array=coordi[(i+1)*m:(i+2)*m]
#    
#    mtx1, mtx2, disparity = procrustes(array0, array)          #mtx2 is rotated
#    
#    print(mtx2)
#    fig = plt.figure()  
#    ax = fig.add_subplot(111, projection='3d')  
#    X = mtx2[:,0] 
#    Y = mtx2[:,1]  
#    Z = mtx2[:,2]  
#    ax.scatter(X, Y, Z)
#    fig.savefig('G:/Anaconda/code/output/scatter3D_5_%i.png' % i, dpi=600)
#
#m=6   
#coordi=np.loadtxt('coordi_6_10^5.txt')
#array0=coordi[0:m]
#fig = plt.figure()  
#ax = fig.add_subplot(111, projection='3d')  
#X = array0[:,0] 
#Y = array0[:,1]  
#Z = array0[:,2]  
#ax.scatter(X, Y, Z)
#fig.savefig('G:/Anaconda/code/output/scatter3D_6_%standard.png', dpi=600)
#for i in range(10**n):                                        #for 3 molecules
#    array=coordi[(i+1)*m:(i+2)*m]
#    
#    mtx1, mtx2, disparity = procrustes(array0, array)          #mtx2 is rotated
#    
#    print(mtx2)
#    fig = plt.figure()  
#    ax = fig.add_subplot(111, projection='3d')  
#    X = mtx2[:,0] 
#    Y = mtx2[:,1]  
#    Z = mtx2[:,2]  
#    ax.scatter(X, Y, Z)
#    fig.savefig('G:/Anaconda/code/output/scatter3D_6_%i.png' % i, dpi=600)    