
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 05:21:25 2017

@author: SKY
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
from scipy.spatial import distance
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
from sklearn import manifold, datasets
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist

def euclideandistances(A, B):                                    #using matrix to compute distance between points 
    BT = B.transpose()
    vecProd = A * BT
    SqA =  A.getA()**2                                          # change it to adarray and square it
    sumSqA = np.sum(np.mat(SqA), axis=1)                        #construct mode^2 of every row of matrix
    sumSqAEx = np.tile(sumSqA, (1, vecProd.shape[1]))    
    SqB = B.getA()**2
    sumSqB = np.sum(SqB, axis=1)
    sumSqBEx = np.tile(sumSqB, (vecProd.shape[0], 1))    
    SqED = sumSqBEx + sumSqAEx - 2*vecProd
    SqED=abs(SqED)
    ED = (SqED.getA())**0.5
    return np.matrix(ED)

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

# Left with 2 dimensions


n=2               #power of number of samples
m=3               #number of molecules

coordi_internal=np.loadtxt('coordi_3_internal.txt')
coordi=np.loadtxt('coordi_3_10^5.txt')
out=np.empty((0,m))

F_lis_1_15=[]
F_lis_2_25=[]
F_lis_16_39=[]
F_lis_51_92=[]
F_lis_61_86=[]
testpoint=0
co_internal = coordi_internal[0:10 **4]
iso = manifold.Isomap(n_neighbors=15, n_components=3)
iso.fit(co_internal)
mani = iso.transform(co_internal)
print("standrdok")
standard = cdist(mani, mani, 'euclidean')
standard=standard[0:100,0:100]
divergence=[]


for i in range(990):
    co_internal = coordi_internal[0:10 ** n+i*10]
    iso = manifold.Isomap(n_neighbors=15, n_components=3)
    iso.fit(co_internal)
    mani = iso.transform(co_internal)
    F=cdist(mani, mani, 'euclidean')
    F_lis_1_15.append(F[1][15])
    F_lis_2_25.append(F[2][25])
    F_lis_16_39.append(F[16][39])
    F_lis_51_92.append(F[51][92])
    F_lis_61_86.append(F[61][86])
    F=F[0:100,0:100]
    di = np.sum(np.absolute(F - standard))
    divergence.append(di)
    if i==10:
        print(10)
    if i==200:
        print(200)
    if i == 500:
        print(500)
    if i==800:
        print(800)
    if i==900:
        print(900)

d100=[x/100 for x in divergence ]
d10000=[x/10000 for x in divergence ]

print("finish part")



fig = plt.figure()
plt.scatter( list(range(100,10000,10)),F_lis_1_15, s=1, alpha=0.5,figure = fig)
plt.savefig("isomap_1_15.png",dpi=600)

fig = plt.figure()
plt.scatter( list(range(100,10000,10)),F_lis_2_25, s=1, alpha=0.5,figure = fig)
plt.savefig("isomap_2_25.png",dpi=600)

fig = plt.figure()
plt.scatter( list(range(100,10000,10)),F_lis_16_39, s=1, alpha=0.5,figure = fig)
plt.savefig("isomap_16_39.png",dpi=600)

fig = plt.figure()
plt.scatter( list(range(100,10000,10)),F_lis_51_92, s=1, alpha=0.5,figure = fig)
plt.savefig("isomap_51_92.png",dpi=600)

fig = plt.figure()
plt.scatter( list(range(100,10000,10)),F_lis_61_86, s=1, alpha=0.5,figure = fig)
plt.savefig("isomap_61_86.png",dpi=600)

fig = plt.figure()
plt.scatter( list(range(100,10000,10)),divergence, s=1, alpha=0.5,figure = fig)
plt.savefig("isomap_convergence.png",dpi=600)

fig = plt.figure()
plt.scatter( list(range(100,10000,10)),d100 , s=1, alpha=0.5,figure = fig)
plt.savefig("isomap_convergence_100.png",dpi=600)

fig = plt.figure()
plt.scatter( list(range(100,10000,10)),d10000 , s=1, alpha=0.5,figure = fig)
plt.savefig("isomap_convergence_10000.png",dpi=600)
