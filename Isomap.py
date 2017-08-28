
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 05:21:25 2017

@author: SKY
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
from sklearn import manifold, datasets
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist



n=2               #power of number of samples
m=3               #number of molecules

# coordi_internal=np.loadtxt('coordi_3_internal.txt')
coordi=np.loadtxt('coordi_3_10^5.txt')

co = coordi[0:3*10**4]
del coordi
co=co.reshape(10**4,9)

points=5000
for j in range(4,30):
    iso = manifold.Isomap(n_neighbors=j, n_components=3)
    iso.fit(co)
    mani = iso.transform(co)
    print("standrdok")
    standard = cdist(mani, mani, 'euclidean')
    standard=standard[0:100,0:100]
    divergence=[]
    F_lis_1_15 = []
    F_lis_2_25 = []
    F_lis_16_39 = []
    F_lis_51_92 = []
    F_lis_61_86 = []
    
    for i in range(490):
        coirdinate = co[0:10 ** n+i*10]
        iso = manifold.Isomap(n_neighbors=j, n_components=3)
        iso.fit(coirdinate)
        mani = iso.transform(coirdinate)
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
    
    print("finish part")
    
    
    
    fig = plt.figure()
    plt.scatter( list(range(100,points,10)),F_lis_1_15, s=1, alpha=0.5,figure = fig)
    plt.savefig("/home/sky8/structure/Cartisian_isomap_1_15_neighbor%d.png" % j,dpi=600)
    
    fig = plt.figure()
    plt.scatter( list(range(100,points,10)),F_lis_2_25, s=1, alpha=0.5,figure = fig)
    plt.savefig("/home/sky8/structure/Cartisian_isomap_2_25_neighbor%d.png" % j,dpi=600)

    fig = plt.figure()
    plt.scatter( list(range(100,points,10)),F_lis_16_39, s=1, alpha=0.5,figure = fig)
    plt.savefig("/home/sky8/structure/Cartisian_isomap_16_39_neighbor%d.png" % j,dpi=600)

    fig = plt.figure()
    plt.scatter( list(range(100,points,10)),F_lis_51_92, s=1, alpha=0.5,figure = fig)
    plt.savefig("/home/sky8/structure/Cartisian_isomap_51_92_neighbor%d.png" % j,dpi=600)

    fig = plt.figure()
    plt.scatter( list(range(100,points,10)),F_lis_61_86, s=1, alpha=0.5,figure = fig)
    plt.savefig("/home/sky8/structure/Cartisian_isomap_61_86_neighbor%d.png" % j,dpi=600)

    fig = plt.figure()
    plt.scatter( list(range(100,points,10)),divergence, s=1, alpha=0.5,figure = fig)
    plt.savefig("/home/sky8/structure/Cartisian_isomap_convergence_neighbor%d.png" % j,dpi=600)
