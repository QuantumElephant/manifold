import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy.spatial.distance import cdist
from sklearn import manifold, datasets

n=4              #power of number of samples
m=3               #number of atomics
mode=1            #mode=0 means cartician coordinates
compo=2
path = '/home/sky8/structure/'
method = 'LLE'
if mode == 0:
    mo = 'Cartisian'
else:
    mo = 'Internal'

if mode==0:
    coordi = np.loadtxt('coordi_3_10^5.txt')
    co = coordi[0:3 * 10 ** n]
    co = co.reshape(10 ** n, 9)
else:
    coordi=np.loadtxt('coordi_3_internal.txt')
    co = coordi[0: 10 ** n]
del coordi
dimen=co.shape[1]
points=5000
fig = plt.figure()

for j in range(4, 90,5):
    LLEoutput=manifold.LocallyLinearEmbedding(n_neighbors=j, n_components=compo).fit_transform(co)
    LLEoutput = LLEoutput[0:100]
    print("standrdok")
    standard=cdist(LLEoutput, LLEoutput, 'euclidean')
    divergence=[]
    F_lis_1_15 = []
    F_lis_2_25 = []
    F_lis_16_39 = []
    F_lis_51_92 = []
    F_lis_61_86 = []
    
    for i in range(490):
        cooridinate = co[0:10 ** n+i*10]
        LLEoutput = manifold.LocallyLinearEmbedding(n_neighbors=j, n_components=compo).fit_transform(cooridinate)
        LLEoutput = LLEoutput[0:100]
        F = cdist(LLEoutput, LLEoutput, 'euclidean')
        F_lis_1_15.append(F[1][15])
        F_lis_2_25.append(F[2][25])
        F_lis_16_39.append(F[16][39])
        F_lis_51_92.append(F[51][92])
        F_lis_61_86.append(F[61][86])
        di = np.sum(np.absolute(F - standard))
        divergence.append(di)
        if i == 10:
            print(10)
        if i == 200:
            print(200)
        if i == 400:
            print(400)
        if i == 800:
            print(800)
        if i == 900:
            print(900)
    
    print("finish part")
    

    plt.scatter( list(range(100,points,10)),F_lis_1_15, s=1, alpha=0.5,figure = fig)
    plt.savefig("%s%s_%s_1_15_neighbor%d_dimension%d_%d.png" % (path,mo,method,j,dimen,compo),dpi=600)
    plt.clf()

    plt.scatter( list(range(100,points,10)),F_lis_2_25, s=1, alpha=0.5,figure = fig)
    plt.savefig("%s%s_%s_2_25_neighbor%d_dimension%d_%d.png" % (path,mo,method,j,dimen,compo),dpi=600)
    plt.clf()

    plt.scatter( list(range(100,points,10)),F_lis_16_39, s=1, alpha=0.5,figure = fig)
    plt.savefig("%s%s_%s_16_39_neighbor%d_dimension%d_%d.png" % (path,mo,method,j,dimen,compo),dpi=600)
    plt.clf()

    plt.scatter( list(range(100,points,10)),F_lis_51_92, s=1, alpha=0.5,figure = fig)
    plt.savefig("%s%s_%s_51_92_neighbor%d_dimension%d_%d.png" % (path,mo,method,j,dimen,compo),dpi=600)
    plt.clf()

    plt.scatter( list(range(100,points,10)),F_lis_61_86, s=1, alpha=0.5,figure = fig)
    plt.savefig("%s%s_%s_61_86_neighbor%d_dimension%d_%d.png" % (path,mo,method,j,dimen,compo),dpi=600)
    plt.clf()

    plt.scatter( list(range(100,points,10)),divergence, s=1, alpha=0.5,figure = fig)
    plt.savefig("%s%s_%s_convergence_neighbor%d_dimension%d_%d.png" % (path,mo,method,j,dimen,compo),dpi=600)
    plt.clf()
