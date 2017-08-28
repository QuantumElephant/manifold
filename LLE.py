import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy.spatial.distance import cdist
from sklearn import manifold, datasets


n=2               #power of number of samples
m=3               #number of molecules

coordi_internal=np.loadtxt('coordi_3_internal.txt')
# coordi=np.loadtxt('coordi_3_10^5.txt')


co_internal = coordi_internal[0:10 **4]
points=5000

for j in range(4,30):
    LLEoutput=manifold.LocallyLinearEmbedding(n_neighbors=j, n_components=3).fit_transform(co_internal)
    print("standrdok")
    standard=cdist(LLEoutput, LLEoutput, 'euclidean')
    standard=standard[0:100,0:100]
    divergence=[]
    F_lis_1_15 = []
    F_lis_2_25 = []
    F_lis_16_39 = []
    F_lis_51_92 = []
    F_lis_61_86 = []
    for i in range(490):
        co_internal = coordi_internal[0:10 ** n+i*10]
        LLEoutput = manifold.LocallyLinearEmbedding(n_neighbors=j, n_components=3).fit_transform(co_internal)
        F = cdist(LLEoutput, LLEoutput, 'euclidean')
        F_lis_1_15.append(F[1][15])
        F_lis_2_25.append(F[2][25])
        F_lis_16_39.append(F[16][39])
        F_lis_51_92.append(F[51][92])
        F_lis_61_86.append(F[61][86])
        F=F[0:100,0:100]
        di = np.sum(np.absolute(F - standard))
        divergence.append(di)
        if i == 10:
            print(10)
        if i == 200:
            print(200)
        if i == 500:
            print(500)
        if i == 800:
            print(800)
        if i == 900:
            print(900)

    print("finish part")

    fig = plt.figure()
    plt.scatter( list(range(100,points,10)),F_lis_1_15, s=1, alpha=0.5,figure=fig)
    plt.savefig("/home/sky8/structure/Internal_LLE_1_15_neighbor%d.png" % j,dpi=600)

    fig = plt.figure()
    plt.scatter( list(range(100,points,10)),F_lis_2_25, s=1, alpha=0.5,figure=fig)
    plt.savefig("/home/sky8/structure/Internal_LLE_2_25_neighbor%d.png" % j,dpi=600)

    fig = plt.figure()
    plt.scatter( list(range(100,points,10)),F_lis_16_39, s=1, alpha=0.5,figure=fig)
    plt.savefig("/home/sky8/structure/Internal_LLE_16_39_neighbor%d.png" % j,dpi=600)

    fig = plt.figure()
    plt.scatter( list(range(100,points,10)),F_lis_51_92, s=1, alpha=0.5,figure=fig)
    plt.savefig("/home/sky8/structure/Internal_LLE_51_92_neighbor%d.png" % j,dpi=600)

    fig = plt.figure()
    plt.scatter( list(range(100,points,10)),F_lis_61_86, s=1, alpha=0.5,figure=fig)
    plt.savefig("/home/sky8/structure/Internal_LLE_61_86_neighbor%d.png" % j,dpi=600)

    fig = plt.figure()
    plt.scatter( list(range(100,points,10)),divergence, s=1, alpha=0.5,figure=fig)
    plt.savefig("/home/sky8/structure/Internal_LLE_convergence_neighbor%d.png" % j, dpi=600)
