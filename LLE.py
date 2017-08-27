import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy import ndimage
from scipy.spatial.distance import cdist
from sklearn import manifold, datasets


n=2               #power of number of samples
m=3               #number of molecules

coordi_internal=np.loadtxt('coordi_3_internal.txt')
# coordi=np.loadtxt('coordi_3_10^5.txt')
out=np.empty((0,m))

F_lis_1_15=[]
F_lis_2_25=[]
F_lis_16_39=[]
F_lis_51_92=[]
F_lis_61_86=[]
testpoint=0
co_internal = coordi_internal[0:10 **4]
LLEoutput=manifold.LocallyLinearEmbedding(n_neighbors=15, n_components=3).fit_transform(co_internal)
print("standrdok")
standard=cdist(LLEoutput, LLEoutput, 'euclidean')
standard=standard[0:100,0:100]
divergence=[]


for i in range(990):
    co_internal = coordi_internal[0:10 ** n+i*10]
    LLEoutput = manifold.LocallyLinearEmbedding(n_neighbors=15, n_components=3).fit_transform(co_internal)
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

d100 = [x / 100 for x in divergence]
d10000 = [x / 10000 for x in divergence]

print("finish part")

fig = plt.figure()
plt.scatter( list(range(100,10000,10)),F_lis_1_15, s=1, alpha=0.5,figure = fig)
plt.savefig("LLE_1_15.png",dpi=600)

fig = plt.figure()
plt.scatter( list(range(100,10000,10)),F_lis_2_25, s=1, alpha=0.5,figure = fig)
plt.savefig("LLE_2_25.png",dpi=600)

fig = plt.figure()
plt.scatter( list(range(100,10000,10)),F_lis_16_39, s=1, alpha=0.5,figure = fig)
plt.savefig("LLE_16_39.png",dpi=600)

fig = plt.figure()
plt.scatter( list(range(100,10000,10)),F_lis_51_92, s=1, alpha=0.5,figure = fig)
plt.savefig("LLE_51_92.png",dpi=600)

fig = plt.figure()
plt.scatter( list(range(100,10000,10)),F_lis_61_86, s=1, alpha=0.5,figure = fig)
plt.savefig("LLE_61_86.png",dpi=600)

fig = plt.figure()
plt.scatter( list(range(100,10000,10)),divergence, s=1, alpha=0.5,figure = fig)
plt.savefig("LLE_convergence.png",dpi=600)


fig = plt.figure()
plt.scatter( list(range(100,10000,10)),d100 , s=1, alpha=0.5,figure = fig)
plt.savefig("LLE_convergence_100.png",dpi=600)

fig = plt.figure()
plt.scatter( list(range(100,10000,10)),d10000 , s=1, alpha=0.5,figure = fig)
plt.savefig("LLE_convergence_10000.png",dpi=600)