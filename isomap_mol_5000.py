import numpy as np
from scipy.spatial import distance
import matplotlib
matplotlib.use('Agg')
import math
from sklearn import manifold, datasets
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist


def euclideandistances(A, B):  # using matrix to compute distance between points
    BT = B.transpose()
    vecProd = A * BT
    SqA = A.getA() ** 2  # change it to adarray and square it
    sumSqA = np.sum(np.mat(SqA), axis=1)  # construct mode^2 of every row of matrix
    sumSqAEx = np.tile(sumSqA, (1, vecProd.shape[1]))
    SqB = B.getA() ** 2
    sumSqB = np.sum(SqB, axis=1)
    sumSqBEx = np.tile(sumSqB, (vecProd.shape[0], 1))
    SqED = sumSqBEx + sumSqAEx - 2 * vecProd
    SqED = abs(SqED)
    ED = (SqED.getA()) ** 0.5
    return np.matrix(ED)


def distance_mol(molecule1, molecule2):
    a = np.square(molecule1 - molecule2)
    dis = math.sqrt(np.sum(a))
    return dis

# Left with 2 dimensions


n = 2  # power of number of samples
m = 3  # number of molecules

coordi_internal = np.loadtxt('coordi_3_internal.txt')
# coordi = np.loadtxt('coordi_3_10^5.txt')
out = np.empty((0, m))

F_lis_1_15 = []
F_lis_2_25 = []
F_lis_16_39 = []
F_lis_51_92 = []
F_lis_61_86 = []
testpoint = 0
co_internal = coordi_internal[0:5000]
iso = manifold.Isomap(n_neighbors=15, n_components=3)
iso.fit(co_internal)
mani = iso.transform(co_internal)
print("standrdok")
standard = cdist(mani, mani, 'euclidean')
standard = standard[0:100, 0:100]
divergence = []

for i in range(490):
    co_internal = coordi_internal[0:10 ** n + i * 10]
    iso = manifold.Isomap(n_neighbors=15, n_components=3)
    iso.fit(co_internal)
    mani = iso.transform(co_internal)
    F = cdist(mani, mani, 'euclidean')
    F_lis_1_15.append(F[1][15])
    F_lis_2_25.append(F[2][25])
    F_lis_16_39.append(F[16][39])
    F_lis_51_92.append(F[51][92])
    F_lis_61_86.append(F[61][86])
    F = F[0:100, 0:100]
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

d100 = [x / 100 for x in divergence]
d10000 = [x / 10000 for x in divergence]

print("finish part")

fig = plt.figure()
plt.scatter(list(range(100, 5000, 10)), F_lis_1_15, s=1, alpha=0.5, figure=fig)
plt.savefig("isomap_1_15_5000.png", dpi=600)

fig = plt.figure()
plt.scatter(list(range(100, 5000, 10)), F_lis_2_25, s=1, alpha=0.5, figure=fig)
plt.savefig("isomap_2_25_5000.png", dpi=600)

fig = plt.figure()
plt.scatter(list(range(100, 5000, 10)), F_lis_16_39, s=1, alpha=0.5, figure=fig)
plt.savefig("isomap_16_39_5000.png", dpi=600)

fig = plt.figure()
plt.scatter(list(range(100, 5000, 10)), F_lis_51_92, s=1, alpha=0.5, figure=fig)
plt.savefig("isomap_51_92_5000.png", dpi=600)

fig = plt.figure()
plt.scatter(list(range(100, 5000, 10)), F_lis_61_86, s=1, alpha=0.5, figure=fig)
plt.savefig("isomap_61_86_5000.png", dpi=600)

fig = plt.figure()
plt.scatter(list(range(100, 5000, 10)), divergence, s=1, alpha=0.5, figure=fig)
plt.savefig("isomap_convergence_5000.png", dpi=600)

fig = plt.figure()
plt.scatter(list(range(100, 5000, 10)), d100, s=1, alpha=0.5, figure=fig)
plt.savefig("isomap_convergence_100_5000.png", dpi=600)

fig = plt.figure()
plt.scatter(list(range(100, 5000, 10)), d10000, s=1, alpha=0.5, figure=fig)
plt.savefig("isomap_convergence_10000_5000.png", dpi=600)