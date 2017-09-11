import numpy as np
from scipy.linalg import orthogonal_procrustes
from scipy.spatial import distance
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
from sklearn import manifold, datasets
import matplotlib.pyplot as plt


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




cordi = np.loadtxt('coordi_3_internal.txt')
# coordi = np.loadtxt('coordi_3_10^5.txt')


co=cordi[:,3:6]
c=np.reciprocal(np.tan(np.arccos(co)))
out=np.concatenate((cordi,c),axis=1)

np.savetxt('coordi_3_internal_cot.txt',out)

