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
    # mtx1 /= norm1
    # mtx2 /= norm2

    # transform mtx2 to minimize disparity
    R, s = orthogonal_procrustes(mtx1, mtx2)  # R=The matrix solution of the orthogonal Procrustes problem.
    mtx2 = np.dot(mtx2, R.T)  # s=Sum of the singular values of ``dot(A.T, B)``
    # measure the dissimilarity between the two datasets
    distance_mol = math.sqrt(np.sum(np.square(mtx1 - mtx2)))

    return mtx1, mtx2, distance_mol


n = 2  # power of number of samples
m = 3  # number of molecules

coordi_internal = np.loadtxt('coordi_3_internal.txt')
coordi = np.loadtxt('coordi_3_10^5.txt')