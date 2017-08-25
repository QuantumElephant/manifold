import numpy as np
from scipy.linalg import orthogonal_procrustes
from scipy.spatial import distance
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
from sklearn import manifold, datasets
import matplotlib.pyplot as plt

def(standard,sample):
    difference=sample-standard

