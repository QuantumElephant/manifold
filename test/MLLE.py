import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from sklearn import manifold, datasets
from scipy.spatial.distance import cdist
from Function_setting import Parameters, Manifold_convergence_test

parame = Parameters(10**3, 3, "MLLE", "Internal")
parame.set_output_path('/home/sky8/structure/')
co=parame.getdata('coordi_3_internal.txt')
parame.set_testParameter(component=2,testingnumber=10**3,point_interval=10,startneighbor=6,endneighbor=8)
Manifold_convergence_test(manifold.LocallyLinearEmbedding, co, parame, method = 'modified')
