
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 05:21:25 2017

@author: SKY
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from sklearn import manifold, datasets
from scipy.spatial.distance import cdist
from Function_setting import Parameters, Manifold_convergence_test, Manifold_residual_variance

# parame = Parameters(10**3, 3, "Isomap", "Internal")
# parame.set_output_path('/home/sky8/structure/')
# co=parame.getdata('coordi_3_internal.txt')
# parame.set_neighbor(component=2,testingnumber=10**3,point_interval=10,startneighbor=6,endneighbor=8)
# Manifold_convergence_test(manifold.Isomap, co, parame)


par = Parameters(5*10**4, 3, "Isomap", "Internal")
par.set_output_path('/home/sky8/structure/')
co=par.getdata('coordi_3_internal.txt')
# co = co [:,3:6]
par.set_component(neighbor=17, smallest=1, biggest=6)

Manifold_residual_variance(manifold.Isomap,co, par)
