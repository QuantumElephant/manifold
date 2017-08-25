# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 09:09:49 2017

@author: SKY
"""

import c

vec, seed = c.i4_sobol(4, 1)
vec
n=5
n

# array([ 0.5,  0.5,  0.5,  0.5])
seed
# 2

# generate the next vector in the sequence:
vec,seed=c.i4_sobol(4, seed)


print("a")