# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 05:39:11 2017

@author: SKY
"""

import itertools
from itertools import cycle, islice, dropwhile
import numpy as np
from scipy.spatial import distance

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def atomic_order(atomic_number, molecule_number, coordinate):
    m=atomic_number
    lis = list(range(1,m))
    permutations=list(itertools.permutations(lis,(m-1)))
    per=np.array(permutations)
    rows, collums= per.shape
    per=np.concatenate((np.zeros(shape=(rows,1),dtype=int),per,np.zeros(shape=(rows,1),dtype=int)),axis=1)
    per=per.astype(int)
    out=np.empty((0,m))
    for i in range(molecule_number):
        array=coordinate[i*m:(i+1)*m]
        dis_arr=np.empty([m,m])
        path_length=[]
        shortest_step=[]
        # everystep=np.empty([rows,collums+1])
        for j in range(rows):
            ring=0
            step=[]
            for k in range(collums+1):
                ring += dis_arr[per[j][k]][per[j][k+1]]
                step.append(dis_arr[per[j][k]][per[j][k+1]])
            # everystep[j][k]=dis_arr[per[j][k]][per[j][k+1]]
            path_length.append(ring)
            step_index, min_value = min(enumerate(step), key=(lambda x: x[1]))
            shortest_step.append(step_index)
        shortest_length=min(path_length)
        index=[i for i, j in enumerate(path_length) if isclose(j, shortest_length)]
        line=index[0]   #the index of chosen path
        first_step=shortest_step[line]   #first step index
        can_first=[per[line][first_step],per[line][first_step+1]]  #breaking point candidate first point
        # trace=per[line]
        rows2, collums2= per.shape
        per_=per[line,0:collums2-1]
        left=dis_arr[per_[(first_step-1) % per_.shape[0]]][per_[first_step]]     #two directions -compare the second step for choosing the smaller one
        right=dis_arr[per_[(first_step+1)% per_.shape[0]]][per_[(first_step+2) % per_.shape[0]]]
        result_order=[]
        if left < right:
            ls_reverse_per=list(per_)
            ls_reverse_per.reverse()
            cycled =  cycle(ls_reverse_per)
            skipped = dropwhile(lambda x: x != can_first[1], cycled)
            sliced = islice(skipped, None, m)
            result_order = list(sliced)  # create a list from iterator
        else:
            # ls_per=list(per_)
            cycled =  cycle(per_)
            skipped = dropwhile(lambda x: x != can_first[0], cycled)
            sliced = islice(skipped, None, m)
            result_order = list(sliced)  # create a list from iterator
        order=np.array(result_order)
        order=order.reshape(1,m).astype(int)
        out=np.append(out,order,axis=0)
    return out


coordi =np.loadtxt('coordi_3_10^5.txt')         #for 4 molecules
out=atomic_order(3, 10**3, coordi)
np.savetxt('coordi_3_order.txt',out,fmt='%d')
