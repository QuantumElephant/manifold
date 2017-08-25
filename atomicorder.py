# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 05:39:11 2017

@author: SKY
"""

import itertools
from itertools import cycle, islice, dropwhile
from itertools import islice
import numpy as np
from scipy.spatial import distance
import operator
import math

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


n=int(0)               #power of number of samples
m=int(3)               #number of molecules

coordi=np.loadtxt('coordi_3_10^5.txt')
out=np.empty((0,m))
for i in range(10**n):                                        #for 3 molecules
    array=coordi[i*m:(i+1)*m]
    dist01 = distance.euclidean(array[0],array[1])
    dist12 = distance.euclidean(array[1],array[2])
    dist20 = distance.euclidean(array[2],array[0])
    #print(i)
    #print(str(dist01)+" "+str(dist12)+" "+str(dist20))
    D={"01":dist01, "12":dist12, "20":dist20 }
    # x = {a: 2, b: 4, c: 3, d: 1, e: 0}
    #sorted_x = sorted(x.items(), key=operator.itemgetter(1))
    m1=min(D, key=D.get)         #smallest
    delete=D.pop(m1)
    m2=min(D, key=D.get)
    l1=list(m1)
    l2=list(m2)
    second=list(set(l1).intersection(l2))   #second order
    first=list(set(l1)-set(second))         #first order
    third=list(set(l2)-set(second))          #third order
    order=np.array([first[0],second[0],third[0]])
    order=order.reshape(1,m).astype(int)
    """
    print(distance.euclidean(array[int(first[0])],array[int(second[0])]))
    print(distance.euclidean(array[int(second[0])],array[int(third[0])]))
    print(distance.euclidean(array[int(third[0])],array[int(first[0])]))
    """
    out=np.append(out,order,axis=0)
np.savetxt('coordi_3_order.txt',out,fmt='%d')
    
m=4    
coordi=np.loadtxt('coordi_4_10^5.txt')         #for 4 molecules
list1 = '123'
iter =itertools.permutations(list1,(m-1))
permutations=list(iter)
per=np.array(permutations)
rows, collums= per.shape
per=np.concatenate((np.zeros(shape=(rows,1),dtype=int),per,np.zeros(shape=(rows,1),dtype=int)),axis=1)
per=per.astype(int)
out=np.empty((0,m))
for i in range(10**n):
    array=coordi[i*m:(i+1)*m]
    dis_arr=np.empty([m,m])
    path=[]
    shor_step=[]
    step_2=np.empty([rows,collums+1])
    for j in range(m):
        for k in range(m):
            dis_arr[j][k]=distance.euclidean(array[j],array[k])
            
    for j in range(rows):
        ring=0
        step=[]
        for k in range(collums+1):
            
            ring+=dis_arr[per[j][k]][per[j][k+1]]
            step.append(dis_arr[per[j][k]][per[j][k+1]])
            step_2[j][k]=dis_arr[per[j][k]][per[j][k+1]]
        path.append(ring)
        step_index, min_value = min(enumerate(step), key=(lambda x: x[1]))
        shor_step.append(step_index)
    #print(i)     
    #min_index, min_value = min(enumerate(shortest), key=(lambda x: x[1]))
    shortest=min(path)
    #print(shortest)
    index=[i for i, j in enumerate(path) if isclose(j, shortest)]
    chosen=shor_step[index[0]]                                      #first step index
    line=index[0]                                                   #the index of chosen path
    can_fir1=[per[line][chosen],per[line][chosen+1]]                #breaking point  candidate first point
    trace=per[line]
    #print(can_fir1)
    #print(step_2[line])
    rows2, collums2= per.shape
    per_=per[line,0:collums2-1]
    #print(per_)
    left=dis_arr[per_[(chosen-1) % per_.shape[0]]][per_[chosen]]     #two directions -compare the second step for choosing the smaller one
    right=dis_arr[per_[(chosen+1)% per_.shape[0]]][per_[(chosen+2) %per_.shape[0]]]    
    result_order=[]
    if left < right:
        ls_re_per=list(per_)
        ls_re_per.reverse()
        cycled =  cycle(ls_re_per)
        skipped = dropwhile(lambda x: x != can_fir1[1], cycled)
        sliced = islice(skipped, None, m)
        result_order = list(sliced)  # create a list from iterator
        #print(result_order)
    else:
        ls_per=list(per_)
        cycled =  cycle(per_)
        skipped = dropwhile(lambda x: x != can_fir1[0], cycled)
        sliced = islice(skipped, None, m)
        result_order = list(sliced)  # create a list from iterator
        #print(result_order)        
    order=np.array(result_order)
    order=order.reshape(1,m).astype(int)   
    out=np.append(out,order,axis=0)
np.savetxt('coordi_4_order.txt',out,fmt='%d')    
        
m=5    
coordi=np.loadtxt('coordi_5_10^5.txt')         #for 5 molecules
list1 = '1234'
iter =itertools.permutations(list1,(m-1))
permutations=list(iter)
per=np.array(permutations)
rows, collums= per.shape
per=np.concatenate((np.zeros(shape=(rows,1),dtype=int),per,np.zeros(shape=(rows,1),dtype=int)),axis=1)
per=per.astype(int)
out=np.empty((0,m))
for i in range(10**n):
    array=coordi[i*m:(i+1)*m]
    dis_arr=np.empty([m,m])
    path=[]
    shor_step=[]
    step_2=np.empty([rows,collums+1])
    for j in range(m):
        for k in range(m):
            dis_arr[j][k]=distance.euclidean(array[j],array[k])
            
    for j in range(rows):
        ring=0
        step=[]
        for k in range(collums+1):
            ring+=dis_arr[per[j][k]][per[j][k+1]]
            step.append(dis_arr[per[j][k]][per[j][k+1]])
            step_2[j][k]=dis_arr[per[j][k]][per[j][k+1]]
        path.append(ring)
        step_index, min_value = min(enumerate(step), key=(lambda x: x[1]))
        shor_step.append(step_index)
    #print(i)     
    #min_index, min_value = min(enumerate(shortest), key=(lambda x: x[1]))
    shortest=min(path)
    #print(shortest)
    index=[i for i, j in enumerate(path) if isclose(j, shortest)]
    chosen=shor_step[index[0]]                                      #first step index
    line=index[0]                                                   #the index of chosen path
    can_fir1=[per[line][chosen],per[line][chosen+1]]                #breaking point  candidate first point
    trace=per[line]
    #print(can_fir1)
    #print(step_2[line])
    rows2, collums2= per.shape
    per_=per[line,0:collums2-1]
    #print(per_)
    left=dis_arr[per_[(chosen-1) % per_.shape[0]]][per_[chosen]]     #two directions -compare the second step for choosing the smaller one
    right=dis_arr[per_[(chosen+1)% per_.shape[0]]][per_[(chosen+2) %per_.shape[0]]]
    result_order=[]
    if left < right:
        ls_re_per=list(per_)
        ls_re_per.reverse()
        cycled =  cycle(ls_re_per)
        skipped = dropwhile(lambda x: x != can_fir1[1], cycled)
        sliced = islice(skipped, None, m)
        result_order = list(sliced)  # create a list from iterator
        #print(result_order)
    else:
        ls_per=list(per_)
        cycled =  cycle(per_)
        skipped = dropwhile(lambda x: x != can_fir1[0], cycled)
        sliced = islice(skipped, None, m)
        result_order = list(sliced)  # create a list from iterator
        #print(result_order)                
    order=np.array(result_order)
    order=order.reshape(1,m).astype(int)   
    out=np.append(out,order,axis=0)
np.savetxt('coordi_5_order.txt',out,fmt='%d')         
        
m=6   
coordi=np.loadtxt('coordi_6_10^5.txt')         #for 6 molecules
list1 = '12345'
iter =itertools.permutations(list1,(m-1))
permutations=list(iter)
per=np.array(permutations)
rows, collums= per.shape
per=np.concatenate((np.zeros(shape=(rows,1),dtype=int),per,np.zeros(shape=(rows,1),dtype=int)),axis=1)
per=per.astype(int)
out=np.empty((0,m))
for i in range(10**n):
    array=coordi[i*m:(i+1)*m]
    dis_arr=np.empty([m,m])
    path=[]
    shor_step=[]
    step_2=np.empty([rows,collums+1])
    for j in range(m):
        for k in range(m):
            dis_arr[j][k]=distance.euclidean(array[j],array[k])
            
    for j in range(rows):
        ring=0
        step=[]
        for k in range(collums+1):
            ring+=dis_arr[per[j][k]][per[j][k+1]]
            step.append(dis_arr[per[j][k]][per[j][k+1]])
            step_2[j][k]=dis_arr[per[j][k]][per[j][k+1]]
        path.append(ring)
        step_index, min_value = min(enumerate(step), key=(lambda x: x[1]))
        shor_step.append(step_index)
    #print(i)     
    #min_index, min_value = min(enumerate(shortest), key=(lambda x: x[1]))
    shortest=min(path)
    #print(shortest)
    index=[i for i, j in enumerate(path) if isclose(j, shortest)]
    chosen=shor_step[index[0]]                                      #first step index
    line=index[0]                                                   #the index of chosen path
    can_fir1=[per[line][chosen],per[line][chosen+1]]                #breaking point  candidate first point
    trace=per[line]
    #print(can_fir1)
    #print(step_2[line])
    rows2, collums2= per.shape
    per_=per[line,0:collums2-1]
    #print(per_)
    left=dis_arr[per_[(chosen-1) % per_.shape[0]]][per_[chosen]]     #two directions -compare the second step for choosing the smaller one
    right=dis_arr[per_[(chosen+1)% per_.shape[0]]][per_[(chosen+2) %per_.shape[0]]]    
    result_order=[]
    if left < right:
        ls_re_per=list(per_)
        ls_re_per.reverse()
        cycled =  cycle(ls_re_per)
        skipped = dropwhile(lambda x: x != can_fir1[1], cycled)
        sliced = islice(skipped, None, m)
        result_order = list(sliced)  # create a list from iterator
        #print(result_order)
    else:
        ls_per=list(per_)
        cycled =  cycle(per_)
        skipped = dropwhile(lambda x: x != can_fir1[0], cycled)
        sliced = islice(skipped, None, m)
        result_order = list(sliced)  # create a list from iterator
        #print(result_order)                 
        
    order=np.array(result_order)
    order=order.reshape(1,m).astype(int)   
    out=np.append(out,order,axis=0)
np.savetxt('coordi_6_order.txt',out,fmt='%d')     
        
        
        