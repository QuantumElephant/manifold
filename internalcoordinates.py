# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 06:50:41 2017

@author: SKY
"""


from saddle.internal import Internal
import numpy as np
import math

def distance(v1,v2):                                             #not used here
    d=math.sqrt((v1-v2)*((v1-v2).T))
    return d


n=int(5)               #power of number of samples
m=int(3)               #number of molecules

coordi=np.loadtxt('coordi_3_10^5.txt')
order=np.loadtxt('coordi_3_order.txt',dtype=int)
out=np.empty((0,2*m))
for i in range(10**n-100):                                        #for 3 molecules
    array=coordi[i*m:(i+1)*m]
    o=order[i]
    #    #    print("\n")
    #    print(o)
    inter=Internal(array,o,0,1)          #construct a internal coordinate object(Derric's function)
    for j in range(m):
        inter.add_bond(o[j],o[(j+1)% m])
        
    for j in range(m):
        inter.add_angle_cos(o[j],o[(j+1)% m],o[(j+2)% m])    
        
        #inter.add_dihedral(o[j],o[(j+1)% m],o[(j+2)% m],o[(j+3)% m])
    #    print(inter.ic)
    #    print(inter.ic_values)
    #    print('\n')
    result= inter.ic_values.reshape(1,len(inter.ic_values))
    out=np.append(out,result,axis=0)
    
#    
np.savetxt('coordi_3_internal.txt',out)
    
    
    
m=4    
coordi=np.loadtxt('coordi_4_10^5.txt')
order=np.loadtxt('coordi_4_order.txt',dtype=int)
out=np.empty((0,3*m))
for i in range(10**n-100):                                        #for 4 molecules
    array=coordi[i*m:(i+1)*m]
    o=order[i]
    #    print("\n")
#    #    print(o)
    inter=Internal(array,o,0,1)          #construct a internal coordinate object(Derric's function)
    for j in range(m):
        inter.add_bond(o[j],o[(j+1)% m])
        
    for j in range(m):
        inter.add_angle_cos(o[j],o[(j+1)% m],o[(j+2)% m])
    for j in range(m):    
        inter.add_dihedral(o[j],o[(j+1)% m],o[(j+2)% m],o[(j+3)% m])
    result= inter.ic_values.reshape(1,len(inter.ic_values))
    out=np.append(out,result,axis=0)
    #    print(inter.ic)
    #    print(inter.ic_values)
    
np.savetxt('coordi_4_internal.txt',out)   
    
m=5    
coordi=np.loadtxt('coordi_5_10^5.txt')
order=np.loadtxt('coordi_5_order.txt',dtype=int)
out=np.empty((0,3*m))
for i in range(10**n-100):                                        #for 5 molecules
    array=coordi[i*m:(i+1)*m]
    o=order[i]
#    #    print("\n")
#    #    print(o)
    inter=Internal(array,o,0,1)          #construct a internal coordinate object(Derric's function)
    for j in range(m):
        inter.add_bond(o[j],o[(j+1)% m])
        
    for j in range(m):
        inter.add_angle_cos(o[j],o[(j+1)% m],o[(j+2)% m])
    for j in range(m):    
        inter.add_dihedral(o[j],o[(j+1)% m],o[(j+2)% m],o[(j+3)% m]) #
    
    result= inter.ic_values.reshape(1,len(inter.ic_values))
    out=np.append(out,result,axis=0)
#    #    print(inter.ic)
#    #    print(inter.ic_values)
    
np.savetxt('coordi_5_internal.txt',out)        #-?
 
m=6    
coordi=np.loadtxt('coordi_6_10^5.txt')
order=np.loadtxt('coordi_6_order.txt',dtype=int)
out=np.empty((0,3*m))
for i in range(10**n-100):                                        #for 6 molecules
    array=coordi[i*m:(i+1)*m]
    o=order[i]
#    #    print("\n")
#    #    print(o)
    inter=Internal(array,o,0,1)          #construct a internal coordinate object(Derric's function)
    for j in range(m):
        inter.add_bond(o[j],o[(j+1)% m])
        
    for j in range(m):
        inter.add_angle_cos(o[j],o[(j+1)% m],o[(j+2)% m])
    for j in range(m):    
        inter.add_dihedral(o[j],o[(j+1)% m],o[(j+2)% m],o[(j+3)% m])
    
    result= inter.ic_values.reshape(1,len(inter.ic_values))
    out=np.append(out,result,axis=0)
#    #    print(inter.ic)
##   #    print(inter.ic_values)
#    result= inter.ic
np.savetxt('coordi_6_internal.txt',out)    