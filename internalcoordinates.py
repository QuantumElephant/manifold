# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 06:50:41 2017

@author: SKY
"""


from saddle.internal import Internal
import numpy as np
import math

def internal_generator(atomic_number, molecule_number, coordi, order):
    m = atomic_number
    if m != 3:
        out=np.empty((0,3*m))
    else:
        out=np.empty((0,2*m))
    for i in range(molecule_number):                                        #for 4 molecules
        array=coordi[i*m:(i+1)*m]
        o=order[i]
        inter=Internal(array,o,0,1)  #construct a internal coordinate object(Derric's function)
        for j in range(m):
            inter.add_bond(o[j],o[(j+1)% m])
        for j in range(m):
            inter.add_angle_cos(o[j],o[(j+1)% m],o[(j+2)% m])
        if m != 3:
            for j in range(m):
                inter.add_dihedral(o[j],o[(j+1)% m],o[(j+2)% m],o[(j+3)% m])
        result= inter.ic_values.reshape(1,len(inter.ic_values))
        out=np.append(out,result,axis=0)
    return out

coordi=np.loadtxt('coordi_3_10^5.txt')
order=np.loadtxt('coordi_3_order.txt',dtype=int)
out = internal_generator(3,10**3, coordi, order)
np.savetxt('coordi_3_internal_test.txt',out)

coordi=np.loadtxt('coordi_4_10^5.txt')
order=np.loadtxt('coordi_4_order.txt',dtype=int)
out = internal_generator(4,10**3, coordi, order)
np.savetxt('coordi_4_internal_test.txt',out)
