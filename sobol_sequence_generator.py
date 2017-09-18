# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 08:37:52 2017

@author: SKY
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 00:57:48 2017

@author: SKY
"""
import numpy as np
import sobol_sequence as sobol_seq
from scipy.spatial.distance import cdist

def generator(atom_number, molecule_number ):
    total=sobol_seq.i4_sobol_generate_std_normal(3,atom_number*molecule_number)
    out=np.empty((0,3))

    for i in range(molecule_number):
        array=total[i*atom_number:(i+1)*atom_number]
        F = cdist(array, array, 'euclidean')
        F.sort()
        F=F[:,1:]
        R_smallest=F.min()
        coordinate=array/R_smallest
        out=np.append(out,coordinate,axis=0)
    return out
 
out=generator(atom_number=3, molecule_number=10**3)
np.savetxt('coordi_3.txt',out)

out=generator(atom_number=4, molecule_number=10**3)
np.savetxt('coordi_4.txt',out)

out=generator(atom_number=5, molecule_number=10**3)
np.savetxt('coordi_5.txt',out)

out=generator(atom_number=6, molecule_number=10**3)
np.savetxt('coordi_6.txt',out)
