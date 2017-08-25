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
import math
import numpy as np
import sobol_sequence as sobol_seq
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

"""def distance(v1,v2):                                             #not used here
    d=math.sqrt((v1-v2)*((v1-v2).T))
    return d
"""
def euclideandistances(A, B):                                    #using matrix to compute distance between points 
    BT = B.transpose()
    vecProd = A * BT
    SqA =  A.getA()**2                                          # change it to adarray and square it
    sumSqA = np.sum(np.mat(SqA), axis=1)                        #construct mode^2 of every row of matrix
    sumSqAEx = np.tile(sumSqA, (1, vecProd.shape[1]))    
    SqB = B.getA()**2
    sumSqB = np.sum(SqB, axis=1)
    sumSqBEx = np.tile(sumSqB, (vecProd.shape[0], 1))    
    SqED = sumSqBEx + sumSqAEx - 2*vecProd
    SqED=abs(SqED)
    ED = (SqED.getA())**0.5
    return np.matrix(ED)
                        
total3=sobol_seq.i4_sobol_generate_std_normal(3, 3*10**2)             #for 3 molecues
out=np.empty((0,3))
for i in range(10**2):                                            
    array=total3[i*3:(i+1)*3]
    
    F=euclideandistances(np.mat(array), np.mat(array))
    #print(F)
    
    Farray=F.getA()
    Farray.sort()
    Farray=Farray[:,1:]
    R_smallest=Farray.min()
    atomicpositions=array/R_smallest
    #print(R_smallest)
    #print(distance(np.mat(array[0]), np.mat(array[3])))
    out=np.append(out,atomicpositions,axis=0)
    """fig = plt.figure()  
    ax = fig.add_subplot(111, projection='3d')  
    X = atomicpositions[:,0] 
    Y = atomicpositions[:,1]  
    Z = atomicpositions[:,2]  
    ax.scatter(X, Y, Z)  
    plt.show()
    fig.savefig('G:/Anaconda/code/output/scatter3D_3_%i.png' % i, dpi=600) """   
np.savetxt('coordi_3.txt',out)
    
'''   f = open('G:/Anaconda/code/output/coordinates_3.txt', "a+")
    f.write(str(i)+"\n"+str(atomicpositions)+"\n")
    f.close()
'''

total4=sobol_seq.i4_sobol_generate_std_normal(3, 4*10**2)
out=np.empty((0,3))
for i in range(10**2):                                       #for 4 molecues
    array=total4[i*4:(i+1)*4]
    
    F=euclideandistances(np.mat(array), np.mat(array))
    #print(F)
    
    Farray=F.getA()
    Farray.sort()
    Farray=Farray[:,1:]
    R_smallest=Farray.min()
    atomicpositions=array/R_smallest
    #print(R_smallest)
    """ax = fig.add_subplot(111, projection='3d')  
    X = atomicpositions[:,0] 
    Y = atomicpositions[:,1]  
    Z = atomicpositions[:,2]  
    ax.scatter(X, Y, Z)  
    plt.show()
    fig.savefig('G:/Anaconda/code/output/scatter3D_4_%i.png' % i, dpi=600) 
    f = open('G:/Anaconda/code/output/coordinates_4.txt', "a+")
    f.write(str(i)+"\n"+str(atomicpositions)+"\n")
    f.close()"""
    out=np.append(out,atomicpositions,axis=0)
np.savetxt('coordi_4.txt',out)
    
total5=sobol_seq.i4_sobol_generate_std_normal(3, 5*10**2)
out=np.empty((0,3))
for i in range(10**2):                                       #for 5 molecues
    array=total5[i*5:(i+1)*5]
    
    F=euclideandistances(np.mat(array), np.mat(array))
    #print(F)
    
    Farray=F.getA()
    Farray.sort()
    Farray=Farray[:,1:]
    R_smallest=Farray.min()
    atomicpositions=array/R_smallest
    #print(R_smallest)
    """
        ax = fig.add_subplot(111, projection='3d')  
        X = atomicpositions[:,0] 
        Y = atomicpositions[:,1]  
        Z = atomicpositions[:,2]  
        ax.scatter(X, Y, Z)  
        plt.show()
        fig.savefig('G:/Anaconda/code/output/scatter3D_5_%i.png' % i, dpi=600)  
    
    f = open('G:/Anaconda/code/output/coordinates_5.txt', "a+")
    f.write(str(i)+"\n"+str(atomicpositions)+"\n")
    f.close()"""
    out=np.append(out,atomicpositions,axis=0)
np.savetxt('coordi_5.txt',out)
    
total6=sobol_seq.i4_sobol_generate_std_normal(3, 6*10**2)
out=np.empty((0,3))
for i in range(10**2):                                       #for 6 molecues
    array=total6[i*6:(i+1)*6]
    
    F=euclideandistances(np.mat(array), np.mat(array))
    #print(F)
    
    Farray=F.getA()
    Farray.sort()
    Farray=Farray[:,1:]
    R_smallest=Farray.min()
    atomicpositions=array/R_smallest
    #print(R_smallest)
    #print(distance(np.mat(array[0]), np.mat(array[3])))
    """   
        fig = plt.figure()  
        ax = fig.add_subplot(111, projection='3d')  
        X = atomicpositions[:,0] 
        Y = atomicpositions[:,1]  
        Z = atomicpositions[:,2]  
        ax.scatter(X, Y, Z)  
    #    plt.show()
        fig.savefig('G:/Anaconda/code/output/scatter3D_6_%i.png' % i, dpi=600)  
        
    f = open('G:/Anaconda/code/output/coordinates_6.txt', "a+")
    f.write(str(i)+"\n"+str(atomicpositions)+"\n")
    f.close()"""
    out=np.append(out,atomicpositions,axis=0)
np.savetxt('coordi_6.txt',out)
