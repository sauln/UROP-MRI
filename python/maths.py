# -*- coding: utf-8 -*-
"""
Created on Wed Feb 04 16:01:09 2015

@author: nathaniel
"""

'''

math and statistical functions


'''

def choose(n,k):
    return np.math.factorial(n)/( np.math.factorial(n-k)*np.math.factorial(k))

def differences(P):
    #this takes the shape object and returns a new matrix
    #that describes each of the pairwise differences
    #A = np.ndarray()
    n = choose(P.x.shape[0], 2)
    
    
    A = np.zeros((n, 3))
    B = np.zeros((n, 1))
    ind = 0
    for j in range(0, P.x.shape[0]-1 ):
        for i in range(j+1, P.x.shape[0]):
            A[ind] =  P.x[j]-P.x[i] 
            B[ind] =  L2metric(P.x[j], P.x[i])
            ind +=1
   
    return A,B
        
def L2metric(A,B):
    return LA.norm(A - B)
    
    
    
    