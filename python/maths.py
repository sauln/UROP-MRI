# -*- coding: utf-8 -*-
"""
Created on Wed Feb 04 16:01:09 2015

@author: nathaniel
"""

'''

math and statistical functions


'''


import numpy as np
from numpy import linalg as LA


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
    
def returnNormsVector(diff):
    '''takes a shape(specifically a shape representing the differences between 
    landmarks) and returns a vector of the norms of each of these differences
    '''
    norm = np.zeros((diff.shape[0],1))
    for each, n  in zip(diff, range(0, norm.shape[0])): 
        norm[n] = LA.norm(each)
        
    return norm 
    
def meanOfNorms(diff):
    '''
    probably an unnecessary function, but it returns the mean of a vector,
    specifically intended to be the mean of the vector returned from returnNormsVector
    '''
    return returnNormsVector(diff).mean()

def doPCA(scatt):
    '''
    attempted to do some pca analysis on the data.  
    '''
    from sklearn.decomposition import PCA
    pca = PCA(n_components=3)
    pca.fit(scatt)
    print(pca.explained_variance_ratio_) 
    print type(pca)
    print pca.components_[0]
    plt.plot([0,0,0], pca.components_[0])
    plt.plot([0,0,0], pca.components_[1])
    
    
    