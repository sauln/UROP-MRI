# -*- coding: utf-8 -*-
"""
Created on Thu Feb 05 23:29:40 2015

@author: nathaniel


This method of analysis has been followed from:


Euclidean Distance Matrix Analysis: A Coordinate-Free Approach for Comparing Biological Shapes Using Landmark Data
by Subhash Lele and Joan Richtmeier



"""

import numpy as np
import matplotlib.pyplot as plt
import random
import scipy


import patientGen
import manager





def EDMA():
    
    drift = 0.2
    control, test = patientGen.patientList(drift)
    
    
    aveCont = euclideanDistanceMatrix( coordinateWiseAverageOfShapes(control) )
    aveTest = euclideanDistanceMatrix( coordinateWiseAverageOfShapes(test) )
    D, Dbar = averageFormDifferenceMatrix(aveCont,aveTest)
    T4 = fourTtest(D, Dbar)
    
    
    T4s = bootstrappingTheNullDistribution(control,test)
    
    print "Our T4: %s" %T4
    print "The random permuation distribution: %s" %np.mean(T4s)
    plt.clf()
    fig = plt.hist(T4s)
    
    ymax = max(np.histogram(T4s)[0])
    plt.plot((T4, T4), (0, ymax*1.15), 'k-')
    plt.title("EDMA with drift %s and sigma: 0.1"%drift)
    plt.show()
    

def coordinateWiseAverageOfShapes(listOfShapes):
    suma = np.zeros(listOfShapes[0].x.shape)
    for each in listOfShapes:
        suma += each.x
        
    return suma/len(listOfShapes)
    
def euclideanDistanceMatrix(shape):
    return scipy.spatial.distance_matrix(shape, shape)
    

   
def averageFormDifferenceMatrix(a,b):
    '''
        takes 2 euclidean distance matrices
        finds the form differene matrix
    
        D_ij (Xhat, Yhat) = F_ij(Xhat)/F_ij(Yhat)
            for i>j=1,2,...,K
    '''
    
    D = np.zeros(a.shape)
    
    for j in range(0, a.shape[0]-1):
        for i in range(j+1, a.shape[0]-1):
            D[i,j] = a[i,j] / b[i,j]
         
         
         
         
    #try to put 1 instead of the mask.
         #look to make sure the masks are only over the diagonal
    
    D = np.ma.masked_equal(D,0.0)
    Dbar = D.flatten().mean()
    return D, Dbar

    
def fourTtest(D, Dbar): 
    D = np.ma.compressed(D)
    T1 = np.sum([(a-Dbar)**2 for a in D])
        #\sum_ij[ Dij(X,Y) - Dbar]^2
    T2 = np.sum([(a-Dbar) for a in D])
        #\sum_ij[ Dij(X,Y) - Dbar]    
    T3 = np.max(D) - np.min(D)
        # max_ij Dij(X,Y) - min_ij D_ij(X,Y)    
    T4 = np.max(D)/np.min(D)
        # max Dij(X,Y)/ min Dij(X,Y)    
    #print T1, T2, T3, T4
    return T4  
    
def bootstrappingTheNullDistribution(X, Y):
    '''
    we want to compare our T4 value for these two groups
    with many T4 values for these 2 groups shuffled up in random ways
    
    '''
    
    n = X+Y
    T4s = []
    
    for x in range(1000):
        random.shuffle(n)
        
        Xstar = n[0:20]
        Ystar = n[20:40]
        #now do all the EDMA stuff with these 2 sets,
        #and track the T4 value.
        aveCont = euclideanDistanceMatrix( coordinateWiseAverageOfShapes(Xstar) )
        aveTest = euclideanDistanceMatrix( coordinateWiseAverageOfShapes(Ystar) )
        D, Dbar = averageFormDifferenceMatrix(aveCont,aveTest)
        T4s.append(fourTtest(D, Dbar))
    
    return T4s

        
    
    
if __name__ == "__main__":
    manager.main()