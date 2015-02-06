# -*- coding: utf-8 -*-
"""
Created on Thu Feb 05 23:29:40 2015

@author: nathaniel


This method of analysis has been followed from:


Euclidean Distance Matrix Analysis: A Coordinate-Free Approach for Comparing Biological Shapes Using Landmark Data
by Subhash Lele and Joan Richtmeier



"""

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt


from scipy.spatial import distance
import scipy

import pylab as P
from numpy import linalg as LA
from numpy import random as ran


np.set_printoptions(2)


import maths
import patientGen
import visualization



def EDMA():
    control, test = patientGen.patientList(0.01)
    aveCont = euclideanDistanceMatrix( coordinateWiseAverageOfShapes(control) )
    aveTest = euclideanDistanceMatrix( coordinateWiseAverageOfShapes(test) )
    D, Dbar = averageFormDifferenceMatrix(aveCont,aveTest)
    fourTtest(D, Dbar)

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
            
    Dbar = calcDbar(D)
    print "woohoo"
    
    D = np.ma.masked_equal(D,0.0)
    return D, Dbar

def calcDbar(D):
    D = np.ma.masked_equal(D,0.0)
    return D.flatten().mean()
    
    
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
    print T1, T2, T3, T4
    