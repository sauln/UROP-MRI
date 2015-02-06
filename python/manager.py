# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
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
import EDMA


def differencesIn2Shapes():
    s1 = Shape()
    s2 = Shape()
    #s1.plot()
    
    
    A1,B1 = differences(s1)
    A2,B2 = differences(s2)
    
    c = np.zeros((B1.size,1))
    for b1, b2, x in zip(B1, B2, range(0, B1.size)):
        c[x] = L2metric(b1,b2)
        
    plt.plot(c)
    plt.show()
    
    
    


def pointwiseOverShapes():
    '''
    now we want to find the difference between individual points of each shape

    '''
    s1 = Shape(0)
    s2 = Shape(1)
    
    c = np.zeros((s1.x.shape))

    cn = np.zeros((s1.x.shape[0],1))
    for b1, b2, x in zip(s1.x, s2.x, range(0, s1.x.size)):
        c[x] = b1-b2
        cn[x] = LA.norm(b1 - b2)
        
    #print cn
    return cn



def twoPointDiff():
    print "geting twopoint"
    mu = 0
    sigma = 0.1
    
    a = ran.normal(mu, sigma, (10000,2))
    print a
    
    plt.scatter(a[:,1], a[:,0])
    heatMap(a[:,1], a[:,0])


def compareEdgeLengthsBetweenShapes(drift):
    control, test = patientList(drift)

    t = 0
    c = 0
    for shape in control:
        c += shape.afl
    for shape in test:
        t += shape.afl
            
    return (t/len(control), c/len(test))
    

#twoPointDiff()

def edgeLengthMonteCarlo():
    testD = []
    contD = []
    
    l = np.arange(-0.5, 0.5, 0.001)
    print len(l)
    for drift in l:
        t,c = compareEdgeLengthsBetweenShapes(drift)
        testD.append(t)
        contD.append(c)
            
    testD = np.array(testD)
    contD = np.array(contD)
    print "Average of control: %s \nAverage of test: %s" %(np.mean(contD), np.mean(testD))
    plt.plot(l,testD)
    plt.plot(l,contD)
    print len(contD)
    plt.show()


EDMA.EDMA()
#edgeLengthMonteCarlo()