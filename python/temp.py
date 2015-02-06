# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt


import pylab as P
from numpy import linalg as LA
from numpy import random as ran

'''

mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')

x = [1,2,3,4]
y = [4,3,2,1]
z = [3,4,3,4]
ax.plot(x, y, z, label='parametric curve')
ax.legend()

plt.show()

'''

class Shape():
    # 8 points, in an array


    def __init__(self):
                
        sigma = 0.1
        
        ran.normal(0,sigma)
        self.points = list()
        self.x = np.zeros((3,8))  
        self.x = np.array([[ran.normal(0,sigma),ran.normal(0,sigma),ran.normal(0,sigma)],\
            [ran.normal(0,sigma),ran.normal(0,sigma),ran.normal(1,sigma)],\
            [ran.normal(0,sigma),ran.normal(1,sigma),ran.normal(1,sigma)], \
            [ran.normal(0,sigma),ran.normal(1,sigma),ran.normal(0,sigma)] ,\
            [ran.normal(1,sigma),ran.normal(0,sigma),ran.normal(0,sigma)],\
            [ran.normal(1,sigma),ran.normal(0,sigma),ran.normal(1,sigma)], \
            [ran.normal(1,sigma),ran.normal(1,sigma),ran.normal(1,sigma)], \
            [ran.normal(1,sigma),ran.normal(1,sigma),ran.normal(0,sigma)] ])
        
           
    def plot(self):
        fig = plt.figure()
        ax  = fig.add_subplot(111, projection = '3d')

        for i,j,k in zip(self.x[:,0], self.x[:,1], self.x[:,2]):
            ax.plot([self.x[0,0],i],[self.x[0,1],j],[self.x[0,2],k], color = 'g')
            
        ax.scatter(self.x[:,0], self.x[:,1], self.x[:,2], color = 'g', marker = "o")
        plt.show()

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

def differencesIn2Shapes():
    s1 = Shape()
    s2 = Shape()
    #s1.plot()
    
    
    A1,B1 = differences(s1)
    A2,B2 = differences(s2)
    
    c = np.zeros((B1.size,1))
    for b1, b2, x in zip(B1, B2, range(0, B1.size)):
        c[x] = L2metric(b1,b2)
        
    #print c
        
    print c.mean()
    plt.plot(c)
    plt.show()
    
    
    
'''
now we want to find the difference between individual points of each shape

'''

def pointwiseOverShapes():
    s1 = Shape()
    s2 = Shape()
    
    c = np.zeros((s1.x.shape))
    #print c.shape
    cn = np.zeros((s1.x.shape[0],1))
    for b1, b2, x in zip(s1.x, s2.x, range(0, s1.x.size)):
        
        
        #print b1-b2
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



print "running"





twoPointDiff()


print "hello world"