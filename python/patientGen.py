# -*- coding: utf-8 -*-
"""
Created on Wed Feb 04 15:59:40 2015

@author: nathaniel
"""

'''

here I want to generate data

'''

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt


import pylab as P
from numpy import linalg as LA
from numpy import random as ran

def patientList(drift):
    control = []
    test = []
    #print "generate patient list"
    for i in range(0,20):
        control.append(Shape(0))
        test.append(Shape(drift))
        
    return (control, test)
    
    

class Shape():
    # 8 points, in an array
    def __init__(self, drift='a'):
                
        sigma = 0.1
        if drift == 'a':
            self.points = list()
            self.x = np.zeros((3,8))  
            self.x = np.array([[0,0,0],[0,0,1], [0,1,1],[0,1,0] ,\
                [1,0,0], [1,0,1],[1,1,1],[1,1,0] ])
            
            
        elif drift  == 0:
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
        else:
            self.points = list()
            self.x = np.zeros((3,8))  
            self.x = np.array([[ran.normal(0,sigma),ran.normal(0,sigma),ran.normal(0,sigma)],\
                [ran.normal(0,sigma),ran.normal(0,sigma),ran.normal(1,sigma)],\
                [ran.normal(0,sigma),ran.normal(1,sigma),ran.normal(1,sigma)], \
                [ran.normal(0,sigma),ran.normal(1,sigma),ran.normal(0,sigma)] ,\
                [ran.normal(1,sigma),ran.normal(0,sigma),ran.normal(0,sigma)],\
                [ran.normal(1,sigma),ran.normal(0,sigma),ran.normal(1,sigma)], \
                [ran.normal(1,sigma),ran.normal(1,sigma),ran.normal(1+drift,sigma)], \
                [ran.normal(1,sigma),ran.normal(1,sigma),ran.normal(0,sigma)] ])
            
        self.edges = [ [0,1], [0,3], [0,4], [1,2], [1,5], [2,3], [2,6], \
            [3,7], [4,5], [4,7], [5,6], [6,7]]
        self.fl = self.frameLength() 
        self.afl = self.fl / len(self.edges)
        
        
        
    def frameLength(self):
        s = 0
        for each in self.edges:
            s +=  LA.norm(self.x[each[0]] - self.x[each[1]])

        return s
           
           
           
    def plot(self):
        fig = plt.figure()
        ax  = fig.add_subplot(111, projection = '3d')

        for i,j,k in zip(self.x[:,0], self.x[:,1], self.x[:,2]):
            ax.plot([self.x[0,0],i],[self.x[0,1],j],[self.x[0,2],k], color = 'g')
            
        ax.scatter(self.x[:,0], self.x[:,1], self.x[:,2], color = 'g', marker = "o")
        plt.show()
        
        
        
        
        
        
        
        
        
