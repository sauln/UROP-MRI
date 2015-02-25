# -*- coding: utf-8 -*-
"""
Created on Wed Feb 04 15:59:40 2015

@author: nathaniel
"""

'''

here I want to generate data

'''


import numpy as np
import matplotlib.pyplot as plt

from numpy import linalg as LA
from numpy import random as ran
import manager

import objectReader


## we need to set up the data
f = r"..\image\original\rawlandmarks.txt"





def newShapeFromFile(f):
    f = open(f, "r")
    sh = Shape(kwarg="read")
    ind = 0
    
    for line in f:
        l = [float(x) for x in line.split(' ')]
        sh.x[ind] = l[1:]
        ind+=1
    return sh
    


def patientList(drift, kwarg=None):
    '''
    This creates two sample sets - 
    one set of shapes that has no drift, just noise samples,
    and one set of shapes that has some drift
    
    we want to adjust this so it is easier to adapt
    ''' 
    size = 50
    #if kwarg == "just_matrix":
    control = genSet(0, size)
    test    = genSet(drift, size, kwarg)
    
    return (control, test)
    
def genSet(drift=None, size = 20, kwarg = None):
    pset = []
    
    if kwarg =="justM":
        for i in range(0,size):
            pset.append(justTheMatrix(drift))
    else:
        for i in range(0,size):
            pset.append(Shape(drift))
        
    return pset

def justTheMatrix(drift=None, kwarg=None):
    return Shape(drift).x


class Shape():
    # 8 points, in an array
    def __init__(self, drift=None,kwarg=None):
        self.sigma = 0.1
    
        if kwarg is not "read":
            self.makeEdges() 
            self.standardSet() 
            self.adjustSet(drift)
        else:
            self.x = np.zeros((8,3))
        
        #self.fillShapeFaux(drift, kwarg)
        
        
        #self.fl = self.frameLength()  #frame length and average frame lenght
        #self.afl = self.fl / len(self.edges) #not used any more
        
        
        
    def frameLength(self):
        s = 0
        for each in self.edges:
            s +=  LA.norm(self.x[each[0]] - self.x[each[1]])

        return s
           

           
    def makeEdges(self):
        
        self.adjMatrix = np.array([[0,1,0,1,0,0,0,1],
                                   [0,0,1,0,0,0,1,0],
                                   [0,0,0,1,0,1,0,0],
                                   [0,0,0,0,1,0,0,0],
                                   [0,0,0,0,0,1,0,1],
                                   [0,0,0,0,0,0,1,0],
                                   [0,0,0,0,0,0,0,1],
                                   [0,0,0,0,0,0,0,0]] )
        edges =  np.where(self.adjMatrix ==1)
        edgeSet = set()
        #print self.adjMatrix
        for a,b in zip(edges[0], edges[1]):
            if a<=b:
                edgeSet.add((a,b))
            else:
                edgeSet.add((b,a))
                
        self.edges = edgeSet

  
    def standardSet(self):
        #create all standard set of points
        global sh
    
        self.x = sh.x
        self.x = self.x.astype(float)
        
    def adjustSet(self, drift):
        #adjust the set of standard points randomly and with drift
        
        if not hasattr(self, 'x'):#this should have some sort of inforce, or try/exept
            self.standardSet()
   
        #this adds random noise to the shape
        self.x = self.x + np.random.normal(0, self.sigma, self.x.shape)
        
        if drift is not None:
            #create a 4x3 matrix random matrix centered at drift
            #add this drift matrix to the matrix
            d = np.random.normal(drift, self.sigma, (self.x.shape[0] - 4, self.x.shape[1]))
            os = np.zeros(self.x.shape)
            os[1:5,:] = d
            self.x += os
           
    def plotMe(self):
        fig = plt.figure()
        ax  = fig.add_subplot(111, projection = '3d',aspect='equal')
        print self.x
        #m=  max(self.x.flatten)
        m = 1
        print type(self.x)
        m = np.max(self.x.flatten())
        print m
        
        ax.set_xlim([-m, m])
        ax.set_ylim([-m, m])
        ax.set_zlim([-m, m])

        for each in self.edges:
        #self.x describes the vertices, self.adjMatrix describes the edges
            ax.plot([self.x[each[0],0],self.x[each[1],0]],
                    [self.x[each[0],1],self.x[each[1],1]],
                    [self.x[each[0],2],self.x[each[1],2]])
            
        ax.scatter(self.x[:,0], self.x[:,1], self.x[:,2], color = 'g', marker = "o")
        plt.show()
            
    def printMe(self):
        for x in range(0, self.x.shape[0]):
            print self.x[x]       


sh = newShapeFromFile(f)

                     
 
if __name__ == "__main__":
    manager.main()      
        
       
    '''
    def fillShapeFaux(self,drift, kwarg):
        
        
        
        #what if instead of all these different cases, we generate random
        #noise by a random permutation to the set of data
        sigma = 0.1
        if drift == 'a':
            self.standardSet()

            
            
        elif drift  == 0:
            self.points = list()
            self.x = np.zeros((3,8))  
            self.x = np.array([
                [ran.normal(0,sigma),ran.normal(0,sigma),ran.normal(0,sigma)],\
                [ran.normal(0,sigma),ran.normal(0,sigma),ran.normal(1,sigma)],\
                [ran.normal(0,sigma),ran.normal(1,sigma),ran.normal(1,sigma)], \
                [ran.normal(0,sigma),ran.normal(1,sigma),ran.normal(0,sigma)] ,\
                [ran.normal(1,sigma),ran.normal(0,sigma),ran.normal(0,sigma)],\
                [ran.normal(1,sigma),ran.normal(0,sigma),ran.normal(1,sigma)], \
                [ran.normal(1,sigma),ran.normal(1,sigma),ran.normal(1,sigma)], \
                [ran.normal(1,sigma),ran.normal(1,sigma),ran.normal(0,sigma)] ])
            self.x = self.x.astype(float)
        elif kwarg == "minus": #this is just some random key that could be changed
                             #to something that makes more sense in the future
            self.points = list()
            self.x = np.zeros((3,8))  
            self.x = np.array([[ran.normal(0,sigma),ran.normal(0,sigma),ran.normal(0,sigma)],\
                [ran.normal(0,sigma),ran.normal(0,sigma),ran.normal(1-drift,sigma)],\
                [ran.normal(0,sigma),ran.normal(1-drift,sigma),ran.normal(1,sigma)], \
                [ran.normal(0,sigma),ran.normal(1-drift,sigma),ran.normal(0,sigma)] ,\
                [ran.normal(1,sigma),ran.normal(0,sigma),ran.normal(0,sigma)],\
                [ran.normal(1,sigma),ran.normal(0,sigma),ran.normal(1,sigma)], \
                [ran.normal(1,sigma),ran.normal(1,sigma),ran.normal(1-drift,sigma)], \
                [ran.normal(1,sigma),ran.normal(1,sigma),ran.normal(0,sigma)] ])
            self.x = self.x.astype(float)  
        elif kwarg == "allShift": #this is just some random key that could be changed
                             #to something that makes more sense in the future
            self.points = list()
            self.x = np.zeros((3,8))  
            self.x = np.array([[ran.normal(0,sigma),ran.normal(0,sigma),ran.normal(0+drift,sigma)],\
                [ran.normal(0,sigma),ran.normal(0,sigma),ran.normal(1+drift,sigma)],\
                [ran.normal(0,sigma),ran.normal(1,sigma),ran.normal(1+drift,sigma)], \
                [ran.normal(0,sigma),ran.normal(1,sigma),ran.normal(0+drift,sigma)] ,\
                [ran.normal(1,sigma),ran.normal(0,sigma),ran.normal(0+drift,sigma)],\
                [ran.normal(1,sigma),ran.normal(0,sigma),ran.normal(1+drift,sigma)], \
                [ran.normal(1,sigma),ran.normal(1,sigma),ran.normal(1+drift,sigma)], \
                [ran.normal(1,sigma),ran.normal(1,sigma),ran.normal(0+drift,sigma)] ])
            self.x = self.x.astype(float)
        
        else:
            self.points = list()
            self.x = np.zeros((3,8))  
            self.x = np.array([[ran.normal(0,sigma),ran.normal(0,sigma),ran.normal(0,sigma)],\
                [ran.normal(0,sigma),ran.normal(0,sigma),ran.normal(1+drift,sigma)],\
                [ran.normal(0,sigma),ran.normal(1+drift,sigma),ran.normal(1,sigma)], \
                [ran.normal(0,sigma),ran.normal(1+drift,sigma),ran.normal(0,sigma)] ,\
                [ran.normal(1,sigma),ran.normal(0,sigma),ran.normal(0,sigma)],\
                [ran.normal(1,sigma),ran.normal(0,sigma),ran.normal(1,sigma)], \
                [ran.normal(1,sigma),ran.normal(1,sigma),ran.normal(1+drift,sigma)], \
                [ran.normal(1,sigma),ran.normal(1,sigma),ran.normal(0,sigma)] ])
            self.x = self.x.astype(float)
                
                       
     '''   
            
