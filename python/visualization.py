# -*- coding: utf-8 -*-
"""
Created on Wed Feb 04 15:54:42 2015

@author: nathaniel
"""

'''

All of the plotting functions 

'''

import pylab as P


import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt


from scipy.spatial import distance
import scipy


from numpy import linalg as LA
from numpy import random as ran
import random

np.set_printoptions(2)



def heatMap(x,y):
 
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=50)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    
    plt.clf()
    plt.imshow(heatmap, extent=extent)
    plt.show()
    
    
def histo(x, c):
    n, bins, patches = P.hist(x, 20, normed=1, histtype='bar')
    P.setp(patches, 'facecolor', c, 'alpha', 0.75)
    plt.show()

def complicatedHisto(x, dists, var, mean):
    
    '''
    this histogram will take a 2d x with each column being a different
    set of data and then plot each one with its own color.. like in one of the
    histogram matplotlib examples
    
    '''
    fig = plt.figure()
    dist = ['{:.2f}'.format(k) for k in dists]    
    colors = ('g', 'b', 'r','c', 'm','y', 'k','w')    

    n, bins, patches = P.hist(x, 20, normed=1, histtype='bar',
                            color=colors[0:len(dist)],
                            label=dist)
                            
    P.legend()
    P.show()
    
    print "woohoo"
    fig2 = plt.figure(2)
    xbins = np.linspace(-0.5,3,200)
    for sigma, mu,d,c  in zip(var, mean, dists,colors):
        
        y = P.normpdf(xbins,  mu, np.sqrt(sigma))
        l = P.plot(xbins, y, c, label=d, linewidth=1.5)
        P.legend()

    
    
  
    # n, bins, patches = P.hist(x, 20, normed=1, histtype='bar')
    # P.setp(patches, 'facecolor', c, 'alpha', 0.75)
    #plt.show()

