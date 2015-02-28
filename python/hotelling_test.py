# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 21:27:19 2015

@author: nathaniel


hotelling analysis and test functions

"""
import numpy as np
import scipy as sp
from visualization import histo
from sanity_check import almost_equal


def analyze_varying_sigma():
    '''
    plots a histogram of T2 values for many iterations of varying sigma
    '''
    mc = 100
    sigma = np.linspace(eps, 1, mc)
    mu_known = 0


    ts = np.zeros(mc)
    
    #change keep track of T2 over varying sigma
    for i, s in zip(range(mc), sigma):
        
        data = make_data(mu_known, s)
        
        ts[i] = return_T2_value(data)
    
    
    from visualization import histo
    histo(ts,'r', 'varying sigma')
    plt.plot(ts)
        
  
def analyze_lots_of_trials():
    ''''
        plots a histogram of T2 values for many iterations
    
    '''
    sigma_known = 0.2
    mu_known = 0

    #faux data
    mc = 100
    ts = np.zeros(mc)
    for i in range(mc):
        data = make_data(mu_known, sigma_known)
        
        ts[i] = return_T2_value(data)
    
    from visualization import histo
    histo(ts,'r', 'lots of trials') 
     
def inv_safe(data):
    '''I was having some strange things happening where matrices that should be
    invertible were having incorrect inverses.
    '''
    
    eps =1e-2
    d_sinv = sp.linalg.inv(data)

    #print data*d_sinv
    return almost_equal(sp.dot(data,d_sinv), np.eye(data.shape[0]))
       
       
def minor(data, i, j):
    m = np.delete(data, (i), axis=0)
    m = np.delete(m, (j), axis=1)
    return m
    
       
def inverse_3x3(data):
    '''
        WARNING: DO NOT USE
        This function does not compute inverse correctly.

    '''    
    
    
    if not data.shape == (3,3):
        print "matrix needs to be 3x3 for this formula to hold"
        return False
    
    
    '''
    This formula has been taken from MathWorld
    
    A^{-1} = 1/det(A) ( det(Mij) ) 
    '''
    d_inv = np.zeros(data.shape)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            d_inv[i,j] = sp.linalg.det(minor(data,i,j))
    d_inv = (1/sp.linalg.det(data)) * d_inv
    

    return d_inv
    
        
        
        
    
def test_inv():
    print "Test that inverse is working:"
    x = True
    #a,b,c,d,e,f,g,h,i=range(1,10)
    #a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g 
    for i in range(5):
        data = np.random.rand(3,3) * 10
        
        #print "next"
        #print np.dot(sp.linalg.invhilbert(data), data) 
        #print np.dot(sp.linalg.inv(data), data)
        #print sp.linalg.inv(data) * data
        #print np.dot(d_inv,  data)
        if not inv_safe(data) and not almost_equal(sp.linalg.det(data),0.0) :
            print "broken"
            x = False
    if x:
        print "Inverse is working"
    
    
    
 
def test_minor():
    data = np.random.rand(3,3) * 10 
    m = minor(data,1,1)
    #print data, m
    
    
    
def test_pvalues_of_F():
    #According to a great online F distribution calculator:
    # If X = F(3,12), P(X<1.75) = .79
    #                 P(X< 0.5) = 0.31
    sigFig = 2
    F = stats.f( dimension, sample_size-dimension)
    print almost_equal(f.cdf(1.75), 0.79, sigFig)
    print almost_equal(f.cdf(0.5), 0.31, sigFig)
 
  
test_minor()
test_inv() 
  
  
  
  