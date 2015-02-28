# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 21:50:50 2015

@author: nathaniel

going to write the hotelling stuffs in python


to do this test we will need a handful of important piecies.

all 





Let x_i, i =1,...,n be p=3 dimensional points that are normally distributed around mu with covariance Sigma

Let W = (1/n-1) sum_{i=1}^n (x_i - xbar)(x_i - xbar)^t
where xbar = (x_1 + ... + x_n)/n 
and (.)^t represents the transpose.

Let t^2 = n(xbar - mu)^t W^{-1}(xbar - mu)

Then we can see that (n-p)/(p(n-1)) t^2 ~ F_{p,n-p}

"""

import numpy as np
import scipy as sp
from scipy import stats
import matplotlib.pyplot as plt
from visualization import histo

from sanity_check import almost_equal

import hotelling_test

eps = 1e-5
sample_size = 25.0
dimension = 3.0




def make_data(mu, sigma):
    data = np.random.normal(mu, sigma, sample_size*dimension)
    data = np.reshape(data,(sample_size,dimension))
    return data

           
def t2_value(data):
    '''
    takes data and returns the T2 value - floating point number
    this assumes that mu=0
    '''

    xbar = np.asarray((np.mean(data[:,0]), 
                       np.mean(data[:,1]), np.mean(data[:,2])) )
    
    #where W is the covariance matrix
    w = np.cov(np.transpose(data))
    
    #where t^2 = n(xbar - mu)^t W^{-1}(xbar - mu)
    t2 =  sample_size * np.dot(xbar,  np.dot( np.linalg.inv(w) , xbar ))
    
    return t2


def t2_test():
    mc = 100
    ts = np.zeros(mc)
    for i in range(mc):
        data = make_data(0,1)
        ts[i] = t2_value(data)
    histo(ts)
    

def t_value(data):
    stdErr = np.std(data) /  np.sqrt(float(len(data)))
    t =  np.mean(data)/ stdErr
    return t
    

def p_value(t2):
    ''' 
    Calculate the p-value of the F distribution at t2    
    '''

    T2 = (sample_size-dimension)/(dimension*(sample_size-1)) * t2
    f = stats.f( dimension, sample_size-dimension)
    return  f.cdf(T2)
    

def hotelling_T2(data):
    '''
       takes multivariant data source and calculates the p value that
       the mean is zero 
    '''
    
    t = t2_value(data)
    p = p_value(t)
    return p
    



    
    
def hotelling_T2_test():
    '''
    This is not working!
    There is no rhyme or reason in the p value - 
    
    '''
    mc = 100
    p = np.zeros(mc)
    for i in range(mc):
        data = make_data(0,1)
        p[i] = hotelling_T2(data)
    histo(p)
    
    
t2_test()
hotelling_T2_test()

#analyze_varying_sigma()
#analyze_lots_of_trials()
#plot_f_distro(15,3)
#test_pvalues()




