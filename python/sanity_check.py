# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 10:39:02 2015

@author: nathaniel




sanity checks

"""

import numpy as np

eps = 1e-5
def dotTest():
    a = np.asarray((1,3))
    b = np.asarray([[1,2],[3,4]])
    c = np.asarray([6,4])
    
    expected_answer = 116
    
    dot_answer_1 = np.dot(a,  np.dot( b ,c ))
    dot_answer_2 = np.dot(a, np.dot(b,np.transpose(c)))
    if not dot_answer_1 == expected_answer:
        print "np.dot didn't work as I expected"
    if dot_answer_1 == dot_answer_2:
        print "np.dot is transpose invariant"
 
def almost_equal(a,b, sigFig = 5):
    eps = 1*10**(-sigFig)
    if (a - b <= eps).all() and (b - a <= eps).all():
        return True
       
        
def test_inverse(a):
    a_inv = np.linalg.inv(a)

    tr = (np.eye(3) - np.dot(a_inv, a) <= eps)
    
    tr[0,0] = False

    return almost_equal(np.eye(3), np.dot(a_inv, a))
    

