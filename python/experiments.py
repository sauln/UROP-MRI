# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 12:54:56 2015

@author: nathaniel


This will be storage for some various tests we have done on the data


most of these are completely disjoint experiments








"""
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

import patientGen
import maths
import visualization



def edgeLengthMonteCarlo():
    '''
    This compares the edge lengths of each shape between shapes that have a varying
    amount of drift and plot those differences over the drift
    '''    
    testD = []
    contD = []
    
    l = np.arange(-0.3, 0.3, 0.002)
    print len(l)
    for drift in l:
        t,c = compareEdgeLengthsBetweenShapes(drift)
        testD.append(t)
        contD.append(c)
            
    testD = np.array(testD)
    contD = np.array(contD)
    #print "Average of control: %s \nAverage of test: %s" %(np.mean(contD), np.mean(testD))
    
    fig = plt.figure(10)
    plt.plot(l,testD)
    plt.plot(l,contD)
    plt.title("Average of control: %s \nAverage of test: %s" %(np.mean(contD), np.mean(testD)))
    plt.show()
    
    
    
def compareEdgeLengthsBetweenShapes(drift):
    '''
    returns the average edge length for the control set and the drift set
    for a variable drift
    '''    
    
    control, test = patientGen.patientList(drift)

    t = 0
    c = 0
    for shape in control:
        c += shape.afl
    for shape in test:
        t += shape.afl
            
    return (t/len(control), c/len(test))
 
def genHistogramsOfData():
    '''
        I want to generate a histogram of the differences
        
        1 histogram of the differences of the same landmarks over time
        and with different levels of drift
    
    '''
    drifts =  (2,1,0.5, 0.25,0.05,0.01,0)
    #drifts =  (0.5, 0.25,0.05,0.01,0)
    kwarg = None#'allShift' #'minus' None
    #test case- to get the sizing right
    diff = differencesBetweenLandmarksOverTime(0)
    x = np.zeros((diff.shape[0], len(drifts)))
    
    var = np.zeros((len(drifts), 1))
    mean = np.zeros((len(drifts), 1))
    #gather our data at with the different drifts
    for d,i in zip(drifts, range(0,len(drifts))):
        tmp = maths.returnNormsVector(differencesBetweenLandmarksOverTime(d, kwarg))
        var[i] = np.var(tmp)
        mean[i] = np.mean(tmp)
        x[:,i] = tmp.reshape(1, tmp.shape[0])
        
    
    #plot them all in a histogram
    visualization.complicatedHisto(x, drifts, var, mean)
    
   
def compareLandmarksOverTime():
    '''
    This function will compute the differences between landmarks over time and 
    over a few different levels of drift.  it will then plot the difference vectors
    and calculate the mean of the norms.


    '''    

    fig = plt.figure()
    ax = []

    ax.append(fig.add_subplot(221, projection='3d'))
    ax.append(fig.add_subplot(222, projection='3d'))
    ax.append(fig.add_subplot(223, projection='3d'))
    ax.append(fig.add_subplot(224, projection='3d'))
    
    for a, d in zip(ax, (0.01, 0.1,0.5 ,1.0)):
        diff1 = differencesBetweenLandmarksOverTime(0)
        diff2 = differencesBetweenLandmarksOverTime(d)
        d1m = maths.meanOfNorms(diff1)
        d2m = maths.meanOfNorms(diff2)
        a.scatter(diff1[:,0], diff1[:,1], diff1[:,2],  c='b', marker='o')
        a.scatter(diff2[:,0], diff2[:,1], diff2[:,2],  c='g', marker='^' )
        a.set_title('difference in means: %.3f\ndrift: %s'%(abs(d1m - d2m), d))
    
    


def differencesBetweenLandmarksOverTime(drift='a', kwarg=None):
    
    '''
    takes a variable drift and generates the control and test patient set
    then calculates the differences between each landmark and returns a difference vector.

    '''    
    
    before, after = patientGen.patientList(drift, kwarg)
    
    diff = []
    x = before[0].x.shape[0]* len(before)
    y =  before[0].x.shape[1]
        
    for eachB, eachA in zip(before, after): 
        diff.append( eachB.x-eachA.x)
   
    diff =  np.array(diff).reshape(x,y)

    return diff
   



   
"""
def pointwiseOverShapes():
    '''
    I don't think this is used any more either
    we want to find the difference between individual points of each shape
    '''
    s1 = Shape(0)
    s2 = Shape(1)
    
    c = np.zeros((s1.x.shape))

    cn = np.zeros((s1.x.shape[0],1))
    for b1, b2, x in zip(s1.x, s2.x, range(0, s1.x.size)):
        c[x] = b1-b2
        cn[x] = LA.norm(b1 - b2)
        
    return cn
   
def twoPointDiff():
    
    '''
    I think this is garbage
    '''
    print "geting twopoint"
    mu = 0
    sigma = 0.1
    
    a = ran.normal(mu, sigma, (10000,2))
    print a
    
    plt.scatter(a[:,1], a[:,0])
    heatMap(a[:,1], a[:,0])

def differencesIn2Shapes():
    '''
    more old stuff that is now garbage

    '''    
    
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
"""