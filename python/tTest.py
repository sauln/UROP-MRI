# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 10:48:34 2015

@author: nathaniel


going to experiment with some hotelling t tests










"""


import patientGen
from scipy import stats

def createSets():
    before = patientGen.genSet(kwarg="justM")
    after = patientGen.genSet(0.1, kwarg="justM")
    return before, after



def someTtestTesting(before = None, after = None):  

    if before == None:
        before, after = createSets()
        
        

    for each in before:
        print each[0,0]
    #print before
    
    
    print "hello world".upper()
    
    print "the ttest_ind is intended for comparing 2 samples from different populations so we want all of the first landmark from each population"
    print stats.ttest_ind(before, after)
    
    
    
    
    print "This is a two-sided test for the null hypothesis that 2 related or repeated samples have identical average (expected) values."
    print stats.ttest_rel(before, after)

    
    #print stats.ttest_ind(after[0].x, before[1].x)
  
