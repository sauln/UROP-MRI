# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np


import EDMA
import ourProcrustes
import experiments
import tTest

import patientGen
    
def example_procrustes(): 
    '''
    example of how we will want to use procrustes
    
    the function 'ourProcrustes' takes as input just the shape data,
    not the Shape class,  ie, just shape.x
    
    '''
    test = patientGen.genSet(0.1, 10)
    testShapes = [s.x for s in test]
    

    mean_shape, aligned_shapes = ourProcrustes.ourProcrustes(np.asarray(testShapes))
    #ourProcrustes.plotOurProcrustes(aligned_shapes, mean_shape)
    sh = patientGen.Shape()
    sh.x = mean_shape
    sh.plotMe()
    
    print mean_shape
    
def main():
    
    
    tTest.someTtestTesting()
    
    #experiments.someTtestTesting()
    #EDMA.EDMA()
    #experiments.edgeLengthMonteCarlo()
    #experiments.compareLandmarksOverTime()
    #experiments.genHistogramsOfData()


    
    
    

if __name__ == "__main__":
    main()
    
    
    
 