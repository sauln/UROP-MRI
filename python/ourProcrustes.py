# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 12:40:50 2015

@author: nathaniel

This will house our script to take a stack of shapes, run general procustes
to find the mean shape and use ordinary procrustes to move all our shapes
as close as possible to the mean shape 


input:   stack of shapes - a list of shapes, each shape a nx2 or nx3 numpy array representing the
         n points in that shape
returns: mean shape - a nx2 or nx3 ndarray
         moved shapes - a list of the input shapes, but translated so they are
         as close to the mean shape as possible

"""


import numpy as np
import procrupy
import matplotlib.pyplot as plt

import manager

def testOurProcrustes():
    ''' woo hoo this is the first official test I've ever written '''
    first_shape = np.array([[0,0], [0, 1], [1,1],[1,0]] ) 
    second_shape = np.array([[0,0],[0, 1], [1,1.5], [1,0.5]])  
    third_shape = np.array([[0.25,-0.54],[0, 1.57], [.97,1.23], [1.11,0.13]]) 
    stack_of_shapes =  np.asarray([first_shape, second_shape, third_shape])

   
    expected_mean_shape = [[-0.35056829, -0.43866594],
                    [-0.28789217,  0.37891625],
                    [ 0.34764185,  0.35111662],
                    [ 0.29081861, -0.29136693]]
  
    mean_shape, aligned_shapes = ourProcrustes(stack_of_shapes)
    np.testing.assert_allclose(mean_shape, expected_mean_shape)
    
    

def ourProcrustes(stack_of_shapes, kwarg = None):
    

    #use GPA to find the mean shape
    mean_shape = procrupy.generalized_procrustes_analysis(stack_of_shapes)


    #use OPA to align each of the shapes as close as possible with the mean shape
    aligned_shapes = []
    for each in stack_of_shapes:
        rotation, scale, translation = procrupy.procrustes_analysis(each, mean_shape)
        aligned_shapes.append ( np.dot(scale * each, rotation) + translation )
        
        
    return mean_shape, aligned_shapes
 

def plotOurProcrustes(stack_of_shapes, mean_shape):  
    #this has gotta be cleaned up a little bit:
    #currently it doesn't put our edges where we want them to be.
    
    
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import matplotlib.pyplot as plt
    
    
    
    fig = plt.figure()
    ax = Axes3D(fig)

    ax.add_collection3d(Poly3DCollection([mean_shape.tolist()]))
    for each in stack_of_shapes:
        ax.add_collection3d(Poly3DCollection([each.tolist()]))
    #plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
if __name__ == "__main__":
    manager.main() 
    
    
    
    
    '''
    ax = plt.axes()

    colors = (['red', 'green', 'blue'])
    
    for each,c in zip(stack_of_shapes, colors): 
        plt.gca().add_line(plt.Polygon(each, fill=None, edgecolor=c))
    
    polygon4 = plt.Polygon(mean_shape , fill=None, edgecolor='k')
    plt.gca().add_line(polygon4)
    
    
    plt.axis('scaled')
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    plt.show()  
    '''