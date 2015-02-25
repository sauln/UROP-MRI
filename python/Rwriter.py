# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 16:32:07 2015

@author: nathaniel

this is still dealing with only landmarks

we want to first generate a bunch of fake images,
then we will pretend that these fake images are our acctual images.



with our set of images - 20 before, 20 after,

we will use procrustes so align and normalize,

then save the images into a format that R can use.

x1 y1 z1 for each 8 landmarks


"""
import numpy as np


def shapesToR(stack_of_shapes, root= "../image"):
    '''
    we should take the stack of shapes and write each one to a file
    '''
    print "starting"
    fn = r'%s/stackToR.txt' %root
    f = open(fn, 'w')
    print len(stack_of_shapes[0])
    names = ["x%s y%s z%s" %(x,x,x) for x in xrange(1,len(stack_of_shapes[0])+1)]
    names = ' '.join(names) + '\n'
    f.write(names)
    print names
    for each in stack_of_shapes:
        data = np.array_str(each.flatten(), precision=12)[1:-1]
        data = data.replace("\n"," ") + '\n'
        f.write(data)
       
        

        
    f.close() 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

