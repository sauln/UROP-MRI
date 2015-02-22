# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 12:44:46 2015

@author: nathaniel


the functions in this file will handle bringing objects in from text files 
and load them into a shape object.


"""
import os
import patientGen


#from OpenGL.GL import *
import objRead


def newShapeFromFile(f):
    f = open(f, "r")
    sh = patientGen.Shape()
    ind = 0
    print type(sh.x[0,0])
    for line in f:
        l = [float(x) for x in line.split(' ')]
        sh.x[ind] = l[1:]
        ind+=1
    return sh

        
        
def main():
    f = r"..\image\rawlandmarks.txt"
    sh = newShapeFromFile(f)   
    sh.plotMe()
    

    verts, norms = objRead.loadOBJ(r"../image/surfaceWOmtl.obj")
    print len(verts)
    
    surface = objRead.OBJ("../image/surface.obj")  
        
    
    
    
    
if __name__ == "__main__":
    main() 
    
    
    
    
    
    
    
    