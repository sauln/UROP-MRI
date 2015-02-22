# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 12:44:46 2015

@author: nathaniel


the functions in this file will handle bringing objects in from text files 
and load them into a shape object.


"""
import os
import patientGen
import numpy as np
import time


#from OpenGL.GL import *
import objRead

import matplotlib.pyplot as plt


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
    #f = r"..\image\rawlandmarks.txt"
    #sh = newShapeFromFile(f)   
    #sh.plotMe()
    
    print "import file"
    verts, norms, faces = objRead.loadOBJ(r"surfaceWOmtl.obj")
    
    print "convert to array"
    verts = np.asarray(verts)
    edges = np.asarray(faces)
    
    print edges.shape
    print verts.shape

    st = 0
    end = st+10
    print "plot the vertices"
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(verts[st:end,0], verts[st:end,1], verts[st:end,2],  c='b', marker='o')
    
    
    
    
    for i in range(100):
        line1 = np.transpose(np.asarray([verts[edges[i][0]] , verts[edges[i][1]]]))
        line2 = np.transpose(np.asarray([verts[edges[i][1]] , verts[edges[i][2]]]))
        line3 = np.transpose(np.asarray([verts[edges[i][2]] , verts[edges[i][0]]]))
        ax.plot(line1[0], line1[1], line1[2])
        ax.plot(line2[0], line2[1], line2[2])
        ax.plot(line3[0], line3[1], line3[2])
    #ax.plot(verts[edges[i][1]] , verts[edges[i][2]])
    #ax.plot(verts[edges[i][2]] , verts[edges[i][1]])
        
        #ax.scatter(verts[,0], verts[0:100,1], verts[0:100,2])
    #plt.scatter(verts)
    #print "begin loading obj"
    #surface = objRead.OBJ("surface.obj") 
    #print type(surface)
    print "finish"   
    
    
    
    
if __name__ == "__main__":
    main() 
    
    
    
    
    
    
    
    