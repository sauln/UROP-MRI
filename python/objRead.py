# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 23:28:34 2015

@author: nathaniel
"""



 	
 	

#This code provides a function for loading a model from a wavefront OBJ file (3d model format) into an OpenGL display list. It additionally loads any referenced material and texture files. It does not support multiple models in the OBJ file, nor any material properties other than diffuse colour (Kd) and diffuse texture (map_Kd).

#TODO:

#    support more of the OBJ format.
#    convert to using vertex/normal/color/texture_coord arrays 

#import pygame
#from OpenGL.GL import *
import time
 
''' loadOBJ was copied from this blog:  http://www.nandnor.net/?p=86'''
def loadOBJ(filename):  
    print "start loading"
    numVerts = 0  
    verts = []  
    norms = []  
    vertsOut = []  
    normsOut = [] 
    edgesOut = []
    for line in open(filename, "r"):  
        vals = line.split()  
        if vals[0] == "v":  
            v = map(float, vals[1:4])  
            verts.append(v)  
        if vals[0] == "vn":  
            n = map(float, vals[1:4])  
            norms.append(n)  
        if vals[0] == "f":  
            #print vals 
            tmp = []
           
            #what is this doing with the faces?
            for f in vals[1:]:  
                w = f.split("/")  
                # OBJ Files are 1-indexed so we must subtract 1 below  
                vertsOut.append(list(verts[int(w[0])-1]))  
                normsOut.append(list(norms[int(w[2])-1])) 
                tmp.append(int(f.split("/")[0]) - 1)
                numVerts += 1 
            edgesOut.append(tmp)
            
    return vertsOut, normsOut, edgesOut
    
    
    
    
    
    
def MTL(filename):
    contents = {}
    mtl = None
    for line in open(filename, "r"):
        if line.startswith('#'): continue
        values = line.split()
        if not values: continue
        if values[0] == 'newmtl':
            mtl = contents[values[1]] = {}
        elif mtl is None:
            raise ValueError, "mtl file doesn't start with newmtl stmt"
        elif values[0] == 'map_Kd':
            # load the texture referred to by this declaration
            mtl[values[0]] = values[1]
            surf = pygame.image.load(mtl['map_Kd'])
            image = pygame.image.tostring(surf, 'RGBA', 1)
            ix, iy = surf.get_rect().size
            texid = mtl['texture_Kd'] = glGenTextures(1)
            glBindTexture(GL_TEXTURE_2D, texid)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                GL_LINEAR)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
                GL_LINEAR)
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, ix, iy, 0, GL_RGBA,
                GL_UNSIGNED_BYTE, image)
        else:
            mtl[values[0]] = map(float, values[1:])
    return contents
 
class OBJ:
    def __init__(self, filename, swapyz=False):
        """Loads a Wavefront OBJ file. """
        self.vertices = []
        self.normals = []
        self.texcoords = []
        self.faces = []
 
        material = None
        for line in open(filename, "r"):
            if line.startswith('#'): continue
            values = line.split()
            if not values: continue
            if values[0] == 'v':
                v = map(float, values[1:4])
                if swapyz:
                    v = v[0], v[2], v[1]
                self.vertices.append(v)
            elif values[0] == 'vn':
                v = map(float, values[1:4])
                if swapyz:
                    v = v[0], v[2], v[1]
                self.normals.append(v)
            elif values[0] == 'vt':
                self.texcoords.append(map(float, values[1:3]))
            elif values[0] in ('usemtl', 'usemat'):
                material = values[1]
            elif values[0] == 'mtllib':
                self.mtl = MTL(values[1])
            elif values[0] == 'f':
                face = []
                texcoords = []
                norms = []
                for v in values[1:]:
                    w = v.split('/')
                    face.append(int(w[0]))
                    if len(w) >= 2 and len(w[1]) > 0:
                        texcoords.append(int(w[1]))
                    else:
                        texcoords.append(0)
                    if len(w) >= 3 and len(w[2]) > 0:
                        norms.append(int(w[2]))
                    else:
                        norms.append(0)
                self.faces.append((face, norms, texcoords, material))
 
        self.gl_list = glGenLists(1)
        glNewList(self.gl_list, GL_COMPILE)
        glEnable(GL_TEXTURE_2D)
        glFrontFace(GL_CCW)
        for face in self.faces:
            vertices, normals, texture_coords, material = face
 
            mtl = self.mtl[material]
            if 'texture_Kd' in mtl:
                # use diffuse texmap
                glBindTexture(GL_TEXTURE_2D, mtl['texture_Kd'])
            else:
                # just use diffuse colour
                glColor(*mtl['Kd'])
 
            glBegin(GL_POLYGON)
            for i in range(len(vertices)):
                if normals[i] > 0:
                    glNormal3fv(self.normals[normals[i] - 1])
                if texture_coords[i] > 0:
                    glTexCoord2fv(self.texcoords[texture_coords[i] - 1])
                glVertex3fv(self.vertices[vertices[i] - 1])
            glEnd()
        glDisable(GL_TEXTURE_2D)
        glEndList()
