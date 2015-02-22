# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 11:59:59 2015

@author: nathaniel
"""

#playing with surfaces and triangulations

import numpy as np
from mayavi import mlab
vertices = np.array([[0, 1, 0, 0],[0, 0, 1, 0],[0, 0, 0, 1]])
faces = np.array([[0, 1, 0, 0],[1, 2, 1, 2],[2, 3, 3, 3]])
mlab.triangular_mesh(vertices[0,:], vertices[1,:], vertices[2,:], faces.T)
mlab.show()