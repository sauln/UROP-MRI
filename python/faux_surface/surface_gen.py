# -*- coding: utf-8 -*-
"""
Created on Wed Mar 04 20:28:25 2015

@author: nathaniel
"""





import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtr


def make_edges(vertices):
    ''' return a set of edges according to an arbitrary method
    I made up late at night - this is deprecated '''

    edges = dict()
    edges_n = dict()
    length = np.sqrt(vertices.shape[0])

    for i in xrange(length*length):
        #add the vertix to the dictionary, and add

        tmp = []
        tmp_n = []  #the same graph but coordinate free.

        if not i%length == 0:
            # add left edge
            tmp_n.append(i-1)
            tmp.append(vertices[i-1])
        if not i%length == length-1:
            #add right edge
            tmp_n.append(i+1)
            tmp.append(vertices[i+1])
        if not i < length:
            #add bottom edge
            tmp_n.append(i-length)
            tmp.append(vertices[i - length])
        if not i > length*length - length-1:
            # add top edge
            tmp_n.append(i+length)
            tmp.append(vertices[i + length])
        if (not i%length == length -1) and (not i > length*length - length - 1):
            #add upper diagaonal
            tmp_n.append(i+length + 1)
            tmp.append( vertices[i + length + 1])
        if (not i%length == 0) and (not i < length):
            #add lower diagonal
            tmp_n.append(i - length - 1)
            tmp.append( vertices[i - length - 1])


        edges[tuple(vertices[i])] = tmp
        edges_n[i] = tmp_n

    return edges, edges_n


def make_triangles(vertices):
    from scipy.spatial import Delaunay
    tri = Delaunay(vertices)
    return tri

def plot_triangles_2d(vertices, triangles):
    plt.plot(vertices[:,0], vertices[:,1], 'o')
    plt.triplot(vertices[:,0], vertices[:,1], triangles.simplices.copy())
    plt.show()

def f(x,y):
    ''' f: R2 -> R3 ; smooth transform'''
    return np.asarray((2*x, y+1, x**3+y))
    
def f_r(x,y):
    ''' f: R2 -> R3 ; z is random'''
    return np.asarray((2*x, y+1, np.random.normal(2, 0.5, len(x)) ))

def plot_edges(vertices, edges_n = None):
    if edges_n == None:
        edges_n = make_edges(vertices)[1]

    xs = [e[0] for e in vertices]
    ys = [e[1] for e in vertices]

    plt.scatter(xs, ys)

    for by in edges_n.keys():
        for to in edges_n[by]:
            #print by, to
            x = (vertices[by][0], vertices[to][0])
            y = (vertices[by][1], vertices[to][1])
            plt.plot(x, y)


def transform_vertices_3d(vertices, f_map):
    vertices_n =  f_map(vertices[:,0], vertices[:,1])
    vertices_n = np.transpose(vertices_n)
    return vertices_n

def plot_triangles_3d(vertices, simplices):
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.plot_trisurf(vertices[:,0], vertices[:,1],
                    vertices[:,2], triangles=simplices, cmap=cm.Spectral)


    plt.show()


def save_as_off(vertices, triangles, filename="surface.off"):
    f = open(filename, 'w')
    f.write("OFF\n")
    f.write("%s %s 0\n"%(vertices.shape[0], len(triangles.simplices)))
    for v in vertices:
        f.write("%s %s %s\n"%(v[0], v[1], v[2]))
    for t in triangles.simplices:
        f.write("3 %s %s %s\n"%(t[0], t[1], t[2]))
    f.close()
    
    
def save_graph(vertices, edges, filename="graph.txt"):      
    f = open(filename, 'w')
    f.write("v\n")
    #f.write("%s %s 0\n"%(vertices.shape[0], len(edges)))
    for v,i in zip(vertices, range(vertices.shape[0])):
        f.write("%s %s %s %s\n"%(i, v[0], v[1], v[2]))
        
    f.write("\n\ne\n")
    for t in edges:
        f.write("%s %s\n"%(t[0], t[1]))
    f.close()

def open_off(filename="surface.off"):
    f = open(filename, 'r')
    f.readline()
    (num_vert, num_tri, num_edges) =  f.readline().rstrip().split(' ')

    num_vert = int(num_vert)
    num_tri = int(num_tri)
    
    vertices = np.zeros((num_vert, 3))
    for i in xrange(num_vert):
        vertices[i] = [float(x) for x in f.readline().split()]
        
    triangles = np.zeros((num_tri, 3))
    for i in xrange(num_tri):
        triangles[i] = [float(x) for x in f.readline().split()][1:]
       
    f.close()
    return vertices, triangles

def edges_from_triangles(triangles):
    ''' take a set of triangles and create a graph of edges '''
    
    edges = set()
    
    for i in triangles:
        #print tuple(sorted(list([i[0], i[1]])))
        edges.add(tuple(sorted(list([i[0], i[1]]))))
        edges.add(tuple(sorted(list([i[2], i[1]]))))
        edges.add(tuple(sorted(list([i[0], i[2]]))))
    
        
    return edges
        

    
''' create very regular vertices '''
length = 5
x = np.linspace(0,1,length)
y = np.linspace(0,1,length)

vertices = []
for i in x:
    for j in y:
        vertices.append([i,j])

vertices =  np.asarray(vertices)


''' create triangles '''
triangles = make_triangles(vertices)
edges = edges_from_triangles(triangles.simplices)
#plot_triangles_2d(vertices, triangles)

''' transform vertices '''
#vertices_n = transform_vertices_3d(vertices, f_r)

#plot_triangles_3d(vertices_n, triangles.simplices)


''' save file and reopen '''
#save_as_off(vertices_n, triangles)
save_graph(vertices_n, edges)
#verts, tris = open_off()
#plot_triangles_3d(verts, tris)














