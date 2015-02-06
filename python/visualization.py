# -*- coding: utf-8 -*-
"""
Created on Wed Feb 04 15:54:42 2015

@author: nathaniel
"""

'''

All of the plotting functions 

'''

def heatMap(x,y):
 
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=50)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    
    plt.clf()
    plt.imshow(heatmap, extent=extent)
    plt.show()
    
    
def histo(x):
    n, bins, patches = P.hist(x, 20, normed=1, histtype='stepfilled')
    P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
    plt.show()
