import numpy as np
from scipy.spatial import distance
from scipy.spatial import Voronoi,Delaunay, voronoi_plot_2d

#----------------------------------
def tetravol(a,b,c,d):
    '''Calculates the volume of a tetrahedron, given vertices a,b,c,d'''
    tetravol=abs(np.dot((a-d),np.cross((b-d),(c-d))))/6
    return tetravol

def vol(vor,p):
    '''Calculate volume of 3d Voronoi cell based on point p.'''
    dpoints=[]
    vol=0
    for v in vor.regions[vor.point_region[p]]:
        dpoints.append(list(vor.vertices[v]))
    tri=Delaunay(np.array(dpoints),qhull_options='Q12 Qs Qz')
    for simplex in tri.simplices:
        vol+=tetravol(np.array(dpoints[simplex[0]]),np.array(dpoints[simplex[1]]),np.array(dpoints[simplex[2]]),np.array(dpoints[simplex[3]]))
    return vol

def voronoi_density(positions):
    galvor=Voronoi(positions)
    galden=np.zeros((len(positions),))
    for i in range(len(positions)):
        galden[i]=1./vol(galvor,i)
    return galden