from itertools import permutations

import numpy as np
from numpy.linalg import det
from scipy.spatial import Delaunay, Voronoi, distance, voronoi_plot_2d
from scipy.spatial import cKDTree


#----------------------------------
def tetravol(a,b,c,d):
    '''Calculates the volume of a tetrahedron, given vertices a,b,c,d'''
    tetravol=abs(np.dot((a-d),np.cross((b-d),(c-d))))/6
    return tetravol

def triangle_area(a,b,c):
    " Ref : https://kr.mathworks.com/matlabcentral/answers/119984-how-to-find-the-area-of-the-triangles-formed-as-a-result-of-delaunay-triangulation"
    #area = abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
    area = abs((b[0]-a[0])*(c[1]-a[1]) - (c[0]-a[0])*(b[1]-a[1]))/2.
    return area

def intersection2(p1,p2,p3,p4):
    " Ref : https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection"
    t = det([[p1[0]-p3[0],p3[0]-p4[0]],[p1[1]-p3[1],p3[1]-p4[1]]])/det([[p1[0]-p2[0],p3[0]-p4[0]],[p1[1]-p2[1],p3[1]-p4[1]]])
    u = det([[p1[0]-p3[0],p1[0]-p2[0]],[p1[1]-p3[1],p1[1]-p2[1]]])/det([[p1[0]-p2[0],p3[0]-p4[0]],[p1[1]-p2[1],p3[1]-p4[1]]])

    return [[p1[0]+t*(p2[0]-p1[0]), p1[1]+t*(p2[1]-p1[1])],[p3[0]+u*(p4[0]-p3[0]),p3[1]+u*(p4[1]-p3[1])]]

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

"""
Below codes are for 2D case. 
"""

def find_boundary(vertex,mid_point,x1,x2,y1,y2,angle):
    if angle > 45. and angle < 135.:
        intersect = intersection2(vertex,mid_point,[x1,y2],[x2,y2]) # upper boundary
    elif angle > 135. or angle < -135.:
        intersect = intersection2(vertex,mid_point,[x1,y1],[x1,y2]) # left boundary
    elif angle > -135. and angle < -45.:
        intersect = intersection2(vertex,mid_point,[x1,y1],[x2,y1]) # bottom boundary
    elif angle > -45. or angle < 45.:
        intersect = intersection2(vertex,mid_point,[x2,y1],[x2,y2]) # right boundary
    elif angle == 45.:
        intersect = np.array([x2,y2])
    elif angle == 135.:
        intersect = np.array([x1,y2])
    elif angle == -135.:
        intersect = np.array([x1,y1])
    elif angle == -45.:
        intersect = np.array([x2,y1])
    else:
        print('Calculation of angle is wrong!')
        intersect = [np.inf,np.inf]
    return intersect


"""
Now code can find far_point. I borrow the code that calculate the far_point from
"scipy : voronoi_plot_2d" function.
"""
def area(vor,p,boundaries):
    '''Calculate area of 2D Voronoi cell based on point p'''
    area = 0
    x1,x2,y1,y2 = boundaries[0],boundaries[1],boundaries[2],boundaries[3]
    ind_vert_in_cell = vor.regions[vor.point_region[p]].copy()
    #print('length of index vertices in a cell : ',len(ind_vert_in_cell))
    center = vor.points.mean(axis=0)
    add_vertices = []

    if  -1 in ind_vert_in_cell:
        test_ridge_vertices = vor.ridge_vertices.copy()
        test_ridge_point = vor.ridge_points.copy()
        for permute in list(permutations(ind_vert_in_cell,2)):
            if -1 in permute: # e.g. [-1, 25 ] or [196, -1]
                try:
                    #print(permute)
                    ind_vert = test_ridge_vertices.index(list(permute))
                    ind_point = test_ridge_point[ind_vert]
                    #print('ridge point : ',ind_point)
                    # Find tangent vector of two ridge points
                    mid1 = vor.points[ind_point[0]]
                    mid2 = vor.points[ind_point[1]]
                    t = mid2-mid1
                    t /= np.linalg.norm(t)
                    n = np.array([-t[1],t[0]]) # Normal vector

                    mid_point = [(mid1[0]+mid2[0])/2. , (mid1[1]+mid2[1])/2.]
                    #print('median point : ',mid_point)

                    # Find direction
                    direction = np.sign(np.dot(mid_point - center, n)) * n
                    # Find boundary
                    linked_index = list(permute)[1]
                    vertex = list(vor.vertices[linked_index])
                    #print('linked vertex : ',vertex)
                    far_point = vertex + direction
                    # check where is the base?
                    test_vector = np.array([1.,0.])
                    inner_product = np.dot(test_vector,direction)
                    if inner_product < 0:
                        angle = -np.arccos(inner_product)*180./np.pi
                    else:
                        angle = np.arccos(inner_product)*180./np.pi
                    # First, check the postion of vertex in quadrant

                    intersect = find_boundary(vertex,far_point,x1,x2,y1,y2,angle)

                    add_vertices.append(intersect[0])
                except:
                    pass
    else:
        pass

    #print(ind_vert_in_cell)

    if -1 in ind_vert_in_cell:
        ind_vert_in_cell.remove(-1)
    origin_vertice = vor.vertices[ind_vert_in_cell].tolist()
    new_vertices = origin_vertice.copy()
    new_vertices = new_vertices + add_vertices

    area = 0
    tri=Delaunay(np.array(new_vertices),qhull_options='Q12 Qs Qz')
    for simplex in tri.simplices:
        area+=triangle_area(new_vertices[simplex[0]],new_vertices[simplex[1]],new_vertices[simplex[2]])
    #print(area)
    return area
    
def voronoi_2d_density(positions,boundaries):
    galvor=Voronoi(positions)
    galden=np.zeros((len(positions),))
    for i in range(len(positions)):
        galden[i]=1./area(galvor,i,boundaries)
    return galden

def voronoi_density(positions):
    galvor=Voronoi(positions)
    galden=np.zeros((len(positions),))
    for i in range(len(positions)):
        galden[i]=1./vol(galvor,i)
    return galden

def spherical_density3d(positions,radius):
    galden=np.zeros((len(positions),))
    volume = 4.*np.pi*radius*radius*radius/3.
    tree = cKDTree(positions)
    for i in range(len(positions)):
       num_gal = tree.query_ball_point(positions[i], r=radius,return_length=True)
       galden[i] = num_gal/volume
    return galden

def spherical_density2d(positions,radius):
    galden=np.zeros((len(positions),))
    area = np.pi*radius*radius
    tree = cKDTree(positions)
    for i in range(len(positions)):
       num_gal = tree.query_ball_point(positions[i], r=radius,return_length=True)
       galden[i] = num_gal/area
    return galden
       
    
