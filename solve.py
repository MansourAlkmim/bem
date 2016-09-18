# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 09:47:52 2016

numerical line integration using linear quadrature quadrature

@author: mansour
"""
import numpy as np
import geometry

def quad(nodes, nodes_bound_coord, coord_font):
    """
    
    """
    jacobian, xsi, weight = gauss_legendre(nodes, nodes_bound_coord)
    coord_node = nodes_bound_coord[nodes[0]]
    rad, norm = geometry.vector(coord_font, coord_node)
    magr = np.linalg.norm(rad)
    
    sumF = 0.0
    sumG = 0.0
    # conductivity
    k = 1
    for i in xsi:
        sumF = sumF + np.log(magr)* jacobian * weight
        sumG = sumG + ((rad[0]*norm[0] + rad[1]*norm[1])/magr**2) * jacobian * weight
    intF = -sumF/(2*np.pi*k)
    intG = sumG/(2*np.pi)

    return intF, intG

# linear relation (x,y) > (-1,1)
def gauss_legendre(nodes, nodes_bound_coord):
    
    x1, y1 = nodes_bound_coord[nodes[0]]    
    x2, y2 = nodes_bound_coord[nodes[1]]
    
    xsi, weight = np.polynomial.legendre.leggauss(1)
    N1 = (1/2)*(1-xsi); dN1_dxsi = np.diff(N1, xsi)
    N2 = (1/2)*(1+xsi); dN2_dxsi = np.diff(N2, xsi)
    dx_dxsi = dN1_dxsi*x1 + dN2_dxsi*x2
    dy_dxsi = dN1_dxsi*y1 + dN2_dxsi*y2
    jacobian = np.sqrt(dx_dxsi**2 + dy_dxsi**2)

    return jacobian, xsi, weight

def int_point(nn, nodes_int_cord, nodes_bound_coord, line_type):
    # interior point
    Fin = np.zeros((len(nodes_int_cord), nn))
    Gin = np.zeros((len(nodes_int_cord), nn))
    for ii, elem in enumerate(nodes_int_cord):
        
        jj = 0
        for line, conn in line_type.items():
            for nodes in conn:
                intF, intG = quad(nodes, nodes_bound_coord, elem)
                Fin[ii, jj] = intF
                Gin[ii, jj] = intG
                jj = jj + 1  
    return Fin, Gin
