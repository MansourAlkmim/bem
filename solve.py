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
    jacobian, xsi, weight,x,y,nx,ny = gauss_legendre(nodes, nodes_bound_coord)
    coord_node = nodes_bound_coord[nodes[0]]
    rad, norm = geometry.vector(coord_font, coord_node)
    magr = np.linalg.norm(rad)
    
    sumF = 0.0
    sumG = 0.0
    comprimento=0;
    # conductivity
    k = 1
    for i in range(len(xsi)):
        rx=x[i]-coord_font[0]
        ry=y[i]-coord_font[1]
        r=np.sqrt(rx**2+ry**2)
        sumG = sumG + np.log(r)* jacobian * weight[i]
        comprimento = comprimento + jacobian * weight[i]
        sumF = sumF + (rx*nx + ry*ny)/r**2 * jacobian * weight[i]
    intF = sumF/(2*np.pi)
    intG = -sumG/(2*np.pi*k)

    return intF, intG

# linear relation (x,y) > (-1,1)
def gauss_legendre(nodes, nodes_bound_coord):
    
    x1, y1 = nodes_bound_coord[nodes[0]]    
    x2, y2 = nodes_bound_coord[nodes[1]]
    
    xsi, weight = np.polynomial.legendre.leggauss(4)
    N1 = (1/2)*(1-xsi); dN1_dxsi = -1./2.
    N2 = (1/2)*(1+xsi); dN2_dxsi = 1./2.
    dx_dxsi = dN1_dxsi*x1 + dN2_dxsi*x2
    dy_dxsi = dN1_dxsi*y1 + dN2_dxsi*y2
    x=N1*x1+N2*x2;
    y=N1*y1+N2*y2;
    jacobian = np.sqrt(dx_dxsi**2 + dy_dxsi**2)
    nx=dy_dxsi/jacobian
    ny=-dx_dxsi/jacobian

    return jacobian, xsi, weight,x,y,nx,ny

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
