# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 22:38:45 2016

@author: mansour
"""
import numpy as np
import solve
import geometry

def mount_matrix(nn, node_med, nodes_bound_coord, line_type):
    F = np.zeros((nn, nn))
    G = np.zeros((nn, nn))
    ii = 0
    for elem in node_med.values():
        for coord_med in elem:
            jj = 0
            
            for line, conn in line_type.items():
                for nodes in conn:
                    if ii == jj:
                       length = geometry.length(nodes, nodes_bound_coord)
                       k = 1
                       F[ii, jj] = (length/(2*np.pi*k))*(1 - np.log(length/2))
                       G[ii, jj] = 0
                       jj = jj + 1
                    else:
                        intF, intG = solve.quad(nodes, nodes_bound_coord, coord_med)
                        F[ii, jj] = intF
                        G[ii, jj] = intG
                        jj = jj + 1
            ii = ii + 1
    return F, G

def mount_linear_system(nn, F, G, bound):
    A = np.zeros((nn, nn))
    b = np.zeros(nn)
    for line in range(len(F)):
        for column, _ in enumerate(F[line, ...]):
            if bound[0, column] == 0: # T is known
                A[line, column] = -F[line, column]
                if line != column:
                    b[column] = -bound[1, column] * G[line, column]
                else: #line == column
                    b[column] = -bound[1, column] * G[line, column] - bound[1, column]/2
            else: #bound[line][column][0] == 1:  q is known
                b[column] = bound[1, column] * F[line, column]
                if line != column:
                    A[line, column] = G[line, column]
                else: #line == column
                     A[line, column] = G[line, column] - 1/2
    return A, b
    
def nodes_coord(line_type, nodes_bound_coord):
    nodes_coord_ord = np.zeros(np.shape(nodes_bound_coord))
    ii = 0
    for nodes in line_type.values():
        for node1, _ in nodes:
            nodes_coord_ord[ii] = nodes_bound_coord[node1]
            ii = ii + 1
    return nodes_coord_ord