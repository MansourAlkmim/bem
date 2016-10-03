# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 20:49:28 2016

@author: mansour

gmsh returns:

        # dictionary XYZ to get the coord of the boundary nodes
        self.nodes_bound_coord = self.XYZ[node_bound]
       
        # dictionary XYZ to get the coord of the interior nodes
        self.nodes_int_cord = self.XYZ[node_int]
        
        # number of elements         
        self.ne = len(CONN1)
       
       # number of nodes in the boundary            
        self.nn = len(self.nodes_bound_coord)   
    
"""
import numpy as np
import matplotlib.pyplot as plt
import gmsh
import geometry
import boundary
import solve 
import index

# boundary type [segm1 segm1 ...]
# types: 0: dirichlet; 1: neumann
bt = [1, 0, 1, 0]
# boundary condition [segm1 segm2 ...]
bc = [0, 1 , 0, 0]

example = gmsh.Parse('squaremesh2')

nodes_coord_ord = index.nodes_coord(example.line_type, example.nodes_bound_coord)

node_med = geometry.node_med(example.line_type, example.nodes_bound_coord)

bound = boundary.set_elem(bt, bc, example.line_type, example.nn)

F, G = index.mount_matrix(example.nn, node_med, example.nodes_bound_coord, example.line_type)

A, b = index.mount_linear_system(example.nn, F, G, bound)

z = np.linalg.solve(A, b)           

T, q = boundary.mount_vector(z, example.nn, example.ne, bound)

Fin, Gin = solve.int_point(example.nn, example.nodes_int_cord, example.nodes_bound_coord, example.line_type)

# Tint = Fin @ T - Gin @ q          
Tint = np.dot(Fin,T) - np.dot(Gin,q)
# for key in node_med:
#     xy=node_med[key]
#     for i in range(len(xy)):
#         xyno=xy[0]
#         x=xyno[0]
#         y=xyno[1]
#         print(x)
x=np.array([]);
y=np.array([]);
for i in range(len(node_med)):
    v=node_med[i]
    for j in range(len(v)):
            xy=v[j]
            x=np.append(x,np.array(xy[0]))
            y=np.append(y,np.array(xy[1]))
fig = plt.figure()
ax = fig.add_axes([.1, .1, .8, .8])
# x = np.append(nodes_coord_ord[:, 0], example.nodes_int_cord[:, 0])
# y = np.append(nodes_coord_ord[:, 1], example.nodes_int_cord[:, 1])
x = np.append(x, example.nodes_int_cord[:, 0])
y = np.append(y, example.nodes_int_cord[:, 1])
Z = np.append(T, Tint)
# C = ax.tricontourf(x, y, Z, 50, cmap='viridis')
C = ax.tricontourf(x, y, Z, 50)
ax.set_aspect('equal')
plt.colorbar(C)

plt.show()

# print(T,Tint)       
# print('boundary nodes coord', example.nodes_bound_coord)
# print('interior nodes', example.nodes_int_cord)
# print('connective matrix', example.CONN1)
# print('number of nodes', example.nn)
# print('number of elements', example.ne)
# print('line type phisical group', example.line_type)
# print('element half node', node_med)
# print('length',length)    
# print('boundary condition ',bound)
# print('radius vector between nodes and font node',rad)
# print('normal vector',norm)