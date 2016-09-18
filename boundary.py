"""

"""
import numpy as np

def set_elem(bt, bc, line_type, nn):
    bound = np.zeros((2, nn))
    jj = 0
    for line, nodes in line_type.items():
        for _ in nodes:
            bound[0, jj] = bt[line]
            bound[1, jj] = bc[line]
            jj = jj + 1
    return bound

def mount_vector(z, nn, ne, bound):
    # mount the boundary T and q
    T = np.zeros(nn)
    q = np.zeros(nn)
    for elem in range(ne):
        if bound[0, elem] == 0:
            T[elem] = bound[1, elem]
            q[elem] = z[elem]
        else:
            T[elem] = z[elem]
            q[elem] = bound[1, elem]
    return T, q
    