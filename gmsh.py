# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 11:16:37 2016

@author: mansour

review object oriented definitions:
1) classes can have objects
2) Object can have multiple instances
3) there are classes atributtes that are for all objects
4) there are object atributtes for each instance
5) the input of a object are called arguments
"""

import numpy as np
import os
import re


def find_num(string):
    """Find all numbers in a string
    """
    num = re.findall(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)', string)
    return num


class Parse:
    """Parse the .msh:
    
    Arguments:
        filename: (string) name of the file with mesh data
        
    Attributes:
        NODE: (array of floats) [nn x 3] array with nodes coordinates
        CONN: (array of int) [ne x 2] array with line connectivity
        INTNODE: (array of floats) [nn x 3] array with interior nodes coordinates
 
    Instructions:   
        In order to create a mesh for the code follow the steps:
        1. Draw the geometry
            1.1 Draw Points.
            1.2 Connect these points with straight lines.
            1.3 Define a plane surface.
        2. Add physical groups where the BC are going to be applied
            2.1 Add lines in order to get the nodes in that line.
            2.2 Add a Surface in order to get the connectivity.
        3. Change the subdivision algorithm to "all quads" ontools-mesh-general tab.
        4. Click on 2D to create the mesh. It should contain quad elements only.
        5. Save the .msh file with the mesh with the same name as the .geo.
  
    Method:
        1. Read a 2.0 ascii .msh file
   """
    def __init__(self, filename):
        msh_path = os.path.join(filename+'.msh')
        msh_file = open(msh_path, 'r')
        
        # dict node [coordx coordy] 
        XYZ = {}
        # list [node1_tag node2_tag]
        CONN1 = []
        # list [node1_tag node2_tag node3_tag node4_tag]
        CONNint = []
        # list boundary nodes [node1_tag node2_tag]
        node_bound = []
        # list line type physical group
        line_type = {}

        for txt_line in msh_file:
            num_list = find_num(txt_line)
                      
            # nodes coordinates xyz
            if len(num_list) == 4:
                n_tag = int(num_list[0]) -1
                XYZ[n_tag] = [float(f) for f in num_list[1:3]]          

            # type 1 = line with two nodes
            if len(num_list) == 7: 
               CONN1.append([int(f) - 1 for f in num_list[5:]]) 
               n_tag = int(num_list[4]) -1
               line_type.setdefault(n_tag,[])
               line_type[n_tag].append([int(f) - 1 for f in num_list[5:]])
               

            # type 2 = line with 3 nodes- triangle elements - interior points
            if len(num_list) == 8: 
               CONNint.append([int(f) - 1 for f in num_list[5:]])             
             
             # type 3 = line with 4 nodes- quad elements - interior points
            if len(num_list) == 9: 
               CONNint.append([int(f) - 1 for f in num_list[5:]])
        
        self.CONN1 = np.array(CONN1)
        self.CONNint = np.array(CONNint)
        
        # unique nodes ignore repeated nodes
        node_bound = np.unique(CONN1)
        node_all = np.unique(CONNint)

        # get nodes from node_all that are not in the boundary
        temp = set(node_bound)
        node_int = [x for x in node_all if x not in temp]
   
        # make it an array
        self.XYZ = np.array(list(XYZ.values()))
        
        # dictionary XYZ to get the coord of the boundary nodes
        self.nodes_bound_coord = self.XYZ[node_bound]
       
       # dictionary XYZ to get the coord of the interior nodes
        self.nodes_int_cord = self.XYZ[node_int]
       
       # number of elements         
        self.ne = len(CONN1)
        
        # number of nodes in the boundary            
        self.nn = len(self.nodes_bound_coord)
       
       # physical group, line type
        self.line_type = line_type
    