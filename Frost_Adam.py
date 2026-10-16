#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 22:57:25 2026

@author: Adam Frost
"""

import math
import numpy as np # does math stuff    
import json # interfaces with .json input file 
import matplotlib.pyplot as plt # does plotting stuff
#%% helper functions (some imported from class documentation)
def load_json(json_filename):
    # extract json dictionary info
    json_string = open(json_filename).read() # read it as txt
    input_vals = json.loads(json_string) # load it as a json dictionary using the json package
    # now the info is easy to extract! 
    geometry = input_vals['geometry']  
    alpha_rad = np.deg2rad(input_vals["alpha[deg]"])
    vel_inf = input_vals["freestream_velocity"]
    return geometry, alpha_rad, vel_inf

def get_nodes_from_geom(geometry): 
    # this get's the nodes from the geometry txt file 
    with open(geometry, "r") as f:
        info = []
        for line in f:
            info.append([float(coordinate) for coordinate in line.split()])
    # now make the info list a numpy array 
    nodes = np.array(info)
    return nodes

def calc_panel_lengths(nodes):
    ## calculate the panel lengths ##
    diffs = np.array(nodes[1:]-nodes[:-1])
    lengths = np.linalg.norm(diffs, axis=1)
    return lengths

def calc_control_point_locations(nodes): 
    ## calculate the control point locations ##
    sums = np.array(nodes[:-1]+nodes[1:])
    control_points = sums*0.5
    return control_points

def populate_a(nodes, N_nodes, N_panels, lengths, cpoints):
    A=np.zeros((N_nodes,N_nodes))
    # loop through control points
    for i in range(N_panels):
        xi, yi = cpoints[i]
        if i == N_panels-1:
            xi1, yi1 = cpoints[0]
        else:
            xi1, yi1 = cpoints[i+1]
        li = lengths[i]
        # get P for each panel
        for j in range(N_panels):
            x, y = cpoints[j]
            xj, yj = nodes[j]
            xj1, yj1 = nodes[j+1]
            lj = lengths[j]
            Pj = calc_P(x, y, lj, xj, yj, xj1, yj1)
            # add panel effect to A (incorporating n vector calculation)
            A[i][j] += ((xi1-xi)/li)*Pj[1][0] - (((yi1-yi)/li)*Pj[0][0])
            A[i][j+1] += ((xi1-xi)/li)*Pj[1][1] - (((yi1-yi)/li)*Pj[0][1])
    return A

def calc_P(x, y, lj, xj, yj, xj1, yj1):
    # get goofy greek letters
    # start with xi and eta
    node_arr = np.array([[xj1-xj, yj1-yj],
                         [-(yj1-yj), xj1-xj]])
    cp_arr = np.array([[x-xj],
                       [y-yj]])
    point_trans = (1/lj)*(node_arr@cp_arr)
    xi, eta = point_trans.flatten()
    # then grab phi
    phi = math.atan2(eta*lj, eta**2 + xi**2 - (xi*lj))
    # and finally psi!
    log_internals = ((xi**2 + eta**2) / ((xi-lj)**2 + eta**2))
    psi = 0.5*(math.log(log_internals))
    # then calculate final P value
    mat1 = np.array([[xj1-xj,   -(yj1-yj)],
                     [(yj1-yj), xj1-xj]])
    mat2 = np.array([[(lj-xi)*phi + eta*psi,        xi*phi - (eta*psi)],
                     [eta*phi - ((lj-xi)*psi) - lj, -eta*phi - (xi*psi) + lj]])
    P = (1/(2*math.pi*(lj**2)))*(mat1@mat2)
    return P

#%% main
if __name__ == "__main__":
    json_filename = "input.json" 
    geometry, alpha_rad, vel_inf = load_json(json_filename)
    nodes = get_nodes_from_geom(geometry)
    N_nodes = len(nodes)
    N_panels = N_nodes - 1
    lengths = calc_panel_lengths(nodes)
    cpoints = calc_control_point_locations(nodes)
    A = populate_a(nodes, N_nodes, N_panels, lengths, cpoints)

#%% final printout
    # print final information
    print('============GEOMETRY INFORMATION================')
    print(f"\nfile:\n{geometry}")
    print(f"\nalpha [deg]:\n{np.rad2deg(alpha_rad):.1f}")
    print(f"\nFreestream velocity:\n{vel_inf:.1f}")
    print(f"\nnodes:\n{nodes}")
    print(f"\nN (nodes):\n{N_nodes}")
    print(f"\nN (panels):\n{N_panels}")
    print('\n============CALCULATION RESULTS=================')
    print(f'\nA:\n{A}')

