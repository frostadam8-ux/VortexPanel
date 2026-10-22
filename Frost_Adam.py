#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 22:57:25 2026

@author: Adam Frost
"""

import math
import numpy as np # does math stuff
import json # interfaces with .json input file 

DEBUG = False

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

def calc_lift_moment_table(geometry, alpha_set, vel_inf):
    ret_mat = [[r'$alpha$ ($^\circ$)', 'C_L', 'C_{m_{LE}}', 'C_{m_{c/4}}']]
    alpha_set = np.array(alpha_set)
    nodes = get_nodes_from_geom(geometry)
    N_nodes = len(nodes)
    N_panels = N_nodes - 1
    lengths = calc_panel_lengths(nodes)
    cpoints = calc_control_point_locations(nodes)
    A = calc_a(nodes, N_nodes, N_panels, lengths, cpoints)
    for alpha in alpha_set:
        B = calc_B(nodes, N_nodes, lengths, alpha)
        gamma = calc_gamma(A, B, vel_inf)
        Cl = calc_cl(N_panels, lengths, gamma, vel_inf)[0]
        Cmle = calc_cmle(lengths, nodes, N_panels, gamma, vel_inf, alpha)[0]
        # simple transformation to cmc4
        Cmc4 = Cmle + (0.25*Cl*math.cos(alpha))
        ret_mat.append([np.rad2deg(alpha), Cl, Cmle, Cmc4])
    return np.array(ret_mat)

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

def calc_a(nodes, N_nodes, N_panels, lengths, cpoints):
    A=np.zeros((N_nodes,N_nodes))
    # loop through control points
    for i in range(N_panels):
        x, y = cpoints[i]
        # panel i for the normal vector
        xi, yi = nodes[i]
        xi1, yi1 = nodes[i+1]
        li = lengths[i]
        # get P for each panel, evaluated at control point i
        for j in range(N_panels):
            xj, yj = nodes[j]
            xj1, yj1 = nodes[j+1]
            lj = lengths[j]
            Pj = calc_P(x, y, lj, xj, yj, xj1, yj1)
            # add panel effect to A (incorporating n vector calculation)
            A[i][j] += ((xi1-xi)/li)*Pj[1][0] - (((yi1-yi)/li)*Pj[0][0])
            A[i][j+1] += ((xi1-xi)/li)*Pj[1][1] - (((yi1-yi)/li)*Pj[0][1])
    # Kutta condition: gamma at the trailing edge sums to zero
    A[N_nodes-1][0] = 1.0
    A[N_nodes-1][N_nodes-1] = 1.0
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

def calc_B(nodes, N_nodes, lengths, alpha_rad):
    # pretty straightforward ngl
    B = np.zeros((N_nodes,1))
    for i in range(N_nodes-1):
        xi, yi = nodes[i]
        xi1, yi1 = nodes[i+1]
        li = lengths[i]
        B[i] = ((yi1-yi)*math.cos(alpha_rad) - ((xi1-xi)*math.sin(alpha_rad)))/li
    return B
        
def calc_gamma(A, B, vel_inf):
    # even more straightforward
    return np.linalg.solve(A, vel_inf*B)

def calc_cl(N_panels, lengths, gamma, vel_inf, c=1.0):
    # coordinates are normalized, but i'll keep a dummy c in here
    cl = 0.0
    for i in range(N_panels):
        # add to cl for each value in the sum
        li = lengths[i]
        gami = gamma[i]
        gami1 = gamma[i+1]
        cl += (li/c)*(gami+gami1)/vel_inf
    return cl

def calc_cmle(lengths, nodes, N_panels, gamma, vel_inf, alpha_rad, c=1.0):
    cmle = 0.0
    for i in range(N_panels):
        li = lengths[i]
        xi, yi = nodes[i]
        xi1, yi1 = nodes[i+1]
        gami = gamma[i]
        gami1 = gamma[i+1]
        # calculate the two main bits individually
        fracx = (2*xi*gami + (xi*gami1) + (xi1*gami) + (2*xi1*gami1))/(c*vel_inf)
        fracy = (2*yi*gami + (yi*gami1) + (yi1*gami) + (2*yi1*gami1))/(c*vel_inf)
        cmle += (li/c) * (fracx*math.cos(alpha_rad) + (fracy*math.sin(alpha_rad)))
    return -(1/3)*cmle

#%% main
if __name__ == "__main__":
    json_filename = "input.json" 
    geometry, alpha_rad, vel_inf = load_json(json_filename)
    nodes = get_nodes_from_geom(geometry)
    N_nodes = len(nodes)
    N_panels = N_nodes - 1
    lengths = calc_panel_lengths(nodes)
    cpoints = calc_control_point_locations(nodes)
    A = calc_a(nodes, N_nodes, N_panels, lengths, cpoints)
    B = calc_B(nodes, N_nodes, lengths, alpha_rad)
    gamma = calc_gamma(A, B, vel_inf)
    Cl = calc_cl(N_panels, lengths, gamma, vel_inf)[0]
    Cmle = calc_cmle(lengths, nodes, N_panels, gamma, vel_inf, alpha_rad)[0]
    # simple transformation to cmc4
    Cmc4 = Cmle + (0.25*Cl*math.cos(alpha_rad))

#%% final printout
    # print final information
    print('============GEOMETRY INFORMATION================')
    print(f"\nfile:\n{geometry}")
    print(f"\nalpha [deg]:\n{np.rad2deg(alpha_rad):.1f}")
    print(f"\nFreestream velocity:\n{vel_inf:.1f}")
    if DEBUG:
        print(f"\nnodes:\n{nodes}")
    print(f"\nN (nodes):\n{N_nodes}")
    print(f"\nN (panels):\n{N_panels}")
    print('\n============CALCULATION RESULTS=================')
    if DEBUG:
        print(f'\nA matrix:\n{A}')
        print(f'\nB matrix:\n{B}')
        print(f'\ngamma matrix:\n{gamma}')
    print(f'\nCl: {Cl}\n')
    print(f'CmLE: {Cmle}\n')
    print(f'Cmc4: {Cmc4}\n')
    print('*** For table/plot results, use table_plot.py ***')