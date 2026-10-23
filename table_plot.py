#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 11:53:21 2026

@author: Adam Frost
"""

import matplotlib.pyplot as plt
import Frost_Adam as vortex
import math
import numpy as np
import sympy as sp

def plot_three(plot1x, plot1y, name1, 
               plot2x, plot2y, name2, 
               plot3x, plot3y, name3, 
               plottitle):
    # plots three lift coefficients on the same plot
    fig, ax = plt.subplots(figsize=(5, 3), layout='constrained')

    # vortex panel
    ax.plot(
        plot1x,
        plot1y,
        linestyle='None',
        marker='o',
        color='black',
        label=name1
    )

    # TAT
    ax.plot(
            plot2x,
            plot2y,
            linestyle='None',
            marker='o',
            color='red',
            label=name2
    )

    # experimental
    ax.plot(
            plot3x,
            plot3y,
            linestyle='None',
            marker='o',
            color='blue',
            label=name3
    )
        
    # set up plot
    ax.set_xlabel(r'Angle of Attack $\alpha$ ($^\circ$)')
    ax.set_ylabel(r'Lift Coefficient $\tilde{C}_L$')
    ax.set_title(plottitle)
    ax.legend(loc="best")
    ax.grid(True)
    plt.show()

#%% a2. determine cl, clme, and cmc4 values for 2412 airfoil from alpha = -12 to 16
vortex_table = []
geometry = '2412_200.txt'
alphas = [np.deg2rad(i) for i in np.linspace(-12,16,num=15)]
alphas_deg = np.rad2deg(alphas)
vel_inf = 10.0
vortex_table = vortex.calc_lift_moment_table(geometry, alphas, vel_inf)
# format data into a latex table for easy report insertion
header = vortex_table[0]
data = vortex_table[1:].astype(np.float64)
data[:,1:]=np.round(data[:,1:],10)
data[:,0]=np.round(data[:,0],1)
print(f'A-2 LaTeX-formatted code:\n{sp.latex(sp.TableForm(data))}')

#%% a3. plot this using matplotlib along with TAT and reference
# get TAT approximations
alphal0 = np.deg2rad(-2.077)
TAT = []
for alpha in alphas:
    TAT.append(2*math.pi*(alpha-alphal0))
    
# import experimental data
alphas_ex = [
    20.437000, 18.246800, 16.312200, 15.796700, 15.281500, 
    14.220000, 12.064400, 10.004900, 7.9769700, 5.8516700, 
    3.7577400, 1.6630400, -0.1747960, -2.3033500, -4.2074000, 
    -6.4343300, -8.4368700
]
cl_ex = [
    1.1337800, 1.2089600, 1.5095600, 1.5650900, 1.5944800, 
    1.5275300, 1.3984900, 1.2432800, 1.0438800, 0.8558550, 
    0.6481280, 0.4566830, 0.2455180, 0.0473698, -0.1803510, 
    -0.3753870, -0.6049370
]

plot_three(alphas_deg, np.array(vortex_table[1:,1],  dtype=float), 'Vortex Panel Method',
           alphas_deg, np.array(TAT, dtype=float), 'Thin Airfoil Theory',
           alphas_ex, cl_ex, 'Experimental Data',
           'NACA 2412 Lift Coefficient Data')

#%% a4. fill out a2 table for NACA 2421 and plot lift coefficient like a3
# table first
vortex_table = []
geometry = '2421_200.txt'
# use the same alphas and vel_inf
vortex_table = vortex.calc_lift_moment_table(geometry, alphas, vel_inf)
# format data into a latex table for easy report insertion
header = vortex_table[0]
data = vortex_table[1:].astype(np.float64)
data[:,1:]=np.round(data[:,1:],10)
data[:,0]=np.round(data[:,0],1)
print(f'\nA-4 LaTeX-formatted code:\n{sp.latex(sp.TableForm(data))}')

# now plot the three lift coefficients
# TAT returns the same result since alphal0 is unchanged
# experimental is different
alphas_ex = [
    24.515300, 22.295800, 20.306400, 18.185900, 16.288500, 
    14.071300, 12.265600, 10.004300, 8.194470, 6.057190, 
    4.082450, 1.911920, -0.031767, -2.170220, -3.987940, 
    -6.160960, -10.310200, -12.256900, -14.396600
]
cl_ex = [
    1.083100, 1.061840, 1.131360, 1.169530, 1.179620, 
    1.229380, 1.173390, 1.039790, 0.925950, 0.726202, 
    0.557802, 0.353030, 0.176306, -0.018620, -0.235073, 
    -0.458164, -0.828423, -1.021860, -1.202110
]
plot_three(alphas_deg, np.array(vortex_table[1:,1],  dtype=float), 'Vortex Panel Method',
           alphas_deg, np.array(TAT, dtype=float), 'Thin Airfoil Theory',
           alphas_ex, cl_ex, 'Experimental Data',
           'NACA 2421 Lift Coefficient Data')

#%% a5. fill out a2 table for NACA 0015 using vortex and XFOIL then plot Cl from vortex, XFOIL, and TAT
# vortex table first
vortex_table = []
geometry = '0015_200.txt'
# use the same alphas and vel_inf
vortex_table = vortex.calc_lift_moment_table(geometry, alphas, vel_inf)
# format data into a latex table for easy report insertion
header = vortex_table[0]
data = vortex_table[1:].astype(np.float64)
data[:,1:]=np.round(data[:,1:],10)
data[:,0]=np.round(data[:,0],1)
print(f'\nA-5 LaTeX-formatted code:\n{sp.latex(sp.TableForm(data))}')

# then the xfoil data
x_alphas = [
    -12.0, -11.0, -10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, 
    -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 
    8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0
]
x_cl = [
    -1.3137, -1.2258, -1.1020, -0.9792, -0.8737, -0.7711, -0.6649, -0.5568, -0.4474, -0.3369, 
    -0.2253, -0.1128, 0.0000, 0.1127, 0.2253, 0.3368, 0.4474, 0.5568, 0.6648, 0.7712, 
    0.8737, 0.9792, 1.1021, 1.2259, 1.3137, 1.3885, 1.4594, 1.5220, 1.5797
]
x_cd = [
    0.01335, 0.01221, 0.01098, 0.00989, 0.00900, 0.00821, 0.00752, 0.00693, 0.00645, 0.00608, 
    0.00585, 0.00570, 0.00565, 0.00571, 0.00585, 0.00609, 0.00645, 0.00693, 0.00752, 0.00820, 
    0.00900, 0.00989, 0.01098, 0.01221, 0.01335, 0.01462, 0.01621, 0.01827, 0.02098
]
x_cmc4 = [
    -0.0039, -0.0001, -0.0037, -0.0067, -0.0056, -0.0036, -0.0024, -0.0015, -0.0008, -0.0004, 
    -0.0001, -0.0000, -0.0000, 0.0001, 0.0001, 0.0004, 0.0008, 0.0015, 0.0024, 0.0036, 
    0.0056, 0.0067, 0.0036, 0.0000, 0.0039, 0.0101, 0.0162, 0.0225, 0.0278
]
# calculate cmle using cn which we get from alpha, cl, and cd
x_cmle = []
for i, alpha in enumerate(np.deg2rad(x_alphas)):
    cn = x_cl[i]*math.cos(alpha)+x_cd[i]*math.sin(alpha)
    x_cmle.append(x_cmc4[i]-((1/4)*cn))
    
# xfoil table time!
x_table = [['alpha (circ)','CL','CmLE','Cmc/4']]
for i, alpha in enumerate(x_alphas):
    if i%2 == 0:
        x_table.append([alpha,x_cl[i],x_cmle[i],x_cmc4[i]])
x_table = np.array(x_table)
# format data into a latex table for easy report insertion
header = x_table[0]
data = x_table[1:].astype(np.float64)
data[:,1:]=np.round(data[:,1:],10)
data[:,0]=np.round(data[:,0],1)
print(f'\nA-5 LaTeX-formatted code for Xfoil:\n{sp.latex(sp.TableForm(data))}')

# TAT stays the same
# now plot each of them

plot_three(alphas_deg, np.array(vortex_table[1:,1],  dtype=float), 'Vortex Panel Method',
           x_alphas, x_cl, 'XFOIL Values',
           alphas_deg, np.array(TAT, dtype=float), 'Thin Airfoil Theory',
           'NACA 0015 Lift Coefficient Data')