#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 12:48:22 2022

@author: cotton

want to save the minima, what type of event occurs and the leaving time if possible

Runs the code through for a range of an initial variable:
    Saves
    1. the IVP.y solution
    2. whether they lose contact (1) in event_type
    3. whether they collide (2) in event_type
"""

import sys
import os

# Calculate the path to other_folder
current_dir = os.path.dirname(__file__)
parent_dir = os.path.dirname(current_dir)
fuc_path_dir = os.path.join(parent_dir, 'tests')

# Add the other_folder to the system path
sys.path.append(fuc_path_dir)

import numpy as np
from scipy.integrate import odeint, solve_ivp, quad
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter
from matplotlib.patches import Arc, RegularPolygon
from numpy import radians as rad
from rotation_solver_fall_off_correct import *
#from rotation_solver_horiz_only import differential_eqn_horiz_only
#from rotation_solver_vert_only import differential_eqn_vert_only
#from plot_test_vert_horiz import *
from plot_rot_simple import *
from plot_multiple import plot_multjple_all_grid, plot_multiple_relative_grid
#from rotation_solver_omega_only import differential_eqn_omega_only


gamma = 1
hmin = 1
#vel = -100
#x_a, x_b, y_a, y_b, theta_a, theta_b, u_a, u_b, v_a, v_b, omega_a, omega_b = state

pathname = '/Users/cotton/Documents/DPhil/Polynas/Code/non_eqm/Rotation/Python 1D case/Data/gamma_' + str(gamma) + '/'
pathname = '/Users/cotton/Documents/DPhil/Polynas/Code/non_eqm/Rotation/Python 1D case/Data/highres/gamma_' + str(gamma) + '/'
thetas = np.load(pathname + 'thetas_gamma_' + str(gamma) + '.npy')
vels = np.load(pathname + 'vels_gamma_' + str(gamma) + '.npy')

vel = vels[48]
theta = thetas[162]

theta2 = thetas[163] #thetas[45]
 

vel = -2.36449
theta = 4.011*np.pi/180
theta2 = theta

y0 = [0, 0, hmin + abs(theta)/2, 0, theta, 0, 0, 0, vel, 0, 0, 0]
y1 = [0, 0, hmin + abs(theta)/2, 0, theta2, 0, 0, 0, vel, 0, 0, 0]
#y0horiz = [0, 0, 1, 0, 5*np.pi/180, 0, 1, 0, 0, 0, 0, 0]
#y0omega = [0, 0, 1, 0, 5*np.pi/180, 0, 0, 0, 0, 0, 0.1, 0]
t_span = (0.0, 100.0)

tmax = 1000000/-vel
t_span = (0.0, tmax)

#result_solve_ivp = solve_ivp(differential_eqn, t_span, y0, events = [collision, no_contact_LHS, no_contact_RHS, move_away, relative_angle_too_large, theta_change], method = 'RK45', args = (gamma, ), rtol = 1e-10, atol = 1e-11)
#result_solve_ivp2 = solve_ivp(differential_eqn, t_span, y1, events = [collision, no_contact_LHS, no_contact_RHS, move_away], method = 'RK45', args = (gamma, ), rtol = 1e-10, atol = 1e-11)



result_solve_ivp = solve_ivp(differential_eqn, t_span, y0, events = [collision, no_contact_LHS, no_contact_RHS, move_away, relative_angle_too_large, theta_change], method = 'RK45', rtol = 1e-13, atol = 1e-14, args = (gamma, ))
result_solve_ivp2 = solve_ivp(differential_eqn, t_span, y1, events = [collision, no_contact_LHS, no_contact_RHS, move_away], method = 'RK45', args = (gamma, ), rtol = 1e-14, atol = 1e-15)

if result_solve_ivp.t_events[5].size > 0:
    ys = np.concatenate((result_solve_ivp.y, result_solve_ivp.y_events[5].T), 1)
    ts = np.concatenate((result_solve_ivp.t, result_solve_ivp.t_events[5]), 0)
    t_sorted_indices = np.argsort(ts, axis=0).flatten()
    result_solve_ivp.t = np.sort(ts)
    result_solve_ivp.y = ys[:,t_sorted_indices]


#result_solve_ivp_vert = solve_ivp(differential_eqn_vert_only, t_span, y0, events = [collision, no_contact_LHS, no_contact_RHS, move_away_RHS, move_away_LHS], method = 'RK45', rtol = 1e-9, atol = 1e-12, args = (gamma, ))#, max_step = 0.001)
#result_solve_ivp_horiz = solve_ivp(differential_eqn_horiz_only, t_span, y0horiz, events = [collision], method = 'RK45', rtol = 1e-9, atol = 1e-12, args = (gamma, ))#, max_step = 0.001)
#result_solve_ivp_omega = solve_ivp(differential_eqn_omega_only, t_span, y0omega, events = [collision], method = 'RK45', rtol = 1e-7, atol = 1e-9, args = (gamma, ))#, max_step = 0.001)
#result_solve_ivp_omega = solve_ivp(differential_eqn_omega_only, t_span, y0omega, events = [collision], method = 'RK45', args = (gamma, ))#, max_step = 0.001)


plot_forces_torques(result_solve_ivp2, gamma)
plot_forces_torques_components(result_solve_ivp2, gamma)
plot_relative_grid(result_solve_ivp2)# plot the relative orientation as a grid
#plot_all_grid(result_solve_ivp) # plot the complete orientation as a grid
#check_no_trans_rot_solution(result_solve_ivp_vert, gamma)
#check_no_vert_rot_solution(result_solve_ivp_horiz, gamma)
#plot_relative_stills(result_solve_ivp, 1, 1, n_squares = 32, aspect = 'equal') # plot a series of stills showing the relative motion
# plot_relative_stills_one_box(result_solve_ivp, a_rad, b_rad, n_squares = 32, aspect = ' not equal') # plot a series of stills showing the relative motion
# plot_stills(result_solve_ivp, a_rad, b_rad, n_squares = 32, aspect = 'not equal') # plot a series of stills showing the complete motion
#plot_pressure_int_result(result_solve_ivp, a_rad)

#plot_om_vels_test(result_solve_ivp_vert, result_solve_ivp_horiz, result_solve_ivp_omega, gamma)

figure_sol, axs_all = plt.subplots(3, 4, constrained_layout=True, sharex = True)
plot_multjple_all_grid(axs_all, result_solve_ivp)
plot_multjple_all_grid(axs_all, result_solve_ivp2)

figure_sol_rel, axs_rel = plt.subplots(2, 4, constrained_layout=True, sharex = True)
plot_multiple_relative_grid(axs_rel, result_solve_ivp)
plot_multiple_relative_grid(axs_rel, result_solve_ivp2)
