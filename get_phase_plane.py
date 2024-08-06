#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 16:40:05 2024

@author: cotton

Gets the phase plane starting off with an initial range of velocities and angles

"""

import numpy as np
from scipy.integrate import odeint, solve_ivp, quad
import matplotlib.pyplot as plt
from rotation_solver import * #all the functions used to solve the ode

def find_scaled_minhmin_S_prime(y):
    theta = abs(y[4,:]-y[5,:])
    theta_b = y[5,:]
    
    z_aflat = y[2,:]*np.cos(theta_b) + y[0,:]*np.sin(theta_b)
    z_bflat = y[3, :]*np.cos(theta_b) + y[1, :]*np.sin(theta_b)
    
    hmin = z_aflat - theta/2 - z_bflat
    minhmin = np.min(hmin)
    
    return minhmin

# initial vertical velocities and angles to iterative over
vels = -np.logspace(-4, 3)
thetas = np.linspace(0.01, 6.99)*np.pi/180
gamma = 1e-5
hmin = 1 # initiate with hmin of 1

thetas = thetas[2:]

pathname = '/Users/cotton/Documents/DPhil reading/Polynas/Code/non_eqm/Rotation/Python 1D case/Data/gamma_' + str(gamma) + '/'

# find which events are triggered
# find the minimum distance reached and associated angle

events = np.zeros((len(vels), len(thetas)))-1

for thet_indx, theta in enumerate(thetas):
    for vel_indx, vel in enumerate(vels):
        
        y0 = [0, 0, hmin + abs(theta)/2, 0, theta, 0, 0, 0, vel, 0, 0, 0]
        
        print(thet_indx, vel_indx)
        
        
        try:
           y =  np.load(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '.npy') 
        except Exception:
            tmax = 10/-vel
            t_span = (0.0, tmax)
            result_solve_ivp = solve_ivp(differential_eqn, t_span, y0, events = [collision, no_contact_LHS, no_contact_RHS, move_away, relative_angle_too_large], method = 'RK45', rtol = 1e-5, atol = 1e-7, args = (gamma, ))#, max_step = 0.001)
            np.save(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '.npy', result_solve_ivp.y, allow_pickle=True)
            np.save(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_t.npy', result_solve_ivp.t, allow_pickle=True)
            
            for event_indx, event in enumerate(result_solve_ivp.t_events):
                if event.size > 0:
                    events[vel_indx, thet_indx] = event_indx
                    np.save(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_event.npy', event_indx, allow_pickle=True)

np.save(pathname + 'thetas_gamma_' + str(gamma) + '.npy', thetas)
np.save(pathname + 'vels_gamma_' + str(gamma) + '.npy', vels)
np.save(pathname + 'events_gamma_' + str(gamma) + '.npy', events)   