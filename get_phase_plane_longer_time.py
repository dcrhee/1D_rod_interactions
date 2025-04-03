#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 11:32:52 2024

@author: cotton

Rund the code for longer with the initial condition if we didn't reach the end of the run'

"""


import numpy as np
from scipy.integrate import odeint, solve_ivp, quad
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter
from matplotlib.patches import Arc, RegularPolygon
from numpy import radians as rad
#from rotation_solver import * #all the functions used to solve the ode
#from rotation_solver_fall_off import * #all the functions used to solve the ode
from rotation_solver_fall_off_correct import *
from plot_rot_simple import *

def find_scaled_minhmin_S_prime(y):
    theta = abs(y[4,:]-y[5,:])
    theta_b = y[5,:]
    
    z_aflat = y[2,:]*np.cos(theta_b) + y[0,:]*np.sin(theta_b)
    z_bflat = y[3, :]*np.cos(theta_b) + y[1, :]*np.sin(theta_b)
    hmin = z_aflat - theta/2 - z_bflat
    minhmin = np.min(hmin)
    
    return minhmin

# initial vertical velocities and angles to iterative over
vels = -np.logspace(-4, 3, 100)

#thetas = np.array([0]) #, np.linspace(0.01, 6.99, 100)])*np.pi/180

thetas = np.concatenate((np.array([0]), np.linspace(0.01, 6.99, 100)), axis =0)*np.pi/180
vels = -np.logspace(-2, 3, 100)
#thetas = np.concatenate((np.array([0]), np.linspace(0.01, 6.99, 500)), axis =0)*np.pi/180

gamma = 0.9
gamma = 1/2
hmin = 1 # initiate with hmin of 1

thetas = thetas[0:]
vels = vels[0:]

#thetas = [5*np.pi/180]
#vels = [-0.01]

pathname = '/Users/cotton/Documents/DPhil/Polynas/Code/non_eqm/Rotation/Python 1D case/Data/highres/gamma_' + str(gamma) + '/'

# find which events are triggered
# find the minimum distance reached and associated angle

vertatol = 1e-11
horizatol = 1e-11
thetaatol = 1e-11
atols = np.array([horizatol, horizatol, vertatol, vertatol, thetaatol, thetaatol, horizatol, horizatol, vertatol, vertatol, thetaatol, thetaatol])


events = np.zeros((len(vels), len(thetas)))-1
for vel_indx, vel in enumerate(vels):
    #if vel_indx > 96:
        for thet_indx, theta in enumerate(thetas):
        
            
            y0 = [0, 0, hmin + abs(theta)/2, 0, theta, 0, 0, 0, vel, 0, 0, 0]
            
            
            
            # first try and see if there is any event
            try:
                event_num = np.load(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_event.npy')
                #y = np.load(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '.npy')
                #t = np.load(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_t.npy')
                #ys= y[:, -1]
                #ts = t[-1]
                #np.save(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '.npy', ys, allow_pickle=True)
                #np.save(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_t.npy', ts, allow_pickle=True)
                
            except:
                # now try and see if there is any data
                try:
                   print(thet_indx, vel_indx) 
                    
                   y =  np.load(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '.npy')
                   t1 =  np.load(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_t.npy')
                   y0 = y[:, -1]
                   
                   tmax = 1000000/-vel
                   t_span = (0.0, tmax)
                   #result_solve_ivp = solve_ivp(differential_eqn, t_span, y0, events = [collision, no_contact_LHS, no_contact_RHS, move_away, relative_angle_too_large, theta_change], method = 'RK45', rtol = 1e-10, atol = atols, args = (gamma, ), max_step = 0.01)
                   result_solve_ivp = solve_ivp(differential_eqn, t_span, y0, events = [collision, no_contact_LHS, no_contact_RHS, move_away, relative_angle_too_large, theta_change], method = 'RK45', rtol = 1e-10, atol = atols, args = (gamma, ))
                   
                   
                   ys = np.concatenate((y, result_solve_ivp.y), 1)
                   ts = np.concatenate((t1, result_solve_ivp.t + t1[-1]))
                   
                   if result_solve_ivp.t_events[5].size > 0:
                       ys = np.concatenate((ys, result_solve_ivp.y_events[5].T), 1)
                       ts = np.concatenate((ts, result_solve_ivp.t_events[5]), 0)
                   
                   np.save(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '.npy', ys, allow_pickle=True)
                   np.save(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_t.npy', ts, allow_pickle=True)
                   
                   for event_indx, event in enumerate(result_solve_ivp.t_events):
                       if event_indx < 4:
                           if event.size > 0:
                               events[vel_indx, thet_indx] = event_indx
                               np.save(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_event.npy', event_indx, allow_pickle=True)
                               print(event_indx)
                       if event_indx == 4:
                           if event.size > 0:
                               np.save(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_too_large_angle.npy', event_indx, allow_pickle=True)
                               print(event_indx)
                           
                   
                except Exception:
                    tmax = 1000000/-vel
                    t_span = (0.0, tmax)
                    #result_solve_ivp = solve_ivp(differential_eqn, t_span, y0, events = [collision, no_contact_LHS, no_contact_RHS, move_away, relative_angle_too_large, theta_change], method = 'RK45', rtol = 1e-10, atol = atols, args = (gamma, ), max_step = 0.01)
                    result_solve_ivp = solve_ivp(differential_eqn, t_span, y0, events = [collision, no_contact_LHS, no_contact_RHS, move_away, relative_angle_too_large, theta_change], method = 'RK45', rtol = 1e-10, atol = atols, args = (gamma, ))
                    
                    
                    if result_solve_ivp.t_events[5].size > 0:
                       ys = np.concatenate((result_solve_ivp.y, result_solve_ivp.y_events[5].T), 1)
                       ts = np.concatenate((result_solve_ivp.t, result_solve_ivp.t_events[5]), 0)
                    else:
                       ys = result_solve_ivp.y
                       ts = result_solve_ivp.t
                    
                    np.save(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '.npy', ys, allow_pickle=True)
                    np.save(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_t.npy', ts, allow_pickle=True)
                    
                    for event_indx, event in enumerate(result_solve_ivp.t_events):
                        if event_indx < 4:
                            if event.size > 0:
                                events[vel_indx, thet_indx] = event_indx
                                np.save(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_event.npy', event_indx, allow_pickle=True)
                                print(event_indx)
                        if event_indx == 4:
                            if event.size > 0:
                                np.save(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_too_large_angle.npy', event_indx, allow_pickle=True)
                                print(event_indx)

np.save(pathname + 'thetas_gamma_' + str(gamma) + '.npy', thetas)
np.save(pathname + 'vels_gamma_' + str(gamma) + '.npy', vels)
#np.save(pathname + 'events_gamma_' + str(gamma) + '.npy', events)   
           
fig, (ax1) = plt.subplots(figsize=(13, 3), ncols=1)
pos = ax1.imshow(events, interpolation='none')

# add the colorbar using the figure's method,
# telling which mappable we're talking about and
# which Axes object it should be near
fig.colorbar(pos, ax=ax1)


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
Y = vels
X = thetas
X, Y = np.meshgrid(X, Y)
Z = events

surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

#ax.elev = 0
#ax.azim = 270  # xz view

#ax.elev = 0
#ax.azim = 0    # yz view

ax.elev = 90
ax.azim = 0  # xy view