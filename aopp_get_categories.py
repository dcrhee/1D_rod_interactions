#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 17:16:03 2024

@author: cotton

Categorises the runs into four different events

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import sys


def find_scaled_minhmin_S_prime(y):
    theta = abs(y[4,:]-y[5,:])
    theta_b = y[5,:]
    
    z_aflat = y[2,:]*np.cos(theta_b) + y[0,:]*np.sin(theta_b)
    z_bflat = y[3, :]*np.cos(theta_b) + y[1, :]*np.sin(theta_b)
    
    hmin = z_aflat - theta/2 - z_bflat
    minhmin = np.min(hmin)
    
    return minhmin

def count_oscillations(y):
    thetas = y[4,:]-y[5,:]
    thetas[thetas >= 0] = 1
    thetas[thetas < 0] = 0
    oscill_count = np.sum(np.abs(np.diff(thetas)))
    return oscill_count


def trim_ys(y, t):
    t_end = np.argmax(t)
    if t_end == np.size(t) - 1:
        y = y
    else:
        y = y[:, :t_end+1]
    return y


fig, _axs = plt.subplots(1, 4, constrained_layout=True, sharex = True, sharey = True)

axs = _axs.flatten()

fig2, _axs2 = plt.subplots(1, 4, constrained_layout=True, sharex = True, sharey = True)

axs2 = _axs2.flatten()

fig3, _axs3 = plt.subplots(1, 4, constrained_layout=True, sharex = True, sharey = True)

axs3 = _axs3.flatten()

fig4, _axs4 = plt.subplots(1, 4, constrained_layout=True, sharex = True, sharey = True)

axs4 = _axs4.flatten()

gamma = 0.75

pathname1 = '/network/group/aopp/oceans/AW006_COTTON_1DDISKS/1D_case/'
pathname2 = '/network/group/aopp/oceans/AW006_COTTON_1DDISKS/1D_case_075/'

thetas = np.load(pathname1 + 'thetas_gamma_' + str(gamma) + '.npy')
vels = np.load(pathname1 + 'vels_gamma_' + str(gamma) + '.npy')

hmins = np.zeros((len(vels), len(thetas)))-1
num_oscillations = np.zeros((len(vels), len(thetas)))-1
events = np.zeros((len(vels), len(thetas)))-1
categories = np.zeros((len(vels), len(thetas)))-1
theta_collision =  np.zeros((len(vels), len(thetas)))-1
theta_final =  np.zeros((len(vels), len(thetas)))-1

try:
    theta_final = np.load(pathname1 + '_theta_final.npy')
except Exception:
    for thet_indx, theta in enumerate(thetas):
        print(thet_indx)
        for vel_indx, vel in enumerate(vels):
            try:
                y = np.load(pathname1 + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '.npy')
                t = np.load(pathname1 + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_t.npy')                        
                y = trim_ys(y, t)
                if np.size(y) == 12:
                    theta_final[vel_indx, thet_indx] = (y[4] - y[5])*180/np.pi
                else:
                    theta_final[vel_indx, thet_indx] = (y[4, -1] - y[5, -1])*180/np.pi
            except Exception:
                try:
                    y = np.load(pathname2 + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '.npy')
                    t = np.load(pathname2 + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_t.npy')                        
                    y = trim_ys(y, t)
                    if np.size(y) == 12:
                        theta_final[vel_indx, thet_indx] = (y[4] - y[5])*180/np.pi
                    else:
                        theta_final[vel_indx, thet_indx] = (y[4, -1] - y[5, -1])*180/np.pi
                except Exception:
                    theta_final[vel_indx, thet_indx] = -1
                    
                
    np.save(pathname1 + '_theta_final', theta_final)


try:
    hmins = np.load(pathname1 + '_hmins.npy')
    categories = np.load(pathname1 + '_categories.npy')
    events = np.load(pathname1 + '_events.npy')
    theta_collision = np.load(pathname1 + '_theta_collision.npy')
except Exception:
    for thet_indx, theta in enumerate(thetas):
        for vel_indx, vel in enumerate(vels):
            
            try:
                event = np.load(pathname1 + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_event.npy')
                events[vel_indx, thet_indx] = event
                if event == 0:
                    y = np.load(pathname1 + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '.npy')
                    t = np.load(pathname1 + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_t.npy')                        
                    y = trim_ys(y, t)
                    
                    theta_collision[vel_indx, thet_indx] = (y[4, -1] - y[5, -1])*180/np.pi
            except Exception:
                try:
                    event = np.load(pathname2 + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_event.npy')
                    events[vel_indx, thet_indx] = event
                    if event == 0:
                        y = np.load(pathname2 + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '.npy')
                        t = np.load(pathname2 + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_t.npy')                        
                        y = trim_ys(y, t)
                        
                        theta_collision[vel_indx, thet_indx] = (y[4, -1] - y[5, -1])*180/np.pi
                except Exception:
                    events[vel_indx, thet_indx] = -1
    
    for thet_indx, theta in enumerate(thetas):
        print(thet_indx)
        for vel_indx, vel in enumerate(vels):
            
            try:
                y = np.load(pathname1 + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '.npy')
                t = np.load(pathname1 + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_t.npy')
                hmins[vel_indx, thet_indx] = find_scaled_minhmin_S_prime(y)
                
                y = trim_ys(y, t)
                
                num_oscillations[vel_indx, thet_indx] = count_oscillations(y)
                
            except Exception:
                try:
                    y = np.load(pathname2 + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '.npy')
                    t = np.load(pathname2 + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta*180/np.pi, 3)) + '_t.npy')
                    hmins[vel_indx, thet_indx] = find_scaled_minhmin_S_prime(y)
                    
                    y = trim_ys(y, t)
                    
                    num_oscillations[vel_indx, thet_indx] = count_oscillations(y)
                    
                except Exception:
                    hmins[vel_indx, thet_indx] = 1
                    num_oscillations[vel_indx, thet_indx] = 1
    
    for thet_indx, theta in enumerate(thetas):
        for vel_indx, vel in enumerate(vels):
            if events[vel_indx, thet_indx] == 0:
                
                if num_oscillations[vel_indx, thet_indx] == 0:
                    categories[vel_indx, thet_indx]= 2
                else:
                    categories[vel_indx, thet_indx]= 3
            elif events[vel_indx, thet_indx] == 3:
                if num_oscillations[vel_indx, thet_indx] == 0:
                    categories[vel_indx, thet_indx]= 0
                else:
                    categories[vel_indx, thet_indx]= 1
                
    np.save(pathname1 + '_hmins', hmins)
    np.save(pathname1 + '_categories', categories)
    np.save(pathname1 + '_events', events)
    np.save(pathname1 + '_theta_collision', theta_collision)

 