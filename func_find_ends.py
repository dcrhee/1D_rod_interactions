#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 15:09:47 2024

@author: cotton
"""

import numpy as np
import sys

sys.path.append('/Users/cotton/Documents/DPhil reading/Polynas/Code/non_eqm/Rotation/Python 1D case/functions/')
from func_get_case import func_get_case

def find_ends(y, gamma, COM = False):
    v_a = y[8,:]
    u_a = y[6,:]
    x_a = y[0, :]
    x_b = y[1, :]
    y_a = y[2, :]
    y_b = y[3, :]
    theta_a = y[4, :]
    theta_b = y[5,:]
    v_b = y[9,:]
    u_b = y[7,:]
    
    
    x_top = np.zeros((2, len(theta_a)))
    y_top = np.zeros((2, len(theta_a)))
    
    x_top[0, :] = x_a - 1/2*np.cos(theta_a)
    x_top[1, :] = x_a + 1/2*np.cos(theta_a)
    y_top[0, :] = y_a + 1/2*np.sin(theta_a)
    y_top[1, :] = y_a - 1/2*np.sin(theta_a)
    
    x_bottom = np.zeros((2, len(theta_a)))
    y_bottom = np.zeros((2, len(theta_a)))
    
    x_bottom[0, :] = x_b - 1/gamma*1/2*np.cos(theta_b)
    x_bottom[1, :] = x_b + 1/gamma*1/2*np.cos(theta_b)
    y_bottom[0, :] = y_b + 1/gamma*1/2*np.sin(theta_b)
    y_bottom[1, :] = y_b - 1/gamma*1/2*np.sin(theta_b)  
    
    if COM: # if instead give the values at the centre of mass
        xCOM_LHS = (x_top[0, :] + x_bottom[0, :])/2
        xCOM_RHS = (x_top[1, :] + x_bottom[1, :])/2
        yCOM_LHS = (y_top[0, :] + y_bottom[0, :])/2
        yCOM_RHS = (y_top[1, :] + y_bottom[1, :])/2
        yCOM = (y_a + y_b)/2
    
        #x_top[0, :] = x_top[0, :] - xCOM_LHS
        #x_top[1, :] = x_top[1, :] - xCOM_RHS
        y_top[0, :] = y_top[0, :] - yCOM
        y_top[1, :] = y_top[1, :] - yCOM
        
        #x_bottom[0, :] = x_bottom[0, :] - xCOM_LHS
        #x_bottom[1, :] = x_bottom[1, :] - xCOM_RHS
        y_bottom[0, :] = y_bottom[0, :] - yCOM
        y_bottom[1, :] = y_bottom[1, :] - yCOM
        
    
    return x_top, y_top, x_bottom, y_bottom

def find_ends_b_flat(y, gamma, COM = False):
    v_a = y[8,:]
    u_a = y[6,:]
    x_a = y[0, :]
    x_b = y[1, :]
    y_a = y[2, :]
    y_b = y[3, :]
    theta_a = y[4, :]
    theta_b = y[5,:]
    v_b = y[9,:]
    u_b = y[7,:]
    
    x_aflat = x_a*np.cos(theta_b) - y_a*np.sin(theta_b)
    x_bflat = x_b*np.cos(theta_b) - y_b*np.sin(theta_b)
    
    z_aflat = y_a*np.cos(theta_b) + x_a*np.sin(theta_b)
    z_bflat = y_b*np.cos(theta_b) + x_b*np.sin(theta_b)
    
    
    x_top = np.zeros((2, len(theta_a)))
    y_top = np.zeros((2, len(theta_a)))
    
    theta = theta_a-theta_b
    
    x_top[0, :] = x_aflat - 1/2*np.cos(theta)
    x_top[1, :] = x_aflat + 1/2*np.cos(theta)
    y_top[0, :] = z_aflat + 1/2*np.sin(theta)
    y_top[1, :] = z_aflat - 1/2*np.sin(theta)
    
    x_bottom = np.zeros((2, len(theta_a)))
    y_bottom = np.zeros((2, len(theta_a)))
    
    x_bottom[0, :] = x_bflat - 1/gamma*1/2
    x_bottom[1, :] = x_bflat + 1/gamma*1/2
    y_bottom[0, :] = z_bflat
    y_bottom[1, :] = z_bflat  
    
    if COM: # if instead give the values at the centre of mass
        xCOM_LHS = (x_top[0, :] + x_bottom[0, :])/2
        xCOM_RHS = (x_top[1, :] + x_bottom[1, :])/2
        yCOM_LHS = (y_top[0, :] + y_bottom[0, :])/2
        yCOM_RHS = (y_top[1, :] + y_bottom[1, :])/2
        yCOM = (y_a + y_b)/2
    
        #x_top[0, :] = x_top[0, :] - xCOM_LHS
        #x_top[1, :] = x_top[1, :] - xCOM_RHS
        y_top[0, :] = y_top[0, :] - yCOM
        y_top[1, :] = y_top[1, :] - yCOM
        
        #x_bottom[0, :] = x_bottom[0, :] - xCOM_LHS
        #x_bottom[1, :] = x_bottom[1, :] - xCOM_RHS
        y_bottom[0, :] = y_bottom[0, :] - yCOM
        y_bottom[1, :] = y_bottom[1, :] - yCOM
        
        
        x_top[0, :] = x_top[0, :] - x_bottom[0, :]
        x_top[1, :] = x_top[1, :] - x_bottom[0, :]
        
        x_bottom[1, :] = x_bottom[1, :] - x_bottom[0, :]
        x_bottom[0, :] = 0
        
    
    return x_top, y_top, x_bottom, y_bottom

def calc_acute_theta_h_inc(theta_a, theta_b):
    if theta_a - theta_b > 0: # clockwise is positive to h will be decreasing
        hinc = False
    else:
        hinc = True
    theta = abs(theta_a-theta_b)
    return theta, hinc

def s_prime_to_s_array(x_a, x_b, y_a, y_b, theta_a, theta_b, u_a, u_b, v_a, v_b, omega_a, omega_b):
    
    theta, hinc = calc_acute_theta_h_inc(theta_a, theta_b)
    
    
    z_aflat = y_a*np.cos(theta_b) + x_a*np.sin(theta_b)
    z_bflat = y_b*np.cos(theta_b) + x_b*np.sin(theta_b)
    
    hmin = z_aflat - theta/2 - z_bflat
    
    x_aflat = x_a*np.cos(theta_b) - y_a*np.sin(theta_b)
    x_bflat = x_b*np.cos(theta_b) - y_b*np.sin(theta_b)
    
    if hinc:
        x_as = 1/2 #np.cos(theta)/2
        x_bs = x_bflat - x_aflat + 1/2 #np.cos(theta)/2 #1/2 #np.cos(theta_b)/2
    else:
        x_as = -1/2#-np.cos(theta)/2
        x_bs = x_bflat - x_aflat - 1/2 #np.cos(theta)/2 #1/2 #np.cos(theta_b)/2    
    
    z_as = hmin + theta/2
    z_bs = 0
    
    theta_as = theta_a - theta_b
    theta_bs = 0
    
    
    return x_as, x_bs, z_as, z_bs, theta_as, theta_bs

def find_ends_according_to_small_angle(y, gamma, COM = False):
    hmins = np.zeros(len(y[5,:]))
    hminsactual = np.zeros(len(y[5,:]))
    x_top  = np.zeros((2, len(y[5,:])))
    y_top  = np.zeros((2, len(y[5,:])))
    
    x_bottom  = np.zeros((2, len(y[5,:])))
    y_bottom  = np.zeros((2, len(y[5,:])))
    
    xedge  = np.zeros(len(y[5,:]))
    xedgeacc  = np.zeros(len(y[5,:]))
    
    for i in range(len(y[5,:])):
        x_as, x_bs, z_as, z_bs, theta_as, theta_bs = s_prime_to_s_array(y[0, i], y[1, i], y[2, i], y[3, i], y[4, i], y[5, i], y[6, i], y[7, i], y[8, i], y[9, i], y[10, i], y[11, i])

    
        theta = theta_as - theta_bs
        x1 = max(x_as - 1/2, x_bs - 1/(2*gamma)) # 
        x2 = min(x_as + 1/2, x_bs + 1/(2*gamma))
        
        x1acc = max(x_as - 1/2*np.cos(theta), x_bs - 1/(2*gamma)) # 
        x2acc = min(x_as + 1/2*np.cos(theta), x_bs + 1/(2*gamma))
        
        x_top[0, i] = x_as - 1/2
        x_top[1, i] = x_as + 1/2
        
        x_bottom[0, i] = x_bs - 1/(2*gamma)
        x_bottom[1, i] = x_bs + 1/(2*gamma)
        
        y_top[0, i] = z_as + 1/2*theta_as
        y_top[1, i] = z_as - 1/2*theta_as
        
        if theta < 0:
            hmins[i] = abs(z_as) + theta*(x_as - x1)
            hminsactual[i] = abs(z_as) + np.tan(theta)*(x_as - x1acc)
            xedge[i] = x1
            xedgeacc[i] = x1acc
        else:
            
            hmins[i] = abs(z_as) + theta*(x_as - x2)
            hminsactual[i] = abs(z_as) + np.tan(theta)*(x_as - x2acc)
            xedge[i] = x2
            xedgeacc[i] = x2acc
            
    if COM: # if instead give the values at the centre of mass
         xCOM_LHS = (x_top[0, :] + x_bottom[0, :])/2
         xCOM_RHS = (x_top[1, :] + x_bottom[1, :])/2
         yCOM_LHS = (y_top[0, :] + y_bottom[0, :])/2
         yCOM_RHS = (y_top[1, :] + y_bottom[1, :])/2
         yCOM = (z_as + z_bs)/2
     
         #x_top[0, :] = x_top[0, :] - xCOM_LHS
         #x_top[1, :] = x_top[1, :] - xCOM_RHS
         y_top[0, :] = y_top[0, :] - yCOM
         y_top[1, :] = y_top[1, :] - yCOM
         
         #x_bottom[0, :] = x_bottom[0, :] - xCOM_LHS
         #x_bottom[1, :] = x_bottom[1, :] - xCOM_RHS
         y_bottom[0, :] = y_bottom[0, :] - yCOM
         y_bottom[1, :] = y_bottom[1, :] - yCOM
         
         xedge = xedge - x_bottom[0, :]
         xedgeacc = xedgeacc - x_bottom[0, :]
         
         x_top[0, :] = x_top[0, :] - x_bottom[0, :]
         x_top[1, :] = x_top[1, :] - x_bottom[0, :]
         
         x_bottom[1, :] = x_bottom[1, :] - x_bottom[0, :]
         x_bottom[0, :] = 0
        
    return x_top, y_top, x_bottom, y_bottom, hmins, xedge, hminsactual, xedgeacc

def get_closest_time(t, num_times):
    # inputs:
        # t = times of the simulation
        # num_times = the number of rows to split the times into
    
    t_indxs = np.zeros(num_times) # stores the indices of the closest match
    calc_times = np.linspace(0, np.max(t), num_times)
    for indx, time in enumerate(calc_times):
        t_indxs[indx] = np.argmin(abs(time - t))    
    return t_indxs


def find_hmin_vrel_theta_rel(y, gamma):
    # find the relative motion as a grid of all the variables
    
    x_a = y[0,:]
    x_b = y[1,:]
    y_a = y[2,:]
    y_b = y[3,:]
    
    v_a = y[8,:]
    u_a = y[6,:]
    
    theta_b = y[5,:]
    theta_a =  y[4,:]
    
    v_b = y[9,:]
    u_b = y[7,:]
    
    omega_a = y[10,:]
    omega_b = y[11,:]
    
    # convert to frame S
    theta = theta_a-theta_b
    
    x_aflat = x_a*np.cos(theta_b) - y_a*np.sin(theta_b)
    x_bflat = x_b*np.cos(theta_b) - y_b*np.sin(theta_b)
    
    z_aflat = y_a*np.cos(theta_b) + x_a*np.sin(theta_b)
    z_bflat = y_b*np.cos(theta_b) + x_b*np.sin(theta_b)
    
    hmin = z_aflat - abs(theta)/2 - z_bflat
    xrel = x_aflat - x_bflat
  
    u_as = u_a*np.cos(theta_b) - v_a*np.sin(theta_b)
    u_bs = u_b*np.cos(theta_b) - v_b*np.sin(theta_b)
    urel = u_as - u_bs
    
    v_as = v_a*np.cos(theta_b) + u_a*np.sin(theta_b)
    v_bs = v_b*np.cos(theta_b) + u_b*np.sin(theta_b)
    vrel = v_as - v_bs
    
    theta_diff_degrees = theta*180/np.pi
    
    omega_diff = omega_a- omega_b
    
    move_away_LHS =  v_a - v_b + (omega_a - omega_b)/2 - omega_b*(x_a - x_b)
    move_away_RHS =  v_a - v_b - (omega_a - omega_b)/2 - omega_b*(x_a - x_b)
    
    
    
    # check if the rods completely overlap or not    
    for i in range(len(x_a)):
        z_as = z_aflat[i]- z_bflat[i]
        if theta[i] < 0:
            x_as = 1/2 #np.cos(theta)/2
            x_bs = x_bflat[i] - x_aflat[i] + 1/2 #np.cos(theta)/2 #1/2 #np.cos(theta_b)/2
        else:
            x_as = -1/2#-np.cos(theta)/2
            x_bs = x_bflat[i] - x_aflat[i] - 1/2 #np.cos(theta)/2 #1/2 #np.cos(theta_b)/2    
        
        x1 = max(x_as - 1/2, x_bs - 1/(2*gamma)) # 
        x2 = min(x_as + 1/2, x_bs + 1/(2*gamma))
    

        
        if theta[i] < 0:
            hmin[i] = abs(z_as) + theta[i]*(x_as - x1)
            #print(hmin[i], x_as - x1)
        else:
            
            hmin[i] = abs(z_as) + theta[i]*(x_as - x2)
            #print(hmin[i], x_as - x2)
            
        # determine which case we are
        case_num = func_get_case(x_as, x_bs, theta[i], gamma) # determine if the disks are fully overlapping or not
        # calculate x1 and P
        x1 = x_bs - 1/(2*gamma) - (x_as - 1/2) #
        P = x_bs + 1/(2*gamma) - (x_as - 1/2)
        
        
        
        if case_num == 1:
            h = z_as - abs(theta[i])*(1/2 - x1)
        elif case_num == 4:
            h = z_as - abs(theta[i])*(P-1/2)
        else:
            h = z_as - abs(theta[i]/2)
        #print(case_num, h)
        
        #print(case_num, h-hmin[i])
        
        
    
    return hmin, theta_diff_degrees, move_away_LHS, move_away_RHS