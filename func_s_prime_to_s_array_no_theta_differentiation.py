#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 16:28:37 2024

@author: cotton
"""

import numpy as np

def s_prime_to_s_array_no_theta_differentiation(x_a, x_b, y_a, y_b, theta_a, theta_b, u_a, u_b, v_a, v_b, omega_a, omega_b):
    # LHS of rod A is always at theta = 0
    
    theta, hinc = calc_acute_theta_h_inc(theta_a, theta_b)
    
    z_aflat = y_a*np.cos(theta_b) + x_a*np.sin(theta_b)
    z_bflat = y_b*np.cos(theta_b) + x_b*np.sin(theta_b)
    
    hmin = z_aflat - theta/2 - z_bflat
    
    x_aflat = x_a*np.cos(theta_b) - y_a*np.sin(theta_b)
    x_bflat = x_b*np.cos(theta_b) - y_b*np.sin(theta_b)
    
    x_as = 1/2 #np.cos(theta)/2
    x_bs = x_bflat - x_aflat + 1/2 #np.cos(theta)/2 #1/2 #np.cos(theta_b)/2
    
    z_as = hmin + theta/2
    z_bs = 0
    
    theta_as = theta_a - theta_b
    theta_bs = 0

    
    u_as = u_a*np.cos(theta_b) - v_a*np.sin(theta_b)
    u_bs = u_b*np.cos(theta_b) - v_b*np.sin(theta_b)
    v_as = v_a*np.cos(theta_b) + u_a*np.sin(theta_b)
    v_bs = v_b*np.cos(theta_b) + u_b*np.sin(theta_b)
    
    omega_as = omega_a
    omega_bs = omega_b
    
    return x_as, x_bs, z_as, z_bs, theta_as, theta_bs, u_as, u_bs, v_as, v_bs, omega_as, omega_bs, hinc

def calc_acute_theta_h_inc(theta_a, theta_b):
    if theta_a - theta_b > 0: # clockwise is positive to h will be decreasing
        hinc = False
    else:
        hinc = True
    theta = abs(theta_a-theta_b)
    return theta, hinc
