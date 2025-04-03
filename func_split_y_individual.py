#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 16:12:03 2024

@author: cotton
"""

def func_split_y(y, indx):
    # splits y at a certain time indx
    x_a = y[0, indx]
    x_b = y[1, indx]
    y_a = y[2, indx]
    y_b = y[3, indx]
    theta_a = y[4, indx]
    theta_b = y[5, indx]
    u_a = y[6, indx]
    u_b = y[7, indx]
    v_a = y[8, indx]
    v_b = y[9, indx]
    omega_a = y[10, indx]
    omega_b = y[11, indx]
    
    return x_a, x_b, y_a, y_b, theta_a, theta_b, u_a, u_b, v_a, v_b, omega_a, omega_b