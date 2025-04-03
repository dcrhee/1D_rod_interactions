#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 16:22:55 2024

@author: cotton
"""


def calc_V_Omega_individual(hinc, theta, x_a, x_b, u_a, u_b, v_a, v_b, omega_a, omega_b, z_a):
    # calcluates the scaled values of U and Omega
    
    
    if hinc: # if the height is increasing with x
        V = (v_a - v_b + omega_a*x_a - omega_b*x_b - theta/2*(u_a - u_b - omega_a*z_a))
        Vztrans = v_a - v_b
        Vrot =  + omega_a*x_a - omega_b*x_b 
        Vxtrans =  - theta/2*(u_a - u_b - omega_a*z_a)
    else:
        V = (v_a - v_b + omega_a*x_a - omega_b*x_b + theta/2*(u_a - u_b- omega_a*z_a))
        Vztrans = v_a - v_b
        Vrot =  + omega_a*x_a - omega_b*x_b 
        Vxtrans =  + theta/2*(u_a - u_b - omega_a*z_a)

    Omega = (omega_a - omega_b)
    
    
    return V, Omega, Vztrans, Vrot, Vxtrans