#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 10:54:35 2024

@author: cotton
"""

import numpy as np

def func_get_theoretical_collision_boundary(gamma, theta, hin):
    # inputs: gamma = La/Lb, theta = initial angles in radians, hin = initial minimum distance
    v_boundary = (1+gamma)*(1-np.log(2))/theta**2# - (1+gamma)/theta**3*(-hin*np.log(1 + theta/hin) + theta*(np.log(2 + theta/hin) - np.log(1 + theta/hin)) )
    #v_boundary = -(1+gamma)/theta**3*(-hin*np.log(1 + theta/hin) + theta*(np.log(2 + theta/hin) - np.log(1 + theta/hin)) )
    theta = theta[5]
    print((1+gamma)*np.log(2)/theta**2)
    print(- (1+gamma)/theta**3*(-hin*np.log(1 + theta/hin) + theta*(np.log(2 + theta/hin) - np.log(1 + theta/hin)) ))
    print(hin, theta)
    #v_boundary = - (1+gamma)/theta**3*(-hin*np.log(1 + theta/hin) + theta*(np.log(2 + theta/hin) - np.log(1 + theta/hin)) )
    return v_boundary