#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 09:48:32 2024

@author: cotton
"""

def weighted_val(theta_min, theta_max, theta, LHS, RHS):
    # find the factor multiplying both of the outputs
    # all thetas given in radians
    if theta <= theta_min:
        combine_factor = 0
    elif theta >= theta_max:
        combine_factor = 1
    else:
        combine_factor = (theta - theta_min)**3/ ( (theta - theta_min)**3 + (theta_max - theta )**3);

    weighted_val = RHS*combine_factor + (1-combine_factor)*LHS;
    
    return weighted_val