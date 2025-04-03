#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 17:39:26 2024

@author: cotton
"""

import numpy as np

def func_get_all_maxes(x, y2_deriv, y2_deriv_p_V, y2_deriv_p_Om):
    # x max is how far along the maximum is
    
    
    max_deriv = np.max(np.abs(y2_deriv))
    max_loc = np.argmax(np.abs(y2_deriv))
    xmax = x[max_loc] - x[0]
    max_pV = y2_deriv_p_V[max_loc]
    max_pOmega = y2_deriv_p_Om[max_loc]
    
    return max_deriv, xmax, max_pV, max_pOmega