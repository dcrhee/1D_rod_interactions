#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 15:36:02 2024

@author: cotton
"""

def func_get_case(x_as, x_bs, theta, gamma):
    # we have to either have case 0 or 2
    if x_as - 1/2 < x_bs - 1/(2*gamma):
        if theta < 0:
            return 1
        else:
            return 3
    elif x_as + 1/2 > x_bs + 1/(2*gamma):
        if theta < 0:
            return 2
        else:
            return 4
    else:
        return 0