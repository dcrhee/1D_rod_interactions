#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 16:01:55 2024

@author: cotton
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import sys

sys.path.append('/Users/cotton/Documents/DPhil/Polynas/Code/non_eqm/Rotation/Python 1D case/functions/')

from func_split_y_individual import func_split_y
from func_s_prime_to_s_array import s_prime_to_s_array
from func_get_case import func_get_case
from func_calc_V_Omega_no_fall_off import calc_V_Omega_individual
from func_s_prime_to_s_array_no_theta_differentiation import s_prime_to_s_array_no_theta_differentiation
from func_sim_curvature_a import func_get_all_maxes_a
from func_get_curvature_case1 import func_get_all_curvatures_a_case1
from func_get_curvature_case2 import func_get_all_curvatures_a_case2
from func_get_curvature_case3 import func_get_all_curvatures_a_case3
from func_get_curvature_case4 import func_get_all_curvatures_a_case4
from func_get_curvature_case0 import func_get_all_curvatures_a_case0
from func_get_all_maxes import func_get_all_maxes

def get_max_curvature_case_a_time(y, t, gamma):
                x = np.linspace(0, 1, 1000)


                # get the relavent components
                
                
                thetas_non_abs = y[4,:] - y[5, :]
                
                individual_curvature = np.zeros(len(t))
                individual_x = np.zeros(len(t))
                individual_pV = np.zeros(len(t))
                individual_pOmega = np.zeros(len(t))
                Vs = np.zeros(len(t))
                Omegas = np.zeros(len(t))
                hmins = np.zeros(len(t))
                thetas = np.zeros(len(t))
                hincs = np.zeros(len(t))
                Vztrans = np.zeros(len(t))
                Vrots = np.zeros(len(t))
                Vxtrans = np.zeros(len(t))
                cases = np.zeros(len(t))
                
                for indx in range(len(t)):
                    x_a, x_b, y_a, y_b, theta_a, theta_b, u_a, u_b, v_a, v_b, omega_a, omega_b = func_split_y(y, indx)
                    
                    # convert to frame where rod A is between 0 and 1 and rod B is at theta = 0
                    x_as, x_bs, z_as, z_bs, theta_as, theta_bs, u_as, u_bs, v_as, v_bs, omega_as, omega_bs, hinc =  s_prime_to_s_array_no_theta_differentiation(x_a, x_b, y_a, y_b, theta_a, theta_b, u_a, u_b, v_a, v_b, omega_a, omega_b)
                    
                    theta = np.abs(theta_as)
                    h = z_as-theta/2 # hmin
                    
                    # determine which case we are
                    case_num = func_get_case(x_as, x_bs, theta_as, gamma) # determine if the disks are fully overlapping or not
                    
                    # calculate Omega and V
                    V, Omega, Vztran, Vrot, Vxtran = calc_V_Omega_individual(hinc, theta_as, x_as, x_bs, u_as, u_bs, v_as, v_bs, omega_as, omega_bs, z_as)# get the V components

                    # calculate x1 and P
                    x1 = x_bs - 1/(2*gamma) - (x_as - 1/2) #
                    P = x_bs + 1/(2*gamma) - (x_as - 1/2)
                    
                    
                    if case_num == 0:
                        # get h and theta (check this one)
                        y2_deriv, y2_deriv_p_V, y2_deriv_p_Om = func_get_all_curvatures_a_case0(x, h, theta, V, Omega, hinc)
                    elif case_num == 1:
                        h = z_as - theta*(1/2 - x1)
                        y2_deriv, y2_deriv_p_V, y2_deriv_p_Om = func_get_all_curvatures_a_case1(x, h, theta, V, Omega, x1)
                    elif case_num == 2:
                        y2_deriv, y2_deriv_p_V, y2_deriv_p_Om = func_get_all_curvatures_a_case2(x, h, theta, V, Omega, P) # h stays the same
                    elif case_num == 3:
                       y2_deriv, y2_deriv_p_V, y2_deriv_p_Om = func_get_all_curvatures_a_case3(x, h, theta, V, Omega, x1)
                    else:
                        h = z_as - theta*(P-1/2)
                        y2_deriv, y2_deriv_p_V, y2_deriv_p_Om = func_get_all_curvatures_a_case4(x, h, theta, V, Omega, P)

                    
                    max_deriv, xmax, max_pV, max_pOmega = func_get_all_maxes(x, y2_deriv, y2_deriv_p_V, y2_deriv_p_Om)
                    
                    Vs[indx] = V
                    Omegas[indx] = Omega
                    hmins[indx] = h
                    thetas[indx] = theta
                    Vztrans[indx] = Vztran
                    Vrots[indx] = Vrot
                    Vxtrans[indx] = Vxtran
                    individual_curvature[indx] = max_deriv
                    individual_x[indx] = xmax
                    individual_pV[indx] = max_pV
                    individual_pOmega[indx] = max_pOmega
                    cases[indx] = case_num
                
                return individual_curvature