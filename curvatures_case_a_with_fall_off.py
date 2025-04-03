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

# check the collision 0 one

# first load in the data by defining the hs and thetas
gamma = 0.5

mask_curvatures = True

x = np.linspace(0, 1, 1000)

pathname = '/Users/cotton/Documents/DPhil/Polynas/Code/non_eqm/Rotation/Python 1D case/Data/gamma_' + str(gamma) + '/'
#pathname = '/Users/cotton/Documents/DPhil reading/Polynas/Code/non_eqm/Rotation/Python 1D case/Data/gamma_' + str(gamma) + '_no_move_away/'
thetas_data = np.load(pathname + 'thetas_gamma_' + str(gamma) + '.npy')
vels = np.load(pathname + 'vels_gamma_' + str(gamma) + '.npy')

#thetas_data = np.array([thetas_data[10], thetas_data[10]])
#vels = np.array([vels[46], vels[46]])

#thetas_data = [0];

pVcurvatures = np.zeros((len(vels), len(thetas_data)))-1
pOmegacurvatures = np.zeros((len(vels), len(thetas_data)))-1
xcurvatures = np.zeros((len(vels), len(thetas_data)))-1
curvatures = np.zeros((len(vels), len(thetas_data)))-1
theta_max = np.zeros((len(vels), len(thetas_data)))-1
theta_abs_max = np.zeros((len(vels), len(thetas_data)))-1
hmin_max = np.zeros((len(vels), len(thetas_data)))-1
Omega_max = np.zeros((len(vels), len(thetas_data)))-1
V_max = np.zeros((len(vels), len(thetas_data)))-1
Vz_max = np.zeros((len(vels), len(thetas_data)))-1
Vx_max = np.zeros((len(vels), len(thetas_data)))-1
Vrot_max = np.zeros((len(vels), len(thetas_data)))-1
Cases = np.zeros((len(vels), len(thetas_data)))-1

theta_min_data = np.zeros((len(vels), len(thetas_data)))-1
h_min_data = np.zeros((len(vels), len(thetas_data)))-1
V_max_data = np.zeros((len(vels), len(thetas_data)))-1
Omega_max_data = np.zeros((len(vels), len(thetas_data)))-1

pvzcurvatures = np.zeros((len(vels), len(thetas_data)))-1
pvxcurvatures = np.zeros((len(vels), len(thetas_data)))-1
pomegacurvatures = np.zeros((len(vels), len(thetas_data)))-1


t_max = np.zeros((len(vels), len(thetas_data)))-1
t_max_type = np.zeros((len(vels), len(thetas_data)))-1

try:
    Cases = np.load(pathname + 'Cases_gamma_' + str(gamma) + 't.npy')
    pVcurvatures = np.load(pathname + 'pVcurvatures_gamma_' + str(gamma) + '.npy')
    pOmegacurvatures = np.load(pathname + 'pOmegacurvatures_gamma_' + str(gamma) + '.npy')
    xcurvatures = np.load(pathname + 'xcurvatures_gamma_' + str(gamma) + '.npy')
    curvatures = np.load(pathname + 'curvatures_gamma_' + str(gamma) + '.npy')
    theta_max = np.load(pathname + 'theta_max_gamma_' + str(gamma) + '.npy')
    theta_abs_max = np.load(pathname + 'theta_abs_max_gamma_' + str(gamma) + '.npy')
    hmin_max = np.load(pathname + 'hmin_max_gamma_' + str(gamma) + '.npy')
    Omega_max = np.load(pathname + 'Omega_max_gamma_' + str(gamma) + '.npy')
    V_max = np.load(pathname + 'V_max_gamma_' + str(gamma) + '.npy')
    Vz_max = np.load(pathname + 'Vz_max_gamma_' + str(gamma) + '.npy')
    Vx_max = np.load(pathname + 'Vx_max_gamma_' + str(gamma) + '.npy')
    Vrot_max = np.load(pathname + 'Vrot_max_gamma_' + str(gamma) + '.npy')

    theta_min_data = np.load(pathname + 'theta_min_data_gamma_' + str(gamma) + '.npy')
    h_min_data = np.load(pathname + 'h_min_data_gamma_' + str(gamma) + '.npy')
    V_max_data = np.load(pathname + 'V_max_data_gamma_' + str(gamma) + '.npy')
    Omega_max_data = np.load(pathname + 'Omega_max_data_gamma_' + str(gamma) + '.npy')
    
    pvzcurvatures = np.load(pathname + 'pvzcurvatures_gamma_' + str(gamma) + '.npy')
    pvxcurvatures = np.load(pathname + 'pvxcurvatures_gamma_' + str(gamma) + '.npy')
    pomegacurvatures = np.load(pathname + 'pomegacurvatures_gamma_' + str(gamma) + '.npy')
    
    t_max = np.load(pathname + 't_max_gamma_' + str(gamma) + '.npy')
    t_max_type = np.load(pathname + 't_max_type_gamma_' + str(gamma) + '.npy')
    
    
    Cases = np.load(pathname + 'Cases_gamma_' + str(gamma) + '.npy')
    

except Exception:

    for thet_indx, theta_data in enumerate(thetas_data):
        print('theta', thet_indx, theta_data*180/np.pi)
        for vel_indx, vel in enumerate(vels):
            #print('vel', vel_indx, vel, theta_data)
            #try:
                y = np.load(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta_data*180/np.pi, 3)) + '.npy')
                t = np.load(pathname + 'v_' + str(round(vel, 5)) + '_theta_' + str(round(theta_data*180/np.pi, 3)) + '_t.npy')
                
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
                
                pVcurvatures[vel_indx, thet_indx] = individual_pV[np.argmax(np.abs(individual_curvature))]
                pOmegacurvatures[vel_indx, thet_indx] = individual_pOmega[np.argmax(np.abs(individual_curvature))]
                xcurvatures[vel_indx, thet_indx] = individual_x[np.argmax(np.abs(individual_curvature))]
                curvatures[vel_indx, thet_indx] = np.max(np.abs(individual_curvature))
                theta_abs_max[vel_indx, thet_indx] = thetas[np.argmax(np.abs(individual_curvature))]
                theta_max[vel_indx, thet_indx] = thetas_non_abs[np.argmax(np.abs(individual_curvature))]
                hmin_max[vel_indx, thet_indx] = hmins[np.argmax(np.abs(individual_curvature))]
                Omega_max[vel_indx, thet_indx] = Omegas[np.argmax(np.abs(individual_curvature))]
                V_max[vel_indx, thet_indx] = Vs[np.argmax(np.abs(individual_curvature))]
                Vz_max[vel_indx, thet_indx] = Vztrans[np.argmax(np.abs(individual_curvature))]
                Vx_max[vel_indx, thet_indx] = Vxtrans[np.argmax(np.abs(individual_curvature))]
                Vrot_max[vel_indx, thet_indx] = Vrots[np.argmax(np.abs(individual_curvature))]
                Cases[vel_indx, thet_indx] = cases[np.argmax(np.abs(individual_curvature))]
                
                theta_min_data[vel_indx, thet_indx] = thetas[np.argmin(np.abs(thetas))]
                h_min_data[vel_indx, thet_indx] = hmins[np.argmin(np.abs(hmins))]
                V_max_data[vel_indx, thet_indx] = np.max(np.abs(Vs))
                Omega_max_data[vel_indx, thet_indx] = np.max(np.abs(Omegas))
                
                # also find the curvatures about the centre of mass
                pvzcurvatures[vel_indx, thet_indx] = individual_pV[np.argmax(np.abs(individual_curvature))]/Vs[np.argmax(np.abs(individual_curvature))]*Vztrans[np.argmax(np.abs(individual_curvature))]
                pvxcurvatures[vel_indx, thet_indx] = individual_pV[np.argmax(np.abs(individual_curvature))]/Vs[np.argmax(np.abs(individual_curvature))]*Vxtrans[np.argmax(np.abs(individual_curvature))]
                pomegacurvatures[vel_indx, thet_indx] = individual_pOmega[np.argmax(np.abs(individual_curvature))] + individual_pV[np.argmax(np.abs(individual_curvature))]/Vs[np.argmax(np.abs(individual_curvature))]*Vrots[np.argmax(np.abs(individual_curvature))]
                
                # also find whether the max curvature occured at the end point of the simulation
                t_max[vel_indx, thet_indx] = t[np.argmax(np.abs(individual_curvature))]
                t_maximum = t[np.argmax(np.abs(individual_curvature))] - np.max(t)
                if t_maximum == 0:
                    t_max_type[vel_indx, thet_indx] = 1
                else:
                    t_max_type[vel_indx, thet_indx] = 0
                    
                
                
                #print(curvatures[vel_indx, thet_indx])
            #except Exception:
                #curvatures[vel_indx, thet_indx] = 0
    
    np.save(pathname + 'pVcurvatures_gamma_' + str(gamma) + '.npy', pVcurvatures)
    np.save(pathname + 'pOmegacurvatures_gamma_' + str(gamma) + '.npy', pOmegacurvatures)
    np.save(pathname + 'xcurvatures_gamma_' + str(gamma) + '.npy', xcurvatures)
    np.save(pathname + 'curvatures_gamma_' + str(gamma) + '.npy', curvatures)
    np.save(pathname + 'theta_max_gamma_' + str(gamma) + '.npy', theta_max)
    np.save(pathname + 'theta_abs_max_gamma_' + str(gamma) + '.npy', theta_abs_max)
    np.save(pathname + 'hmin_max_gamma_' + str(gamma) + '.npy', hmin_max)
    np.save(pathname + 'Omega_max_gamma_' + str(gamma) + '.npy', Omega_max)
    np.save(pathname + 'V_max_gamma_' + str(gamma) + '.npy', V_max)
    
    np.save(pathname + 'Vx_max_gamma_' + str(gamma) + '.npy', Vx_max)
    np.save(pathname + 'Vz_max_gamma_' + str(gamma) + '.npy', Vz_max)
    np.save(pathname + 'Vrot_max_gamma_' + str(gamma) + '.npy', Vrot_max)
    
    np.save(pathname + 'theta_min_data_gamma_' + str(gamma) + '.npy', theta_min_data)
    np.save(pathname + 'h_min_data_gamma_' + str(gamma) + '.npy', h_min_data)
    np.save(pathname + 'V_max_data_gamma_' + str(gamma) + '.npy', V_max_data)
    np.save(pathname + 'Omega_max_data_gamma_' + str(gamma) + '.npy', Omega_max_data)

    np.save(pathname + 'pvzcurvatures_gamma_' + str(gamma) + '.npy', pvzcurvatures)
    np.save(pathname + 'pvxcurvatures_gamma_' + str(gamma) + '.npy', pvxcurvatures)
    np.save(pathname + 'pomegacurvatures_gamma_' + str(gamma) + '.npy', pomegacurvatures)

    np.save(pathname + 't_max_gamma_' + str(gamma) + '.npy', t_max)
    np.save(pathname + 't_max_type_gamma_' + str(gamma) + '.npy', t_max_type)
    np.save(pathname + 'Cases_gamma_' + str(gamma) + '.npy', Cases)
    
theta_original = np.zeros((len(vels), len(thetas_data)))-1
for thet_indx, theta_data in enumerate(thetas_data):
    theta_original[:, thet_indx] = theta_data
np.save(pathname + 'theta_original_gamma_' + str(gamma) + '.npy', theta_original)

