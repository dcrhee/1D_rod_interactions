#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 13:48:26 2024

@author: cotton
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

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

def calc_integrals_small_theta(theta, h):
    
    pV = -1/(12*h**3)*(1 - 3*theta/(2*h) + 33/20*(theta/h)**2)
    pomega = 1/(24*h**3)*(1 - 17*theta/(10*h) + 41/20*(theta/h)**2)    
    return pV, pomega

def calc_integrals(theta, h, V, Omega):
    # calculates the scaled integrals
    # inputs: theta = relative angle between the rods, h = minimum distance/length of the top rod
    theta_min = min(h/100, 1e-3)
    theta_max = min(h/10, 1e-2)
    
    if theta == 0:
        pV, pomega = calc_integrals_small_theta(theta, h)
        
    else:
        pVLHS, pomegaLHS = calc_integrals_small_theta(theta, h)
        
        pVRHS = 1/theta**3*(np.log(h/(h+theta)) + 2*theta/(2*h+theta))
        pomegaRHS = 1/(2*theta**3)*(1/(2*h+theta))*(2*h*(3*h/theta+2)*np.log(h/(h+theta)) + 6*h+theta)
        
        pV = weighted_val(theta_min, theta_max, theta, pVLHS, pVRHS)
        pomega = weighted_val(theta_min, theta_max, theta, pomegaLHS, pomegaRHS)
    
    return V*pV, Omega*pomega

# point force derivations
def func_force_V_comp(V, theta, h):
    if theta == 0:
        forceV = V/h**3
    else:
        # the integral of p_V dx
        forceV = V/theta**3*(np.log(h/(h+theta)) + 2*theta/(2*h+theta))
    return forceV

def func_force_Om_comp(omega, theta, h):
    if theta == 0:
        forceOm = omega/(24*h**3)
    else:        
        # the integral of p_Om dx
        forceOm = omega/(2*theta**3)*(1/(2*h+theta))*(2*h*(3*h/theta+2)*np.log(h/(h+theta)) + 6*h+theta)
    return forceOm

def func_x0V(h, theta):
    x0 = h/(2*h+theta)
    return x0

def func_x0Om(h, theta):
    def mysqrt(x): return np.sqrt(x)
    
    if theta == 0:    
        x0 = 1/mysqrt(3)
    else:
        aux0=(-(theta**-2.)*((theta*((2.*h)+(3.*theta)))+(2.*((((h+theta)**2))*((np.log(h))-(np.log((h+theta))))))))
        x0=(((((h**-2.)*((2.*h)+theta))/theta)**-0.5)*(mysqrt(aux0)))/theta
    return x0

def fun_x0VCM(h, theta):
    if theta == 0:
        x0 = 1/2
    else:
        aux0=(0.5*((2.*(h*(((3.*h)+(2.*theta))*(np.log(((h+theta)/h))))))-(theta*((6.*h)+theta))))/((2.*theta)+(((2.*h)+theta)*(\
                                                                                                                                np.log((h/(h+theta))))));
        x0=aux0/theta;
    return x0

def fun_x0OmCM(h, theta):
    if theta == 0:
        x0 = 8/15
    else:
        aux0=(np.log((h/(h+theta))))*((-2.*(theta*((4.*h)+(3.*theta))\
                                            ))+((((h+theta)**2))*(np.log(((h+theta)/h)))));
        x0=((((-7.*(h**2))/theta)+((0.25*theta)+((h**2)*((theta**-\
                                                          3.)*aux0))))-h)/((6.*h)+(theta+(2.*(h*((2.+((3.*h)/theta))*(np.\
                                                                                                                      log((h/(h+theta)))))))));
    return x0

def max_deriv_point_force(F, x0):
    # point force curvaure to the RHS of where the force is applied
    if x0 < 0.17:
        y2max = abs(F)*(1/3-x0)**3/(1/2-x0)**2
        xmax = 1/(6*(1/2-x0))
    elif x0 > 0.83:
        y2max = abs(F)*(x0-2/3)**3/(x0-1/2)**2
        xmax = 1-1/(6*(x0-1/2))
    else:
        y2max = 2*abs(F)*x0**2*(1-x0)**2
        xmax = x0
            
    return y2max, xmax

def func_get_all_maxes_a(x, hmin, theta, V, omega):
    # get the forces
    #forceV = func_force_V_comp(V, theta, hmin)
    #forceOm = func_force_Om_comp(omega, theta, hmin)

    forceV, forceOm = calc_integrals(theta, hmin, V, omega)

    # get the position of the peak
    x0V = func_x0V(hmin, theta)
    x0Om = func_x0Om(hmin, theta)
    
    x0VCM = fun_x0VCM(hmin, theta)
    x0OmCM = fun_x0OmCM(hmin, theta)
    
    
    [y2maxFV, xmaxFV] = max_deriv_point_force(forceV, x0V)
    [y2maxFOm, xmaxFOm] = max_deriv_point_force(forceOm, x0Om)
    
    [y2maxFVCM, xmaxFVCM] = max_deriv_point_force(forceV, x0VCM)
    [y2maxFOmCM, xmaxFOmCM] = max_deriv_point_force(forceOm, x0OmCM)
    
    return y2maxFV, xmaxFV, y2maxFOm, xmaxFOm, y2maxFVCM, xmaxFVCM, y2maxFOmCM, xmaxFOmCM
V = 1
omega = 1
x = np.linspace(0, 1, 1000)
hmin = 0.1
theta = 0
func_get_all_maxes_a(x, hmin, theta, V, omega)
