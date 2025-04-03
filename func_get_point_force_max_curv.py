#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 11:42:00 2024

@author: cotton

gets the points force curvatures

"""

import numpy as np

# point force derivations
def func_force_V_comp(V, theta, h):
    # the integral of p_V dx
    forceV = V/theta**3*(np.log(h/(h+theta)) + 2*theta/(2*h+theta))
    return forceV

def func_force_Om_comp(omega, theta, h):
    # the integral of p_Om dx
    forceOm = omega/(2*theta**3)*(1/(2*h+theta))*(2*h*(3*h/theta+2)*np.log(h/(h+theta)) + 6*h+theta)
    return forceOm

def func_x0V(h, theta):
    x0 = h/(2*h+theta)
    return x0

def func_x0Om(h, theta):
    def mysqrt(x): return np.sqrt(x)
    aux0=(-(theta**-2.)*((theta*((2.*h)+(3.*theta)))+(2.*((((h+theta)**2))*((np.log(h))-(np.log((h+theta))))))))
    x0=(((((h**-2.)*((2.*h)+theta))/theta)**-0.5)*(mysqrt(aux0)))/theta
    return x0

def fun_x0VCM(h, theta):
    aux0=(0.5*((2.*(h*(((3.*h)+(2.*theta))*(np.log(((h+theta)/h))))))-(theta*((6.*h)+theta))))/((2.*theta)+(((2.*h)+theta)*(\
                                                                                                                            np.log((h/(h+theta))))));
    x0=aux0/theta;
    return x0

def fun_x0OmCM(h, theta):
    aux0=(np.log((h/(h+theta))))*((-2.*(theta*((4.*h)+(3.*theta))\
                                        ))+((((h+theta)**2))*(np.log(((h+theta)/h)))));
    x0=((((-7.*(h**2))/theta)+((0.25*theta)+((h**2)*((theta**-\
                                                      3.)*aux0))))-h)/((6.*h)+(theta+(2.*(h*((2.+((3.*h)/theta))*(np.\
                                                                                                                  log((h/(h+theta)))))))));
    return x0

def func_y2_deriv_point_force_case_b(xLHS, xRHS, F, x0, M, P):
    LHS = -((F*(M - xLHS)**2 *(-2* P**2 + M * xLHS - (M + 2 * xLHS) * x0 + P * (xLHS + 3 * x0)))/(M - P)**3)
    RHS = (F * (P - xRHS)**2 *(2* M**2 + 2* xRHS* x0 + P* (-xRHS + x0) - M* (xRHS + 3* x0)))/(M - P)**3
    return LHS, RHS


def get_maxes_b_point_force(x, hmin, theta, V, omega, L, m):
    # inputs
    # x = M --> P, the points along rod B
    # L = length of rod B
    
    # get the forces
    forceV = func_force_V_comp(V, theta, hmin)
    forceOm = func_force_Om_comp(omega, theta, hmin)

    # get the position of the peak
    x0V = func_x0V(hmin, theta)
    x0Om = func_x0Om(hmin, theta)
    
    x0VCM = fun_x0VCM(hmin, theta)
    x0OmCM = fun_x0OmCM(hmin, theta)
    
    # get the delta distribution integrated forces
    xV_LHS = x[x < x0V]
    xV_RHS = x[x >= x0V]
    xOm_LHS = x[x < x0Om]
    xOm_RHS = x[x >= x0Om]
    
    xVCM_LHS = x[x < x0VCM]
    xVCM_RHS = x[x >= x0VCM]
    xOmCM_LHS = x[x < x0OmCM]
    xOmCM_RHS = x[x >= x0OmCM]
    
    y2_deriv_F_V_LHS, y2_deriv_F_V_RHS = func_y2_deriv_point_force_case_b(xV_LHS, xV_RHS, forceV, x0V, x[0], x[-1])    
    y2_deriv_F_V = np.concatenate((y2_deriv_F_V_LHS, y2_deriv_F_V_RHS))
    
    y2_deriv_F_omega_LHS, y2_deriv_F_omega_RHS = func_y2_deriv_point_force_case_b(xOm_LHS, xOm_RHS, forceOm, x0Om, x[0], x[-1])
    y2_deriv_F_omega = np.concatenate((y2_deriv_F_omega_LHS, y2_deriv_F_omega_RHS))
    
    y2maxFV = np.max(np.abs(y2_deriv_F_V))
    xmaxFV = x[np.argmax(np.abs(y2_deriv_F_V))]
    
    y2maxFOm = np.max(np.abs(y2_deriv_F_omega))
    xmaxFOm = x[np.argmax(np.abs(y2_deriv_F_omega))]    
    
    y2_deriv_F_V_LHSCM, y2_deriv_F_V_RHSCM = func_y2_deriv_point_force_case_b(xVCM_LHS, xVCM_RHS, forceV, x0VCM, x[0], x[-1])    
    y2_deriv_F_VCM = np.concatenate((y2_deriv_F_V_LHSCM, y2_deriv_F_V_RHSCM))
    
    y2_deriv_F_omega_LHSCM, y2_deriv_F_omega_RHSCM = func_y2_deriv_point_force_case_b(xOmCM_LHS, xOmCM_RHS, forceOm, x0OmCM, x[0], x[-1])
    y2_deriv_F_omegaCM = np.concatenate((y2_deriv_F_omega_LHSCM, y2_deriv_F_omega_RHSCM))
    
    y2maxFVCM = np.max(np.abs(y2_deriv_F_VCM))
    xmaxFVCM = x[np.argmax(np.abs(y2_deriv_F_VCM))]
    
    y2maxFOmCM = np.max(np.abs(y2_deriv_F_omegaCM))
    xmaxFOmCM = x[np.argmax(np.abs(y2_deriv_F_omegaCM))]    
    
    
    return y2maxFV, xmaxFV, y2maxFOm, xmaxFOm, y2maxFVCM, xmaxFVCM, y2maxFOmCM, xmaxFOmCM


