#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 17:41:15 2024

@author: cotton

Calculates the maximum curvature for each simulation

"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 11:20:10 2024

@author: cotton

Plots the  max curvature for a certain h and theta for the different rods

Find the value of theta, h, V, Omega, vrel that this occurs at
Check if this is the minimum value of theta or h

"""

# have a new euqation for alpha, beta, C and D for the middle derivative that needs to be inputted, the other ones are fine for case B and then test it

import numpy as np

import sys

sys.path.append('/Users/cotton/Documents/DPhil reading/Polynas/Code/non_eqm/Rotation/Python 1D case/functions/')
from func_s_prime_to_s_array import s_prime_to_s_array

# continous pressure derivations
def func_y2_deriv_px_omega_0_case_a(h, theta, x, V):
    # the solution for the top rod assuming the pressure only comes from the V term
    # the end positions of the rod are scaled such that they are at 0 and 1
    aux0=(((h**2)*(3.+(6.*x)))+((x*(theta**2))+(2.*(h*(theta+(3.*(x*\
    theta)))))))*(np.log(h));
    aux1=((((h**2)*(9.+(-6.*x)))+(2.*(h*((5.+(-3.*x))*theta))))-((-2.+\
    x)*(theta**2)))*(np.log((h+theta)));
    aux2=((((-1.+x)**2))*aux0)+(x*(((-1.+x)*(theta*((h*(-3.+(6.*x)))+((\
    -2.+(3.*x))*theta))))+(x*aux1)));
    aux3=aux2-(((3.*(h**2))+((2.*(h*((1.+x)*theta)))+(x*(theta**2)))\
    )*(np.log((h+(x*theta)))));
    output=(V*((theta**-4.)*aux3))/((2.*h)+theta)
    return output

def func_y2_deriv_px_V_0_case_a(h, theta, x, omega):
    # the solution for the top rod assuming the pressure only comes from the Omega term
    # the end positions of the rod are scaled such that they are at 0 and 1

    aux0=(h**2)*((((-1.+x)**2))*((1.+(2.*x))*((((h+theta)**2))*(\
    omega*(((np.log(h))**2))))));
    aux1=(h**2)*((x**2)*((-3.+(2.*x))*((((h+theta)**2))*(omega*(((\
    np.log((h+theta)))**2))))));
    aux2=((h**2)*((-4.+((61.+(-36.*x))*x))*theta))+((-2.*(h*((1.+(x*(-\
    9.+(4.*x))))*(theta**2))))+(x*(theta**3.)));
    aux3=(x*(theta*((((h**3.)*(-2.+(16.*((3.+(-2.*x))*x))))+aux2)*\
    omega)))+(2.*((h**2)*((((h+theta)**2))*(omega*(np.log((h+(x*\
    theta))))))));
    aux4=x*(theta*(((14.*((h**2)*(-1.+(2.*x))))+((h*((-9.+(16.*x))*\
    theta))+(x*(theta**2))))*omega));
    aux5=(3.*((h**2)*((3.+(4.*x))*theta)))+((2.*(h*(x*((3.+x)*(theta\
    **2)))))+((x**2)*(theta**3.)));
    aux6=(2.*((h**2)*(-7.+(8.*(x*(-1.+(2.*x)))))))+((h*((-9.+(x*(-25.+(36.\
    *x))))*theta))+(8.*((-1.+x)*(x*(theta**2)))));
    aux7=(((h+theta)**2))*(((1.+((-6.*(x**2))+(4.*(x**3.))))*(np.log((\
    h+theta))))+(np.log((h+(x*theta)))));
    aux8=(theta*(((-1.+x)*aux4)-(((14.*(h**3.))+aux5)*(omega*(np.\
    log((h+(x*theta))))))))+(h*(omega*((np.log(h))*(((-1.+x)*(\
    theta*aux6))+(-2.*(h*aux7))))));
    output=(0.25*((theta**-6.)*((2.*aux0)+((2.*aux1)+(((np.log((h+\
    theta)))*aux3)+aux8)))))/((2.*h)+theta);   
        
    return output

def func_y2_deriv_px_omega_0_case_a_small_theta(h, theta, x, V):
    zeroth_term = V*(x**2 - 2*x**3 + x**4)/(24*h**3)
    first_term = -V*(3*x**2 - 2*x**3 - 5*x**4 + 4*x**5)/(80*h**4)*theta
    
    aux0=(3.*(x**2))+((2.*(x**3.))+((-5.*(x**4.))+((-8.*(x**5.))+(8.*(\
    x**6.)))));
    second_term=0.00625*((h**-5.)*(V*(aux0*(theta**2))));
    
    aux0=(234.*(x**3.))+((-105.*(x**4.))+((-168.*(x**5.))+((-280.*(x**6.))\
    +(320.*(x**7.)))));
    third_term=-0.00014881*((h**-6.)*(V*((aux0-(x**2))*(theta**3.))));
    
    aux0=(-105.*(x**4.))+((-168.*(x**5.))+((-280.*(x**6.))+((-480.*(x**7.)\
    )+(600.*(x**8.)))));
    fourth_term=0.0000744048*((h**-7.)*(V*(((-201.*(x**2))+((634.*(x**3.))+\
    aux0))*(theta**4.))));

    aux0=(-168.*(x**5.))+((-280.*(x**6.))+((-480.*(x**7.))+((-840.*(x**8.)\
    )+(1120.*(x**9.)))));
    aux1=(h**-8.)*(V*(((-681.*(x**2))+((1434.*(x**3.))+((-105.*(x**4.))+\
    aux0)))*(theta**5.)));
    fifth_term=0.0000372024*aux1;
    
    output = zeroth_term + first_term + second_term + third_term + fourth_term + fifth_term
    
    return output

def func_y2_deriv_px_V_0_case_a_small_theta(h, theta, x, omega):
    zeroth_term = omega*-((-1 + x)**2*x**2*(2+x))/(120*h**3)
    first_term = omega*(4*x**2 - 2*x**3 - 5*x**4 + 3*x**6)/(240*h**4)*theta
    
    aux0=(-66.*((x**3.)*omega))+((105.*((x**4.)*omega))+((112.*((\
    x**5.)*omega))+(-96.*((x**7.)*omega))));
    second_term=0.00014881*((h**-5.)*((theta**2)*((-55.*((x**2)*omega))+\
    aux0)))
        
    
    aux0=(-399.*((x**4.)*omega))+((-504.*((x**5.)*omega))+((-560.*(\
    (x**6.)*omega))+(600.*((x**8.)*omega))));
    third_term=0.0000248016*((h**-6.)*((theta**3.)*((-87.*((x**2)*omega)\
    )+((950.*((x**3.)*omega))+aux0))))

    aux0=(1064.*((x**5.)*omega))+((1400.*((x**6.)*omega))+((1600.*((\
    x**7.)*omega))+(-2000.*((x**9.)*omega))));
    aux1=(theta**4.)*((1533.*((x**2)*omega))+((-4402.*((x**3.)*\
    omega))+((805.*((x**4.)*omega))+aux0)));
    fourth_term = 1/134400*((h**-7.)*aux1)
    
    aux0=(-5320.*((x**6.)*omega))+((-7200.*((x**7.)*omega))+((-8400.\
    *((x**8.)*omega))+(11760.*((x**10.)*omega))));
    aux1=(30982.*((x**3.)*omega))+((-2895.*((x**4.)*omega))+((-3864.\
    *((x**5.)*omega))+aux0))
    fifth_term = 1/806400*((h**-8.)*((theta**5.)*((-15063.*((x**2)*omega))+\
    aux1)))

    output = zeroth_term + first_term + second_term + third_term + fourth_term + fifth_term
    
    return output

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


def func_get_curvature_case_a(x, h, theta, V, omega, hinc):
    # transition function locations
    theta_min = h/10
    theta_max = h/2
    
    if theta == 0:

        y2_deriv_p_V = func_y2_deriv_px_omega_0_case_a_small_theta(h, theta, x, V)
        y2_deriv_p_omega = func_y2_deriv_px_V_0_case_a_small_theta(h, theta, x, omega)
        
    else:
        
        y2_deriv_p_V_LHS = func_y2_deriv_px_omega_0_case_a_small_theta(h, theta, x, V)
        y2_deriv_p_omega_LHS = func_y2_deriv_px_V_0_case_a_small_theta(h, theta, x, omega)
        
        y2_deriv_p_V_RHS = func_y2_deriv_px_omega_0_case_a(h, theta, x, V)
        y2_deriv_p_omega_RHS = func_y2_deriv_px_V_0_case_a(h, theta, x, omega)
        
        
        
    
        y2_deriv_p_V = weighted_val(theta_min, theta_max, theta, y2_deriv_p_V_LHS, y2_deriv_p_V_RHS)
        y2_deriv_p_omega = weighted_val(theta_min, theta_max, theta, y2_deriv_p_omega_LHS, y2_deriv_p_omega_RHS)        
    
    if not hinc:
        # reverse the direction of the arrays and switch the sign of the y''VOmega array, jusitification in curvature sign change.nb
        y2_deriv_p_V = y2_deriv_p_V[::-1]
        y2_deriv_p_omega = -y2_deriv_p_omega[::-1]
    
    y2deriv = y2_deriv_p_V + y2_deriv_p_omega
    
    return y2deriv, y2_deriv_p_V, y2_deriv_p_omega

def func_get_all_maxes_a(x, h, theta, V, omega, hinc):
    
    # transition function locations
    theta_min = h/10
    theta_max = h/2
    
    if theta == 0:

        y2_deriv_p_V = func_y2_deriv_px_omega_0_case_a_small_theta(h, theta, x, V)
        y2_deriv_p_omega = func_y2_deriv_px_V_0_case_a_small_theta(h, theta, x, omega)
        
    else:
        
        y2_deriv_p_V_LHS = func_y2_deriv_px_omega_0_case_a_small_theta(h, theta, x, V)
        y2_deriv_p_omega_LHS = func_y2_deriv_px_V_0_case_a_small_theta(h, theta, x, omega)
        
        y2_deriv_p_V_RHS = func_y2_deriv_px_omega_0_case_a(h, theta, x, V)
        y2_deriv_p_omega_RHS = func_y2_deriv_px_V_0_case_a(h, theta, x, omega)
        
        
        
    
        y2_deriv_p_V = weighted_val(theta_min, theta_max, theta, y2_deriv_p_V_LHS, y2_deriv_p_V_RHS)
        y2_deriv_p_omega = weighted_val(theta_min, theta_max, theta, y2_deriv_p_omega_LHS, y2_deriv_p_omega_RHS)        
    
    if not hinc:
        # reverse the direction of the arrays and switch the sign of the y''VOmega array, jusitification in curvature sign change.nb
        y2_deriv_p_V = y2_deriv_p_V[::-1]
        y2_deriv_p_omega = -y2_deriv_p_omega[::-1]
    
    y2deriv = y2_deriv_p_V + y2_deriv_p_omega
    max_deriv = np.max(np.abs(y2deriv))
    arg_max_deriv = np.argmax(np.abs(y2deriv))
    
    xmax = x[arg_max_deriv] # the location of the maximum curvature
    max_pV = y2_deriv_p_V[arg_max_deriv]
    max_pOmega = y2_deriv_p_omega[arg_max_deriv]

    return max_deriv, xmax, max_pV, max_pOmega

def calc_acute_theta_h_inc(theta_a, theta_b):
    if theta_a - theta_b > 0: # clockwise is positive to h will be decreasing
        hinc = False
    else:
        hinc = True
    theta = abs(theta_a-theta_b)
    return theta, hinc

def calc_V_Omega(hinc, theta, x_a, x_b, u_a, u_b, v_a, v_b, omega_a, omega_b, z_a):
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

def get_V_omega_hmin_theta(y):
    
    Vs = np.zeros(len(y[0,:]))
    Omegas = np.zeros(len(y[0,:]))
    hmins = np.zeros(len(y[0,:]))
    thetas = np.zeros(len(y[0,:]))
    hincs = np.zeros(len(y[0,:]))
    Vztrans = np.zeros(len(y[0,:]))
    Vrot = np.zeros(len(y[0,:]))
    Vxtrans = np.zeros(len(y[0,:]))
    
    for indx in range(len(y[0,:])):
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
    
        x_as, x_bs, z_as, z_bs, theta_as, theta_bs, u_as, u_bs, v_as, v_bs, omega_as, omega_bs, hinc = s_prime_to_s_array(x_a, x_b, y_a, y_b, theta_a, theta_b, u_a, u_b, v_a, v_b, omega_a, omega_b)# first convert to frame S
        theta = np.abs(theta_as)
        hmin = z_as - theta/2
        
        V, Omega,  Vztransi, Vroti, Vxtransi = calc_V_Omega(hinc, theta, x_as, x_bs, u_as, u_bs, v_as, v_bs, omega_as, omega_bs, z_as) # then calculate V and Omega
        
        Vs[indx] = V
        Omegas[indx] = Omega
        hmins[indx] = hmin
        thetas[indx] = theta
        hincs[indx] = hinc
        Vztrans[indx] = Vztransi
        Vrot[indx] = Vroti
        Vxtrans[indx] = Vxtransi
        
    return Vs, Omegas, hmins, thetas, hincs ,  Vztrans, Vrot, Vxtrans

# Custom formatter to display labels as 10^x
def fmt(x):
    exponent = np.log10(x)
    return r'$10^{{{:.0f}}}$'.format(exponent)

def func_mask_curvature(t_max, array_to_mask):
    array_to_mask = np.where(t_max == 0, array_to_mask, 0)
    masked_array = np.ma.masked_equal(array_to_mask, 0)
    return masked_array
