#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 13:46:36 2024

@author: cotton
"""

import numpy as np

# applies the forces and torques and stops when the rods start to move away from each other

# 1. start in frame S'
# 2. Convert to frame S
# 3. Calculate forces and torques in frame S
# 4. Convert forces and torques to frame S'
# 5. Apply forces and torques
# 6. Check if end criteria has been met
# 7. Repeat

# definitions (clockwise is positive)

# to do:
    # actually update the ODE
    # run is and hope for the best!

def differential_eqn(t, state, gamma):
    # x_a = x CM of disk A
    # x_b = x CM of disk B
    # y_a = y CM of disk A
    # y_b = y CM of disk B
    # theta_a = angle of disk A (clockwise is positive)
    
    #print(t)
    
    x_a, x_b, y_a, y_b, theta_a, theta_b, u_a, u_b, v_a, v_b, omega_a, omega_b = state
    #print(state)
    
    # convert to frame S (x_a, CM = L*cos(theta)/2, zaCM = hmin + L*theta/2, zbCM = 0)
    x_as, x_bs, z_as, z_bs, theta_as, theta_bs, u_as, u_bs, v_as, v_bs, omega_as, omega_bs =  s_prime_to_s_array(x_a, x_b, y_a, y_b, theta_a, theta_b, u_a, u_b, v_a, v_b, omega_a, omega_b)
    
    # calculate the forces and torques in frame S
    Fxas, Fzas, taua, taub = calc_forces_and_torques_s(theta_as, theta_bs, x_as, x_bs, z_as, u_as, u_bs, v_as, v_bs, omega_as, omega_bs)

    # convert the forces and torques to frame S'
    Fxas_prime, Fzas_prime = calc_forces_s_prime(Fxas, Fzas, theta_b)

    #apply the forces    

    d_ua = Fxas_prime
    d_ub = -gamma*Fxas_prime
    
    d_va = Fzas_prime
    d_vb = gamma*-Fzas_prime

    d_omega_a = 12*taua
    d_omega_b = 12*gamma**3*taub
    
    d_theta_a = omega_a
    d_theta_b = omega_b

    d_xa = u_a
    d_xb = u_b
    
    d_ya = v_a
    d_yb = v_b
        

    return [d_xa, d_xb, d_ya, d_yb, d_theta_a, d_theta_b, d_ua, d_ub, d_va, d_vb, d_omega_a, d_omega_b]

def calc_integrals(theta, h):
    # calculates the scaled integrals
    # inputs: theta = relative angle between the rods, h = minimum distance/length of the top rod
    
    aux_pxV = (theta**-4.)*((2.*(h*(((3.*h)+(2.*theta))*(np.log(((h+theta)/h))))))-(theta*((6.*h)+theta))) 
    
    aux0_pxOm=(np.log((h/(h+theta))))*((-2.*(theta*((4.*h)+(3.*theta))\
    ))+((((h+theta)**2))*(np.log(((h+theta)/h)))));
    aux1_pxOm=0.5*((theta**-3.)*((((-7.*(h**2))/theta)+((0.25*theta)+(\
    (h**2)*((theta**-3.)*aux0_pxOm))))-h));
    
    
    
    one_over_h = 1/theta*np.log((h + theta)/h)
    pV = 1/theta**3*(np.log(h/(h+theta)) + 2*theta/(2*h+theta))
    pomega = 1/(2*theta**3)*(1/(2*h+theta))*(2*h*(3*h/theta+2)*np.log(h/(h+theta)) + 6*h+theta)
    xpV = (0.5*(aux_pxV))/((2.*h)+theta);
    xpomega = aux1_pxOm/((2.*h)+theta);
    
    return one_over_h, pV, pomega, xpV, xpomega

def  find_scaled_hmin(za, theta):
    # works out hmin
    # inputs
    # variables in frame S, acute angle theta

    hmin = (za - 1/2*theta)
        
    return hmin

def find_scaled_hmin_S_prime(y):
    theta = abs(y[4]-y[5])
    theta_b = y[5]
    
    z_aflat = y[2]*np.cos(theta_b) + y[0]*np.sin(theta_b)
    z_bflat = y[3]*np.cos(theta_b) + y[1]*np.sin(theta_b)
    
    hmin = z_aflat - theta/2 - z_bflat
    return hmin

def calc_V_Omega(hinc, theta, x_a, x_b, u_a, u_b, v_a, v_b, omega_a, omega_b):
    # calcluates the scaled values of U and Omega
    
    if hinc: # if the height is increasing with x
        V = (v_a - v_b + omega_a*x_a - omega_b*x_b - theta/2*(u_a - u_b))
    else:
        V = (v_a - v_b + omega_a*x_a - omega_b*x_b + theta/2*(u_a - u_b))
            
    Omega = (omega_a - omega_b)
    
    return V, Omega

def calc_forces_and_torques_s(theta_a, theta_b, x_a, x_b, z_a, u_a, u_b, v_a, v_b, omega_a, omega_b):
    # calculates the scaled forces in frame S, scaled forces = F/(12*eta*W*L) and scaled torques = tau/(12*eta*W*L^2) 
    
    # inputs:
        # acute relative angle theta
        # minimum distance
        # position of xaCM etc.
    
    theta, hinc = calc_acute_theta_h_inc(theta_a, theta_b) # returns the absolute value of delta theta and which orientation we are in
    h = find_scaled_hmin(z_a, theta)
    one_over_h, pV, pomega, xpV, xpomega = calc_integrals(theta, h)
    V, Omega = calc_V_Omega(hinc, theta, x_a, x_b, u_a, u_b, v_a, v_b, omega_a, omega_b)
    
    if hinc: # if the height is increasing with x
        Fza = V*pV + Omega*pomega
        Fxa = -(1/2*theta*Fza + 1/12*(u_a-u_b)*one_over_h + 1/12*omega_a)
        taua = -(V*xpV + Omega*xpomega - x_a*Fza)
        taub = (V*xpV + Omega*xpomega - x_b*Fza)
    else:
        # if the height is decreasing with 
        Fza = V*pV - Omega*pomega
        Fxa = -(-1/2*theta*Fza + 1/12*(u_a-u_b)*one_over_h + 1/12*omega_a)
        taua = -(-V*xpV + Omega*xpomega - x_a*Fza)
        taub = (-V*xpV + Omega*xpomega - x_b*Fza)    
    
    return Fxa, Fza, taua, taub

def calc_forces_s_prime(Fxas, Fzas, thetab):
    # inputs:
       # Fxas: horizontal force on disk A in frame S
       # Fzas: vertical force on disk A in frame S
       # thetab: the angle of disk B is frame S'

   # outputs:
       # Fxas_prime: horizontal force on disk A in frame S'
       # Fzas_prime: vertical force on disk A in frame S'
   
   Fxas_prime = Fxas*np.cos(thetab) + Fzas*np.sin(thetab);
   Fzas_prime = Fzas*np.cos(thetab) - Fxas*np.sin(thetab); 
   return Fxas_prime, Fzas_prime

def calc_acute_theta_h_inc(theta_a, theta_b):
    if theta_a - theta_b > 0: # clockwise is positive to h will be decreasing
        hinc = False
    else:
        hinc = True
    theta = abs(theta_a-theta_b)
    return theta, hinc

def s_prime_to_s_array(x_a, x_b, y_a, y_b, theta_a, theta_b, u_a, u_b, v_a, v_b, omega_a, omega_b):
    
    theta, hinc = calc_acute_theta_h_inc(theta_a, theta_b)
    
    
    z_aflat = y_a*np.cos(theta_b) + x_a*np.sin(theta_b)
    z_bflat = y_b*np.cos(theta_b) + x_b*np.sin(theta_b)
    
    hmin = z_aflat - theta/2 - z_bflat
    
    x_aflat = x_a*np.cos(theta_b) - y_a*np.sin(theta_b)
    x_bflat = x_b*np.cos(theta_b) - y_b*np.sin(theta_b)
    
    if hinc:
        x_as = np.cos(theta)/2
        x_bs = x_bflat - x_aflat + np.cos(theta_b)/2
    else:
        x_as = -np.cos(theta)/2
        x_bs = x_bflat - x_aflat - np.cos(theta_b)/2    
    
    z_as = hmin + theta/2
    z_bs = 0
    
    theta_as = theta_a - theta_b
    theta_bs = 0
    
    u_as = u_a*np.cos(theta_b) - v_a*np.sin(theta_b)
    u_bs = u_b*np.cos(theta_b) - v_b*np.sin(theta_b)
    v_as = v_a*np.cos(theta_b) + u_a*np.sin(theta_b)
    v_bs = v_b*np.cos(theta_b) + u_b*np.sin(theta_b)
    
    omega_as = omega_a
    omega_bs = omega_b
    
    return x_as, x_bs, z_as, z_bs, theta_as, theta_bs, u_as, u_bs, v_as, v_bs, omega_as, omega_bs

def collision(t, state, gamma):
    '''terminates the run if the disks are within collision_dist of each other.
    '''
    
    collision_dist = 10**-6
    
    # find the minumum distance between the rods
    hmin = find_scaled_hmin_S_prime(state) - collision_dist
    return hmin
collision.terminal = True

def no_contact_LHS(t, state, gamma):
    
    theta_b = state[5]              
    x_aflat = state[0]*np.cos(theta_b) - state[2]*np.sin(theta_b)
    x_bflat = state[1]*np.cos(theta_b) - state[3]*np.sin(theta_b)
    
    
    return x_bflat - x_aflat - 1/2*(1/gamma-1)
no_contact_LHS.terminal = True

def no_contact_RHS(t, state, gamma):
    
    theta_b = state[5]              
    x_aflat = state[0]*np.cos(theta_b) - state[2]*np.sin(theta_b)
    x_bflat = state[1]*np.cos(theta_b) - state[3]*np.sin(theta_b)
        
    return x_bflat - x_aflat + 1/2*(1/gamma + 1)
no_contact_LHS.terminal = True

def move_away_RHS(t, state, gamma):
    
    v_a = state[8]
    v_b = state[9]

    omega_a = state[10]
    omega_b = state[11]    

    return v_a - v_b - (omega_a - omega_b)/2 - omega_b*(state[0] - state[1])
move_away_RHS.terminal = False

def move_away_LHS(t, state, gamma):
    
    v_a = state[8]
    v_b = state[9]


    omega_a = state[10]
    omega_b = state[11]    

    return v_a - v_b + (omega_a - omega_b)/2 - omega_b*(state[0] - state[1])
move_away_LHS.terminal = False

def move_away(t, state, gamma):
    v_a = state[8]
    v_b = state[9]


    omega_a = state[10]
    omega_b = state[11]    
    move_away_LHS =  v_a - v_b + (omega_a - omega_b)/2 - omega_b*(state[0] - state[1])
    move_away_RHS =  v_a - v_b - (omega_a - omega_b)/2 - omega_b*(state[0] - state[1])
    
    #print(move_away_LHS, move_away_RHS)
    
    move = 0
    if move_away_LHS > 0:
        move = move + 1
    if move_away_RHS > 0:
        move = move + 1
    return move - 1.5
move_away.terminal = True

def relative_angle_too_large(t, state, gamma):
    return abs(state[4]-state[5]) - 7*np.pi/180
relative_angle_too_large.terminal = True