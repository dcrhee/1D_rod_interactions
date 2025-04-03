#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 15:52:37 2024

@author: cotton

Done case 2 for both a and b
"""

import numpy as np
import sys
sys.path.append('/Users/cotton/Documents/DPhil reading/Polynas/Code/non_eqm/Rotation/Python 1D case/functions/')
from func_weighted_value import weighted_val

def get_curvature_case_2_rod_a_LHS_V_small_theta(x, V, h, P, theta):
    aux0=(-2.*(P*(x**3.)))+((-2.*((P**3.)*(x**3.)))+((2.*((P**4.)*(x**3.))\
    )+(x**4.)));
    aux1=(h**-3.)*((P**-4.)*(V*((4.*((P**3.)*(x**2)))+((-3.*((P**4.)*(\
    x**2)))+aux0))));
    aux2=(-30.*((P**4.)*(x**3.)))+((34.*((P**5.)*(x**3.)))+((-15.*(P*(\
    x**4.)))+(12.*(x**5.))));
    aux3=(60.*((P**4.)*(x**2)))+((-51.*((P**5.)*(x**2)))+((-10.*((P**2)*(\
    x**3.)))+aux2));
    aux4=(82.*((P**6.)*(x**3.)))+((-15.*((P**2)*(x**4.)))+((-24.*(P*(x**5.\
    )))+(24.*(x**6.))));
    aux5=(-123.*((P**6.)*(x**2)))+((-10.*((P**3.)*(x**3.)))+((-66.*((P**5.\
    )*(x**3.)))+aux4));
    aux6=0.00208333*((h**-5.)*((P**-6.)*(V*(((132.*((P**5.)*(x**2)))+aux5)\
    *(theta**2)))));
    output=(0.0416667*aux1)+((-0.00416667*((h**-4.)*((P**-5.)*(V*(aux3*\
    theta)))))+aux6);
    return output

def get_curvature_case_2_rod_a_LHS_Om_small_theta(x, omega, h, P, theta):
    aux0=(-10.*((P**2)*(x**3.)))+((-15.*((P**4.)*(x**3.)))+((16.*((P**5.)*\
    (x**3.)))+(3.*(x**5.))));
    aux1=(h**-3.)*((P**-5.)*(((30.*((P**4.)*(x**2)))+((-24.*((P**5.)*(\
    x**2)))+aux0))*omega));
    aux2=(-17.*((P**5.)*(x**3.)))+((20.*((P**6.)*(x**3.)))+((-5.*((P**2)*(\
    x**4.)))+(3.*(x**6.))));
    aux3=(34.*((P**5.)*(x**2)))+((-30.*((P**6.)*(x**2)))+((-5.*((P**3.)*(\
    x**3.)))+aux2));
    aux4=(2186.*((P**7.)*(x**3.)))+((-315.*((P**3.)*(x**4.)))+((-336.*((\
    P**2)*(x**5.)))+(288.*(x**7.))));
    aux5=(-3279.*((P**7.)*(x**2)))+((-266.*((P**4.)*(x**3.)))+((-1722.*((\
    P**6.)*(x**3.)))+aux4));
    aux6=-0.0000496032*((h**-5.)*((P**-7.)*(((3444.*((P**6.)*(x**2)))+\
    aux5)*((theta**2)*omega))));
    output=(-0.00277778*aux1)+((0.00416667*((h**-4.)*((P**-6.)*(aux3*(\
    theta*omega)))))+aux6);
    return output

def get_curvature_case_2_rod_a_RHS_V_small_theta(x, V, h, P, theta):
    aux0=-0.00416667*((h**-4.)*(V*((((-1.+x)**2))*(((17.*P)+((-30.*x)+(34.\
    *(P*x))))*theta))));
    aux1=(h**-5.)*(V*((((-1.+x)**2))*(((41.*P)+((-66.*x)+(82.*(P*x))))*(\
    theta**2))));
    output=((0.0416667*((h**-3.)*(V*((((-1.+x)**2))*(P+((-2.*x)+(2.*(P*x))\
    ))))))/P)+((aux0/P)+((0.00208333*aux1)/P));
    return output

def get_curvature_case_2_rod_a_RHS_Om_small_theta(x, omega, h, P, theta):
    aux0=-0.00277778*((h**-3.)*((((-1.+x)**2))*(((8.*P)+((-15.*x)+(16.*(P*\
    x))))*omega)));
    aux1=0.00416667*((h**-4.)*((((-1.+x)**2))*(((10.*P)+((-17.*x)+(20.*(P*\
    x))))*(theta*omega))));
    aux2=(h**-5.)*((((-1.+x)**2))*(((1093.*P)+((-1722.*x)+(2186.*(P*x))))*\
    ((theta**2)*omega)));
    output=(aux0/P)+((aux1/P)+((-0.0000496032*aux2)/P));
    return output


def get_curvature_case_2_rod_b_LHS_V_small_theta(x, V, h, x1, P, theta):
    aux0=(h**-3.)*(V*(((P-x1)**-3.)*((((x-x1)**2))*((P**2)+((P*x1)+(-2.*(\
    x*x1)))))));
    aux1=(((x1-x)**2))*(((9.*(P**2))+((4.*(P*x))+((17.*(P*x1))+(-30.*(x*\
    x1)))))*theta);
    aux2=(((x1-x)**2))*(((9.*(P**2))+((16.*(P*x))+((41.*(P*x1))+(-66.*(x*\
    x1)))))*(theta**2));
    aux3=((-0.00416667*((h**-4.)*(V*(((P-x1)**-3.)*aux1))))/P)+((0.00208333*((h**-5.)*(V*(((P-x1)**-3.)*aux2))))/P);
    output=((0.0416667*aux0)/P)+aux3;
    return output

def get_curvature_case_2_rod_b_LHS_Om_small_theta(x, omega, h, x1, P, theta):
    aux0=((P-x1)**-3.)*((((x-x1)**2))*(((6.*(P**2))+((P*x)+((8.*(P*x1))+(-\
    15.*(x*x1)))))*omega));
    aux1=(((x-x1)**2))*(((4.*(P**2))+((3.*(P*x))+((10.*(P*x1))+(-17.*(x*\
    x1)))))*(theta*omega));
    aux2=((165.*(P**2))+((464.*(P*x))+((1093.*(P*x1))+(-1722.*(x*x1)))))*(\
    (theta**2)*omega);
    aux3=((0.00416667*((h**-4.)*(((P-x1)**-3.)*aux1)))/P)+((-0.0000496032*\
    ((h**-5.)*(((P-x1)**-3.)*((((x-x1)**2))*aux2))))/P);
    output=((-0.00277778*((h**-3.)*aux0))/P)+aux3;
    return output

def get_curvature_case_2_rod_b_RHS_V_small_theta(x, V, h, x1, P, theta):
    aux0=(-0.75*((h**-2.)*((P**-3.)*(x**4.))))+(0.666667*((h**-2.)*((P**-\
    3.)*((x**3.)*(P+x)))));
    aux1=(-0.5*((h**-2.)*((-2.*(P**2))+((P*x)+(x*x1)))))+(0.666667*((h**-\
    2.)*((-5.*(P**2))+((3.*(P*x))+((P*x1)+(x*x1))))));
    aux2=(((-0.5*((h**-2.)*((P**-2.)*(x**3.))))+aux0)*((P-x1)**3.))+((((x-\
    x1)**2))*((-0.75*((h**-2.)*(P*((-3.*P)+((2.*x)+x1)))))+aux1));
    aux3=(0.5*((P**-3.)*(V*(((P-x1)**-2.)*((P+x1)*aux2)))))/(-1.+((P**-2.)\
    *(x1**2)));
    aux4=(10.*((P**2)*((x**3.)*(x1**3.))))+((15.*(P*((x**4.)*(x1**3.))))+(\
    -12.*((x**5.)*(x1**3.))));
    aux5=(36.*(P*((x**5.)*(x1**2))))+((17.*((P**5.)*(x1**3.)))+((-30.*((\
    P**4.)*(x*(x1**3.))))+aux4));
    aux6=(-30.*((P**3.)*((x**3.)*(x1**2))))+((-45.*((P**2)*((x**4.)*(\
    x1**2))))+aux5);
    aux7=(9.*((P**6.)*(x1**2)))+((-30.*((P**5.)*(x*(x1**2))))+((60.*((\
    P**4.)*((x**2)*(x1**2))))+aux6));
    aux8=(9.*((P**5.)*((x**2)*x1)))+((45.*((P**3.)*((x**4.)*x1)))+((-36.*(\
    (P**2)*((x**5.)*x1)))+aux7));
    aux9=(-15.*((P**4.)*(x**4.)))+((12.*((P**3.)*(x**5.)))+((-18.*((P**6.)\
    *(x*x1)))+aux8));
    aux10=((P-x1)**-3.)*(((9.*((P**6.)*(x**2)))+((-6.*((P**5.)*(x**3.)))+\
    aux9))*theta);
    aux11=(15.*((P**2)*((x**4.)*(x1**3.))))+((24.*(P*((x**5.)*(x1**3.))))+\
    (-24.*((x**6.)*(x1**3.))));
    aux12=(-66.*((P**5.)*(x*(x1**3.))))+((10.*((P**3.)*((x**3.)*(x1**3.)))\
    )+aux11);
    aux13=(-72.*((P**2)*((x**5.)*(x1**2))))+((72.*(P*((x**6.)*(x1**2))))+(\
    (41.*((P**6.)*(x1**3.)))+aux12));
    aux14=(-30.*((P**4.)*((x**3.)*(x1**2))))+((-45.*((P**3.)*((x**4.)*(\
    x1**2))))+aux13);
    aux15=(9.*((P**7.)*(x1**2)))+((-66.*((P**6.)*(x*(x1**2))))+((132.*((\
    P**5.)*((x**2)*(x1**2))))+aux14));
    aux16=(45.*((P**4.)*((x**4.)*x1)))+((72.*((P**3.)*((x**5.)*x1)))+((-\
    72.*((P**2)*((x**6.)*x1)))+aux15));
    aux17=(-18.*((P**7.)*(x*x1)))+((9.*((P**6.)*((x**2)*x1)))+((-36.*((\
    P**5.)*((x**3.)*x1)))+aux16));
    aux18=(-15.*((P**5.)*(x**4.)))+((-24.*((P**4.)*(x**5.)))+((24.*((P**3.\
    )*(x**6.)))+aux17));
    aux19=((P-x1)**-3.)*(((9.*((P**7.)*(x**2)))+((6.*((P**6.)*(x**3.)))+\
    aux18))*(theta**2));
    aux20=(-0.00416667*((h**-4.)*((P**-5.)*(V*aux10))))+(0.00208333*((h**-\
    5.)*((P**-6.)*(V*aux19))));
    output=(aux3/h)+aux20;
    return output

def get_curvature_case_2_rod_b_RHS_Om_small_theta(x, omega, h, x1, P, theta):
    aux0=(15.*((P**4.)*(x*((x1**3.)*omega))))+((-10.*((P**2)*((x**3.)*(\
    (x1**3.)*omega))))+(3.*((x**5.)*((x1**3.)*omega))));
    aux1=(30.*((P**3.)*((x**3.)*((x1**2)*omega))))+((-9.*(P*((x**5.)*((\
    x1**2)*omega))))+((-8.*((P**5.)*((x1**3.)*omega)))+aux0));
    aux2=(15.*((P**5.)*(x*((x1**2)*omega))))+((-30.*((P**4.)*((x**2)*((\
    x1**2)*omega))))+aux1);
    aux3=(-15.*((P**4.)*((x**3.)*(x1*omega))))+((9.*((P**2)*((x**5.)*(\
    x1*omega))))+((-6.*((P**6.)*((x1**2)*omega)))+aux2));
    aux4=(-3.*((P**3.)*((x**5.)*omega)))+((12.*((P**6.)*(x*(x1*omega\
    ))))+((-6.*((P**5.)*((x**2)*(x1*omega))))+aux3));
    aux5=((P-x1)**-3.)*((-6.*((P**6.)*((x**2)*omega)))+((9.*((P**5.)*((\
    x**3.)*omega)))+aux4));
    aux6=(5.*((P**3.)*((x**3.)*((x1**3.)*omega))))+((5.*((P**2)*((x**4.\
    )*((x1**3.)*omega))))+(-3.*((x**6.)*((x1**3.)*omega))));
    aux7=(9.*(P*((x**6.)*((x1**2)*omega))))+((10.*((P**6.)*((x1**3.)*\
    omega)))+((-17.*((P**5.)*(x*((x1**3.)*omega))))+aux6));
    aux8=(-15.*((P**4.)*((x**3.)*((x1**2)*omega))))+((-15.*((P**3.)*((\
    x**4.)*((x1**2)*omega))))+aux7);
    aux9=(-17.*((P**6.)*(x*((x1**2)*omega))))+((34.*((P**5.)*((x**2)*((\
    x1**2)*omega))))+aux8);
    aux10=(15.*((P**4.)*((x**4.)*(x1*omega))))+((-9.*((P**2)*((x**6.)*(\
    x1*omega))))+((4.*((P**7.)*((x1**2)*omega)))+aux9));
    aux11=(4.*((P**6.)*((x**2)*(x1*omega))))+((-2.*((P**5.)*((x**3.)*(\
    x1*omega))))+aux10);
    aux12=(-5.*((P**5.)*((x**4.)*omega)))+((3.*((P**3.)*((x**6.)*\
    omega)))+((-8.*((P**7.)*(x*(x1*omega))))+aux11));
    aux13=((P-x1)**-3.)*(theta*((4.*((P**7.)*((x**2)*omega)))+((-2.*\
    ((P**6.)*((x**3.)*omega)))+aux12)));
    aux14=(-315.*((P**3.)*((x**4.)*((x1**3.)*omega))))+((-336.*((P**2)*\
    ((x**5.)*((x1**3.)*omega))))+(288.*((x**7.)*((x1**3.)*omega))));\
    
    aux15=(1722.*((P**6.)*(x*((x1**3.)*omega))))+((-266.*((P**4.)*((\
    x**3.)*((x1**3.)*omega))))+aux14);
    aux16=(-864.*(P*((x**7.)*((x1**2)*omega))))+((-1093.*((P**7.)*((\
    x1**3.)*omega)))+aux15);
    aux17=(945.*((P**4.)*((x**4.)*((x1**2)*omega))))+((1008.*((P**3.)*(\
    (x**5.)*((x1**2)*omega))))+aux16);
    aux18=(-3444.*((P**6.)*((x**2)*((x1**2)*omega))))+((798.*((P**5.)*(\
    (x**3.)*((x1**2)*omega))))+aux17);
    aux19=(-165.*((P**8.)*((x1**2)*omega)))+((1722.*((P**7.)*(x*((\
    x1**2)*omega))))+aux18);
    aux20=(-1008.*((P**4.)*((x**5.)*(x1*omega))))+((864.*((P**2)*((\
    x**7.)*(x1*omega))))+aux19);
    aux21=(924.*((P**6.)*((x**3.)*(x1*omega))))+((-945.*((P**5.)*((\
    x**4.)*(x1*omega))))+aux20);
    aux22=(-288.*((P**3.)*((x**7.)*omega)))+((330.*((P**8.)*(x*(x1*\
    omega))))+((-165.*((P**7.)*((x**2)*(x1*omega))))+aux21));
    aux23=(-198.*((P**7.)*((x**3.)*omega)))+((315.*((P**6.)*((x**4.)*\
    omega)))+((336.*((P**5.)*((x**5.)*omega)))+aux22));
    aux24=(h**-5.)*((P**-7.)*(((P-x1)**-3.)*((theta**2)*((-165.*((P**8.\
    )*((x**2)*omega)))+aux23))));
    output=(0.00277778*((h**-3.)*((P**-5.)*aux5)))+((0.00416667*((h**-4.)*\
    ((P**-6.)*aux13)))+(0.0000496032*aux24));
    return output

def get_curvature_case_2_rod_a_LHS_V(x, V, h, P, theta):
    aux0=(6.*(h*(P*(1.+(P*(x*(-3.+(2.*x))))))))+((x+(P*(4.+(x*(-8.+((-3.*\
    P)+(2.*((2.+P)*x))))))))*theta);
    aux1=(3.*((h**2)*(P+(2.*(P*x)))))+((2.*(h*((P+(x+(2.*(P*x))))*theta\
    )))+(x*(theta**2)));
    aux2=(3.*((h**2)*(P*(-3.+(2.*x)))))+((2.*(h*((-2.+((-3.*P)+(x+(2.*(P*\
    x)))))*theta)))+((-2.+x)*(theta**2)));
    aux3=((3.*((h**2)*P))+((2.*(h*((P+x)*theta)))+(x*(theta**2))))*(\
    np.log((h+((x*theta)/P))));
    aux4=(2.*(P*((((-1.+x)**2))*(aux1*(np.log(h))))))+((-2.*(P*((x**2)*(\
    aux2*(np.log((h+theta)))))))+(-2.*(P*aux3)));
    output=(0.5*((P**-2.)*(V*((theta**-4.)*((x*(theta*aux0))+aux4)))\
    ))/((2.*h)+theta);
    return output

def get_curvature_case_2_rod_a_LHS_Om(x, omega, h, P, theta):
    aux0=(P**2)*((((-1.+x)**2))*((1.+(2.*x))*((((h+theta)**2))*(((np.\
    log(h))**2)))));
    aux1=(h**2)*((P**2)*((x**2)*((-3.+(2.*x))*((((h+theta)**2))*(((np.\
    log((h+theta)))**2))))));
    aux2=((h**2)*((x+(-4.*(P*(1.+(3.*(x*(-2.+((-3.*P)+(x+(2.*(P*x)))))))))\
    ))*theta))+((2.*(h*((x+(P*(-1.+(-4.*((-2.+x)*x)))))*(theta**2)))\
    )+(x*(theta**3.)));
    aux3=(x*(theta*((-2.*((h**3.)*(P*(1.+(8.*(P*(x*(-3.+(2.*x)))))))))+\
    aux2)))+(2.*((h**2)*((P**2)*((((h+theta)**2))*(np.log((h+((x*\
    theta)/P))))))));
    aux4=(2.*(h*(((5.*x)+(P*((-3.+(2.*x))*(-3.+(2.*((3.+P)*x))))))*\
    theta)))+(x*((3.+(P*(-8.+((3.*P)+((4.*x)+(-2.*(P*x)))))))*(theta\
    **2)));
    aux5=(3.*((h**2)*(P*(((3.*P)+(4.*x))*theta))))+((2.*(h*(x*(((3.*P)+\
    x)*(theta**2)))))+((x**2)*(theta**3.)));
    aux6=(x*(theta*((28.*((h**2)*(P*(1.+(P*(x*(-3.+(2.*x))))))))+aux4))\
    )+(-2.*(((14.*((h**3.)*(P**2)))+aux5)*(np.log((h+((x*theta)/P))))));
    aux7=(-4.*(P*(x*(4.+(3.*((-2.+x)*x))))))+(-3.*((P**2)*(3.+(4.*((x**2)*\
    (-3.+(2.*x)))))));
    aux8=(-2.*((h**2)*(P*(x+(P*(7.+(8.*((x**2)*(-3.+(2.*x))))))))))+((h*((\
    (x**2)+aux7)*theta))+(-8.*(P*((((-1.+x)**2))*(x*(theta**2))))));\
    
    aux9=((1.+((-6.*(x**2))+(4.*(x**3.))))*(np.log((h+theta))))+(np.\
    log((h+((x*theta)/P))));
    aux10=(theta*aux6)+(-2.*(h*((np.log(h))*((theta*aux8)+(2.*(h*((\
    P**2)*((((h+theta)**2))*aux9))))))));
    aux11=(theta**-6.)*(omega*((4.*((h**2)*aux0))+((4.*aux1)+((2.*((\
    np.log((h+theta)))*aux3))+aux10))));
    output=(0.125*((P**-2.)*aux11))/((2.*h)+theta);
    
        
    return output

def get_curvature_case_2_rod_a_RHS_V(x, V, h, P, theta):
    aux0=(3.*((h**2)*(P+(2.*(P*x)))))+((2.*(h*((P+(x+(2.*(P*x))))*\
    theta)))+(x*(theta**2)));
    aux1=(theta*((6.*(h*(P+(2.*(P*x)))))+((P+(2.*((2.+P)*x)))*theta)\
    ))+(2.*(aux0*((np.log(h))-(np.log((h+theta))))));
    output=((0.5*(V*((((-1.+x)**2))*((theta**-4.)*aux1))))/((2.*h)+\
    theta))/P;

    return output

def get_curvature_case_2_rod_a_RHS_Om(x, omega, h, P, theta):
    aux0=((28.*((h**2)*(P+(2.*(P*x)))))+(4.*(h*((P+(2.*((3.+P)*x)))*\
    theta))))-((P+(2.*((-2.+P)*x)))*(theta**2));
    aux1=(4.*(x*(theta**3.)))+(h*(P*((1.+(2.*x))*((((h+theta)**2))*(\
    (np.log(h))-(np.log((h+theta))))))));
    aux2=(8.*((h**2)*(P*((1.+(2.*x))*theta))))+((6.*(h*((P+(x+(2.*(P*x)\
    )))*(theta**2))))+aux1);
    aux3=(theta**-6.)*(((theta**2)*(aux0*omega))+(4.*(h*(omega\
    *(aux2*((np.log(h))-(np.log((h+theta)))))))));
    output=((0.125*((((-1.+x)**2))*aux3))/((2.*h)+theta))/P;

    return output

def get_curvature_case_2_rod_b_LHS_V(x, V, h, x1, P, theta):
    # P = the RHS of rod B
    # x1 = the LHS of rod B
    aux0=(6.*(h*(P*((-3.*P)+((2.*x)+x1)))))+(((-11.*(P**2))+((4.*(x*x1))+(\
    P*((6.*x)+x1))))*theta);
    aux1=(2.*(h*(((P*((5.*P)+(-3.*x)))-((P+x)*x1))*theta)))+((((2.*(\
    P**2))-(x*x1))-(P*x))*(theta**2));
    aux2=((3.*((h**2)*(P*(((3.*P)+(-2.*x))-x1))))+aux1)*((np.log((h+\
    theta)))-(np.log(h)));
    aux3=0.5*(V*(((P-x1)**-3.)*((((x-x1)**2))*((theta**-4.)*((theta*\
    aux0)+(2.*aux2))))));
    output=(aux3/((2.*h)+theta))/P;
    return output

def get_curvature_case_2_rod_b_LHS_Om(x, omega, h, x1, P, theta):
    aux0=(4.*(h*((((15.*(P**2))+(-6.*(x*x1)))-(P*((8.*x)+x1)))*theta)))\
    +(((5.*(P**2))+((-4.*(x*x1))+(P*((-2.*x)+x1))))*(theta**2));
    aux1=h*(P*((((3.*P)+(-2.*x))-x1)*((((h+theta)**2))*((np.log(h))-(\
    np.log((h+theta)))))));
    aux2=(6.*(h*(((P*((5.*P)+(-3.*x)))-((P+x)*x1))*(theta**2))))+((4.*(\
    (((2.*(P**2))-(x*x1))-(P*x))*(theta**3.)))+aux1);
    aux3=((8.*((h**2)*(P*((((3.*P)+(-2.*x))-x1)*theta))))+aux2)*((np.\
    log(h))-(np.log((h+theta))));
    aux4=((theta**2)*(((28.*((h**2)*(P*(((3.*P)+(-2.*x))-x1))))+aux0)*\
    omega))+(4.*(h*(omega*aux3)));
    output=(0.125*((P**-4.)*((((x-x1)**2))*(((-1.+(x1/P))**-3.)*((theta\
    **-6.)*aux4)))))/((2.*h)+theta);
    return output

def get_curvature_case_2_rod_b_RHS_V(x, V, h, x1, P, theta):
    aux0=h*(P*((P*((P+(-2.*x))*x))+((3.*(P*(x*x1)))+((-3.*(P*(x1**2)))+(\
    x1**3.)))));
    aux1=(2.*(P*(((5.*P)+(-2.*x))*(x*x1))))+((P*(((-11.*P)+(5.*x))*(x1**2)\
    ))+((P+x)*(x1**3.)));
    aux2=(2.*(h*(((-5.*(P**2))+((3.*(P*x))+((P*x1)+(x*x1))))*theta)))+(\
    ((-2.*(P**2))+((P*x)+(x*x1)))*(theta**2));
    aux3=(((x-x1)**2))*(((3.*((h**2)*(P*((-3.*P)+((2.*x)+x1)))))+aux2)*(\
    np.log((h+theta))));
    aux4=((3.*((h**2)*P))+((2.*(h*((P+x)*theta)))+(x*(theta**2))))*(\
    np.log((h+((x*theta)/P))));
    aux5=(2.*(h*(((P**2)+((3.*(P*(x-x1)))+((x+(-2.*x1))*x1)))*theta)))+\
    (((P*x)+((x+(-2.*x1))*x1))*(theta**2));
    aux6=(((P-x)**2))*(((3.*((h**2)*(P*(P+((2.*x)+(-3.*x1))))))+aux5)*(np.\
    log(h)));
    aux7=((P-x)*(theta*((-6.*aux0)-(((2.*((P**2)*(((2.*P)+(-3.*x))*x)))\
    +aux1)*theta))))+(2.*(P*((aux3+(((P-x1)**3.)*aux4))-aux6)));
    aux8=(0.5*((P**-4.)*(V*(((P-x1)**-2.)*((P+x1)*((theta**-4.)*aux7)))\
    )))/((2.*h)+theta);
    output=aux8/(-1.+((P**-2.)*(x1**2)));
    return output

def get_curvature_case_2_rod_b_RHS_Om(x, omega, h, x1, P, theta):
    aux0=(h**2)*(P*((P*((P+(-2.*x))*x))+((3.*(P*(x*x1)))+((-3.*(P*(x1**2))\
    )+(x1**3.)))));
    aux1=(3.*(P*(((11.*P)+(-4.*x))*(x*x1))))+((3.*(P*(((-10.*P)+(3.*x))*(\
    x1**2))))+(((2.*P)+(5.*x))*(x1**3.)));
    aux2=(2.*(P*(x*(((-5.*P)+(2.*x))*x1))))+((P*(((5.*P)+x)*(x1**2)))+((P+\
    (-3.*x))*(x1**3.)));
    aux3=(-2.*(h*((((P**2)*(((9.*P)+(-16.*x))*x))+aux1)*theta)))+(((2.*\
    ((P**2)*(x**2)))+aux2)*(theta**2));
    aux4=(8.*(((P**2)+((P*x1)+(x1**2)))*(theta**3.)))+(3.*(h*(P*((P+x1)\
    *((((h+theta)**2))*((np.log(h))-(np.log((h+theta)))))))));
    aux5=(24.*((h**2)*(P*((P+x1)*theta))))+((6.*(h*(((5.*(P**2))+((5.*(\
    P*x1))+(2.*(x1**2))))*(theta**2))))+aux4);
    aux6=h*(P*((P-x)*((P+((2.*x)+(-3.*x1)))*(aux5*((np.log(h))-(np.log((h+\
    theta))))))));
    aux7=((P-x1)**-2.)*(omega*((3.*((P+x1)*((theta**2)*((-28.*aux0)+\
    aux3))))+(-4.*aux6)));
    aux8=(h**2)*((P**2)*((P-x1)*((P+x1)*((((h+theta)**2))*(omega*(((\
    np.log((h+theta)))**2)))))));
    aux9=(11.*(P**3.))+((-5.*((P**2)*x))+((-4.*(P*(x**2)))+((-3.*(P*(\
    x1**2)))+(x*(x1**2)))));
    aux10=(6.*((h**2)*(P*((P-x1)*(P+x1)))))+((3.*(h*(aux9*theta)))+(8.*\
    (P*((P-x)*(((2.*P)+x)*(theta**2))))));
    aux11=(P-x1)*((P+x1)*((((h+theta)**2))*((np.log((h+theta)))-(np.\
    log((h+((x*theta)/P)))))));
    aux12=(3.*((h**2)*(P*(((3.*P)+(4.*x))*theta))))+((2.*(h*(x*(((3.*P)\
    +x)*(theta**2)))))+((x**2)*(theta**3.)));
    aux13=(P+x1)*(theta*(((14.*((h**3.)*(P**2)))+aux12)*(omega*(np.\
    log((h+((x*theta)/P)))))));
    aux14=((2.*P)+x)*((((10.*(P**3.))+(P*(x*((-7.*P)+(4.*x)))))-(((6.*P)+\
    x)*(x1**2)))*theta);
    aux15=(-3.*((P**3.)*x))+((P*((x**2)*((3.*P)+(4.*x))))+(-3.*(x*(((3.*P)\
    +x)*(x1**2)))));
    aux16=(2.*(h*(((8.*(P**4.))+aux15)*(theta**2))))+(3.*((x**2)*((P-\
    x1)*((P+x1)*(theta**3.)))));
    aux17=(6.*((h**3.)*(P*(((8.*P)-x)*((P-x1)*(P+x1))))))+((3.*((h**2)*\
    aux14))+aux16);
    aux18=(h**2)*((P**2)*((P-x1)*((P+x1)*((((h+theta)**2))*(omega*(\
    np.log((h+((x*theta)/P)))))))));
    aux19=(2.*(h*(omega*((np.log(h))*(((P-x)*(theta*aux10))+(-6.*(h*\
    ((P**2)*aux11))))))))+((6.*((P-x1)*aux13))+(-2.*((np.log((h+theta))\
    )*((theta*(aux17*omega))+(6.*aux18)))));
    aux20=(0.0416667*((P**-4.)*((theta**-6.)*(((P-x)*aux7)+((12.*aux8)+\
    aux19)))))/((2.*h)+theta);
    output=aux20/(-1.+((P**-2.)*(x1**2)));
    return output


def func_get_all_curvatures_a_case2(x, h, theta, V, Omega, x1):
    xLHS = x[x<=x1]
    xRHS = x[x>x1]
    
    # transition function locations
    theta_min = h/20
    theta_max = h/10

    if theta == 0:

        y2_derivLHSV = get_curvature_case_2_rod_a_LHS_V_small_theta(xLHS, V, h, x1, theta)
        y2_derivLHSOm = get_curvature_case_2_rod_a_LHS_Om_small_theta(xLHS, Omega, h, x1, theta)
        
        y2_derivRHSV = get_curvature_case_2_rod_a_RHS_V_small_theta(xRHS, V, h, x1, theta)
        y2_derivRHSOm = get_curvature_case_2_rod_a_RHS_Om_small_theta(xRHS, Omega, h, x1, theta)
        
    else:
        
        y2_derivLHSVa = get_curvature_case_2_rod_a_LHS_V_small_theta(xLHS, V, h, x1, theta)
        y2_derivLHSOma = get_curvature_case_2_rod_a_LHS_Om_small_theta(xLHS, Omega, h, x1, theta)
        
        y2_derivRHSVa = get_curvature_case_2_rod_a_RHS_V_small_theta(xRHS, V, h, x1, theta)
        y2_derivRHSOma = get_curvature_case_2_rod_a_RHS_Om_small_theta(xRHS, Omega, h, x1, theta)
        
        y2_derivLHSVb = get_curvature_case_2_rod_a_LHS_V(xLHS, V, h, x1, theta)
        y2_derivLHSOmb = get_curvature_case_2_rod_a_LHS_Om(xLHS, Omega, h, x1, theta)
        
        y2_derivRHSVb = get_curvature_case_2_rod_a_RHS_V(xRHS, V, h, x1, theta)
        y2_derivRHSOmb = get_curvature_case_2_rod_a_RHS_Om(xRHS, Omega, h, x1, theta)

        y2_derivLHSV = weighted_val(theta_min, theta_max, theta, y2_derivLHSVa, y2_derivLHSVb)
        y2_derivLHSOm = weighted_val(theta_min, theta_max, theta, y2_derivLHSOma, y2_derivLHSOmb)      
    
        y2_derivRHSV = weighted_val(theta_min, theta_max, theta, y2_derivRHSVa, y2_derivRHSVb)
        y2_derivRHSOm = weighted_val(theta_min, theta_max, theta, y2_derivRHSOma, y2_derivRHSOmb)    
        
    
    y2_deriv_p_V = np.concatenate((y2_derivLHSV, y2_derivRHSV))
    y2_deriv_p_Om = np.concatenate((y2_derivLHSOm, y2_derivRHSOm))
    
    
    y2_deriv = y2_deriv_p_V + y2_deriv_p_Om
    
    return y2_deriv, y2_deriv_p_V, y2_deriv_p_Om

def func_get_all_curvatures_b_case2(x, V, Omega, h, theta):
    xLHS = x[x<=0]
    xRHS = x[x>0]
    M = x[0]
    P = x[-1]
    
    # transition function locations
    theta_min = h/20
    theta_max = h/10

    if theta == 0:

        y2_derivLHSV = get_curvature_case_2_rod_b_LHS_V_small_theta(xLHS, V, h, M, P, theta)
        y2_derivLHSOm = get_curvature_case_2_rod_b_LHS_Om_small_theta(xLHS, Omega, h, M, P, theta)
        
        y2_derivRHSV = get_curvature_case_2_rod_b_RHS_V_small_theta(xRHS, V, h, M, P, theta)
        y2_derivRHSOm = get_curvature_case_2_rod_b_RHS_Om_small_theta(xRHS, Omega, h, M, P, theta)
        
    else:
        
        y2_derivLHSVa = get_curvature_case_2_rod_b_LHS_V_small_theta(xLHS, V, h, M, P, theta)
        y2_derivLHSOma = get_curvature_case_2_rod_b_LHS_Om_small_theta(xLHS, Omega, h, M, P, theta)
        
        y2_derivRHSVa = get_curvature_case_2_rod_b_RHS_V_small_theta(xRHS, V, h, M, P, theta)
        y2_derivRHSOma = get_curvature_case_2_rod_b_RHS_Om_small_theta(xRHS, Omega, h, M, P, theta)
        
        y2_derivLHSVb = get_curvature_case_2_rod_b_LHS_V(xLHS, V, h, M, P, theta)
        y2_derivLHSOmb = get_curvature_case_2_rod_b_LHS_Om(xLHS, Omega, h, M, P, theta)
        
        y2_derivRHSVb = get_curvature_case_2_rod_b_RHS_V(xRHS, V, h, M, P, theta)
        y2_derivRHSOmb = get_curvature_case_2_rod_b_RHS_Om(xRHS, Omega, h, M, P, theta)

        y2_derivLHSV = weighted_val(theta_min, theta_max, theta, y2_derivLHSVa, y2_derivLHSVb)
        y2_derivLHSOm = weighted_val(theta_min, theta_max, theta, y2_derivLHSOma, y2_derivLHSOmb)      
    
        y2_derivRHSV = weighted_val(theta_min, theta_max, theta, y2_derivRHSVa, y2_derivRHSVb)
        y2_derivRHSOm = weighted_val(theta_min, theta_max, theta, y2_derivRHSOma, y2_derivRHSOmb)    
        
    
    y2_deriv_p_V = np.concatenate((y2_derivLHSV, y2_derivRHSV))
    y2_deriv_p_Om = np.concatenate((y2_derivLHSOm, y2_derivRHSOm))
    
    y2_deriv = y2_deriv_p_V + y2_deriv_p_Om
    
    return y2_deriv, y2_deriv_p_V, y2_deriv_p_Om