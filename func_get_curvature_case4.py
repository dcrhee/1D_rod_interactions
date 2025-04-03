#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 14:52:04 2024

@author: cotton

Calculates the curvatures when the rods have fallen off the right hand side i.e. x_a + 1/2 > x_b + 1/(2*gamma)

In both cases we say the left hand end of the rod is located at x = 0, such that when converting to frame S' we don't need to change x_aCM

Done case 4

"""

import numpy as np
import sys
sys.path.append('/Users/cotton/Documents/DPhil reading/Polynas/Code/non_eqm/Rotation/Python 1D case/functions/')
from func_weighted_value import weighted_val

def get_curvature_case_4_rod_a_LHS_V_small_theta(x, V, h, P, theta):
    aux0=(-2.*(P*(x**3.)))+((-2.*((P**3.)*(x**3.)))+((2.*((P**4.)*(x**3.))\
    )+(x**4.)));
    aux1=(h**-3.)*((P**-4.)*(V*((4.*((P**3.)*(x**2)))+((-3.*((P**4.)*(\
    x**2)))+aux0))));
    aux2=(-30.*((P**4.)*(x**3.)))+((26.*((P**5.)*(x**3.)))+((45.*(P*(x**4.\
    )))+(-12.*(x**5.))));
    aux3=(60.*((P**4.)*(x**2)))+((-39.*((P**5.)*(x**2)))+((-50.*((P**2)*(\
    x**3.)))+aux2));
    aux4=(50.*((P**6.)*(x**3.)))+((225.*((P**2)*(x**4.)))+((-120.*(P*(\
    x**5.)))+(24.*(x**6.))));
    aux5=(-75.*((P**6.)*(x**2)))+((-170.*((P**3.)*(x**3.)))+((-66.*((P**5.\
    )*(x**3.)))+aux4));
    aux6=0.00208333*((h**-5.)*((P**-6.)*(V*(((132.*((P**5.)*(x**2)))+aux5)\
    *(theta**2)))));
    output=(0.0416667*aux1)+((-0.00416667*((h**-4.)*((P**-5.)*(V*(aux3*\
    theta)))))+aux6);
    return output

def get_curvature_case_4_rod_a_LHS_Om_small_theta(x, omega, h, P, theta):
    aux0=(-10.*((P**2)*(x**3.)))+((-15.*((P**4.)*(x**3.)))+((16.*((P**5.)*\
    (x**3.)))+(3.*(x**5.))));
    aux1=(h**-3.)*((P**-5.)*(((30.*((P**4.)*(x**2)))+((-24.*((P**5.)*(\
    x**2)))+aux0))*omega));
    aux2=(12.*((P**6.)*(x**3.)))+((5.*((P**2)*(x**4.)))+((6.*(P*(x**5.)))+\
    (-3.*(x**6.))));
    aux3=(-18.*((P**6.)*(x**2)))+((-15.*((P**3.)*(x**3.)))+((-13.*((P**5.)\
    *(x**3.)))+aux2));
    aux4=(1365.*((P**3.)*(x**4.)))+((672.*((P**2)*(x**5.)))+((-1008.*(P*(\
    x**6.)))+(288.*(x**7.))));
    aux5=(-1946.*((P**4.)*(x**3.)))+((-1050.*((P**6.)*(x**3.)))+((842.*((\
    P**7.)*(x**3.)))+aux4));
    aux6=(P**-7.)*(((2100.*((P**6.)*(x**2)))+((-1263.*((P**7.)*(x**2)))+\
    aux5))*((theta**2)*omega));
    aux7=(0.00416667*((h**-4.)*((P**-6.)*(((26.*((P**5.)*(x**2)))+aux3)*(\
    theta*omega)))))+(-0.0000496032*((h**-5.)*aux6));
    output=(-0.00277778*aux1)+aux7;
    return output

def get_curvature_case_4_rod_a_RHS_V_small_theta(x, V, h, P, theta):
    aux0=-0.00416667*((h**-4.)*(V*((((-1.+x)**2))*(((13.*P)+((-30.*x)+(26.\
    *(P*x))))*theta))));
    aux1=(h**-5.)*(V*((((-1.+x)**2))*(((25.*P)+((-66.*x)+(50.*(P*x))))*(\
    theta**2))));
    output=((0.0416667*((h**-3.)*(V*((((-1.+x)**2))*(P+((-2.*x)+(2.*(P*x))\
    ))))))/P)+((aux0/P)+((0.00208333*aux1)/P));

    return output

def get_curvature_case_4_rod_a_RHS_Om_small_theta(x, omega, h, P, theta):
    aux0=-0.00277778*((h**-3.)*((((-1.+x)**2))*(((8.*P)+((-15.*x)+(16.*(P*\
    x))))*omega)));
    aux1=0.00416667*((h**-4.)*((((-1.+x)**2))*(((6.*P)+((-13.*x)+(12.*(P*\
    x))))*(theta*omega))));
    aux2=(h**-5.)*((((-1.+x)**2))*(((421.*P)+((-1050.*x)+(842.*(P*x))))*((\
    theta**2)*omega)));
    output=(aux0/P)+((aux1/P)+((-0.0000496032*aux2)/P));
    return output


def get_curvature_case_4_rod_b_LHS_V_small_theta(x, V, h, x1, P, theta):
    aux0=(h**-3.)*(V*(((P-x1)**-3.)*((((x-x1)**2))*((P**2)+((P*x1)+(-2.*(\
    x*x1)))))));
    aux1=(((x-x1)**2))*(((21.*(P**2))+((-4.*(P*x))+((13.*(P*x1))+(-30.*(x*\
    x1)))))*theta);
    aux2=(((x-x1)**2))*(((57.*(P**2))+((-16.*(P*x))+((25.*(P*x1))+(-66.*(\
    x*x1)))))*(theta**2));
    aux3=((-0.00416667*((h**-4.)*(V*(((P-x1)**-3.)*aux1))))/P)+((0.00208333*((h**-5.)*(V*(((P-x1)**-3.)*aux2))))/P);
    output=((0.0416667*aux0)/P)+aux3;
    return output

def get_curvature_case_4_rod_b_LHS_Om_small_theta(x, omega, h, x1, P, theta):
    aux0=((P-x1)**-3.)*((((x-x1)**2))*(((6.*(P**2))+((P*x)+((8.*(P*x1))+(-\
    15.*(x*x1)))))*omega));
    aux1=(((x-x1)**2))*((((8.*(P**2))+((6.*(P*x1))+(-13.*(x*x1))))-(P*x))*\
    (theta*omega));
    aux2=((837.*(P**2))+((-208.*(P*x))+((421.*(P*x1))+(-1050.*(x*x1)))))*(\
    (theta**2)*omega);
    aux3=((0.00416667*((h**-4.)*(((P-x1)**-3.)*aux1)))/P)+((-0.0000496032*\
    ((h**-5.)*(((P-x1)**-3.)*((((x-x1)**2))*aux2))))/P);
    output=((-0.00277778*((h**-3.)*aux0))/P)+aux3;
    return output

def get_curvature_case_4_rod_b_RHS_V_small_theta(x, V, h, x1, P, theta):
    aux0=(-3.*((P**2)*((x**2)*x1)))+(((P**3.)*(x1**2))+((3.*(P*((x**2)*(\
    x1**2))))+((P**2)*(x1**3.))));
    aux1=((P-x1)**-3.)*((((P**3.)*(x**2))+((-2.*((P**3.)*(x*x1)))+aux0))-(\
    (x**2)*(x1**3.)));
    aux2=(-4.*((P**2)*(x*(x1**3.))))+((-21.*(P*((x**2)*(x1**3.))))+(12.*((\
    x**3.)*(x1**3.))));
    aux3=(63.*((P**2)*((x**2)*(x1**2))))+((-36.*(P*((x**3.)*(x1**2))))+((\
    13.*((P**3.)*(x1**3.)))+aux2));
    aux4=(36.*((P**2)*((x**3.)*x1)))+((21.*((P**4.)*(x1**2)))+((12.*((\
    P**3.)*(x*(x1**2))))+aux3));
    aux5=(-12.*((P**3.)*(x**3.)))+((-42.*((P**4.)*(x*x1)))+((-63.*((P**3.)\
    *((x**2)*x1)))+aux4));
    aux6=(P**-5.)*(V*((((P-x)**2))*(((P-x1)**-3.)*(((21.*((P**4.)*(x**2)))\
    +aux5)*theta))));
    aux7=(-57.*((P**2)*((x**2)*(x1**3.))))+((72.*(P*((x**3.)*(x1**3.))))+(\
    -24.*((x**4.)*(x1**3.))));
    aux8=(72.*(P*((x**4.)*(x1**2))))+((25.*((P**4.)*(x1**3.)))+((-16.*((\
    P**3.)*(x*(x1**3.))))+aux7));
    aux9=(171.*((P**3.)*((x**2)*(x1**2))))+((-216.*((P**2)*((x**3.)*(\
    x1**2))))+aux8);
    aux10=(-72.*((P**2)*((x**4.)*x1)))+((57.*((P**5.)*(x1**2)))+((48.*((\
    P**4.)*(x*(x1**2))))+aux9));
    aux11=(-114.*((P**5.)*(x*x1)))+((-171.*((P**4.)*((x**2)*x1)))+((216.*(\
    (P**3.)*((x**3.)*x1)))+aux10));
    aux12=(57.*((P**5.)*(x**2)))+((-72.*((P**4.)*(x**3.)))+((24.*((P**3.)*\
    (x**4.)))+aux11));
    aux13=0.00208333*((h**-5.)*((P**-6.)*(V*((((P-x)**2))*(((P-x1)**-3.)*(\
    aux12*(theta**2)))))));
    output=(0.0416667*((h**-3.)*((P**-4.)*(V*((((P-x)**2))*aux1)))))+((-0.00416667*((h**-4.)*aux6))+aux13);
    return output

def get_curvature_case_4_rod_b_RHS_Om_small_theta(x, omega, h, x1, P, theta):
    aux0=(-15.*((P**4.)*(x*(x1**3.))))+((10.*((P**2)*((x**3.)*(x1**3.))))\
    +(-3.*((x**5.)*(x1**3.))));
    aux1=(-30.*((P**3.)*((x**3.)*(x1**2))))+((9.*(P*((x**5.)*(x1**2))))+((\
    8.*((P**5.)*(x1**3.)))+aux0));
    aux2=(6.*((P**6.)*(x1**2)))+((-15.*((P**5.)*(x*(x1**2))))+((30.*((\
    P**4.)*((x**2)*(x1**2))))+aux1));
    aux3=(6.*((P**5.)*((x**2)*x1)))+((15.*((P**4.)*((x**3.)*x1)))+((-9.*((\
    P**2)*((x**5.)*x1)))+aux2));
    aux4=(-9.*((P**5.)*(x**3.)))+((3.*((P**3.)*(x**5.)))+((-12.*((P**6.)*(\
    x*x1)))+aux3));
    aux5=-0.00277778*((h**-3.)*((P**-5.)*(((P-x1)**-3.)*(((6.*((P**6.)*(\
    x**2)))+aux4)*omega))));
    aux6=(-5.*((P**2)*((x**4.)*(x1**3.))))+((-6.*(P*((x**5.)*(x1**3.))))+(\
    3.*((x**6.)*(x1**3.))));
    aux7=(-13.*((P**5.)*(x*(x1**3.))))+((15.*((P**3.)*((x**3.)*(x1**3.))))\
    +aux6);
    aux8=(18.*((P**2)*((x**5.)*(x1**2))))+((-9.*(P*((x**6.)*(x1**2))))+((\
    6.*((P**6.)*(x1**3.)))+aux7));
    aux9=(-45.*((P**4.)*((x**3.)*(x1**2))))+((15.*((P**3.)*((x**4.)*(\
    x1**2))))+aux8);
    aux10=(8.*((P**7.)*(x1**2)))+((-13.*((P**6.)*(x*(x1**2))))+((26.*((\
    P**5.)*((x**2)*(x1**2))))+aux9));
    aux11=(-15.*((P**4.)*((x**4.)*x1)))+((-18.*((P**3.)*((x**5.)*x1)))+((\
    9.*((P**2)*((x**6.)*x1)))+aux10));
    aux12=(-16.*((P**7.)*(x*x1)))+((8.*((P**6.)*((x**2)*x1)))+((32.*((\
    P**5.)*((x**3.)*x1)))+aux11));
    aux13=(5.*((P**5.)*(x**4.)))+((6.*((P**4.)*(x**5.)))+((-3.*((P**3.)*(\
    x**6.)))+aux12));
    aux14=((P-x1)**-3.)*(((8.*((P**7.)*(x**2)))+((-16.*((P**6.)*(x**3.)))+\
    aux13))*(theta*omega));
    aux15=(-672.*((P**2)*((x**5.)*(x1**3.))))+((1008.*(P*((x**6.)*(x1**3.)\
    )))+(-288.*((x**7.)*(x1**3.))));
    aux16=(1946.*((P**4.)*((x**3.)*(x1**3.))))+((-1365.*((P**3.)*((x**4.)*\
    (x1**3.))))+aux15);
    aux17=(864.*(P*((x**7.)*(x1**2))))+((421.*((P**7.)*(x1**3.)))+((-1050.\
    *((P**6.)*(x*(x1**3.))))+aux16));
    aux18=(2016.*((P**3.)*((x**5.)*(x1**2))))+((-3024.*((P**2)*((x**6.)*(\
    x1**2))))+aux17);
    aux19=(-5838.*((P**5.)*((x**3.)*(x1**2))))+((4095.*((P**4.)*((x**4.)*(\
    x1**2))))+aux18);
    aux20=(-1050.*((P**7.)*(x*(x1**2))))+((2100.*((P**6.)*((x**2)*(x1**2))\
    ))+aux19);
    aux21=(3024.*((P**3.)*((x**6.)*x1)))+((-864.*((P**2)*((x**7.)*x1)))+((\
    837.*((P**8.)*(x1**2)))+aux20));
    aux22=(-4095.*((P**5.)*((x**4.)*x1)))+((-2016.*((P**4.)*((x**5.)*x1)))\
    +aux21);
    aux23=(-1674.*((P**8.)*(x*x1)))+((837.*((P**7.)*((x**2)*x1)))+((4788.*\
    ((P**6.)*((x**3.)*x1)))+aux22));
    aux24=(672.*((P**5.)*(x**5.)))+((-1008.*((P**4.)*(x**6.)))+((288.*((\
    P**3.)*(x**7.)))+aux23));
    aux25=(837.*((P**8.)*(x**2)))+((-2154.*((P**7.)*(x**3.)))+((1365.*((\
    P**6.)*(x**4.)))+aux24));
    aux26=(0.00416667*((h**-4.)*((P**-6.)*aux14)))+(-0.0000496032*((h**-5.\
    )*((P**-7.)*(((P-x1)**-3.)*(aux25*((theta**2)*omega))))));
    output=aux5+aux26;
    return output



def get_curvature_case_4_rod_a_LHS_V(x, V, h, P, theta):
    aux0=(-6.*(h*(P*(1.+(P*(x*(-3.+(2.*x))))))))+((x+(P*(-2.+((5.*(P*((3.+\
    (-2.*x))*x)))+(4.*((-2.+x)*x))))))*theta);
    aux1=(2.*(h*((-2.+((6.*P)+(x+(-4.*(P*x)))))*theta)))+((-2.+((3.*P)+\
    (x+(-2.*(P*x)))))*(theta**2));
    aux2=(3.*((h**2)*(P+(2.*(P*x)))))+((2.*(h*((((2.*P)+(4.*(P*x)))-x)*\
    theta)))+(((P+(2.*(P*x)))-x)*(theta**2)));
    aux3=((-3.*((h**2)*P))+((2.*(h*(((-2.*P)+x)*theta)))+((x-P)*(\
    theta**2))))*(np.log(((h+theta)-((x*theta)/P))));
    aux4=((x**2)*(((3.*((h**2)*(P*(3.+(-2.*x)))))+aux1)*(np.log(h))))+((((\
    (-1.+x)**2))*(aux2*(np.log((h+theta)))))+aux3);
    output=(0.5*((P**-2.)*(V*((theta**-4.)*((x*(theta*aux0))+(2.*(P*\
    aux4)))))))/((2.*h)+theta);
    return output

def get_curvature_case_4_rod_a_LHS_Om(x, omega, h, P, theta):
    aux0=(h**2)*((P**2)*((x**2)*((-3.+(2.*x))*((((h+theta)**2))*(((np.log(h))**2))))));
    aux1=(P**2)*((((-1.+x)**2))*((1.+(2.*x))*((((h+theta)**2))*(((np.\
    log((h+theta)))**2)))));
    aux2=(-12.*(P*((((-1.+x)**2))*x)))+((x**2)+((P**2)*(19.+(20.*((x**2)*(\
    -3.+(2.*x)))))));
    aux3=(-2.*(P*(x*(3.+(2.*((-2.+x)*x))))))+((P**2)*(5.+(4.*((x**2)*(-3.+\
    (2.*x))))));
    aux4=(2.*((h**2)*(P*(x+(P*(7.+(8.*((x**2)*(-3.+(2.*x))))))))))+((h*(\
    aux2*theta))+(((x**2)+aux3)*(theta**2)));
    aux5=(np.log((h+theta)))*((theta*aux4)+(2.*((h**2)*((P**2)*((h+\
    theta)*(np.log(((h+theta)-((x*theta)/P)))))))));
    aux6=h*(((-5.*x)+(P*(19.+(2.*(x*(12.+((-39.*P)+((-6.*x)+(26.*(P*x)))))\
    )))))*theta);
    aux7=(2.*aux6)+(((-7.*x)+(P*(10.+(x*(40.+((-69.*P)+((-20.*x)+(46.*(P*\
    x)))))))))*(theta**2));
    aux8=(2.*(h*(((12.*(P**2))+((-9.*(P*x))+(x**2)))*(theta**2))))+(((\
    5.*(P**2))+((-6.*(P*x))+(x**2)))*(theta**3.));
    aux9=((14.*((h**3.)*(P**2)))+((3.*((h**2)*(P*(((11.*P)+(-4.*x))*\
    theta))))+aux8))*(np.log(((h+theta)-((x*theta)/P))));
    aux10=theta*((x*(theta*((28.*((h**2)*(P*(1.+(P*(x*(-3.+(2.*x))))\
    ))))+aux7)))+(2.*aux9));
    aux11=(8.*(h*(P*(x*((4.+((-9.*P)+((-2.*x)+(6.*(P*x)))))*(theta**2))\
    ))))+(4.*(P*(x*(((2.+((-3.*P)+(2.*(P*x))))-x)*(theta**3.)))));
    aux12=((h**2)*((x+(2.*(P*(1.+(6.*(x*((2.+((-9.*P)+(6.*(P*x))))-x))))))\
    )*theta))+aux11;
    aux13=((-1.+((6.*(x**2))+(-4.*(x**3.))))*(np.log((h+theta))))+(np.\
    log(((h+theta)-((x*theta)/P))));
    aux14=(x*(theta*((2.*((h**3.)*(P*(1.+(8.*(P*(x*(-3.+(2.*x)))))))))+\
    aux12)))+(2.*((h**2)*((P**2)*((((h+theta)**2))*aux13))));
    aux15=(4.*((h**2)*aux1))+((-2.*((h+theta)*aux5))+(aux10+(2.*((np.\
    log(h))*aux14))));
    output=(0.125*((P**-2.)*((theta**-6.)*(omega*((4.*aux0)+aux15)))\
    ))/((2.*h)+theta);
    return output

def get_curvature_case_4_rod_a_RHS_V(x, V, h, P, theta):
    aux0=(3.*((h**2)*(P+(2.*(P*x)))))+((2.*(h*(((P*(2.+(4.*x)))-x)*\
    theta)))+(((P+(2.*(P*x)))-x)*(theta**2)));
    aux1=(theta*((6.*(h*(P+(2.*(P*x)))))+((-4.*(x*theta))+(5.*(P*((\
    1.+(2.*x))*theta))))))+(2.*(aux0*((np.log(h))-(np.log((h+theta))\
    ))));
    output=((-0.5*(V*((((-1.+x)**2))*((theta**-4.)*aux1))))/((2.*h)+\
    theta))/P;
    return output

def get_curvature_case_4_rod_a_RHS_Om(x, omega, h, P, theta):
    aux0=(4.*(h*(((-6.*x)+(13.*(P*(1.+(2.*x)))))*theta)))+(((-20.*x)+(\
    23.*(P*(1.+(2.*x)))))*(theta**2));
    aux1=(4.*((h**2)*(P+(2.*(P*x)))))+((h*(((-3.*x)+(5.*(P*(1.+(2.*x)))))*\
    theta))+(((P+(2.*(P*x)))-x)*(theta**2)));
    aux2=(2.*(theta*aux1))+((h**2)*(P*((1.+(2.*x))*((h+theta)*((np.\
    log(h))-(np.log((h+theta))))))));
    aux3=((theta**2)*(((28.*((h**2)*(P+(2.*(P*x)))))+aux0)*omega))+(\
    4.*((h+theta)*(omega*(aux2*((np.log(h))-(np.log((h+theta)))))\
    )));
    output=((0.125*((((-1.+x)**2))*((theta**-6.)*aux3)))/((2.*h)+\
    theta))/P;
    return output

def get_curvature_case_4_rod_b_LHS_V(x, V, h, x1, P, theta):
    aux0=(2.*(h*(((4.*(P**2))+((-3.*(P*x))+((-2.*(P*x1))+(x*x1))))*\
    theta)))+((P-x)*((P-x1)*(theta**2)));
    aux1=((3.*((h**2)*(P*(((3.*P)+(-2.*x))-x1))))+aux0)*((np.log(h))-(np.\
    log((h+theta))));
    aux2=(((7.*(P**2))+((-6.*(P*x))+((-5.*(P*x1))+(4.*(x*x1)))))*(theta\
    **2))+(2.*aux1);
    aux3=((-1.+(x1/P))**-3.)*((theta**-4.)*((6.*(h*(P*((((3.*P)+(-2.*x)\
    )-x1)*theta))))+aux2));
    output=(-0.5*((P**-4.)*(V*((((x-x1)**2))*aux3))))/((2.*h)+theta);

    return output

def get_curvature_case_4_rod_b_LHS_Om(x, omega, h, x1, P, theta):
    aux0=(4.*(h*(((27.*(P**2))+((-20.*(P*x))+((-13.*(P*x1))+(6.*(x*x1)))))\
    *theta)))+(((29.*(P**2))+((-26.*(P*x))+((-23.*(P*x1))+(20.*(x*x1)))\
    ))*(theta**2));
    aux1=(h*(((9.*(P**2))+((-7.*(P*x))+((-5.*(P*x1))+(3.*(x*x1)))))*\
    theta))+((P-x)*((P-x1)*(theta**2)));
    aux2=(h**2)*(P*((((3.*P)+(-2.*x))-x1)*((h+theta)*((np.log(h))-(np.\
    log((h+theta)))))));
    aux3=((2.*(theta*((4.*((h**2)*(P*(((3.*P)+(-2.*x))-x1))))+aux1)))+\
    aux2)*((np.log(h))-(np.log((h+theta))));
    aux4=((theta**2)*(((28.*((h**2)*(P*(((3.*P)+(-2.*x))-x1))))+aux0)*\
    omega))+(4.*((h+theta)*(omega*aux3)));
    output=(0.125*((P**-4.)*((((x-x1)**2))*(((-1.+(x1/P))**-3.)*((theta\
    **-6.)*aux4)))))/((2.*h)+theta);
    return output

def get_curvature_case_4_rod_b_RHS_V(x, V, h, x1, P, theta):
    aux0=h*(P*((P*((P+(-2.*x))*x))+((3.*(P*(x*x1)))+((-3.*(P*(x1**2)))+(\
    x1**3.)))));
    aux1=(-4.*(P*(x*(((2.*P)+x)*x1))))+((P*(((7.*P)+(5.*x))*(x1**2)))+(((-\
    5.*P)+x)*(x1**3.)));
    aux2=(2.*(h*(((4.*(P**2))+((-3.*(P*x))+((-2.*(P*x1))+(x*x1))))*\
    theta)))+((P-x)*((P-x1)*(theta**2)));
    aux3=(((x-x1)**2))*(((3.*((h**2)*(P*(((3.*P)+(-2.*x))-x1))))+aux2)*(\
    np.log(h)));
    aux4=(2.*(h*((((P*((2.*P)+(3.*x)))+(2.*(x1**2)))-(((6.*P)+x)*x1))*\
    theta)))+((P+(x+(-2.*x1)))*((P-x1)*(theta**2)));
    aux5=(((P-x)**2))*(((3.*((h**2)*(P*(P+((2.*x)+(-3.*x1))))))+aux4)*(np.\
    log((h+theta))));
    aux6=((3.*((h**2)*P))+((2.*(h*(((2.*P)-x)*theta)))+((P-x)*(theta\
    **2))))*(np.log(((h+theta)-((x*theta)/P))));
    aux7=((P-x)*(theta*((-6.*aux0)+(((-2.*((P**2)*((P+(-3.*x))*x)))+\
    aux1)*theta))))+(2.*(P*((aux3+aux5)-(((P-x1)**3.)*aux6))));
    output=(0.5*((P**-2.)*(V*(((P-x1)**-3.)*((theta**-4.)*aux7)))))/((\
    2.*h)+theta);
    return output

def get_curvature_case_4_rod_b_RHS_Om(x, omega, h, x1, P, theta):
    aux0=(((3.*P)+(-2.*x))-x1)*((((x-x1)**2))*((((h+theta)**2))*(((np.\
    log(h))**2))));
    aux1=(((P-x)**2))*((P+((2.*x)+(-3.*x1)))*((((h+theta)**2))*(((np.\
    log((h+theta)))**2))));
    aux2=(P*((7.*(P**2))+((8.*(P*x))+(-16.*(x**2)))))+((3.*(P*(((-7.*P)+(\
    8.*x))*x1)))+((-3.*(P*(x1**2)))+(x1**3.)));
    aux3=(3.*(P*(((-19.*(P**2))+((17.*(P*x))+(4.*(x**2))))*x1)))+((3.*(P*(\
    ((7.*P)+(-9.*x))*(x1**2))))+((P+x)*(x1**3.)));
    aux4=(P-x)*((P-x1)*(((5.*(P**2))+((4.*(P*x))+((-10.*(P*x1))+(x1**2))))\
    *(theta**2)));
    aux5=(2.*((h**2)*(P*aux2)))+((h*((((P**2)*((19.*(P**2))+((7.*(P*x))+(-\
    28.*(x**2)))))+aux3)*theta))+aux4);
    aux6=2.*((h**2)*((P**2)*(((P-x1)**3.)*((h+theta)*(np.log(((h+\
    theta)-((x*theta)/P))))))));
    aux7=(h**2)*(P*((P*((P+(-2.*x))*x))+((3.*(P*(x*x1)))+((-3.*(P*(x1**2))\
    )+(x1**3.)))));
    aux8=(3.*(P*(x*(((17.*P)+(4.*x))*x1))))+((-9.*(P*(((6.*P)+x)*(x1**2)))\
    )+(((26.*P)+(-5.*x))*(x1**3.)));
    aux9=(-4.*(P*(x*(((7.*P)+(5.*x))*x1))))+((P*(((29.*P)+(19.*x))*(x1**2)\
    ))+(((-23.*P)+(7.*x))*(x1**3.)));
    aux10=(-2.*(h*((((P**2)*(((19.*P)+(-40.*x))*x))+aux8)*theta)))+(((\
    2.*((P**2)*(x*((-5.*P)+(13.*x)))))+aux9)*(theta**2));
    aux11=(2.*(h*(((12.*(P**2))+((-9.*(P*x))+(x**2)))*(theta**2))))+(((\
    5.*(P**2))+((-6.*(P*x))+(x**2)))*(theta**3.));
    aux12=((14.*((h**3.)*(P**2)))+((3.*((h**2)*(P*(((11.*P)+(-4.*x))*\
    theta))))+aux11))*(np.log(((h+theta)-((x*theta)/P))));
    aux13=(3.*(P*(((15.*P)+(-8.*x))*(x*x1))))+((3.*(P*(((-8.*P)+x)*(x1**2)\
    )))+(((8.*P)-x)*(x1**3.)));
    aux14=(3.*(P*(((28.*(P**2))+((2.*(P*x))+(-9.*(x**2))))*(x1**2))))+(((-\
    36.*(P**2))+((14.*(P*x))+(x**2)))*(x1**3.));
    aux15=((P**2)*(x*((-2.*(P**2))+((83.*(P*x))+(-60.*(x**2))))))+((3.*(P*\
    (x*(((-54.*(P**2))+((29.*(P*x))+(4.*(x**2))))*x1))))+aux14);
    aux16=(((x-x1)**2))*(((5.*(P**2))+((-4.*(P*x))+((-3.*(P*x1))+(2.*(x*\
    x1)))))*(theta**2));
    aux17=((h**2)*(aux15*theta))+((8.*(h*(P*aux16)))+(4.*(P*((P-x)*((P-\
    x1)*((((x-x1)**2))*(theta**3.)))))));
    aux18=(-2.*((h**3.)*(P*((P*(x*((P**2)+((-24.*(P*x))+(16.*(x**2))))))+\
    aux13))))+aux17;
    aux19=(P**2)+((2.*(P*x))+((-2.*(x**2))+((-4.*(P*x1))+((2.*(x*x1))+(x1**\
    2)))));
    aux20=((P+((-2.*x)+x1))*(aux19*(np.log((h+theta)))))-(((P-x1)**3.)*\
    (np.log(((h+theta)-((x*theta)/P)))));
    aux21=2.*((np.log(h))*((theta*aux18)+(2.*((h**2)*((P**2)*((((h+\
    theta)**2))*aux20))))));
    aux22=(2.*((h+theta)*((np.log((h+theta)))*(((P-x)*(theta*\
    aux5))+aux6))))+((theta*(((P-x)*(theta*((-28.*aux7)+aux10)))+(-\
    2.*(((P-x1)**3.)*aux12))))+aux21);
    aux23=(theta**-6.)*(omega*((4.*((h**2)*((P**2)*aux0)))+((-4.*((\
    h**2)*((P**2)*aux1)))+aux22)));
    output=(-0.125*((P**-2.)*(((P-x1)**-3.)*aux23)))/((2.*h)+theta);
    return output


def func_get_all_curvatures_a_case4(x, h, theta, V, Omega, x1):
    # x1 = P here
    xLHS = x[x<=x1]
    xRHS = x[x>x1]
    
    # transition function locations
    theta_min = h/20
    theta_max = h/10

    if theta == 0:

        y2_derivLHSV = get_curvature_case_4_rod_a_LHS_V_small_theta(xLHS, V, h, x1, theta)
        y2_derivLHSOm = get_curvature_case_4_rod_a_LHS_Om_small_theta(xLHS, Omega, h, x1, theta)
        
        y2_derivRHSV = get_curvature_case_4_rod_a_RHS_V_small_theta(xRHS, V, h, x1, theta)
        y2_derivRHSOm = get_curvature_case_4_rod_a_RHS_Om_small_theta(xRHS, Omega, h, x1, theta)
        
    else:
        
        y2_derivLHSVa = get_curvature_case_4_rod_a_LHS_V_small_theta(xLHS, V, h, x1, theta)
        y2_derivLHSOma = get_curvature_case_4_rod_a_LHS_Om_small_theta(xLHS, Omega, h, x1, theta)
        
        y2_derivRHSVa = get_curvature_case_4_rod_a_RHS_V_small_theta(xRHS, V, h, x1, theta)
        y2_derivRHSOma = get_curvature_case_4_rod_a_RHS_Om_small_theta(xRHS, Omega, h, x1, theta)
        
        y2_derivLHSVb = get_curvature_case_4_rod_a_LHS_V(xLHS, V, h, x1, theta)
        y2_derivLHSOmb = get_curvature_case_4_rod_a_LHS_Om(xLHS, Omega, h, x1, theta)
        
        y2_derivRHSVb = get_curvature_case_4_rod_a_RHS_V(xRHS, V, h, x1, theta)
        y2_derivRHSOmb = get_curvature_case_4_rod_a_RHS_Om(xRHS, Omega, h, x1, theta)

        y2_derivLHSV = weighted_val(theta_min, theta_max, theta, y2_derivLHSVa, y2_derivLHSVb)
        y2_derivLHSOm = weighted_val(theta_min, theta_max, theta, y2_derivLHSOma, y2_derivLHSOmb)      
    
        y2_derivRHSV = weighted_val(theta_min, theta_max, theta, y2_derivRHSVa, y2_derivRHSVb)
        y2_derivRHSOm = weighted_val(theta_min, theta_max, theta, y2_derivRHSOma, y2_derivRHSOmb)    
        
    
    y2_deriv_p_V = np.concatenate((y2_derivLHSV, y2_derivRHSV))
    y2_deriv_p_Om = np.concatenate((y2_derivLHSOm, y2_derivRHSOm))
    
    
    y2_deriv = y2_deriv_p_V + y2_deriv_p_Om
    
    return y2_deriv, y2_deriv_p_V, y2_deriv_p_Om

def func_get_all_curvatures_b_case4(x, V, Omega, h, theta):
    xLHS = x[x<=0]
    xRHS = x[x>0]
    M = x[0]
    P = x[-1]
    
    # transition function locations
    theta_min = h/20
    theta_max = h/10

    if theta == 0:

        y2_derivLHSV = get_curvature_case_4_rod_b_LHS_V_small_theta(xLHS, V, h, M, P, theta)
        y2_derivLHSOm = get_curvature_case_4_rod_b_LHS_Om_small_theta(xLHS, Omega, h, M, P, theta)
        
        y2_derivRHSV = get_curvature_case_4_rod_b_RHS_V_small_theta(xRHS, V, h, M, P, theta)
        y2_derivRHSOm = get_curvature_case_4_rod_b_RHS_Om_small_theta(xRHS, Omega, h, M, P, theta)
        
    else:
        
        y2_derivLHSVa = get_curvature_case_4_rod_b_LHS_V_small_theta(xLHS, V, h, M, P, theta)
        y2_derivLHSOma = get_curvature_case_4_rod_b_LHS_Om_small_theta(xLHS, Omega, h, M, P, theta)
        
        y2_derivRHSVa = get_curvature_case_4_rod_b_RHS_V_small_theta(xRHS, V, h, M, P, theta)
        y2_derivRHSOma = get_curvature_case_4_rod_b_RHS_Om_small_theta(xRHS, Omega, h, M, P, theta)
        
        y2_derivLHSVb = get_curvature_case_4_rod_b_LHS_V(xLHS, V, h, M, P, theta)
        y2_derivLHSOmb = get_curvature_case_4_rod_b_LHS_Om(xLHS, Omega, h, M, P, theta)
        
        y2_derivRHSVb = get_curvature_case_4_rod_b_RHS_V(xRHS, V, h, M, P, theta)
        y2_derivRHSOmb = get_curvature_case_4_rod_b_RHS_Om(xRHS, Omega, h, M, P, theta)

        y2_derivLHSV = weighted_val(theta_min, theta_max, theta, y2_derivLHSVa, y2_derivLHSVb)
        y2_derivLHSOm = weighted_val(theta_min, theta_max, theta, y2_derivLHSOma, y2_derivLHSOmb)      
    
        y2_derivRHSV = weighted_val(theta_min, theta_max, theta, y2_derivRHSVa, y2_derivRHSVb)
        y2_derivRHSOm = weighted_val(theta_min, theta_max, theta, y2_derivRHSOma, y2_derivRHSOmb)    
        
    
    y2_deriv_p_V = np.concatenate((y2_derivLHSV, y2_derivRHSV))
    y2_deriv_p_Om = np.concatenate((y2_derivLHSOm, y2_derivRHSOm))
    
    y2_deriv = y2_deriv_p_V + y2_deriv_p_Om
    
    return y2_deriv, y2_deriv_p_V, y2_deriv_p_Om