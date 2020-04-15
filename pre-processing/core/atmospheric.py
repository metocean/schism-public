#!/usr/bin/env python2.7
from math import log10,exp
import numpy as np
from numpy import pi
#saturation pressure

Mw=18.0160 # molecular weight of water
Md=28.9660 # molecular weight of dry air
DEG2RADS    = 0.01745329252

def get_convergence(Lat,Lon,inputEPSG):
    if inputEPSG==2193:
        
        lat = Lat * DEG2RADS
        lon = Lon * DEG2RADS
        alpha=lon
        alpha0=173*DEG2RADS
        teta=lat
        FLAT=0.0034
        a= 6378137
        e2  = 2*FLAT - FLAT*FLAT
        p=a*(1-e2)/(1-e2*np.sin(teta)**2)**(3/2)
        v=a/(1-e2*np.sin(teta)**2)**0.5
        F=v/p;

        W=alpha-alpha0
        t=np.tan(teta)

        Term1=-W*np.sin(teta)
        Term2=-(W**3/3)*np.sin(teta)*np.cos(teta)**2*(2*F**2-F)
        Term3=-(W**5/15)*np.sin(teta)*np.cos(teta)**4*(F**4*(11-24*t**2)-F**3*(11-36*t**2)+2*F**2*(1-7*t**2)+F*t**2)
        Term4=-(W**7/315)*np.sin(teta)*np.cos(teta)**6*(17-26*t**2+2*t**4)

        A=Term1+Term2+Term3+Term4;

    elif inputEPSG>32600 and inputEPSG<32761:
        lat = Lat * DEG2RADS
        lon = Lon * DEG2RADS
        zone=np.arange(1,61,1)
        lon0=np.arange(-177,177+6,6)
        if inputEPSG>32700: # South hemishpere
            utm=inputEPSG-32700
        else:
            utm=inputEPSG-32600

        lon0=lon0[utm-1]


        A=np.arctan(np.tan(lon-lon0)*np.sin(lat))

    else:
        A=Lon*0.






    return A*180/pi

def esat(T):
    ''' get sateration pressure (units [Pa]) for a given air temperature (units [K])'''
    from numpy import log10
    TK = 273.15
    e1 = 101325.0
    logTTK = log10(T/TK)
    esat =  e1*10**(10.79586*(1-TK/T)-5.02808*logTTK+ 1.50474*1e-4*(1.-10**(-8.29692*(T/TK-1)))+ 0.42873*1e-3*(10**(4.76955*(1-TK/T))-1)-2.2195983) 
    return esat

def mixr2sh(W):
    '''conversion from mixing ratio (units [kg/kg]) to specific humidity (units also [kg/kg])
    '''
    return W/(1.+W)

def rh2mixr(RH,p,T):
    '''purpose: conversion relative humidity (unitless) to mixing ratio [kg/kg]'''
    es = esat(T)
    return Mw/Md*RH*es/(p-RH*es)

def rh2sh(RH,p,T):
    '''conversion relative humidity to specific humidity (kg Water per m^3 Air)'''
    mixr=rh2mixr(RH, p,T)
    sh=mixr2sh(mixr)
    return sh

if __name__ == "__main__":
    RH=90./100.
    p=100180
    T=290

    sh=rh2sh(RH,p,T)
