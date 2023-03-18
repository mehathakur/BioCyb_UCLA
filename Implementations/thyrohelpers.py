"""
Independent package to store all of our plotting functions to avoid making the project notebook too messy
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def plot_thyroid(tspan, sol, lims, title = "Figure", Methimazole = False, subs = False, data = None):
    if not subs:
        plt.figure(figsize=(8,6))
    plt.title(title)
    plt.plot(tspan, sol[:,0], label="MMI (mg/L)", color="r")
    plt.plot(tspan, sol[:,1], label="FT4 (pg/mL)", color="b")
    plt.plot(tspan, sol[:,2], label="Thyroid Size (mL)", color="g")
    plt.plot(tspan, sol[:,3], label="TRAb (U/mL)", color="#c926b9")
    if Methimazole:
        plt.plot(tspan, sol[:,4], label="Oral MMI (mg/L/day)",color="k")
    plt.legend()
    plt.ylim(0,lims[0])
    plt.xlim(0,lims[1])
    plt.xlabel("Time (days)")
    plt.ylabel("State Variables")

    # TODO: add scatter function to plot emperical data

def plot_dosing(sol, dose, title = "Figure", subs = False, data = None, showdose=False):
    tspan = np.linspace(0,dose[1][-1], dose[1][-1]+1)
    lims=[40,dose[1][-1]]
    if not subs:
        plt.figure(figsize=(8,6))
    plt.title(title)
    plt.plot(tspan, sol[:,0], label="MMI (mg/L)", color="r")
    plt.plot(tspan, sol[:,1], label="FT4 (pg/mL)", color="b")
    plt.plot(tspan, sol[:,2], label="Thyroid Size (mL)", color="g")
    plt.plot(tspan, sol[:,3], label="TRAb (U/mL)", color="#c926b9")
    if showdose:
        plt.plot(tspan, sol[:,4], label="Oral MMI (mg/L/day)",color="k")
    plt.legend(loc='upper right')
    plt.ylim(0,lims[0])
    plt.xlim(0,lims[1])
    plt.xlabel("Time (days)")
    plt.ylabel("State Variables")

def initialize_pandiyan(patient = None, optim = False):
    """
    Initializes the parameter values set out in the Pandian paper.
    """
    #TODO: will most likey need to normalize this. Check here first if plots look wonky
    k_1 =  8.374e-3
    k_2 = 3.3271   
    k_3 = 0.119  
    k_4 = 0.099021 
    k_5 = 1e6   
    k_6 = 0.001 
    k_7 = 0.875
    k_8 = 0.035
    k_a = 0.358068
    k_b = 1.5
    k_d = 0.05
    N = 0.833

    # k3, kd, N, k7, and kb vary
    if patient == 20:
        k_3 = 0.08535  
        k_7 = 0.2625
        k_b = 4.95
        k_d = 0.067
        N = 0.250
    elif patient == 31: 
        k_3 = 0.09075 
        k_7 = 0.0609
        k_b = 11.8
        k_d = 0.07
        N = 0.058
    elif patient == 55:
        k_3 = 0.09031 
        k_7 = 0.2177
        k_b = 4.09
        k_d = 0.081
        N = 0.207733
    elif patient == 70:
        k_3 = 0.11784 
        k_7 = 0.308
        k_b = 3.15
        k_d = 0.075
        N = 0.29333 

    p = [k_1,k_2,k_3,k_4,k_5,k_6,k_7,k_8,k_a,k_b,k_d,N]

    if optim:
        fixed = [k_1,k_2,k_4,k_5,k_6,k_8,k_a]
        unfixed = [k_3,k_7,k_b,k_d,N]
        return fixed, unfixed
    else:
        return p


def thyroid_sim(x, t, k1, k2, k3, k4, k5, k6, k7, k8, ka, kb, kd, N ):
    """Takes x vector and passes it through one step of differentials."""
    dx = lambda t,x,z,s,k1,k2,ka: s-(k1*z*x)/(ka+x)-k2*x
    dy = lambda t,y,z,w,k3,k4,kd: (k3*z*w)/(kd+w)-k4*y
    dz = lambda t,x,z,w,k5,k6,N: k5*(w/z-N)-k6*x*z
    dw = lambda t,x,w,k7,k8,kb: k7-(k7*x)/(kb+x)-k8*w

    dxdt = dx(t,x[0],x[2],x[4],k1,k2,ka)
    dydt = dy(t,x[1],x[2],x[3],k3,k4,kd)
    dzdt = dz(t,x[0],x[2],x[3],k5,k6,N)
    dwdt = dw(t,x[0],x[3],k7,k8,kb)
    dsdt = 0 # no MMI treatment
    x = np.array([dxdt,dydt,dzdt,dwdt,dsdt])

    return x

def use_odeint(x, t, p):
    k1, k2, k3, k4, k5, k6, k7, k8, ka, kb, kd, N = p
    return odeint(thyroid_sim, x, t, args=(k1, k2, k3, k4, k5, k6, k7, k8, ka, kb, kd, N))

def dose_sim(E0,p,dose):
    days = dose[1]
    dose = dose[0]
    sol = np.zeros((days[-1]+1,5))
    for i, day in enumerate(days):
        if i == 0:
            dayrange = np.linspace(0, day+1, day+1)
            sol[:day+1,:]  = use_odeint(E0, dayrange, p)
            EK = sol[day,:] # need to get the last values from plots
        else:
            dayrange = np.linspace(days[i-1], day+1, (day-days[i-1])+1)
            EK[-1] = (0.93 * dose[i] * (day-days[i-1]))/59.71
            sol[days[i-1]:day+1, :] = use_odeint(EK, dayrange, p)
            EK = sol[day,:]
    return sol