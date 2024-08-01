# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 11:30:24 2023

@author: Veliy
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import xlwings as xw
import scipy.optimize as scopt


wb = xw.Book("PFR model.xlsx")
sheet = wb.sheets("Sheet1")

def optimization(parameters):
    [L,D,Tj] = parameters
    sheet.range("B8").value=L
    sheet.range("B9").value=D
    sheet.range("B19").value=Tj

    def odes(x,z):
        Cb = x[0]
        Cp = x[1]
        Ci = x[2]
        Cd = x[3]
        T = x[4]
        
        k01 = sheet.range("B3").value
        k02 = sheet.range("C3").value
        Ea1 = sheet.range("B4").value
        Ea2 = sheet.range("C4").value
        delH1 = sheet.range("B5").value
        delH2 = sheet.range("C5").value
        r1 = k01*np.exp(-Ea1/(8.314556*T))*pow(Cb,0.96)*pow(Cp,0.87)
        r2 = k02*np.exp(-Ea2/(8.314556*T))*pow(Ci,0.61)*pow(Cp,0.92)
        L = sheet.range("B8").value
        D = sheet.range("B9").value
        Ax = sheet.range("B10").value
        Ah = sheet.range("B11").value
        Q = sheet.range("B12").value
        v = sheet.range("B13").value
        Cpb = sheet.range("B14").value
        Cpp = sheet.range("B15").value
        Cpi = sheet.range("B16").value
        Cpd = sheet.range("B17").value
        U = sheet.range("B18").value
        Tj = sheet.range("B19").value
        
        dCbdz = -r1/v
        dCpdz = -(r1+r2)/v
        dCidz = (r1-r2)/v
        dCddz = r2/v
        dTdz = (Ax/((Cb*Q*Cpb)+(Cp*Q*Cpp)+(Ci*Q*Cpi)+(Cd*Q*Cpd)))*((r1*delH1)+(r2*delH2)-((4*U*(T-Tj))/D))
    
        return[dCbdz, dCpdz, dCidz, dCddz, dTdz]  
    Cb0 = sheet.range("B22").value
    Cp0 = sheet.range("B23").value
    Ci0 = sheet.range("B24").value
    Cd0 = sheet.range("B25").value
    T0 = sheet.range("B26").value
    
    x0 = [Cb0, Cp0, Ci0, Cd0, T0]
    zstart = sheet.range("J3").value
    zend = sheet.range("J303").value
    zpoints = len(sheet.range("J3:J303").value)
    z = np.linspace(zstart,zend,zpoints)
    sheet.range("J3:J303").value=z.reshape((zpoints,1))
    
    x = odeint(odes, x0, z)
    Cb = x[:,0]
    Cp = x[:,1]
    Ci = x[:,2]
    Cd = x[:,3]
    T = x[:,4]
    
    sheet.range("E3").value=Cb.reshape((zpoints,1))
    sheet.range("F3").value=Cp.reshape((zpoints,1))
    sheet.range("G3").value=Ci.reshape((zpoints,1))
    sheet.range("H3").value=Cd.reshape((zpoints,1))
    sheet.range("I3").value=T.reshape((zpoints,1))
    
    #plt.plot(z, Cb, 'r')
    #plt.plot(z, Cp, 'b')
    #plt.plot(z, Ci, 'o')
    #plt.plot(z, Cd, 'g')
    #plt.plot(z, T, 'y')
    
    Yield=Ci[-1]/(Cp0)
    sheet.range("M11").value = Yield
    return -Yield

scopt.minimize(optimization, [15,0.035,370])



    