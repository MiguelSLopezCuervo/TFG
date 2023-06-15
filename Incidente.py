# -*- coding: utf-8 -*-
"""
Created on Fri May 26 14:03:51 2023

@author: Miguel Sánchez

CALCULO DE LA ONDA INCIDENTE
"""

import numpy as np

def dftNI_op(Field, frec, t_max, dt):
    n = len(Field)
    t = np.arange(n) * dt
    omega = 2 * np.pi * frec

    exponent = np.exp(-1j * omega * t)
    res = np.sum(Field * exponent) * dt

    return res


def Incidente():
    
    e0 = 8.8541878176e-12
    m0 = 1.2566370614359173e-6
    c0 = 1/np.sqrt(e0*m0)
    Z0 = 376.73031346177  #=m0*c0  EVITAR ESA CUENTA
    
    grid   = 87
    L      = 0.0043
    ubic_mat = 40 # OJO ARBITRARIAMENTE HE ELEGIDO QUE EL MATERIAL COMIENCE AHÍ
    dx     = L/(grid-1)
    Sc     = 1 #impongo que Sc en el vacío sea uno y con eso calculo el dt
    dt     = Sc*dx/c0 
    tsteps = 3500
    
    #Sonda
    Pos_S2=grid-1
    E_s2=[] 
    
    #SIMULO EL VACÍO
    thickness  =  25   
    er         =  1    
    mr         =  1
    sigma_e    =  0  
    sigma_m    =  0
    
    #Esto no es necesario hacerlo en esta versión, pero lo dejo así porq es más general y solo se llama una vez esta funcion
    mr_vect=np.zeros(grid+2)
    er_vect=np.zeros(grid+2)
    sigma_e_vect=np.zeros(grid+2)
    sigma_m_vect=np.zeros(grid+2) #sobrará un número pero dejarlo así para hacer 1 solo for
    for i in range(0 , grid+2):
        if i >= ubic_mat and i < thickness+ubic_mat: 
            mr_vect[i]       =  mr
            er_vect[i]       =  er
            sigma_e_vect[i]  =  sigma_e
            sigma_m_vect[i]  =  sigma_m
        else:
            mr_vect[i] =  1  
            er_vect[i] =  1
            sigma_e_vect[i]  =  0
            sigma_m_vect[i]  =  0
    
    #Por simplicidad de código se crean los siguientes factores Se usarán en las ecuaciones de ev 
    fact1_E=np.zeros(grid+2)
    fact2_E=np.zeros(grid+2)
    fact1_H=np.zeros(grid+2) #sobrará un número al final pero dejarlo así por hacer solo 1 for
    fact2_H=np.zeros(grid+2)
    for i in range(0, grid+2):
        fact1_E[i] = sigma_e_vect[i]*dt/(2*er_vect[i]*e0)
        fact2_E[i] = dt/(dx*er_vect[i]*e0)
        fact1_H[i] = sigma_m_vect[i]*dt/(2*mr_vect[i]*m0)
        fact2_H[i] = dt/(dx*mr_vect[i]*m0)
    
    E = np.zeros(grid+2) #Ojo: crea grid+2 ceros en total Vector dim grid+1 empezando en 0
    H = np.zeros(grid+1) #Vector dim grid empezando en 0
    #el E tiene un cero más porque por las condiciones ABC empiezo y acababo por E
    
    for i in range(0, tsteps):
        #Grabo los campos en la sonda
        E_s2.append(E[Pos_S2])
    
        #Puntos interiores campo H
        for j in range(0, grid+1):
            H[j]=(1-fact1_H[j])/(1+fact1_H[j])*H[j]+fact2_H[j]/(1+fact1_H[j])*(E[j+1]-E[j]) 
        #Le añado la condición onda plana
        H[9] -= np.exp(-((i-30)/10)**2) / Z0
            
        #Condiciones de contorno tipo ABC, las más sencillas porque Sc=1 (vacío)
    
        E[0]=E[1]
        E[grid+1]=E[grid]
        
        #Puntos interiores campo E
        for j in range(1, grid+1):
            E[j]=(1-fact1_E[j])/(1+fact1_E[j])*E[j]+fact2_E[j]/(1+fact1_E[j])*(H[j]-H[j-1])
        #Le añado la condición onda plana
        E[10] += np.exp(-(i+ 0.5-(-0.5)-30)**2 / 100)
        
    #Cálculo de la frecuencia de las ondas
    frec_max = 1/(30*dt)
    frec_min = frec_max*10**-3
    frec_med = (frec_max + frec_min) / 2
    DFT_I = dftNI_op(E_s2, frec_med, tsteps, dt)
    
    return DFT_I
    
