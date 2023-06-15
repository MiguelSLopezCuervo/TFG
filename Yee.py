# -*- coding: utf-8 -*-
"""
Created on Fri May 26 14:03:51 2023

@author: Miguel Sánchez

ALGORITMO DE YEE. SIMULACIÓN ELECTROMAGNÉTICA.

"""

import numpy as np
from Incidente import Incidente

e0 = 8.8541878176e-12
m0 = 1.2566370614359173e-6
c0 = 1/np.sqrt(e0*m0)
Z0 = 376.73031346177  #=m0*c0  EVITAR ESA CUENTA

grid   = 87
L      = 0.0043
ubic_mat = 40 # ARBITRARIAMENE ELIJO QUE EL PRIMER MATERIAL COMIENCE AHÍ
dx     = L/(grid-1)
Sc     = 1 #impongo que Sc en el vacío sea uno y con eso calculo el dt
dt     = Sc*dx/c0 
tsteps = 3500

#Cálculo de la frecuencia de las ondas
frec_max=Sc/(30*dt) #usar ppw puntos por lambda
frec_min=frec_max*10**-3
frec_med=(frec_max+frec_min)/2

# Creo los vectores mr, er, sigma_e y sigma_m
mr_vect=np.zeros(grid+2)
er_vect=np.zeros(grid+2)
sigma_e_vect=np.zeros(grid+2)
sigma_m_vect=np.zeros(grid+2)

# Sonda
Pos_S1=1

# Características de la pintura que no optimizaremos
thickness  =  15  # 0.75mm CADA CAPA DE PINTURA
sigma_m    =  0.0

# Función que servirá para calcular la DFT y pasar al dominio de la frecuencia
def dftNI_op(Field, frec, t_max, dt):
    n = len(Field)
    t = np.arange(n) * dt
    omega = 2 * np.pi * frec

    exponent = np.exp(-1j * omega * t)
    res = np.sum(Field * exponent) * dt

    return res


#Por simplicidad de código se crean los siguientes factores Se usarán en las ecuaciones de ev 
fact1_E=np.zeros(grid+2)
fact2_E=np.zeros(grid+2)
fact1_H=np.zeros(grid+2) 
fact2_H=np.zeros(grid+2)
  

#Lectura de la DFT_I
DFT_I=Incidente()

def Calculo_R(x):
    if x[0]<1 or x[1]<1 or x[2]<1 or x[3]<1 or x[4]<1 or x[5]<1 or x[6]<0 or x[6]>3 or x[7]<0 or x[7]>3 or x[8]<0 or x[8]>3:
        raise Exception("Error, las er y mr han de ser >=1 y las sigma_e entre 0 y 3")
        
    E_s1=[]
    E = np.zeros(grid+2)
    H = np.zeros(grid+1)  
    
    for i in range(0 , grid+2):
        if i >= ubic_mat and i < thickness+ubic_mat: 
            mr_vect[i]       =  x[1]
            er_vect[i]       =  x[0]
            sigma_e_vect[i]  =  x[6]
            sigma_m_vect[i]  =  sigma_m
        elif i >= thickness+ubic_mat and i < 2*thickness+ubic_mat:
            mr_vect[i]       =  x[3]
            er_vect[i]       =  x[2]
            sigma_e_vect[i]  =  x[7]
            sigma_m_vect[i]  =  sigma_m
        elif i >= 2*thickness+ubic_mat and i < 3*thickness+ubic_mat:
            mr_vect[i]       =  x[5]
            er_vect[i]       =  x[4]
            sigma_e_vect[i]  =  x[8]
            sigma_m_vect[i]  =  sigma_m
        else:
            mr_vect[i]       =  1
            er_vect[i]       =  1
            sigma_e_vect[i]  =  0
            sigma_m_vect[i]  =  0
    
    # Utilizo la multiplicación vectorial de numpy para ahorrar tiempo
    fact1_E = sigma_e_vect * dt / (2 * er_vect * e0)
    fact2_E = dt / (dx * er_vect * e0)
    fact1_H = sigma_m_vect * dt / (2 * mr_vect * m0)
    fact2_H = dt / (dx * mr_vect * m0)
    
    for i in range(0, tsteps):
        #Grabo los campos en las sondas
        E_s1.append(E[Pos_S1])
    
        #Puntos interiores campo H
        for j in range(0, grid+1):
            H[j]=(1-fact1_H[j])/(1+fact1_H[j])*H[j]+fact2_H[j]/(1+fact1_H[j])*(E[j+1]-E[j]) 
        #Le añado la condición onda plana
        H[2] -= np.exp(-((i-30)/10)**2) / Z0
            
        #Condiciones de contorno tipo ABC, las más sencillas porque Sc=1 (vacío)
        E[0]=E[1]
        E[grid+1]=E[grid]
        
        #Puntos interiores campo E
        for j in range(1, grid+1):
            E[j]=(1-fact1_E[j])/(1+fact1_E[j])*E[j]+fact2_E[j]/(1+fact1_E[j])*(H[j]-H[j-1])
        #Le añado la condición onda plana
        E[3] += np.exp(-(i+ 0.5-(-0.5)-30)**2 / 100)
        
        #Condiciones tipo PEC después de las capas de pintura
        E[ubic_mat+3*thickness]=0
    
    DFT_R = dftNI_op(E_s1, frec_med, tsteps, dt)
    
    #Calculo R y T para cada frecuencia:
    R = np.abs(DFT_R) / np.abs(DFT_I)
        
    return R