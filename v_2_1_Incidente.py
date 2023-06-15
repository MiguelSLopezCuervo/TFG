# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 17:35:02 2023

@author: Miguel Sánchez

programa 1D vacío, condiciones ABC, ONDA PLANA, VACÍO

        ABC--S(no influye)---PW--VACÍO---S---ABC
        
VERSIÓN 2_1: CALCULA LA DFT DE LA ONDA INCIDENTE PARA UNAS CUANTAS FRECUENCIAS

PROGRAMA DEFINITIVO
"""

import numpy as np
import matplotlib.pyplot as plt
from DFT_Numeric_Integral_op import dftNI_op

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
Pos_S1=1; Pos_S2=grid-1
E_s1=[]; E_s2=[] #la s1 no se utiliza

#SIMULO EL VACÍO

thickness  =  25   
er         =  1    
mr         =  1
sigma_e    =  0  
sigma_m    =  0

#Esto no es necesario hacerlo con ifs en esta versión, pero lo dejo así porq es más general
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
#el E tiene un cero más porque por las condiciones MUR empiezo y acababo por E

for i in range(0, tsteps):
    #Grabo los campos en las sondas
    E_s1.append(E[Pos_S1])
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
    
    """
    #Representar la onda:
    if i%2==0:
        fig, ax = plt.subplots(2, 1)
        
        #gráfica en sí
        ax[0].plot(range(0, grid+2), E)
        ax[1].plot(range(0, grid+1), H)
        
        #Etiquetas eje y
        ax[0].set_ylabel('Campo E')
        ax[1].set_ylabel('Campo H')
        
        ax[0].set_ylim(bottom=-1.3, top=1.3)
        ax[1].set_ylim(bottom=-0.005, top=0.005)
        
        plt.show()
    """
#Cálculo de la frecuencia de las ondas
frec_max=1/(30*dt)
frec_min=frec_max*10**-3
df=(frec_max-frec_min)/300 #300 es aprox el nº de frec en las q elijo evaluar DFT
DFT_I=[]

for f in range(0, 301): #Estoy haciendo el cálculo en 300 frecs, por eso pongo 301
    DFT_I.append(dftNI_op(E_s2, frec_min+f*df, tsteps, dt))

#Escribir dichas frecuencias en fichero:
# Nombre del archivo donde se guardará el vector
archivo = "DFT_Incidente.txt"

# Abrir el archivo en modo escritura
with open(archivo, "w") as file:
    # Escribir cada número complejo del vector en una línea del archivo
    for numero in DFT_I:
        file.write(str(numero) + "\n")

