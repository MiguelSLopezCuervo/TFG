# -*- coding: utf-8 -*-
"""
Created on Sun May 14 12:43:11 2023

@author: Miguel Sánchez

v_3_1

    R en 3 capas de pintura.
    
        ABC--S--PW---111122223333PEC
"""

"""
EJECUTAR PRIMERO EL PROGRAMA INCIDENTE_2_1 PARA GENERAR EL ARCHIVO DE
TEXTO
"""

import numpy as np
import matplotlib.pyplot as plt
from Incidente_Definitivo import Incidente_def, dftNI_op

e0 = 8.8541878176e-12
m0 = 1.2566370614359173e-6
c0 = 1/np.sqrt(e0*m0)
Z0 = 376.73031346177

grid   = 87
L      = 0.0043
ubic_mat = 40 # ARBITRARIAMENE ELIJO QUE EL PRIMER MATERIAL COMIENCE AHÍ
dx     = L/(grid-1)
Sc     = 1 #impongo que Sc en el vacío sea uno y con eso calculo el dt
dt     = Sc*dx/c0 
tsteps = 3500

# Sonda
Pos_S1=1
E_s1=[]

# Características del material que sí optimizaremos
q = [5.50195312 ,1.53710938 ,4.19335938 ,5.32617188 ,5.36523438 ,4.89648438] #añadirle las sigmas si se desea (modificar entonces los ifs del for de la 60)
er_1         =  q[0]      #Esto y mr tendrán que ser mayor que 1
mr_1         =  q[1]
er_2         =  q[2]  #Esto y mr tendrán que ser mayor que 1
mr_2         =  q[3]
er_3         =  q[4]  #Esto y mr tendrán que ser mayor que 1
mr_3         =  q[5] 

# Características de la pintura que no optimizaremos
thickness  =  15 # 0.75mm CADA CAPA DE PINTURA
sigma_e    =  0.1  #Esto y sigma_m tendrán que ser mayor que 0
sigma_m    =  0.0

# Creo los vectores mr, er, sigma_e y sigma_m
mr_vect=np.zeros(grid+2)
er_vect=np.zeros(grid+2)
sigma_e_vect=np.zeros(grid+2)
sigma_m_vect=np.zeros(grid+2) #sobrará un número pero dejarlo así para hacer 1 solo for

for i in range(0 , grid+2):
    if i >= ubic_mat and i < thickness+ubic_mat: 
        mr_vect[i]       =  mr_1
        er_vect[i]       =  er_1
        sigma_e_vect[i]  =  sigma_e
        sigma_m_vect[i]  =  sigma_m
    elif i >= thickness+ubic_mat and i < 2*thickness+ubic_mat:
        mr_vect[i]       =  mr_2
        er_vect[i]       =  er_2
        sigma_e_vect[i]  =  sigma_e
        sigma_m_vect[i]  =  sigma_m
    elif i >= 2*thickness+ubic_mat and i < 3*thickness+ubic_mat:
        mr_vect[i]       =  mr_3
        er_vect[i]       =  er_3
        sigma_e_vect[i]  =  sigma_e
        sigma_m_vect[i]  =  sigma_m
    else:
        mr_vect[i]       =  1
        er_vect[i]       =  1
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

    #Representar la onda:
    if i%1==0 and i<500: #cambiar eso según si se quiere más o menos representacion
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

#Cálculo de la frecuencia de las ondas

frec_max=Sc/(30*dt) #usar ppw puntos por lambda
frec_min=frec_max*10**-3
frec_med=(frec_max+frec_min)/2
df=(frec_max-frec_min)/300

DFT_R=[]
#Leo la DFT_I y calculo las DFT_T y DFT_R
archivo = "DFT_Incidente.txt"
# Abrir el archivo en modo lectura
with open(archivo, "r") as file:
    # Leer cada línea del archivo y convertirla en un número complejo
    DFT_I = [complex(line.strip()) for line in file]

for f in range(0, 301): 
    DFT_R.append(dftNI_op(E_s1, frec_min+f*df, tsteps, dt))

#Calculo R para cada frecuencia:
R=[]
for i in range(0, 301): 
    R.append( np.abs(DFT_R[i]) / np.abs(DFT_I[i]) )

#Representar en no logaritmica
# Crear la figura y los ejes
fig, ax = plt.subplots()

# Graficar los vectores R y T en el mismo eje
ax.plot(np.linspace(frec_min, frec_max, num=len(R)), R, label='R')

# Establecer las etiquetas de los ejes y el título de la figura
ax.set_xlabel('Frecuencia')
ax.set_ylabel('Valor de R')
ax.set_title('Coeficiente R según frecuencia')

# Mostrar la leyenda
ax.legend()

# Mostrar la gráfica
plt.show()

#Calcular R en frecuencia media
DFT_I_fmed = Incidente_def()
DFT_R_fmed = dftNI_op(E_s1, frec_med, tsteps, dt)
R_med = np.abs(DFT_R_fmed) / np.abs(DFT_I_fmed)
R_med_dB = 20*np.log10(R_med)

print("R en f_med:", R_med)
print("R_dB en f_med:", R_med_dB)

# Representar R en log
# Crear la figura y los ejes
fig, ax = plt.subplots()

# Graficar los vectores R y T en el mismo eje
ax.plot(np.linspace(frec_min, frec_max, num=len(R)), 20 * np.log10(R), label='R')

# Establecer las etiquetas de los ejes y el título de la figura
ax.set_xscale('log')  # Cambiar el eje X a escala logarítmica
ax.set_xlabel('Frecuencia (Escala log)')
ax.set_ylabel('20 * log(R)')  # Cambiar el eje Y a 20*log
ax.set_title('Coeficiente R según frecuencia')

# Mostrar la leyenda
ax.legend()

# Mostrar la gráfica
plt.show()