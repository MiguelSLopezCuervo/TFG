# -*- coding: utf-8 -*-
"""
Created on Mon May 22 13:19:11 2023

@author: Miguel Sánchez
"""

import numpy as np
import random

def Decodificar(Cromosoma, m, x_max, x_min):
    # Obtener el número de variables
    n_variables = len(m)
    
    #Pasar el cromosoma de vector a string
    Cromosoma_s=''.join([str(bit) for bit in Cromosoma])
    
    # Inicializar el vector de variables decimales
    variables_dec = np.zeros(n_variables)
    
    # Convertir cada subcadena binaria en un número decimal
    inicio = 0
    for i in range(n_variables):
        fin = inicio + m[i]
        variable_bin = Cromosoma_s[inicio:fin]
        variables_dec[i] = int(variable_bin, 2)/(2**m[i])
        inicio = fin
        
    #Tomar el valor medio, no el del redondeo hacia abajo
    for i in range(len(m)):
        variables_dec[i]+=1/(2*(2**m[i]))
        
    #Transformar de la variable normalizada a la no normalizada
    for i in range(len(m)):
        variables_dec[i]=variables_dec[i]*(x_max[i]-x_min[i])+x_min[i]
        
    return variables_dec

def Mating(Cromosoma_1, Cromosoma_2, m):
    Desc_1=[]
    Desc_2=[]
    
    contador_gen=0
    for i in range(0, len(m)):
        ran = random.randint(1, m[i]-1)
        
        for j in range(0,m[i]):
            if j<ran:
                Desc_1.append( Cromosoma_1[contador_gen + j] )
                Desc_2.append( Cromosoma_2[contador_gen + j] )
            else:
                Desc_1.append( Cromosoma_2[contador_gen + j] )
                Desc_2.append( Cromosoma_1[contador_gen + j] )
        contador_gen += m[i]
    
    return Desc_1, Desc_2

def Mutar(Poblacion, in_mutar, m_rate):
    pop_cromosomas2 = np.copy(Poblacion)
    
    for i in range(in_mutar[0]-1, in_mutar[1]-1): #No altero a los in_mutar[0]-1 mejores individuos
        for j in range(0, pop_cromosomas2.shape[1]):
            al = random.random()
            if al < m_rate[0]:
                if pop_cromosomas2[i][j] == 1: 
                    pop_cromosomas2[i][j] = 0
                if pop_cromosomas2[i][j] == 0: 
                    pop_cromosomas2[i][j] = 1
    
    for i in range(in_mutar[1]-1, in_mutar[2]-1): # Altero desde el in_mutar[1]-1
        for j in range(0, pop_cromosomas2.shape[1]):
            al = random.random()
            if al < m_rate[1]:
                if pop_cromosomas2[i][j] == 1: 
                    pop_cromosomas2[i][j] = 0
                if pop_cromosomas2[i][j] == 0: 
                    pop_cromosomas2[i][j] = 1
                    
    for i in range(in_mutar[2]-1, in_mutar[3]-1): 
        for j in range(0, pop_cromosomas2.shape[1]):
            al = random.random()
            if al < m_rate[2]:
                if pop_cromosomas2[i][j] == 1: 
                    pop_cromosomas2[i][j] = 0
                if pop_cromosomas2[i][j] == 0: 
                    pop_cromosomas2[i][j] = 1
    
    for i in range(in_mutar[3]-1, pop_cromosomas2.shape[0]): 
        for j in range(0, pop_cromosomas2.shape[1]):
            al = random.random()
            if al < m_rate[3]:
                if pop_cromosomas2[i][j] == 1: 
                    pop_cromosomas2[i][j] = 0
                if pop_cromosomas2[i][j] == 0: 
                    pop_cromosomas2[i][j] = 1

    return pop_cromosomas2

def Ordenar_Poblacion(costes, pop_cromosomas):
    # Obtenemos los índices que ordenan el vector costes de menor a mayor
    indices_ordenados = np.argsort(costes)

    # Ordenamos el vector costes y la matriz pop_cromosomas según los índices obtenidos
    costes_ordenados = costes[indices_ordenados]
    pop_cromosomas_ordenada = pop_cromosomas[indices_ordenados, :]

    # Devolvemos el vector costes y la matriz pop_cromosomas ordenados
    return costes_ordenados, pop_cromosomas_ordenada

def Nuevos_Desc_Y_Padre(x1, x2, y, m):
    bit_random1   =   random.randint(0, len(x1)-1)
    bit_random2   =   random.randint(0, len(x2)-1)
    
    if x1[bit_random1] == 1:
        x1[bit_random1]=0
    else:
        x1[bit_random1]=1
        
    if x2[bit_random2] == 1:
        x2[bit_random2]=0
    else:
        x2[bit_random2]=1
    
    genes         =   len(m)
    gen_random_1  =   random.randint(0, genes-1)
    Variacion     =   random.randint(0, 1)

    if Variacion==0:
        gen_random_2=gen_random_1
        while gen_random_2==gen_random_1:
            gen_random_2=random.randint(0, genes-1)
        
        bits=0
        for i in range(0, gen_random_1+1):
            bits+=m[i]
        
        if y[bits-1]==1:
            y[bits-1]=0
        else:
            y[bits-1]=1
            
        bits=0
        for i in range(0, gen_random_2+1):
            bits+=m[i]
        
        if y[bits-1]==1:
            y[bits-1]=0
        else:
            y[bits-1]=1
    else:
        bits=0
        for i in range(0, gen_random_1+1):
            bits+=m[i]
        
        if y[bits-2]==1:
            y[bits-2]=0
        else:
            y[bits-2]=1
        
    return x1, x2, y

