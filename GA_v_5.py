# -*- coding: utf-8 -*-
"""
Created on Sat May 27 14:28:36 2023

@author: Miguel Sánchez

GA v_5
    
    ALGORITMO GENÉTICO QUE OPTIMIZA LA SIMULACIÓN ELECTROMAGNÉTICA
"""
import numpy as np
from Funciones_GA import Decodificar, Ordenar_Poblacion, Mating, Mutar, Nuevos_Desc_Y_Padre
import Yee
import time

#Medimos tiempo:
inicio = time.time()

#Parámetros del GA
n = 9 # Variables a optimizar
m = [8, 8, 8, 8, 8, 8, 7, 7, 7] # PRECISIÓN 0.02 cada variable Con esto suficienete
x_max = [6, 6, 6, 6, 6, 6, 3, 3, 3]
x_min = [1, 1, 1, 1, 1, 1, 0, 0, 0]
N_pop = 10000
N_Keep= 210
N_Desc = N_pop-N_Keep

bits=0 
for i in range(n):
    bits+=m[i]

m_rate=[1./bits, 2./bits, 3./bits, 20./bits]
in_mutar=[20, 200, 5000, 9500]  #Individuo a partir del cual queremos empezar a mutar

#Se calcula el número de apareos
if    N_Keep%2==0:  N_Ap = N_Keep**2 * 0.25 #N_Ap es el número de apareos
else:               N_Ap = (N_Keep**2/2 - 0.5)/2

#Si no son suficientes para reponer la población no se ejecuta el programa
if (2*N_Ap < N_Desc):
    raise Exception("Error, la descendencia tiene que llegar a reponer la población")
    
#Si son tantos que el programa irá lento, se hace que solo se reproduzcan los primeros
N_Rep=N_Keep
while (2*N_Ap > 8*N_pop):
    N_Rep=N_Rep // 2
    if     N_Rep%2==0:  N_Ap = N_Rep**2 * 0.25 
    else:               N_Ap = (N_Rep**2/2 - 0.5)/2

if len(m_rate) !=4 or len(in_mutar) !=4:
    raise Exception("Error: El modelo trabaja con 4 m_rate y 4 in_mutar")
for i in range(len(in_mutar)):
    if in_mutar[i]>(N_pop-1):
        raise Exception("Error, in_mutar ha de ser mayor que el tamaño de la población")

#Generar población inicial aleatoria
pop_cromosomas = np.random.randint(2, size=(N_pop, bits))
pop_variables = np.zeros((N_pop, n))
costes = np.zeros(N_pop)

# Archivo para guardar los resultados
file_name = "Coste_it.txt"
file_name_2 = "Mejores_individuos_it.txt"

# Bucle principal
it = 150
with open(file_name, 'w') as file, open(file_name_2, 'w') as file_2:
    for i in range(0, it):
        # Contamos el tiempo:
        t_it_in = time.time()
        # Decodificar los cromosomas
        for j in range(0, N_pop):
            pop_variables[j] = Decodificar(pop_cromosomas[j], m, x_max, x_min)
            
        # Evaluar los costes
        for j in range(0, N_pop):
            costes[j] = Yee.Calculo_R(pop_variables[j])
        
        # Ordenar la población según los costes (el menor coste el primero)
        costes, pop_cromosomas = Ordenar_Poblacion(costes, pop_cromosomas)
        
        # Descartar los que estén por debajo de N_Keep
        # Cruzar, a los N_Rep
        # Sólo se dará con exito la descendencia que reemplaza la población, ni más ni menos
        Ap_Exitosos = np.zeros(int(2*N_Ap), dtype=int)
        indices = np.random.choice(int(2*N_Ap), N_Desc, replace=False) 
        Ap_Exitosos[indices] = 1 # Los apareos exitosos se seleccionan aleatoriamente
        
        # Apareos
        contador_hijos=0
        contador_intentos=0
        for j in range(0, N_Rep): 
            for h in range(j+1, N_Rep-j):
                desc1, desc2 = Mating(pop_cromosomas[j], pop_cromosomas[h], m)
                if desc1==desc2: 
                    desc1, desc2, pop_cromosomas[h] = Nuevos_Desc_Y_Padre(desc1, desc2, pop_cromosomas[h], m)
                if Ap_Exitosos[contador_intentos] == 1:
                    pop_cromosomas[N_Keep + contador_hijos] = desc1
                    contador_hijos += 1
                contador_intentos += 1
                if Ap_Exitosos[contador_intentos] == 1:
                    pop_cromosomas[N_Keep + contador_hijos] = desc2
                    contador_hijos += 1
                contador_intentos += 1
        
        # Mutaciones, No se muta en la última iteración
        if i!=it-1:
            pop_cromosomas=Mutar(pop_cromosomas, in_mutar, m_rate)
    
        # Contamos el tiempo:
        t_it_fin = time.time()
        t_it=t_it_fin-t_it_in
        
        # Mejor individuo, guardarlo en fichero e imprimirlo por pantalla
        Mejor_ind = Decodificar(pop_cromosomas[0], m, x_max, x_min)
        Mejor_coste = costes[0]
        print("It", i+1,  Mejor_ind, "Cost", round(Mejor_coste, 7), "Time:", round(t_it, 1))
        file.write(f"{i+1} {Mejor_coste}\n"), file_2.write(f"{i+1} {Mejor_ind}\n")
    
    file.close(), file_2.close()  
        
fin = time.time()
tiempo_transcurrido = fin - inicio
print("Tiempo de ejecución:", round(tiempo_transcurrido, 1), "segundos")

