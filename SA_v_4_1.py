# -*- coding: utf-8 -*-
"""
Created on Sat May 27 14:40:56 2023

@author: Miguel Sánchez

Versión v_4_1:
    Igual que la v_4 pero con fichero Y CON DESCENSO DE SIGMA GEOMÉTRICO
    
    parámetros a calibrar: no copiar test raginin.
        Unas 800 it principales, 100 it por it princ
        T_ini = 0.075         
        T_fin = 1.21e-4
        alpha = 0.992
        alpha_s = 0.996
            
        num_iter 
        el parámetro suma de dism_sigma
            ese parámetro se fija en 15 porque la curva que representa
            es más o menos como queremos disminuir
        
Con estas características se llama a la función:
    alpha^(it bucle principal) * T_in = T_fin
    it bucle principal * it = it totales
    
    it bucle principal = ln(T_fin/T_in) / ln(alpha)
    alpha=e^(it_bucle pequeño / it tot *ln(T_fin/T_in))
    
"""

import numpy as np
import random
from Yee import Calculo_R
import time

#Medimos tiempo:
inicio = time.time()

# Definir el límite inferior y superior del espacio de búsqueda (queda definida la varianza inicial)
x_max = np.array([6., 6., 6., 6., 6., 6., 3., 3., 3.])
x_min = np.array([1., 1., 1., 1., 1., 1., 0., 0., 0.])
sigma = np.array(x_max - x_min) 

# Definir la solución inicial
sol_in = np.zeros(len(x_max))
for i in range(len(sol_in)):
    sol_in[i]=random.uniform( x_max[i], x_min[i] )

# Definir la temperatura inicial
T_ini = 0.05         
T_fin = 1e-10 #NO HACE FALTA QUE SEA TAN BAJA


# Definir la tasa de enfriamiento y tasa de reducción de sigma
alpha = 0.994
alpha_s = 0.995

# Definir el número máximo de iteraciones en bucle pequeño
num_iter = 500 
"""
Este parámetro sirve para salir de minimos locales, a mayor sea, mayor
facilidad para salir de mínimos locales. Más tiempo de ejecución también
"""

# Función de aceptación de soluciones peores
def accept_prob(deltaE, T):
    if deltaE < 0:
        return 1.0
    else:
        return np.exp(-deltaE/T)

# Inicializar la solución y temperatura
best_sol = sol_in
Coste_Best = Calculo_R(sol_in)
T = T_ini

# Archivo para guardar los resultados
file_name = "mejor_solucion.txt"

# Bucle de enfriamiento
a=1
with open(file_name, 'w') as file:
    while T>T_fin:
        for i in range(num_iter):
            # Contamos el tiempo:
            t_it_in = time.time()
            
            # Generar una solución vecina siguiendo gaussiana
            new_sol = np.zeros(len(x_max))
            for j in range(len(sol_in)): 
                new = x_min[j] - 1.0
                while new < x_min[j] or new > x_max[j]: #La sol debe estar en el intervalo dado
                    new = random.gauss(best_sol[j], sigma[j])
                new_sol[j] = new
        
            # Calcular el cambio en la función objetivo
            Coste_new  = Calculo_R(new_sol)
            deltaE     = Coste_new - Coste_Best
            
            # Evaluar la aceptación de la solución
            if accept_prob(deltaE, T) > random.random():
                best_sol   = new_sol
                Coste_Best = Coste_new
    
        # Contamos el tiempo:
        t_it_fin = time.time()
        t_it=t_it_fin-t_it_in
        
        # Imprimir resultados 
        if a%20==0:
            print(a, ":", best_sol, "Coste", round(Coste_Best, 4), "Sig", round(sigma[0], 1), "T", round(T, 4), "time", round(t_it))
            file.write(f"{a} {best_sol} {Coste_Best} {sigma[0]} {T}\n")   
        
        # Enfriar la temperatura y disminuir sigma
        T *= alpha 
        sigma *= alpha_s
        a+=1
        

    # Imprimir la solución encontrada
    print(a, ":", best_sol, "Coste", round(Coste_Best, 4), "Sig", round(sigma[0], 1), "T", round(T, 5))
    file.write(f"{a} {best_sol} {Coste_Best} {sigma[0]} {T}\n")
    file.close()
    
fin = time.time()
tiempo_transcurrido = fin - inicio
print("Tiempo: ", tiempo_transcurrido)
