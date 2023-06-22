# -*- coding: utf-8 -*-
"""
Created on Sun Jun  4 18:16:13 2023

@author: Miguel Sánchez
"""

import numpy as np
import random
from Test_2 import f2


# Definir el límite inferior y superior del espacio de búsqueda (queda definida la varianza inicial)
x_max = np.array([5.12, 5.12, 5.12, 5.12, 5.12, 5.12])
x_min = np.array([-5.12, -5.12, -5.12, -5.12, -5.12, -5.12])
sigma = np.array(x_max - x_min) / 50.

# Definir la solución inicial
sol_in = np.zeros(len(x_max))
for i in range(len(sol_in)):
    sol_in[i]=random.uniform( x_max[i], x_min[i] )

# Definir la temperatura inicial
T_ini = 0.1          
T_fin = 0

# Definir la tasa de enfriamiento y tasa de reducción de sigma
alpha = 0.9998
alpha_s = 0.996

# Definir el número máximo de iteraciones en bucle pequeño
num_iter = 130000 

# Función de aceptación de soluciones peores
def accept_prob(deltaE, T):
    if deltaE < 0:
        return 1.0
    else:
        return np.exp(-deltaE/T)

# Inicializar la solución y temperatura
best_sol = sol_in
Coste_Best = f2(sol_in)
T = T_ini

# Archivo para guardar los resultados
file_name = "Coste_it.txt"
file_name2 = "Indiv_it.txt"
# Bucle de enfriamiento

with open(file_name, 'w') as file, open(file_name2, 'w') as file2:
    for i in range(num_iter):
        
        # Generar una solución vecina siguiendo gaussiana
        new_sol = np.zeros(len(x_max))
        for j in range(len(sol_in)): 
            new = x_min[j] - 1.0
            while new < x_min[j] or new > x_max[j]: #La sol debe estar en el intervalo dado
                new = random.gauss(best_sol[j], sigma[j])
            new_sol[j] = new
    
        # Calcular el cambio en la función objetivo
        Coste_new  = f2(new_sol)
        deltaE     = Coste_new - Coste_Best
        
        # Evaluar la aceptación de la solución
        if deltaE < 0:
            best_sol   = new_sol
            Coste_Best = Coste_new
    
        # Imprimir resultados 
        if i%10000==0:
            print(i, ":", best_sol, "Coste", Coste_Best, "Sig", sigma[0], 7)
            file.write(f"{i} {Coste_Best}\n")
            file2.write(f"{i} {best_sol} {sigma[0]} {T}\n")
        
        # Enfriar la temperatura y disminuir sigma
        sigma *= alpha_s
            

    # Imprimir la solución encontrada
    print(i, ":", best_sol, "Coste", Coste_Best, "Sig", sigma[0])
    file.write(f"{i} {Coste_Best}\n")
    file2.write(f"{i} {best_sol} {sigma[0]} {T}\n")
    file.close(), file2.close()
    