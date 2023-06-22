# -*- coding: utf-8 -*-
"""
Created on Sat May  6 12:14:42 2023

@author: Miguel Sánchez

Test: Rastrigin Function

Función test 1 para algoritmo genético
"""

import numpy as np

A = 10

def Rastrigin(x):

    res = A*len(x)
    for i in range(0, len(x)):
        if (x[i] < -5.12 or x[i] > 5.12):
            raise Exception("Error. Fuera del dominio en el que se quiere encontrar el mínimo")
        res += x[i]**2 - A*np.cos( 2 * np.pi * x[i] )
    
    return res
        