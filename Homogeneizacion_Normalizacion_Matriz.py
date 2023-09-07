#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 12:43:12 2023

@author: violetapascuallaborda
"""

import numpy as np

CM = "/Users/violetapascuallaborda/Desktop/TFG/matricesContactos/SriLankaMC.txt"
Population = "/Users/violetapascuallaborda/Desktop/TFG/piramidesPoblacionales/poblacionSriLanka.txt"
matrix = np.loadtxt(CM)
N = np.loadtxt(Population)

Ngrupos = 16
suma = 0
kmedio = 0
Ntotal = 0

#Homogeneización 
for i in np.arange(0, Ngrupos, 1):
    for j in np.arange(0, Ngrupos, 1):
        if(i!=j):
            matrix[i][j] = (matrix[i][j]*N[i]+matrix[j][i]*N[j])/(N[i]*2.0)
  
#Fuente: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006638

#Normalización
for i in np.arange(0, Ngrupos, 1):
    for j in np.arange(0, Ngrupos, 1):
        if(i<=j):
           suma = suma + matrix[i][j]

for i in np.arange(0, Ngrupos, 1):
    for j in np.arange(0, Ngrupos, 1):
           matrix[i][j] = matrix[i][j]/float(suma)
           

file = open('/Users/violetapascuallaborda/Desktop/TFG/MatricesContactos/SriLankaMC_Homogeneizada_Normalizada.txt', 'w')
for i in np.arange(0, Ngrupos, 1):
    for j in np.arange(0, Ngrupos, 1):
        file.write("%lf " %matrix[i][j])
    file.write("\n")
        
file.close()
