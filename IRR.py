#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 15:32:29 2023

@author: violetapascuallaborda
"""
import numpy as np
import matplotlib.pyplot as plt

IRR_fuerza = "/Users/violetapascuallaborda/Desktop/TFG/SriLanka/SriLanka_IRR_fuerzaInfeccion.txt"
IRR_incidencia = "/Users/violetapascuallaborda/Desktop/TFG/SriLanka/SriLanka_IRR_incidencia.txt"

fuerza = np.loadtxt(IRR_fuerza)
incidencia = np.loadtxt(IRR_incidencia)

xvalues = list(range(1, 17))


plt.figure(figsize=(8, 5))  # Crear una figura independiente
# Crear el gráfico de líneas
plt.plot(xvalues, incidencia, marker='o', color= '#b01518ff', label='Incidencia')
plt.plot(xvalues, fuerza, marker='o', label='Fuerza de infección')

# Personalizar el eje x
plt.xticks(xvalues)

# Personalizar el gráfico
#plt.title('IRR (Incidence Reduction Rate)')
plt.xlabel('Número de grupos vacunados')
#plt.ylabel('Reducción de la incidencia')
plt.ylabel('IRR (Incidence Reduction Rate)')

# Agregar leyendas y mostrar la gráfica
plt.legend(loc='lower right')
plt.grid(True, linestyle=':', alpha=0.7)  # Agregar estilo de grid
plt.show()


IRR_fuerza = "/Users/violetapascuallaborda/Desktop/TFG/SriLanka/SriLanka_IRR_fuerzaInfeccion1.txt"
IRR_incidencia = "/Users/violetapascuallaborda/Desktop/TFG/SriLanka/SriLanka_IRR_incidencia1.txt"

fuerza = np.loadtxt(IRR_fuerza)
incidencia = np.loadtxt(IRR_incidencia)

xvalues = list(range(1, 17))


plt.figure(figsize=(8, 5))  # Crear una figura independiente
# Crear el gráfico de líneas
plt.plot(xvalues, incidencia, marker='o', color= '#b01518ff', label='Incidencia')
plt.plot(xvalues, fuerza, marker='o', label='Fuerza de infección')

# Personalizar el eje x
plt.xticks(xvalues)

# Personalizar el gráfico
#plt.title('IRR (Incidence Reduction Rate)')
plt.xlabel('Número de grupos vacunados')
#plt.ylabel('Reducción de la incidencia')
plt.ylabel('IRR (Incidence Reduction Rate)')

# Agregar leyendas y mostrar la gráfica
plt.legend(loc='lower right')
plt.grid(True, linestyle=':', alpha=0.7)  # Agregar estilo de grid
plt.show()
