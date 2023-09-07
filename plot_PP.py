#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 15:45:46 2023

@author: violetapascuallaborda
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

poblacion = "/Users/violetapascuallaborda/Desktop/TFG/piramidesPoblacionales/poblacionSriLanka.txt"

# Datos de ejemplo para las pirámides poblacionales
xvalues = ['0-4', '5-9', '10-14', '15-19', '20-24', '25-29', '30-34',' 35-39', 
          '40-44', '45-49', '50-54', '55-59', '60-64', '65-69', '70-74', '75+']

N = np.loadtxt(poblacion)

Ngrupos = 16
Ntotal = np.sum(N[:Ngrupos])  # Sumar solo los primeros Ngrupos valores de N

yvalues = (N[:Ngrupos] / Ntotal)*100

print(Ntotal)

# Crear la instancia de ScalarFormatter para el eje y
formatter = ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-3, 3))  # Ajustar estos límites según tus necesidades

# Crear la figura y los ejes
fig, ax = plt.subplots()
plt.grid(True, linestyle=':', alpha=0.7)  # Agregar estilo de grid

# Graficar la pirámide poblacional
ax.barh(xvalues, yvalues, align='center', color='#b01518ff')

# Configurar el eje y
ax.set_yticks(np.arange(len(xvalues)))
ax.set_yticklabels(xvalues)

# Agregar una línea vertical en el centro
ax.axvline(0, color='gray', linewidth=0.8)

# Aplicar el formateo personalizado a los números en el eje x
ax.xaxis.set_major_formatter(formatter)

# Etiquetas de los ejes y título del gráfico
ax.set_xlabel('Porcentaje de personas (%)')
ax.set_ylabel('Grupo de edad')

# Mostrar el gráfico
plt.tight_layout()
plt.show()
