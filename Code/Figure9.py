# -*- coding: utf-8 -*-
"""
@author: Maxime Clenet

This file is associated to the Figure 9.

Representation of the co-feasibility domain depending on 
the fixed intra-community interaction.

In particular, we are interested in the case where (kappa11,kappa22) > (kappa12,kappa21).
"""

#Importation of the main packages and functions:
import numpy as np
import matplotlib.pyplot as plt


def f(kappa_12, kappa_11=2):
    a = 1-2/kappa_12**2
    b = 2/kappa_11**2-2/kappa_12**2
    if a/b > 0 and a/b < 1:
        return a/b
    elif a/b > 1:
        return 1
    else:
        return 0
        


def g(kappa_21, kappa_22=2):

    a = 1-2/kappa_22**2
    b = 2/kappa_21**2-2/kappa_22**2
    if a/b > 0 and a/b < 1:
        return a/b
    elif a/b > 1:
        return 1
    else:
        return 0  
        

kappa_11 = 2
kappa_22 = 2

# Part dedicated to the plot of the figures:

x_inf = np.linspace(0.7, kappa_11-0.0001, 100)
x_sup = np.linspace(kappa_11+0.00001,3,100)
y1_inf = np.linspace(0, 1, 100)
y2_inf = np.linspace(0, 1, 100)
y1_sup = np.linspace(0, 1, 100)
y2_sup = np.linspace(0, 1, 100)
for i in range(len(x_inf)):


    y1_inf[i] = f(x_inf[i],kappa_11)
    y2_inf[i] = g(x_inf[i],kappa_22)

    y1_sup[i] = f(x_sup[i],kappa_11)
    y2_sup[i] = g(x_sup[i],kappa_22)



fig = plt.figure(1, figsize=(10, 6))
plt.plot(x_sup, y1_sup, label='Upper bound', color="#0000ff")
plt.plot(x_sup, y2_sup, label='Lower bound', color='red')
plt.plot(x_inf, y1_inf, color='red')
plt.plot(x_inf, y2_inf, color="#0000ff")
plt.fill_between(x_sup, y1_sup, y2_sup, where=y1_sup > y2_sup, color='#4300c84d')
plt.fill_between(x_inf, y1_inf, y2_inf, where=y1_inf < y2_inf, color='#4300c84d')
plt.vlines(kappa_11,-0.1,1.1,linestyles = 'dashed',color = 'black')
plt.xlabel(r"$\kappa_{12},\kappa_{21}$", fontsize=15)
plt.ylabel(r"$\beta_1$", fontsize=15)
plt.legend(loc='upper right')
plt.ylim(-0.1,1.1)
plt.show()