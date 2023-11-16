# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 15:55:47 2021

@author: Maxime
"""


import seaborn as sns
from scipy import linspace, meshgrid, arange, empty, concatenate, newaxis, shape
from numpy.random import randn, shuffle
import numpy
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
import numpy as np
import matplotlib.pyplot as plt

N = 100

Beta = [1/2, 1/2]


def f(Kappa_11, Kappa_22, Beta):

    a = np.sqrt(Beta[1]/(1/2-Beta[0]/Kappa_11**2))
    b = np.sqrt(Beta[0]/(1/2-Beta[1]/Kappa_22**2))
    if (Kappa_11 > Kappa_22):
        return(a)
    else:
        return(b)


# =========================
# generating ordered data:
N = 100
x = np.linspace(1.5, 3, 100)
y = np.linspace(1.5, 3, 100)

X, Y = np.meshgrid(x, y)


Z = np.zeros((N, N))

for i in range(N):
    for j in range(N):
        Z[i, j] = f(x[i], y[j], Beta)


# ======================================
# reference picture (X, Y and Z in 2D):

fig = plt.figure(1, figsize=(20, 6))
ax = fig.add_subplot(projection='3d')

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet, linewidth=0)
fig.colorbar(surf)

#title = ax.set_title("Feasibility phase diagram")
# title.set_y(1.01)

ax.set_xlabel(r"Community 1 $(\kappa_{11})$", fontsize=20)
ax.set_ylabel(r"Community 2 $(\kappa_{22})$", fontsize=20)
ax.set_zlabel(r"Inter-communities $ (\kappa_{12},\kappa_{21}) $", fontsize=20)

ax.xaxis.set_major_locator(MaxNLocator(5))
ax.yaxis.set_major_locator(MaxNLocator(6))
ax.zaxis.set_major_locator(MaxNLocator(5))
ax.set_facecolor('w')
ax.view_init(25, 30)
fig.tight_layout()
fig.savefig('Feasibility_Phase_Diagram.pdf')

# %%

# Heatmap format:

sns.set_theme()
ax = sns.heatmap(Z, cbar=True)
ax.xaxis.set_major_locator(MaxNLocator(10))
ax.set(xlabel='X-Axis', ylabel='Y-Axis')
