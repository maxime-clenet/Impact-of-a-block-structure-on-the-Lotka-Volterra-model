# -*- coding: utf-8 -*-
"""
@author: Maxime Clenet
"""

import numpy as np
import matplotlib.pyplot as plt


N = 10000

print(np.sqrt(np.log(N)))

print(np.sqrt(np.log(N)))

x = np.linspace(0, 5, 100)
y = np.linspace(0, 5, 100)
for i in range(100):
    y[i] = np.sqrt(1/(1-1/x[i]**2))

fig = plt.figure(1, figsize=(10, 6))
plt.plot(x, y)
plt.xlabel(r"$\kappa_1$")
plt.ylabel(r"$\kappa_2$")
plt.title("Phase diagram")
plt.show

# plt.savefig('Phase_diagram_v0.pdf')


# %%

# Phase diagram dans le cas plus complexe.


N = 10000

print(np.sqrt(np.log(N)))

print(np.sqrt(np.log(N)))


beta_1 = 4/5

beta_2 = 1/5

x = np.linspace(0, 5, 100)
y = np.linspace(0, 5, 100)
y1 = np.linspace(0, 5, 100)
y2 = np.linspace(0, 5, 100)
for i in range(100):

    y1[i] = np.sqrt(beta_1/(1/2-beta_2/x[i]**2))

    y2[i] = np.sqrt(beta_2/(1/2-beta_1/x[i]**2))

    if(y1[i] > x[i]):
        y[i] = y1[i]
    else:
        y[i] = y2[i]

fig = plt.figure(1, figsize=(10, 6))
plt.plot(x, y)

plt.plot(x, x)
# plt.plot(x,y)
plt.xlabel(r"$\kappa_1$")
plt.ylabel(r"$\kappa_2$")
plt.title("Phase diagram")
plt.show

# plt.savefig('Phase_diagram_beta_v0.pdf')
