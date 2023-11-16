# -*- coding: utf-8 -*-
"""
@author: Maxime Clenet
"""

import numpy as np
import matplotlib.pyplot as plt
from Sensibility_LCP.functions import block_function
import pandas as pd
import seaborn

# %%


# Sensitivity between two communities: part 1.

x = np.linspace(1, 7, 50)
y1 = np.linspace(1, 3, 50)
y2 = np.linspace(1, 3, 50)
y3 = np.linspace(1, 3, 50)


beta = [1/2, 1/2]
mu1 = np.array([[0, -0.4], [0.4, 0]])
mu2 = np.array([[0, 0.4], [0.4, 0]])
mu3 = np.array([[0, -0.4], [-0.4, 0]])

t = 1

for i in range(50):
    alpha1 = np.array([[2, x[i]], [x[i], 2]])
    y1[i] = block_function(mu1, alpha1, beta)[t]
    alpha2 = np.array([[2, x[i]], [x[i], 2]])
    y2[i] = block_function(mu2, alpha2, beta)[t]
    alpha3 = np.array([[2, x[i]], [x[i], 2]])
    y3[i] = block_function(mu3, alpha3, beta)[t]


fig = plt.figure(1, figsize=(10, 6))

plt.plot(x, y1, label='Antagonist')

plt.plot(x, y2, label='Mutualism')

plt.plot(x, y3, 'k*', label='Competitive')

plt.xlabel(
    r"Interaction between community 1 and community 2 ($\alpha_{21}=\alpha_{12}$)", fontsize=15)
plt.ylabel(r"Persisting species in the 1st community ($p_1$)", fontsize=15)

plt.legend(loc='upper right')
#plt.title('Sensibility curve')
plt.show()

# plt.close()


# %%

# Example 2:

prec = 30

x = np.linspace(0.9, 2, prec)
y = np.linspace(0.4, -0.4, prec)
Beta = [1/2, 1/2]
mu = np.array([[0, 0], [0, 0]])

mat = np.zeros((prec, prec))

for i, v1 in enumerate(x):
    for j, v2 in enumerate(y):
        alpha = np.array([[2, 10**3], [v1, 2]])
        mu = np.array([[v2, 0], [0, 0]])
        mat[j, i] = block_function(mu, alpha, Beta)[1]

df = pd.DataFrame(mat, columns=np.around(x, decimals=2),
                  index=np.around(y, decimals=2))
fig = plt.figure(1, figsize=(10, 6))
seaborn.heatmap(df, cbar=True, cmap='plasma')
plt.tight_layout()
plt.xlabel(
    r"Impact of community 1 on community 2 ($ \alpha_{21}$)", fontsize=15)
plt.ylabel(r"Type of community 1 ($ \mu_{11}$)", fontsize=15)


# plt.close()


# %%

# Example 3:

prec = 30

x = np.linspace(0.9, 2, prec)
y = np.linspace(0.4, -0.4, prec)
Beta = [1/2, 1/2]
mu = np.array([[0, 0], [0, 0]])

mat = np.zeros((prec, prec))

for i, v1 in enumerate(x):
    for j, v2 in enumerate(y):
        alpha = np.array([[2, 10**3], [v1, 2]])
        mu = np.array([[-0.4, 0], [0, v2]])
        mat[j, i] = block_function(mu, alpha, Beta)[1]

df = pd.DataFrame(mat, columns=np.around(x, decimals=2),
                  index=np.around(y, decimals=2))
fig = plt.figure(1, figsize=(10, 6))
seaborn.heatmap(df, cbar=True, cmap='plasma')
plt.tight_layout()
plt.xlabel(
    r"Impact of community 1 on community 2 ($ \alpha_{21}$)", fontsize=15)
plt.ylabel(r"Type of community 2 ($ \mu_{22}$)", fontsize=15)

# plt.close()

# %%


# Example 4

x = np.linspace(1, 7, 50)
y1 = np.linspace(1, 3, 50)
y2 = np.linspace(1, 3, 50)
y3 = np.linspace(1, 3, 50)
y4 = np.linspace(1, 3, 50)


beta = [1/2, 1/2]
c1 = -0.4
c2 = - 0.4

mu1 = np.array([[0.4, c1], [c2, 0.4]])
mu2 = np.array([[0.4, c1], [c2, -0.4]])
mu3 = np.array([[-0.4, c1], [c2, 0.4]])
mu4 = np.array([[-0.4, c1], [c2, -0.4]])

t = 1
for i in range(50):
    alpha1 = np.array([[2, x[i]], [x[i], 2]])
    y1[i] = block_function(mu1, alpha1, beta)[t]
    alpha2 = np.array([[2, x[i]], [x[i], 2]])
    y2[i] = block_function(mu2, alpha2, beta)[t]
    alpha3 = np.array([[2, x[i]], [x[i], 2]])
    y3[i] = block_function(mu3, alpha3, beta)[t]
    alpha4 = np.array([[2, x[i]], [x[i], 2]])
    y4[i] = block_function(mu4, alpha4, beta)[t]


fig = plt.figure(1, figsize=(10, 6))

plt.plot(x, y1, label='M-M')

plt.plot(x, y2, label='M-C')

plt.plot(x, y3, label='C-M')

plt.plot(x, y4, label='C-C')
plt.ylim(0.65, 1)
plt.xlabel(
    r"Interaction between community 1 and community 2 ($\alpha_{21},\alpha_{12}$)", fontsize=15)
plt.ylabel(r"Surviving species in the 2nd community ($p_2$)", fontsize=15)

plt.legend(loc='upper right')
#plt.title('Sensibility curve')
plt.show()

# plt.close()
