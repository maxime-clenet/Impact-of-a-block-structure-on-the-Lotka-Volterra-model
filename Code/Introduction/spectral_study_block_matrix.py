# -*- coding: utf-8 -*-
"""
@author: Maxime Clenet

This file is used to analyse the spectrum of the interaction matrix
"""

# Importation of the packages and implementation of the main
# function.

import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from Sensibility_LCP.functions import block_matrix_normal

# %%

N_SIZE = 20
beta = [3/4, 1/4]
sigma = np.array([[1, 1/3], [1/10, 1]])/np.sqrt(N_SIZE)

mu = np.array([[0, -10], [0, 0]])/N_SIZE

BMN = block_matrix_normal(N_SIZE, beta, mu, sigma)

# plt.colorbar(plt.pcolor(BMN))

sns.set_theme()
ax = sns.heatmap(BMN, vmin=-1.5, vmax=1.5, cmap="seismic",
                 xticklabels=np.arange(1, 21, 1), yticklabels=np.arange(1, 21, 1))

# %%
# Observation of the spectrum of the symmetric of the BMN matrix.

N_SIZE = 1000
beta = [1/2, 1/2]
sigma = np.array([[1/2, 1/10], [1/10, 1]])/np.sqrt(N_SIZE)

mu = np.array([[0, 0], [0, 0]])/N_SIZE

BMN = block_matrix_normal(N_SIZE, beta, mu, sigma)


eig_B = np.linalg.eigvals(BMN)

max(abs(eig_B))

sig = np.array([[1/2, 1/3], [1/10, 1]])

eig = (np.linalg.eigvals(np.diag(beta)**(1/2) @
                         (sig**2)@np.diag(beta)**(1/2)))**(1/2)

fig = plt.figure(1, figsize=(10, 6))
for alpha in eig:

    radius = alpha  # radius of the circle.

    t = np.linspace(0, 2*np.pi, 100)

    plt.plot(radius*np.cos(t), radius*np.sin(t), color='red', linewidth=3)

plt.plot(eig_B.real, eig_B.imag, '.', color='k')
plt.grid(color='lightgray', linestyle='--')
plt.axis("equal")

plt.xlabel(r"Real axis", fontsize=15)
plt.ylabel("Imaginary axis", fontsize=15)
plt.show()

plt.close()


# %%


N_SIZE = 1000
beta = [1/2, 1/2]
sigma = np.array([[1/np.sqrt(2), 1/np.sqrt(2)],
                 [1/np.sqrt(2), 1/np.sqrt(2)]])/np.sqrt(N_SIZE)

mu = np.array([[0, 0], [0, 0]])/N_SIZE

BMN = block_matrix_normal(N_SIZE, beta, mu, sigma)


eig_B = np.linalg.eigvals(BMN)
#eig_B = np.linalg.eigvals(BMN+np.eye(N_SIZE))
print(max(abs(eig_B)))

sing_1 = np.linalg.eigvals(np.linalg.inv(-2*np.linalg.inv(BMN)+np.eye(N_SIZE)))
max(abs(sing_1))

fig = plt.figure(1, figsize=(10, 6))

plt.plot(eig_B.real, eig_B.imag, '.', color='k')
plt.grid(color='lightgray', linestyle='--')

plt.xlabel(r"Real axis", fontsize=15)
plt.ylabel("Imaginary axis", fontsize=15)
plt.xlim(-2, 2)
plt.ylim(-2, 2)
plt.show()


plt.close()


# %%

N_SIZE = 1000
beta = [1/2, 1/2]
sigma = np.array([[1, 1], [1, 1]])/np.sqrt(N_SIZE)

mu = np.array([[0, 0], [0, 0]])/N_SIZE

BMN = block_matrix_normal(N_SIZE, beta, mu, sigma)


eig_B = np.linalg.eigvals(BMN)

max(abs(eig_B))

sig = np.array([[1, 1], [1, 1]])

eig = (np.linalg.eigvals(np.diag(beta)**(1/2) @
                         (sig**2)@np.diag(beta)**(1/2)))**(1/2)
