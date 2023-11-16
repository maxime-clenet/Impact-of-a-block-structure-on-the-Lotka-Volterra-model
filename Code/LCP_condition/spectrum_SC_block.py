# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 09:40:44 2022

@author: Maxime Clenet
"""


import numpy as np
import matplotlib.pyplot as plt
from functions import block_matrix_normal


def Gamma(u, v, z, beta, sigma):
    a1 = beta[0]*2*sigma[0, 0]**2
    a2 = beta[1]*(sigma[0, 1]**2+sigma[1, 0]**2)
    a3 = beta[0]*(sigma[1, 0]**2+sigma[0, 1]**2)
    a4 = beta[1]*2*sigma[1, 1]**2
    return(z+a1*u+a2*v, z+a3*u+a4*v)


# Renvoie la transformée de Stieljes.
def g(x, beta, sigma):
    z = complex(x, 10**(-3))
    u = -1/z
    v = -1/z

    for i in range(10000):

        (u, v) = Gamma(-1/u, -1/v, z, beta, sigma)
    return(beta[0]*(-1/u)+beta[1]*(-1/v))


def plot_block_spectrum(n, beta, sigma):

    mu = np.array([[0, 0], [0, 0]])
    sigma_norm = sigma/np.sqrt(n)
    B = block_matrix_normal(n, beta, mu, sigma_norm)

    eig_B = np.linalg.eigvals(B+B.T)

    # The rest of the function is dedicated to the plot.

    x = np.linspace(-3, 3, 100)
    y = np.linspace(-2, 2, 100)

    for i in range(100):
        y[i] = 1/np.pi*g(x[i], beta, sigma).imag

    # bound = 2*(np.linalg.norm(np.diag(beta)**(1/2) @
    #                           (sigma**2+sigma.T**2)@np.diag(beta)**(1/2),ord = np.inf))**(1/2)
    bound = 2*(np.linalg.norm(np.diag(beta) @
                              (sigma**2+sigma.T**2),ord = np.inf))**(1/2)

    fig = plt.figure(1, figsize=(10, 6))

    plt.vlines(bound, 0, 0.9,
               linestyles='--', color='k', linewidth=2)
    plt.plot(x, y, color='k', linewidth=2)
    plt.hist(eig_B.real, density=True, bins=40,
             edgecolor='black', color='#0000003d')
    plt.grid(color='lightgray', linestyle='--')

    plt.xlim(-bound-0.5, bound+0.5)
    plt.ylim(0, 0.9)

    plt.xlabel(r"Spectrum", fontsize=15)
    plt.ylabel("Density", fontsize=15)
    plt.show()

    plt.close()

    return fig
# plt.savefig('Spectre_Ginibre_Presentation_Lille.pdf')

#Case 1:

beta = [1/2, 1/2]
sigma = np.array([[1/np.sqrt(2), 1/np.sqrt(2)], [1/np.sqrt(2), 1/np.sqrt(2)]])

plot_block_spectrum(n=2000, beta=beta, sigma=sigma)

beta = [1/2, 1/2]
sigma = np.array([[1/2, 1], [1, 1/8]])

plot_block_spectrum(n=2000, beta=beta, sigma=sigma)

beta = [1/2, 1/2]
sigma = np.array([[1, 1/2], [1/8, 1]])

plot_block_spectrum(n=2000, beta=beta, sigma=sigma)

beta = [3/4, 1/4]
sigma = np.array([[1/3, 1/5], [1, 1/2]])


# An exemple of the first scenario:
plot_block_spectrum(n=2000, beta=beta, sigma=sigma)

# An exemple of the second scenario:
