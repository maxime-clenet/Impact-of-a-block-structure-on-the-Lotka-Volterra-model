# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 15:02:15 2022

@author: Maxime



"""
import numpy as np
import matplotlib.pyplot as plt
from Sensibility_LCP.functions import block_matrix_normal

N = 100

beta = [3/4, 1/4]
sig_1 = np.array([[0, 0], [0, 0]])
sigma_norm_1 = sig_1/np.sqrt(N)

mu = np.array([[-2, -10], [-3, -5]])
mu_norm = mu/N


BMN_1 = block_matrix_normal(N, beta, mu_norm, sigma_norm_1)

eig_BMN_1 = np.linalg.eigvals(BMN_1+BMN_1.transpose())

eig = eig_BMN_1[abs(eig_BMN_1) > 0.01]

print(np.linalg.eigvals(np.diag(beta)**(1/2)@(mu+mu.T)@np.diag(beta)**(1/2)))
print(np.linalg.eigvals(np.diag(beta)@(mu+mu.T)))
print(eig)

print((mu[0, 0]+mu[1, 1]+np.sqrt((mu[0, 1]+mu[1, 0])**2+(mu[0, 0]-mu[1, 1])**2))/2)
print((mu[0, 0]+mu[1, 1]-np.sqrt((mu[0, 1]+mu[1, 0])**2+(mu[0, 0]-mu[1, 1])**2))/2)


# %%


def plot_capitaine_spectrum(n, beta, mu, alpha):

    mu_norm = mu/n
    sigma = (1/alpha)*np.ones((2, 2))/np.sqrt(n)
    B = block_matrix_normal(n, beta, mu_norm, sigma)

    eig_B = np.linalg.eigvals(B+B.T)

    # The rest of the function is dedicated to the plot.

    t = np.linspace(-2*np.sqrt(2/alpha**2), 2*np.sqrt(2/alpha**2), 100000)

    fig = plt.figure(1, figsize=(10, 6))
    eig = np.linalg.eigvals(np.diag(beta)**(1/2) @
                            (mu+mu.T)@np.diag(beta)**(1/2))
    print(eig)
    plt.vlines(eig[0]+2/(eig[0]*alpha**2), 0, 0.35,
               linestyles='--', color='k', linewidth=2)
    plt.vlines(eig[1]+2/(eig[1]*alpha**2), 0, 0.35,
               linestyles='--', color='k', linewidth=2)

    plt.plot(t, np.sqrt(4*(2/(alpha**2))-t**2) /
             (2*np.pi*(2/alpha**2)), color='k', linewidth=3)
    plt.hist(eig_B.real, density=True, bins=40,
             edgecolor='black', color='#0000003d')
    plt.grid(color='lightgray', linestyle='--')

    plt.xlim(-7, 7)
    plt.ylim(0, 0.35)

    plt.xlabel(r"Spectrum", fontsize=15)
    plt.ylabel("Density", fontsize=15)
    plt.show()

    plt.close()

    return fig


# plt.savefig('Spectre_Ginibre_Presentation_Lille.pdf')
n = 1000
beta = [1/2, 1/2]
alpha = np.sqrt(2)
mu = np.array([[-2, -4], [-2, -2]])


# An exemple of the first scenario:
plot_capitaine_spectrum(n, beta, mu, alpha)


# %%


def psi(mu):
    return 0.5*(mu[0, 0]+mu[1, 1]+np.sqrt((mu[0, 1]+mu[1, 0])**2+(mu[0, 0]-mu[1, 1])**2))


x = np.linspace(-2, 6, 100)
y = np.linspace(-2, 2, 100)

for i in range(100):
    mu = np.array([[2, x[i]], [-2, -4]])
    y[i] = psi(mu)

fig = plt.figure(1, figsize=(10, 6))
plt.plot(x, y, 'k')

plt.grid(color='lightgray', linestyle='--')

plt.xlim(-2, 6)


plt.xlabel(r"$\mu_{12}$", fontsize=15)
plt.ylabel("$\psi(\mu)$", fontsize=15)
plt.show()
