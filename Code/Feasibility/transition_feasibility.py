# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 13:45:26 2021

@author: Maxime

Transition toward feasibility.
"""

# Importation of the packages
import numpy as np
import matplotlib.pyplot as plt


def feasibility_solution(kappa_1, kappa_2, N, sigma):

    beta_1 = 5/10
    beta_2 = 5/10
    # Initialization
    A_11 = np.random.randn(int(N*beta_1), int(N*beta_1))*sigma/kappa_1
    A_12 = np.random.randn(int(N*beta_1), int(N*beta_2))*sigma/kappa_2
    A_21 = np.random.randn(int(N*beta_2), int(N*beta_1))*sigma/kappa_2
    A_22 = np.random.randn(int(N*beta_2), int(N*beta_2))*sigma/kappa_1

    A_bis1 = np.concatenate([A_11, A_12], axis=1)
    A_bis2 = np.concatenate([A_21, A_22], axis=1)
    A = np.concatenate([A_bis1, A_bis2], axis=0)

    sol = np.dot(np.linalg.inv(
        np.eye(N, N)-(1/(np.sqrt(np.log(N))*np.sqrt(N)))*A), np.ones(N))

    return(sol)


# Set the parameters:
MC_PREC = 500  # number of MC experiments
SIZE_SAMPLE = 40  # definition of the size of the sample
# We introduce the interval of Kappa we want to study.
ind = np.linspace(0.6, 1.85, SIZE_SAMPLE)

# Dimension of the matrix:
n = 5000

sigma = 1

compt = 0  # initial counter of the progress of simulation

# Part of the code dedicated to compute the value
# associated to the transition towards feasibility.

beta_1 = 1/2
beta_2 = 1/2
kappa_2 = 2

res = np.zeros(SIZE_SAMPLE)
compt_res = 0

for kappa in ind:
    compt += 1
    print('Progress:', compt, '/', SIZE_SAMPLE, end='\n')
    count_sol = 0
    kappa_1 = kappa
    for i in range(MC_PREC):

        # Compute the solution
        sol = feasibility_solution(kappa_1, kappa_2, n, sigma)

        if sol[sol < 0].shape[0] == 0:
            count_sol = count_sol+1  # If positive solution.

    res[compt_res] = count_sol/MC_PREC  # Proportion of positive solution
    compt_res = compt_res+1


threshold = np.sqrt(beta_1/(1/2-beta_2/kappa_2**2))


# Part dedicated to the display of the figure:

fig = plt.figure(1, figsize=(10, 6))

plt.plot(ind, res, linestyle='solid',
         color='k')

plt.axvline(threshold, color='black',
            linestyle='dashdot', label='Threshold')
axes = plt.gca()
axes.xaxis.set_ticks([1, threshold, 1.5, 1.7])
axes.xaxis.set_ticklabels(
    ('1', r'$\sqrt{\frac{\beta_1}{\frac{1}{2}-\frac{\beta_2}{\kappa_2^2}}}$', '1.5', '1.7'), color='black', fontsize=10)
plt.xlabel(r"$\nu_1$", fontsize=15)
plt.ylabel("P(Feasibility)", fontsize=15)
plt.legend(loc='upper left')
plt.show()
