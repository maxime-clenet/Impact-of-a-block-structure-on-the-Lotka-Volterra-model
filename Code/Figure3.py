# -*- coding: utf-8 -*-
"""
@author: Maxime Clenet

This code is associated with figure 3.

The purpose of this file is to verify that the heuristics are correct.

"""


# Importation of the important packages and main functions:
import numpy as np
import matplotlib.pyplot as plt
from base_function import block_function, empirical_prop


x = np.linspace(1/5, 1, 20)
y_theo = np.zeros((6, 20))
y_emp = np.zeros((6, 20))

# Parameters of the Block model:
beta = [1/2, 1/2]
mu = np.array([[0, 0], [0, 0]])

# Comparison between the theoretical solutions 
# $(p_1^*,p_2^*,\sigma_1^*,\sigma_2^*)$ of the heuristics
#  and their empirical Monte Carlo counterpart as functions of the off-diagonal 
# block interaction strength $(s_{12},s_{21})$. 

for i in range(20):
    print(i)
    s = np.array([[1/2, x[i]], [x[i], 1/1.4]])
    # Theoretical solution of the heuristics :
    y_theo[:, i] = block_function(mu, 1/s, beta)
    # Empirical Monte Carlo counterpart :
    y_emp[:, i] = empirical_prop(500, 1/s, mu, beta)



# This part of the file is dedicated for the plot of the Figure :

for t in range(6):
    # Proportion of the surviving species in each community :
    if t in {0, 1}:
        fig = plt.figure(1, figsize=(10, 6))

        plt.plot(x, y_theo[t, :], 'k')
        plt.plot(x, y_emp[t, :], 'k*')
        plt.xlabel(
            r"Interaction between community 1 and community 2 ($s_{21}=s_{12}$)", fontsize=15)
        if t == 0:
            plt.ylabel(
                r"Proportion of the surviving species ($p_1^*$)", fontsize=15)
        if t == 1:
            plt.ylabel(
                r"Proportion of the surviving species ($p_2^*$)", fontsize=15)
        #plt.legend(loc='upper right')
        plt.show()
        plt.close()
    # Mean of the surviving species in each community :
    if t in {2, 3}:
        fig = plt.figure(1, figsize=(10, 6))
        plt.plot(x, y_theo[t, :], 'k')
        plt.plot(x, y_emp[t, :], 'k*')
        plt.xlabel(
            r"Interaction between community 1 and community 2 ($s_{21}=s_{12}$)", fontsize=15)
        if t == 2:
            plt.ylabel(r"Mean of the surviving species ($m_1^*$)", fontsize=15)

        if t == 3:
            plt.ylabel(r"Mean of the surviving species ($m_2^*$)", fontsize=15)

        #plt.legend(loc='upper right')
        plt.show()
        plt.close()
        
    # Root mean square of the surviving species in each community :
    if t in {4, 5}:
        fig = plt.figure(1, figsize=(10, 6))
        plt.plot(x, y_theo[t, :], 'k')
        plt.plot(x, y_emp[t, :], 'k*')
        plt.xlabel(
            r"Interaction between community 1 and community 2 ($s_{21}=s_{12}$)", fontsize=15)
        if t == 4:
            plt.ylabel(
                r"Root mean square of the surviving species ($\sigma_1^*$)", fontsize=14)
        if t == 5:
            plt.ylabel(
                r"Root mean square of the surviving species ($\sigma_2^*$)", fontsize=14)

        #plt.legend(loc='upper right')
        plt.show()
        plt.close()


