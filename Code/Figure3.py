# -*- coding: utf-8 -*-
"""
@author: Maxime Clenet
"""

# %%

import numpy as np
import matplotlib.pyplot as plt
from Code.base_function import block_function, empirical_prop


# Test des sensi dans le cas LCP:

beta = [1/2, 1/2]
mu = np.array([[0, 0], [0, 0]])
s = np.array([[1, 1], [1, 1]])/2

#print(empirical_prop(300, alpha, mu, beta, mc_prec=500))
print(block_function(mu, 1/s, beta))



x = np.linspace(1/5, 1, 20)
y_theo = np.zeros((6, 20))
y_emp = np.zeros((6, 20))


beta = [1/2, 1/2]
mu = np.array([[0, 0], [0, 0]])

for i in range(20):
    print(i)
    s = np.array([[1/2, x[i]], [x[i], 1/1.4]])
    y_theo[:, i] = block_function(mu, 1/s, beta)
    y_emp[:, i] = empirical_prop(50, 1/s, mu, beta)


for t in range(6):
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


