# -*- coding: utf-8 -*-
"""
@author: Maxime Clenet

This code is associated with figures 1 and 2.

The aim is to display a dynamics of the LV system in
the case of a block matrix. We also add the matrix of interaction associated
with the dynamics.
"""

# Importation of the important packages and main functions:
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from base_function import block_matrix_normal, dynamics_LV_2B

# alpha type of Figure 1:

def alpha_type_normal(t, N):
    return 1

# alpha type of Figure 2:

def alpha_type_linear(t, N):
    if 15 > t > 5:
        return (t-5)*N
    elif t >= 15:
        return 100
    else:
        return 1
    

# Number of iterations:
nbr_it = 10000
# Selection of the time step:
tau = 0.002
# Number of species:
N = 10
# Initial condition:
x_init = np.random.random(N)*2


# Parameters of the block matrix:
Beta = [1/2, 1/2]
Sigma = np.array([[1/2, 1/80], [1/80, 1/2]])
Mu = np.array([[0, 0], [0, 0]])

A = block_matrix_normal(N, Beta, Mu, Sigma)
A_base = A.copy()


# Part dedicated to the plot of the dynamics:

S = dynamics_LV_2B(alpha_type_normal, A=A, x_init=x_init,
                   N=N, nbr_it=nbr_it, tau=tau)


x = np.linspace(0, nbr_it*tau, nbr_it)

fig = plt.figure(1, figsize=(10, 6))

for i in range(S.shape[0]):
    #lab = '$x_'+str(i)+'$'
    if (i < S.shape[0]//2 and i == 0):
        plt.plot(x, S[i, :], label='Community 1', color='red')
    elif(i < S.shape[0]//2 and i != 0):
        plt.plot(x, S[i, :], color='red')
    elif(i == 5):
        plt.plot(x, S[i, :], label='Community 2', color='blue')
    else:
        plt.plot(x, S[i, :], color='blue')


plt.xlabel("Time (t)", fontsize=15)
plt.ylabel("Abundance ($x_i$)", fontsize=15)
plt.legend(loc='upper left', fontsize=15)
plt.show()
plt.close()



# Part dedicated to the plot of the interaction matrix:

BMN = (A_base/np.sqrt(N)).copy()
t = 15
# BMN[5:, :5] = BMN[5:, :5]*alpha_type_linear(t, 10)
# BMN[:5, 5:] = BMN[:5, 5:]*alpha_type_linear(t, 10)

BMN[5:, :5] = BMN[5:, :5]*alpha_type_normal(t, 10)
BMN[:5, 5:] = BMN[:5, 5:]*alpha_type_normal(t, 10)

# plt.colorbar(plt.pcolor(BMN))

#sns.set_theme()

ax = sns.heatmap(BMN, vmin=-1.5, vmax=1.5, square=True, cmap="seismic",
                 xticklabels=np.arange(1, 11, 1), yticklabels=np.arange(1, 11, 1))

