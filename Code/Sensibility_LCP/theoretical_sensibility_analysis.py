# -*- coding: utf-8 -*-
"""
@author: Maxime Clenet
"""

# %%
# Les packages prérequis:

import numpy as np
import matplotlib.pyplot as plt
from functions import block_function
import pandas as pd
import seaborn

# %%


# Sensitivity between two communities: part 1.

x = np.linspace(0.9, 7, 50)
y_f = np.linspace(1, 3, 50)
y_nonf = np.linspace(1, 3, 50)


Beta = [1/2, 1/2]
mu = np.array([[0, 0], [0, 0]])

for i in range(50):
    Alpha_f = np.array([[10**3, 10**3], [x[i], 0.9]])
    y_f[i] = block_function(mu, Alpha_f, Beta)[1]
    Alpha_nonf = np.array([[1, 10**3], [x[i], 0.9]])
    mu_nonf = np.array([[0, 0], [0, 0]])
    y_nonf[i] = block_function(mu_nonf, Alpha_nonf, Beta)[1]


fig = plt.figure(1, figsize=(10, 6))

plt.plot(x, y_f, label='$p_1 = 1$')

plt.plot(x, y_nonf, label='$p_1 = 0.8$')


plt.xlabel(
    r"Impact of community 1 on community 2 ($ \alpha_{21}$)", fontsize=15)
plt.ylabel(r"Surviving species in the 2nd community ($p_2$)", fontsize=15)

plt.legend(loc='upper right')
#plt.title('Sensibility curve')
plt.show

# plt.close()


# %%
prec = 30

x = np.linspace(0.9, 2, prec)
y = np.linspace(0.9, 2, prec)
Beta = [1/2, 1/2]
mu = np.array([[0, 0], [0, 0]])

mat = np.zeros((prec, prec))

for i, v1 in enumerate(x):
    for j, v2 in enumerate(y):
        alpha = np.array([[v2, 10**3], [v1, 1.4]])
        mat[prec-j-1, i] = block_function(mu, alpha, Beta)[1]

df = pd.DataFrame(mat, columns=np.around(x, decimals=2),
                  index=np.around(y[::-1], decimals=2))
fig = plt.figure(1, figsize=(10, 6))
seaborn.heatmap(df, cbar=True, cmap='plasma')
plt.tight_layout()
plt.xlabel(
    r"Impact of community 1 on community 2 ($ \alpha_{21}$)", fontsize=15)
plt.ylabel(r"Self-regulation of community 1 ($ \alpha_{11}$)", fontsize=15)

# %%


# Sensitivity between two communities: part 2.

x = np.linspace(0.9, 3, 50)
p1 = np.linspace(1, 3, 50)
p2 = np.linspace(1, 3, 50)


Beta = [1/2, 1/2]
mu = np.array([[0, 0], [0, 0]])

for i in range(50):
    Alpha_f = np.array([[x[i], 10**3], [0.9, 3]])
    p2[i] = block_function(mu, Alpha_f, Beta)[1]
    p1[i] = block_function(mu, Alpha_f, Beta)[0]

fig = plt.figure(1, figsize=(10, 6))

#plt.plot(x, p1_f)
plt.plot(p1, p2)

#plt.plot(x, y_nonf, label='$p_1 = 0.8$')


plt.xlabel(
    r"Surviving species in the 1st community ($p_1$)", fontsize=15)
plt.ylabel(r"Surviving species in the 1st community ($p_2$)", fontsize=15)

#plt.legend(loc='upper right')
#plt.title('Sensibility curve')
plt.show

# %%

# Feedback effect

x = np.linspace(1, 3, 50)
p1 = np.linspace(1, 3, 50)


Beta = [1/2, 1/2]
mu = np.array([[0, 0], [0, 0]])

for i in range(50):
    Alpha = np.array([[3, 1], [x[i], 3]])
    p1[i] = block_function(mu, Alpha, Beta)[0]

fig = plt.figure(1, figsize=(10, 6))

#plt.plot(x, p1_f)
plt.plot(x, p1)

#plt.plot(x, y_nonf, label='$p_1 = 0.8$')


plt.xlabel(
    r"Interaction strength of community 1 on community 2 ($\alpha_{21})$", fontsize=15)
plt.ylabel(r"Persisting species in the 1st community ($p_1$)", fontsize=15)

#plt.legend(loc='upper right')
#plt.title('Sensibility curve')
plt.show()


# %%

def f_12(x, a_11, a_21, a_22):
    return x*a_11*a_22/(a_21)


def f_21(x, a_11, a_12, a_22):
    return x*a_11*a_22/a_12


x = np.linspace(0.1, 1, 50)


y_21 = np.array([])

y_12 = np.array([])


Beta = [1/2, 1/2]
mu = np.array([[0, 0], [0, 0]])

for i in range(50):
    Alpha_f_12 = np.array([[3, f_12(x[i], 3, 1, 3)], [1, 3]])
    print(f_21(x[i], 3, 1, 3))
    y_12 = np.append(y_12, block_function(mu, Alpha_f_12, Beta)[0])

for i in range(50):
    Alpha_f_21 = np.array([[3, 1], [f_21(x[i], 3, 1, 3), 3]])
    print(f_21(x[i], 3, 1, 3))
    y_21 = np.append(y_21, block_function(mu, Alpha_f_21, Beta)[0])


fig = plt.figure(1, figsize=(10, 6))

#plt.plot(x, p1_f)
plt.plot(x, y_21, label=r"variation of $\alpha_{21}$")
plt.plot(x, y_12, label=r"variation of $\alpha_{12}$")


#plt.plot(x, y_nonf, label='$p_1 = 0.8$')


plt.xlabel(
    r"Proportion of the interaction ($\alpha_{21}\alpha_{12}/(\alpha_{11}\alpha_{22}))$", fontsize=15)
plt.ylabel(r"Surviving species in the 1st community ($p_1$)", fontsize=15)

plt.legend(loc='lower right')
#plt.title('Sensibility curve')
plt.show
