import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0,2,1000)


beta1 = 1/2
beta2 = 1-beta1

kappa = 0.8



def f1(k,kappa,beta1,beta2):
    return np.sqrt((kappa-beta2*k**2-1/2)/beta1)

def f2(k,beta1,beta2):
    return np.sqrt((1-2*beta2*k**2)/(2*beta1))

y1 = np.zeros(1000)
y2 = np.zeros(1000)

for i in range(1000):
    y1[i] = f1(x[i],kappa,beta1,beta2)
    y2[i] = f2(x[i],beta1,beta2)

plt.plot(x,y1,linestyle = 'dashed')
plt.plot(x,y2)

plt.xlim(0, 1)
plt.ylim(0, 1)


