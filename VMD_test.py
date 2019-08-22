
# coding: utf-8


from VMD import VMD
import numpy as np
import math
import matplotlib.pyplot as plt
f_1 = 2.0
f_2 = 24.0
f_3 = 288.0
T=1000
t=np.linspace(1/float(T),1,T,endpoint=True)
pi=np.pi
# print(t)
v_1 =np.cos(2*pi*f_1*t)
v_2 = 1/4.0*np.cos(2*pi*f_2*t)
v_3 = 1/16.0*np.cos(2*pi*f_3*t)
# print(v_1)
# print(v_2)
# print(v_3)
alpha = 2000.0     
tau = 0         
K = 3           
DC = 0         
init = 1         
tol = 1e-7
signal=v_1+v_2+v_3
(u,u_hat,omega)=VMD(signal, alpha, tau, K, DC, init, tol)
print(u.shape)
print(u_hat.shape)
print(omega.shape)

