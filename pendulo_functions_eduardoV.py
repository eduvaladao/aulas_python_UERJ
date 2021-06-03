"""
Funções importantes do pêndulo duplo: posição x e y de cada massa,
energia cinética, potencial e total 

Autor: Eduardo da Costa Valadão - eduardovaladao98@gmail.com
""" 
import numpy as np

def pendulum_x1(theta1, l1, i):
    x1 = l1*np.sin(theta1[i])
    return x1 

def pendulum_y1(theta1, l1, i):
    y1 = - l1*np.cos(theta1[i])
    return y1 

def pendulum_x2(theta2, l2, x1, i):
    x2 = x1[i] + l2*np.sin(theta2[i])
    return x2

def pendulum_y2(theta2, l2, y1, i):
    y2 = y1[i] - l2*np.cos(theta2[i])
    return y2 

def energia_cinetica(theta1, theta2, alfa1, alfa2, m1, m2, l1, l2, i):
    k = 0.5*m1*(l1**2)*(alfa1[i]**2) + 0.5*m2*((l1**2)*(alfa1[i]**2) + (l2**2)*(alfa2[i]**2) + 2*l1*l2*alfa1[i]*alfa2[i]*np.cos(theta1[i] - theta2[i]))
    return k 

def energia_potencial(theta1, theta2, m1, m2, l1, l2, g, i):
    u = -(m1 + m2)*g*l1*np.cos(theta1[i]) - m2*g*l2*np.cos(theta2[i])
    return u
