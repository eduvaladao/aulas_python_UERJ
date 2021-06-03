"""
Método numérico de Runge-Kutta para a resolução de 
equações diferenciais ordinárias de primeira, segunda 
e quarta ordem

Autor: Eduardo da Costa Valadão - eduardovaladao98@gmail.com
"""
import numpy as np

def runge_kutta_1order(y_prime, t, y_0):
    ''' 
    Dada uma função y(t) e uma edo de primeira ordem 
    y' = y_prime(y, t), com condição de contorno 
    y(0) = y_0 e dada uma lista t de intervalos iguais 
    entre um item e outro, definida pelo usuário,
    o método numérico de Runge-Kutta, é dado por:
    '''
    y = [y_0]      
    dt = t[1] - t[0]
    for i in range(len(t)-1):
        k1 = y_prime(t[i], y[i])
        k2 = y_prime(t[i] + 0.5*dt, y[i] + dt*(0.5*k1))
        k3 = y_prime(t[i] + 0.5*dt, y[i] + dt*(0.5*k2))
        k4 = y_prime(t[i] + dt, y[i] + dt*k3)
        y_iplus1 = y[i] + (k1 + 2*k2 + 2*k3 + k4)*(dt/6)
        y.append(y_iplus1)
    return y 


def runge_kutta_2order(y1_prime, y2_prime, t, y1_0, y2_0):
    ''' 
    Dada duas funções y1(t) e y2(t) com duas edos de primeira ordem
    acopladas y1' = y1_prime(y1, y2, t) e y2' = y2_prime(y1, y2, t), 
    com condições de contorno y1(0) = y1_0 e y2(0) = y2_0, e dada 
    uma lista t de intervalos iguais entre um item e outro, definida 
    pelo usuário, o método numérico de Runge-Kutta, é dado por:
    '''
    y1 = [y1_0]   
    y2 = [y2_0]    
    dt = t[1] - t[0]
    for i in range(len(t)-1):
        k1 = y1_prime(t[i], y1[i], y2[i])
        k2 = y1_prime(t[i] + 0.5*dt, y1[i] + dt*(0.5*k1), y2[i] + dt*(0.5*k1))
        k3 = y1_prime(t[i] + 0.5*dt, y1[i] + dt*(0.5*k2), y2[i] + dt*(0.5*k2))
        k4 = y1_prime(t[i] + dt, y1[i] + dt*k3, y2[i] + dt*k3)
        j1 = y2_prime(t[i], y1[i], y2[i])
        j2 = y2_prime(t[i] + 0.5*dt, y1[i] + dt*(0.5*j1), y2[i] + dt*(0.5*j1))
        j3 = y2_prime(t[i] + 0.5*dt, y1[i] + dt*(0.5*j2), y2[i] + dt*(0.5*j2))
        j4 = y2_prime(t[i] + dt, y1[i] + dt*j3, y2[i] + dt*j3)
        y1_iplus1 = y1[i] + (k1 + 2*k2 + 2*k3 + k4)*(dt/6)
        y2_iplus1 = y2[i] + (j1 + 2*j2 + 2*j3 + j4)*(dt/6)
        y1.append(y1_iplus1)
        y2.append(y2_iplus1)
    return y1, y2 


def inclinations(func, t, y, dt, i):
    k1 = func(t[i], y[0, i], y[1, i], y[2, i], y[3, i])
    k2 = func(t[i] + 0.5*dt, y[0, i] + dt*(0.5*k1), y[1, i] + dt*(0.5*k1), y[2, i] + dt*(0.5*k1), y[3, i] + dt*(0.5*k1))
    k3 = func(t[i] + 0.5*dt, y[0, i] + dt*(0.5*k2), y[1, i] + dt*(0.5*k2), y[2, i] + dt*(0.5*k2), y[3, i] + dt*(0.5*k2))
    k4 = func(t[i] + dt, y[0, i] + dt*k3, y[1, i] + dt*k3, y[2, i] + dt*k3, y[3, i] + dt*k3)
    return k1, k2, k3, k4


def runge_kutta_4order(y1_prime, y2_prime, y3_prime, y4_prime, t, y1_0, y2_0, y3_0, y4_0):
    ''' 
    Dada 4 funções y1(t), y2(t), y3(t) e y4(t) com 4 edos de primeira ordem acopladas 
    y1' = y1_prime(y1, y2, y3, y4, t), y2' = y2_prime(y1, y2, y3, y4, t), 
    y3' = y3_prime(y1, y2, y3, y4, t) e y4' = y4_prime(y1, y2, y3, y4, t) com condições 
    de contorno y1(0) = y1_0, y2(0) = y2_0, y3(0) = y3_0 e y4(0) = y4_0 num intervalo 
    de t0 a tf com dt = t[i+1] - t[i], o método numérico de Runge-Kutta, é dado por:
    '''
    dt = t[1] - t[0]
    y = np.array([len(t)*[y1_0], len(t)*[y2_0], len(t)*[y3_0], len(t)*[y4_0]])    
    for i in range(len(t)-1):
        j = 0
        for func in [y1_prime, y2_prime, y3_prime, y4_prime]:
            k1, k2, k3, k4 = inclinations(func, t, y, dt, i)
            y[j, i+1] = y[j, i] + (k1 + 2*k2 + 2*k3 + k4)*(dt/6) 
            j = j + 1
    return y[0].tolist(), y[1].tolist(), y[2].tolist(), y[3].tolist()


# Comentário: primeiro fiz uma função para resolver uma edo de primeira ordem, depois fiz uma para resolver uma de segunda e percebi que poderia simplificar e muito usando os numpy array mas deixei da forma como está para mostrar o meu passo a passo, pois na última função para resolver uma edo de quarta ordem ou duas acopladas de segunda ordem eu implementei as array, o que simplificou muito. Porém, como não possuo muito conhecimento de uso das array, não consegui implementar uma função mais geral para resolver uma edo de ordem n.


