"""
Simulação de um pêndulo duplo utilizando-se o método RK4 com input pelo terminal, gera uma lista do tempo inicial ao final 
com um intervalo entre cada item definido pelo usuário e então a partir das condições de contorno iniciais, massas e tamanhos
das hastes do pêndulo gera listas: ângulos, velocidades angulares, posições x e y, energia cinética, potencial e total de cada
uma das 2 massas. Além disso, gera gráficos: ângulo x tempo, velocidade angular x tempo, ângulo de uma massa x ângulo de outra,  
velocidade angular de uma massa x velocidade angular de outra, velocidade angular x ângulo, 'x' x 'y', energia (cinética, 
potencial e total) x tempo para cada uma das massas.

Módulos utilizados: matplotlib e numpy;
Módulos criados: runge_kutta_eduardoV onde se encontra o método RK4 e pendulo_functions_eduardoV onde estão funções relacionadas
ao pêndulo duplo, como posições cartesianas e energias. 

Animação do pêndulo duplo: tirado de 'https://matplotlib.org/stable/gallery/animation/double_pendulum.html' 

Autor: Eduardo da Costa Valadão - eduardovaladao98@gmail.com
"""
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from collections import deque
from runge_kutta_eduardoV import * 
from pendulo_functions_eduardoV import *

print ("*---------------------------------------------------------------------------*")
print ("| Simulação de um pêndulo duplo utilizando-se o método de Runge-Kutta (RK4) |")
print ("*---------------------------------------------------------------------------*")

a = float(input('Insira o comprimento da haste ligado a massa 1: '))
b = float(input('Insira o valor da massa 1: '))
c = float(input('Insira o comprimento da segunda haste: '))
d = float(input('Insira o valor da massa 2: '))
t0 = float(input('Escolha um tempo inicial: '))
tf = float(input('Escolha um tempo final: '))
x = float(input('Escolha o intervalo entre cada item da lista de tempo: '))
gr = float(input('Escolha a gravidade local g: '))

print ("*-----------------------*")
print ("| Condições de contorno |")
print ("*-----------------------*")

theta1_0 = float(input('Insira o valor de theta1 em t = 0: '))
alfa1_0 = float(input('Insira o valor da derivada de theta1 em t = 0: '))
theta2_0 = float(input('Insira o valor de theta2 em t = 0: '))
alfa2_0 = float(input('Insira o valor da derivada de theta2 em t = 0: '))

#-------------------------------------------------------------------------------------------------------------------------------------------
# Resolução das duas edos de segunda ordem acopladas do pêndulo duplo
#-------------------------------------------------------------------------------------------------------------------------------------------
def theta1_prime(t, theta1, theta2, alfa1, alfa2):
    return alfa1

def theta2_prime(t, theta1, theta2, alfa1, alfa2):
    return alfa2

def alfa1_prime(t, theta1, theta2, alfa1, alfa2, l1 = a, l2 = c, m1 = b, m2 = d, g = gr):
    return (m2*(g*np.cos(theta1 - theta2)*np.sin(theta2) - np.sin(theta1 - theta2)*(l1*(alfa1**2)*np.cos(theta1 - theta2) + l2*(alfa2**2))) - (m1 + m2)*g*np.sin(theta1))/(l1*((m1 + m2) - m2*(np.cos(theta1 - theta2)**2)))

def alfa2_prime(t, theta1, theta2, alfa1, alfa2, l1 = a, l2 = c, m1 = b, m2 = d, g = gr):
    return ((m1 + m2)/(l2*((m1 + m2) - m2*(np.cos(theta1 - theta2)**2))))*(g*(np.sin(theta1)*np.cos(theta1 - theta2) - np.sin(theta2)) + l1*(alfa1**2)*np.sin(theta1 - theta2) + (m2*l2*(alfa2**2)*np.sin(theta1 - theta2)*np.cos(theta1 - theta2))/(m1 + m2))

t = np.arange(t0, tf+1, x)
theta1, theta2, alfa1, alfa2 = runge_kutta_4order(theta1_prime, theta2_prime, alfa1_prime, alfa2_prime, t, theta1_0, theta2_0, alfa1_0, alfa2_0)
#-------------------------------------------------------------------------------------------------------------------------------------------

# listas das posições cartesianas de cada massa
x1 = []
y1 = []
x2 = []
y2 = []
for i in range(len(theta1)):
    a1 = pendulum_x1(theta1, a, i)
    b1 = pendulum_y1(theta1, a, i)
    x1.append(a1)
    y1.append(b1) 
    a2 = pendulum_x2(theta2, c, x1, i)
    b2 = pendulum_y2(theta2, c, y1, i)
    x2.append(a2)
    y2.append(b2)

# listas das diferentes energias do pêndulo duplo
cinetica = []
potencial = []
total = []
for j in range(len(t)):
    a3 = energia_cinetica(theta1, theta2, alfa1, alfa2, b, d, a, c, j)
    b3 = energia_potencial(theta1, theta2, b, d, a, c, gr, j)
    ab3 = a3 + b3
    cinetica.append(a3)
    potencial.append(b3)
    total.append(ab3)

# função para salvar os gráficos
def save_or_not(nome):
    save = input('Deseja salvar o gráfico a seguir? [sim/nao]: ')
    if save == 'sim':
        plt.savefig(nome, dpi=300, bbox_inches='tight')
        print('O gráfico foi salvo.')
    elif save == 'nao':
        print('O gráfico não será salvo.')
    else:
        print('Erro: valor não aceito.')
    

# gráfico de theta1 x t
plt.plot(t, theta1, 'r')
plt.xlabel(r'$t (s)$')
plt.ylabel(r'$\theta_{1} (radianos)$')
save_or_not('plot1.png')
plt.show()

# gráfico de theta2 x t
plt.plot(t, theta2, 'b')
plt.xlabel(r'$t (s)$')
plt.ylabel(r'$\theta_{2} (radianos)$')
save_or_not('plot2.png')
plt.show()

# gráfico de alfa1 x t
plt.plot(t, alfa1, 'r')
plt.xlabel(r'$t (s)$')
plt.ylabel(r'$\alpha_{1} (radianos)$')
save_or_not('plot3.png')
plt.show()

# gráfico de alfa2 x t
plt.plot(t, alfa2, 'b')
plt.xlabel(r'$t (s)$')
plt.ylabel(r'$\alpha_{2} (radianos)$')
save_or_not('plot4.png')
plt.show()

# gráfico de theta2 x theta1
plt.plot(theta1, theta2, 'purple')
plt.xlabel(r'$\theta_{1} (radianos)$')
plt.ylabel(r'$\theta_{2} (radianos)$')
save_or_not('plot5.png')
plt.show()

# gráfico de alfa2 x alfa1
plt.plot(alfa1, alfa2, 'black')
plt.xlabel(r'$\alpha_{1} (radianos)$')
plt.ylabel(r'$\alpha_{2} (radianos)$')
save_or_not('plot6.png')
plt.show() 

# gráfico de alfa1 x theta1
plt.plot(theta1, alfa1, 'r')
plt.xlabel(r'$\theta_{1} (radianos)$')
plt.ylabel(r'$\alpha_{1} (radianos/s)$')
save_or_not('plot7.png')
plt.show()

# gráfico de alfa2 x theta2
plt.plot(theta2, alfa2, 'b')
plt.xlabel(r'$\theta_{2} (radianos)$')
plt.ylabel(r'$\alpha_{2} (radianos/s)$')
save_or_not('plot8.png')
plt.show()

# gráfico de y por x das duas massas
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(x1, y1, 'r')
plt.plot(x2, y2, 'b')
plt.axis([-(a + c + 0.25), (a + c + 0.25), -(a + c + 0.25), (a + c + 0.25)])
ax.set_aspect('equal', adjustable='box')   # imagem quadrada
plt.xlabel(r'$x (m)$')
plt.ylabel(r'$y (m)$')
save_or_not('plot9.png')
plt.show()

# gráfico da energia cinética, potencial e total x t
plt.plot(t, cinetica, 'g', label = 'Energia Cinética')
plt.plot(t, potencial, 'purple', label = 'Energia Potencial')
plt.plot(t, total, 'b', label = 'Energia Total')
plt.xlabel(r'$t (s)$')
plt.ylabel(r'$E (J)$')
plt.legend()
save_or_not('plot10.png')
plt.show()


#-------------------------------------------------------------------------------------------------------------------------------------------
# Animação do pêndulo duplo baseado num código do próprio site do matplotlib
#-------------------------------------------------------------------------------------------------------------------------------------------
history_len = 500
fig = plt.figure(figsize=(5, 4))
ax1 = fig.add_subplot(autoscale_on=False, xlim=(-(a + c + 0.25), (a + c + 0.25)), ylim=(-(a + c + 0.25), (a + c + 0.25)))
ax1.set_aspect('equal')
ax1.grid()

line, = ax1.plot([], [], 'o-', lw=2)
trace, = ax1.plot([], [], ',-', lw=1)
time_template = 'time = %.1fs'
time_text = ax1.text(0.05, 0.9, '', transform=ax1.transAxes)
history_x, history_y = deque(maxlen=history_len), deque(maxlen=history_len)


def animate(i):
    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]

    if i == 0:
        history_x.clear()
        history_y.clear()

    history_x.appendleft(thisx[2])
    history_y.appendleft(thisy[2])

    line.set_data(thisx, thisy)
    trace.set_data(history_x, history_y)
    time_text.set_text(time_template % (i*x))
    return line, trace, time_text


ani = animation.FuncAnimation(
    fig, animate, len(theta1), interval=1, blit=True)
plt.show()
#ani.save('anima_pendulo_simples.gif', writer='pillow', fps=60)
#-------------------------------------------------------------------------------------------------------------------------------------------
