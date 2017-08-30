
# coding: utf-8

# In[21]:

from sympy import *
from sympy import init_printing;  init_printing(use_latex='mathjax')
import numpy as np


# In[22]:

var('x l m hbar W')
n = int(input('Ingrese el numero de funciones: '))
funcion = []
for i in range(n):
    a = input('Ingrese la funcion: ')
    funcion.append(a)
funciones = sympify(funcion)
H = zeros(n,n)
S = zeros(n,n)
for i in range(n):
    for j in range(n):
        H[i,j] = integrate((-hbar**2/(2*m))*funciones[i]*diff(funciones[j], x, 2),(x,0,l))
        S[i,j] = integrate(funciones[i]*funciones[j],(x,0,l))
#determinante = (H-S*W).det()
Soluciones = solve((H-S*W).det(),W)
t = len(Soluciones)
for i in range(t):
    Soluciones[i] = N(Soluciones[i]*m*l**2/hbar**2)/(4*np.pi**2)
Soluciones.sort()
print('(ml^2/h^2)W = ')
Soluciones

