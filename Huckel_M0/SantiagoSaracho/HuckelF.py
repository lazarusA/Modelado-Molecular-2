from sympy import *
from sympy import init_printing; init_printing(use_latex = 'mathjax')
import numpy as np

#Definicion de variables y funciones
var('x l m hbar W k') #Variables a utilizar

def funciones():
    n = int(input('Indica la cantidad de funciones de prueba: '))
    F=[]
    for i in range (0,n):
        F.append(str(input('Escribre tu funcion %d: '  %(i+1))))
    
    F=sympify(F)
    F=Matrix(F)
    return F

def hami(F): #Arreglo de funciones y cantidad de funciones como argumentos
    H='(-hbar**2/(2*m))*diff(F, x, 2)+'
    H= H + input(str('Escribe el potencial aplicado a F, la parte cinética del Hamiltoniano ya está incluida: '))
    H=eval(H)
    return H

def HW(F,hami):
    ff=F*F.T
    fHf=F*hami.T
    H=Integral(fHf, (x,0,l)).doit()
    S=Integral(ff, (x,0,l)).doit()
    H_SW=H-S*W
    Wn=solve(H_SW.det(),W)
    
    A=N(Matrix(Wn)/Wn[0])
    B=np.array(A)
    C=sorted(B)
    WnSort=zeros(1,len(F))
    for i in range (0,len(F)):
        for j in range (0,len(F)):
            if B[i]==C[j]:
                WnSort[i]=Wn[j]
    
    
    return H, S, H_SW, Wn, WnSort

def FuncPhi(F,H,S,H_SW,Wn):
    #Generalizar para todas las Wn

    n=len(F)
    H_SWn = np.zeros((n,n,n), dtype=object) #Matriz (Hij - Sij*Wn)cj

    c_gral=zeros(n,n) #Matriz de coeficientes c
    for i in range (0,n):
        for j in range (0,n):
            c_gral[i,j]=sympify((str('c_' +str(i+1) + '_' + str(j+1))))

    for z in range (0,n):
        for i in range (0,n):
            for j in range (0,n):
                H_SWn[z,i,j]=sympify(((H[i,j]-S[i,j]*Wn[z])*c_gral[z,j])) 
    
    eq_gral=zeros(n,n) #Matriz de ecuaciones para obtener las soluciones a las c
    for z in range (0,n):
        for i in range (0,n):
            eqn=0
            for j in range (0,n):
                eqn=eqn+H_SWn[z,i,j]
            
            eq_gral[z,i]=eqn
    eq_gral=sympify(eq_gral)

    sol_gral=[]
    for i in range (0,n):
        sol_gral.append(solve(eq_gral[i,:], c_gral[i,:], dict=true))

    soluciones=zeros(n,n)
    if n != 1:
        for i in range(0, n):
            for j in range (0,n):
                if sympify('c_'+str(i+1)+'_'+str(j+1)) in sol_gral[i][0]:
                    soluciones[i,j]= sol_gral[i][0][c_gral[i,j]]
                else:
                    #soluciones[i,j] = c_gral[i,j]
                    soluciones[i,j] = k
                    soluciones[i,:]=soluciones[i,:].subs(c_gral[i,j],k)
    else:
        soluciones[0,0]=F*k
    #Normalizar
    Phi_n=zeros(n,1)
    Norm_gral=zeros(n,1)
    for j in range (0,n):
        for i in range (0,n):
            Phi_n[j,0]=Phi_n[j,0] + soluciones[j,i]*F[i] 
    
    for i in range (0,n):
        Norm_gral[i,0]=Integral(Phi_n[i,0]**2, (x,0,l)).doit()-1
        K=(solve(Norm_gral[i,0],k)[0]).evalf() #Valor positivo de k 
        Phi_n[i,0]=(Phi_n[i,0].subs(k,K)).evalf()

    return Phi_n


def variacional():
    F=funciones()
    ham=hami(F)
    H,S,H_SW,Wn,WnSort=HW(F,ham)
    f=FuncPhi(F,H,S,H_SW,WnSort)
    return WnSort,f