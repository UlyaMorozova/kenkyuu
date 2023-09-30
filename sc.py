import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import cmath
from mpmath import nsum, exp, inf
import math 
from scipy.integrate import quad
import sympy as sym


K = 5
myu = 1.25

def eq_R(p):    # ������� ��� ���������� ������ ������������ ���������
    x = sym.Symbol('x', real = True)
    sol = sym.solve(x**(K+1) - (p+1) * x + p, x)
    return sol

def F(t, lm, a):    # ��������������� �������
    return lm**a * t**(a-1)/(math.factorial(a-1))*np.exp(-lm*t)

sol = []
R = []
ro = np.linspace(0, 5, 20*K)  # ����������� �������� �������

for i in range(len(ro)):   # ����� ������ �� �����, ������� �� ������ ������ 1
    sol = eq_R(ro[i])
    for k in range(len(sol)):
        if abs(sol[k]) < 1:
            R.append(sol[k])
R = [sym.N(solution) for solution in R]

q = []
q_K1 = ((ro[K-1]/(1-R[K-1]) + nsum(lambda i: nsum(lambda k: R[K-1]**k, [0,i]), [0,K-2]))/(ro[K-1] + 1 - R[K-1]**K) + 1)**(-1)
Q_sum = 0
j = 0

while True:
    if Q_sum >= 0.999:
        break
    if 0 <= j <= K-2:
        q_j = q_K1 * nsum(lambda i: R[j]**i, [0,j]) / (ro[j] + 1 - R[j]**K)
        q.insert(j, q_j)
    elif j == K-1:
        q_j = q_K1
        q.insert(j, q_j)
    elif j == K:
        q_j = q_K1 / (ro[j] + 1 - R[j]**K)
        q.insert(j, q_j)
    else:
        q_j = q[K] * R[j]**(j - K)
        q.insert(j, q_j)
    j += 1
    Q_sum += q_j

x = np.linspace(0, 5, len(ro)) # ����������� �
lya = ro / myu
dx = 0.001

def W_function(x):  # ������� W(x)
    j = 0
    W_j = 0
    if x <= 0:
        return 0 
    t = np.append(np.arange(0,x,dx),x)
    bj = j // K
    aj = K * (bj + 1) - j - 1

    while j <= len(q):
        while aj > 0 and bj > 0:
            Ia = (dx/3)*(F(t[0],lya[j],aj)+2*np.sum(F(t[2:-2:2],lya[j],aj))+4*np.sum(F(t[1:-1:2],lya[j],aj))+F(t[-1], lya[j],aj))
            Ib = (dx/3)*(F(t[0],myu,bj)+2*np.sum(F(t[2:-2:2],myu,bj))+4*np.sum(F(t[1:-1:2],myu,bj))+F(t[-1],myu,bj))
            W_j += q[j]*Ia*Ib
            j += 1
            bj = j // K + 1
            aj = K * (bj + 1) - j - 1
            if j > len(q):
                break
                
        while aj == 0 and bj > 0:
            Ib = (dx/3)*(F(t[0],myu,bj)+2*np.sum(F(t[2:-2:2],myu,bj))+4*np.sum(F(t[1:-1:2],myu,bj))+F(t[-1],myu,bj))
            W_j += q[j]*Ib
            j += 1
            bj = j // K + 1
            aj = K * (bj + 1) - j - 1
            if j > len(q):
                break

        while aj > 0 and bj == 0:
            Ia = (dx/3)*(F(t[0],lya[j],aj)+2*np.sum(F(t[2:-2:2],lya[j],aj))+4*np.sum(F(t[1:-1:2],lya[j],aj))+F(t[-1], lya[j],aj))
            W_j += q[j]*Ia
            j += 1
            bj = j // K + 1
            aj = K * (bj + 1) - j - 1
            if j > len(q):
                break
            
        if aj == bj == 0 or j > len(q):
            W_j = 1
            break
    
    return W_j

W_j = []
for j in range(len(x)):
    w = W_function(x[j])
    W_j.append(w)
print(W_j)