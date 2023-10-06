```python
import numpy as np
import cmath
from mpmath import nsum, exp, inf
import math 
from scipy.integrate import quad
import sympy as sym


K = 5
myu = 1.25
ro = 1
RO = 5
```


```python
def eq_R(p):    # Функция для нахождения корней комплексного уравнения
    x = sym.Symbol('x', real = True)
    sol = sym.solve(x**(K+1) - (p+1) * x + p, x)
    return sol
```


```python
# Для RO = 1

sol1 = []
R1 = []

# Берем только те корни, которые по модулю меньше 1
sol1 = eq_R(ro)
for k in range(len(sol)):
    if abs(sol1[k]) < 1:
        R1.append(sol1[k])
R1 = [sym.N(solution) for solution in R1]

```


```python
R1
```




    [0.508660391642004]




```python
q = []
q_K1 = ((ro/(1-R1[0]) + nsum(lambda i: nsum(lambda k: R1[0]**k, [0,i]), [0,K-2]))/(ro + 1 - R1[0]**K) + 1)**(-1)
Q_sum = 0

for j in range(11):
    if Q_sum >= 0.999:
        break
    if 0 <= j <= K-2:
        q_j = q_K1 * nsum(lambda i: R1[0]**i, [0,j]) / (ro + 1 - R1[0]**K)
        q.insert(j, q_j)
    elif j == K-1:
        q_j = q_K1
        q.insert(j, q_j)
    elif j == K:
        q_j = q_K1 / (ro + 1 - R1[0]**K)
        q.insert(j, q_j)
    elif j > K:
        q_j = q[K] * R1[0]**(j - K)
        q.insert(j, q_j)
    Q_sum += q_j
print(Q_sum)
q
```

    0.996535843343198
    




    [0.0982679216715992,
     0.148252921194921,
     0.173678310628679,
     0.186611199175705,
     0.193189647329097,
     0.0982679216715992,
     0.0499849995233214,
     0.0254253894337581,
     0.0129328885470259,
     0.00657844815339256,
     0.00334619601410128]




```python
# Для RO = 5

sol2 = []
R2 = []

# Берем только те корни, которые по модулю меньше 1
sol2 = eq_R(RO)
for k in range(len(sol)):
    if abs(sol2[k]) < 1:
        R2.append(sol2[k])
R2 = [sym.N(solution) for solution in R2]

```


```python
sol2
```




    [1,
     -1/2 - sqrt(-1 + 5/(2*(-5/8 + 5*I/4)**(1/3)) + 2*(-5/8 + 5*I/4)**(1/3))/2 - sqrt(-2 - 2*(-5/8 + 5*I/4)**(1/3) + 4/sqrt(-1 + 5/(2*(-5/8 + 5*I/4)**(1/3)) + 2*(-5/8 + 5*I/4)**(1/3)) - 5/(2*(-5/8 + 5*I/4)**(1/3)))/2,
     -1/2 + sqrt(-1 + 5/(2*(-5/8 + 5*I/4)**(1/3)) + 2*(-5/8 + 5*I/4)**(1/3))/2 - sqrt(-2 - 2*(-5/8 + 5*I/4)**(1/3) - 4/sqrt(-1 + 5/(2*(-5/8 + 5*I/4)**(1/3)) + 2*(-5/8 + 5*I/4)**(1/3)) - 5/(2*(-5/8 + 5*I/4)**(1/3)))/2,
     -1/2 + sqrt(-2 - 2*(-5/8 + 5*I/4)**(1/3) - 4/sqrt(-1 + 5/(2*(-5/8 + 5*I/4)**(1/3)) + 2*(-5/8 + 5*I/4)**(1/3)) - 5/(2*(-5/8 + 5*I/4)**(1/3)))/2 + sqrt(-1 + 5/(2*(-5/8 + 5*I/4)**(1/3)) + 2*(-5/8 + 5*I/4)**(1/3))/2,
     -1/2 + sqrt(-2 - 2*(-5/8 + 5*I/4)**(1/3) + 4/sqrt(-1 + 5/(2*(-5/8 + 5*I/4)**(1/3)) + 2*(-5/8 + 5*I/4)**(1/3)) - 5/(2*(-5/8 + 5*I/4)**(1/3)))/2 - sqrt(-1 + 5/(2*(-5/8 + 5*I/4)**(1/3)) + 2*(-5/8 + 5*I/4)**(1/3))/2]




```python
R2
```




    []


