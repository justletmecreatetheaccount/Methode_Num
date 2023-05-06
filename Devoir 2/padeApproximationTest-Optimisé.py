#
# PYTHON for DUMMIES 22-23
# Problème 2
#
# Approximant de Padé
# ou comment estimer mieux le cosinus qu'avec un développement
# de Mac-laurin : Taylor à l'origine :-)
#
# Vincent Legat 
# Nathan Coppin
#
# -------------------------------------------------------------------------
# 

from math import factorial
import numpy as np
#from numpy.linalg import solve
from perfs import perfs

np.set_printoptions(suppress=True,linewidth=np.nan)

# ============================================================
# FONCTION A MODIFIER [begin]
#
#
@perfs
def padeApproximationCompute(n,dU) :
# 
# A COMPLETER / MODIFIER
#
  fact = 1
  dU[0] = - dU[0]
  for i in range(2*n):
    fact *= i + 1
    dU[i+1] = - dU[i+1] / fact
  for i in range(n):
    dU = np.insert(dU, 0, 0)
  first_half = np.eye(2*n + 1, n + 1)
  second_half = np.array([np.flip(dU[i:n+i]) for i in range(2*n + 1)], dtype=np.float64)

  dU = dU.astype("float64")
  matrix = np.append(first_half, second_half, axis=1).astype("float64")
  return np.linalg.solve(matrix, - dU[n:])


#
# FONCTION A MODIFIER [end]
# ============================================================
#
#     Calcul de développement de Mac-Laurin (Taylor à l'origine)
#     en utilisant le calcul symbolique
#

from sympy import *

def macLaurinCompute(u,x,n,X) :
  ut = 0
  Ut = 0
  dU = np.zeros(n+1)
  for i in range(0,n+1):
    dudx=diff(u,x,i) 
    dUdx=dudx.subs(x,0)  
    dU[i] = dUdx   
    if dUdx != 0: 
       term = dUdx*(x**i)/factorial(i)
       Ut += term.subs(x,X)
       ut += term
  return [ut,Ut,dU] 
  
#
#     Evaluation d'un quotient de Padé 
#
 

def padeEval(a,x) :
  n = len(a) // 2
  chi = np.array([x**i for i in range(1,n+1)]).T
  return (a[0] + chi@a[1:n+1])/(1 + chi@a[n+1:]) 

def main() :

  x = symbols('x') 
  u = cos(x)
  n = 8
  X  = 1.5 
  U = u.subs(x,X)
  
  #
  # -2- Calcul du développement de MacLaurin
  #

  print(" ==== Obtaining MacLaurin expansion of order %d for %s" % (n,u))

  [ut,Ut,dU]  = macLaurinCompute(u,x,n,X)
     
  print("  MacLaurin expansion of %s : %s" % (u,str(ut)))
  print(" ======== Derivatives of %s for x = 0 : " % u, end=''); print(dU)
  print(" ======== Value of %s for x = %4.2f : %.16f" % (u,X,U))
  print(" ======== MacLaurin expansion of order %d  for x = %4.2f : %.16f" % (n,X,Ut))
  print("  Approximation error : %14.7e" % (U-Ut))

  #
  # -2- Calcul de l'approximation de Pade
  #

  n = int(n/2)
  a = padeApproximationCompute(n,dU)
  Up = padeEval(a,X)


  print(" ======== Coefficients of Pade approximation of order %d for %s : " % (n,u))
  print("  ",a[:n+1]);
  print("  ",a[n+1:]);
  #print(" ======== Pade approximation for x = %4.2f : %.16f" % (X,Up))
  #print("  Approximation error : %14.7e" % (U-Up))


  #
  # -3- Le Taylor en compétition :-)
  #

  [ut,Ut,dU]  = macLaurinCompute(u,x,6,X)
     
  print(" ======== MacLaurin expansion of order %d  for x = %4.2f : %.16f" % (6,X,Ut))
  #print("  Approximation error : %14.7e" % (U-Ut))



  #
  # -4- Et un joli plot :-)
  #


  from sympy.utilities.lambdify import lambdify
  from matplotlib import pyplot as plt
  plt.rcParams['toolbar'] = 'None'
  plt.rcParams['figure.facecolor'] = 'lavender'
  plt.rcParams['axes.facecolor'] = 'lavender'


  u  = lambdify(x,u,'numpy') 
  ut = lambdify(x,ut,'numpy') 

  m = 100;
  x = np.linspace(-5,5,m)

  #
  #   Solution de référence :-)
  #
     
  uReference = (15120 - 6900*x*x + 313*x**4)/(15120 +660*x**2 + 13*x**4)
 
  plt.figure("Approximant de Padé : n=%d" % n)
  plt.plot(x,uReference,'--k',label='Référence')
  plt.plot(x,ut(x),'-r',label='Taylor')
  plt.plot(x,u(x),'-b',label='u')
  plt.xlim((-5,5)); plt.ylim((-1.2,1.2))
  plt.plot([0],[1],'ob')
  plt.plot(x,padeEval(a,x),'-k',label='Padé')

  plt.legend(loc='upper right')
  plt.show()


main()
