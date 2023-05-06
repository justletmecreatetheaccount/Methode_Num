#
# PYTHON for DUMMIES 22-23
# Problème 1
#
# Interpolation de Padé
# Interpolation avec des quotients de polynômes
#
# Vincent Legat 
# Nathan Coppin
#
# -------------------------------------------------------------------------
# 

from numpy import *
from numpy.linalg import solve

# ============================================================
# FONCTIONS A MODIFIER [begin]
#
#

def padeInterpolationCompute(X,U):

# 
# A COMPLETER / MODIFIER
#
  matrice = []
  line = []
  for j in range(len(X)):
      line = []
      matrice.append(line)
      for i in range(len(X)//2 + 1):
          matrice[-1].append(X[j] ** i)
      for k in range(1, len(X)//2 + 1):
          matrice[-1].append(- U[j] * X[j] ** k)
  print(matrice)
  a = solve(matrice, U)
  return a

  
def padeEval(a,x) :

# 
# A COMPLETER / MODIFIER
#
  somme_num = 0
  somme_deno = 1
  it = 0
  
  for i in range(len(a)//2 + 1):
      somme_num += a[i] * x ** i
  for i in range(len(a)//2 + 1, len(a)):
      it += 1 
      somme_deno += a[i] * x ** it
  
  return somme_num / somme_deno
  
  
    
#
# FONCTIONS A MODIFIER [end]
# ============================================================

#
# -1- Test de la fonction interpolation
#     On considère un jeu des 3 fonctions u(x)
#

def main() :

  n = 2
  u =  lambda x : cos(x)   
  X = linspace(-2,2,(2*n+1)) 
  U = u(X)

#
# -1- Calcul des coefficients de l'interpolation
#
  print("==== Computing the Padé approximation :-)") 
  a = padeInterpolationCompute(X,U)
  print(" a = ",list(a))

  
#
# -2- Evaluation l'interpolation de Padé
#     et de l'interpolation polynomiale de Lagrange
#
  x = linspace(-5,5,100)
  upade = padeEval(a,x)
  uh = polyval(polyfit(X,U,len(X)-1),x)

#
# -3- Et un joli plot :-)
#

  from matplotlib import pyplot as plt
  plt.rcParams['toolbar'] = 'None'
  plt.rcParams['figure.facecolor'] = 'silver'

  plt.figure('Interpolation de Padé n = %d ' % n)
  plt.plot(x,u(x),'-b',label='cosinus')
  plt.plot(x,uh,'-g',label='Lagrangian interpolation') 
  plt.plot(x,upade,'-r',label='Padé interpolation')
  plt.plot(X,U,'ob')
  plt.xlim((-5,5)); plt.ylim((-1.2,1.2))
  plt.legend(loc='upper right') 
  plt.show()

main()

