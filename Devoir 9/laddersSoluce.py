# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 22-23
# Problème 9
#
# Script de test
#  Vincent Legat
#
# -------------------------------------------------------------------------
#

from scipy.linalg import norm,solve

# ============================================================
# FONCTIONS A MODIFIER [begin]
#
# xn+1 = xn + fxn/dfxn

def f(u):
  global a,b,c
  x = u[0]; y = u[1]
  return [a*a*x*x -(x*x+c*c)*((x+y)**2),
          b*b*y*y -(y*y+c*c)*((x+y)**2)]

def dfdx(u):
  global a,b,c
  x = u[0]; y = u[1]
  return [[2*a*a*x - 2*x*((x+y)**2) - 2*(x*x+c*c)*(x+y) ,       -2*(x*x+c*c)*(x+y) ],
          [-2*(y*y+c*c)*(x+y) ,      2*b*b*y - 2*y*((x+y)**2) - 2*(y*y+c*c)*(x+y) ]]


# -------------------------------------------------------------------------

def laddersSolve(geometry,tol,nmax):
  global a,b,c
  a = geometry[0]
  b = geometry[1]
  c = geometry[2]

#
# Astuce éventuelle pour le cas a=b qui rend la jacobienne de Newton-Raphson singulière
# Dans ce cas précis, on peut obtenir directement une solution analytique
# 
# On utilise la tolérance pour le test du cas particulier
# Dans le cas d'une itération simple, il a été nécessaire 
# d'intéger une tolérance absolue (beurck :-)
#
# 
#   if (abs(a-b) < tol) :
#     x = ((b*b - 4*c*c)**(1/2))/2
#     return solve([[1,0],[0,1]],[x,x])
# 
# En pratique, solve fournira la bonne réponse même sans y ajouter l'astuce
#

#
# Heuristique de Nathan pour obtenir un bon candidat initial
# Cela fonctionne remarquablement bien :-)
#

  r = (min(a,b) - 0.1)/2
  x = [r,r]  
  
  n = 0; delta = tol+1
  while (norm(delta) > tol and n < nmax):
    delta = - solve(dfdx(x),f(x))
    x = x + delta
    n = n+1
    print(" Estimated error %9.2e at iteration %d : " % (norm(delta),n),x)
  return x
  
# -------------------------------------------------------------------------

def laddersIterate(geometry,x):
  global a,b,c
  a = geometry[0]
  b = geometry[1]
  c = geometry[2]

  
  delta = - solve(dfdx(x),f(x))
  x = x + delta
  return x



    
#
# FONCTIONS A MODIFIER [end]
# ============================================================


def main():
  
# ------------------------------------------------------------------------------------ 
#
# Script de test
#
#
# ------------------------------------------------------------------------------------
#
# -1- Calcul de l'écart entre les deux murs :-)
#

  import numpy as np
  geometry = [3,4,1]
  print(" ========= my Newton-Raphson scheme with your proposed step :-)")

  x = np.array([1.0,1.5]); tol = 10e-12; nmax = 50
  n = 0; delta = tol+1
  while (norm(delta) > tol and n < nmax):
    xold = x
    x = laddersIterate(geometry,xold)
    delta = [x[0]-xold[0],x[1]-xold[1]]; n = n+1
    print(" Estimated error %9.2e at iteration %d : " % (norm(delta),n),x)
  print(" Computed distance is : %13.6f " % sum(x))


  print(" ========= your full computation :-)")
  sol = laddersSolve(geometry,1e-14,50)
  print(" Computed distance is : %13.6f " % sum(sol))

  a = geometry[0]
  b = geometry[1]
  c = geometry[2]
  ab = max(a,b)

#
# -2- Et un joli dessin
#

  import matplotlib.pyplot as plt
  import matplotlib
 
  matplotlib.rcParams['toolbar'] = 'None'
  plt.rcParams['figure.facecolor'] = 'lavender'

  plt.figure("Ladders geometry")
  x = sol[0]; y = sol[1]; d = x + y
  hx = np.sqrt(b*b - d*d); hy = np.sqrt(a*a - d*d)
  plt.plot([-x,y],[hx,0],'-r')
  plt.plot([-x,y],[0,hy],'-b')
  plt.plot([-x,-x,y,y],[ab,0,0,ab],'k')
  plt.axis('equal')
  ax = plt.gca()
  ax.yaxis.grid(color='gray',linestyle='dashed')
  ax.xaxis.grid(color='gray',linestyle='dashed')
  plt.xticks(np.arange(-ab,ab+1,1))
  plt.yticks(np.arange(0,ab+1,1))
  plt.show()
  
  
  
main()  
