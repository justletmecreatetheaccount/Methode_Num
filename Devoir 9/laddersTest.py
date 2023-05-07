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

def f(data, x):
  return [(data[0]**2 - (x[0] + x[1])**2) * x[1]**2 - data[2]**2 * (x[0] + x[1])**2, (data[1]**2 - (x[0] + x[1])**2) * x[0]**2 - data[2]**2 * (x[0] + x[1])**2]
  
def df(data, x):
  return [[-2 * x[0] * x[1]**2 - 2 * x[1]**3 - 2 * x[0] * data[2]**2 - 2 * x[1] * data[2]**2, 2 * x[1] * data[0]**2 - 2 * x[1] * x[0]**2 - 6 * x[0] * x[1]**2 - 4 * x[1]**3 - 2 * x[0] * data[2]**2 - 2 * x[1] * data[2]**2],
  [2 * x[0] * data[1]**2 - 2 * x[0] * x[1]**2 - 6 * x[1] * x[0]**2 - 4 * x[0]**3 - 2 * x[0] * data[2]**2 - 2 * x[1] * data[2]**2, -2 * x[1] * x[0]**2 - 2 * x[0]**3 - 2 * x[0] * data[2]**2 - 2 * x[1] * data[2]**2]]

def laddersIterate(geometry,x):
  return [x[0] - f(geometry, x)[0]/norm([df(geometry, x)[0][0], df(geometry,x)[0][1]]), x[1] - f(geometry, x)[1]/norm([df(geometry, x)[1][0], df(geometry,x)[1][1]])]

# ============================================================
  
def laddersSolve(geometry,tol,nmax):
  a = geometry[0]
  b = geometry[1]
  c = geometry[2]
  
  j = (min(a,b) - 0.1)/2
  x = [j,j]  

  for i in range(nmax):
    newx = laddersIterate(geometry, x)
    if abs(norm([x[0] - newx[0], x[1] - newx[1]])) < tol:
      return newx
    x = newx
  
  return [-1, -1]



    
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
