# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 22-23
# Problème 4
#
# Script de test
#  Vincent Legat
#
# -------------------------------------------------------------------------
# 

from numpy import *

# ============================================================
# FONCTIONS A MODIFIER [begin]
#

mem = {}
def b(t,T,i,p):
  if (i,p) in mem.keys():
    return mem[(i,p)]
  if p == 0:
    u = (T[i] <= t)*(t < T[i+1])
    return u
  else:
    u  = 0.0 if T[i+p ]  == T[i]   else (t-T[i])/(T[i+p]- T[i]) * b(t,T,i,p-1)
    u += 0.0 if T[i+p+1] == T[i+1] else (T[i+p+1]-t)/(T[i+p+1]-T[i+1]) * b(t,T,i+1,p-1)
    return u
  
  
# ============================= mainProgram ===============================
 
def bspline(X,Y,t):
    print(X)
    X = [*X, *X[:3]]
    Y = [*Y, *Y[:3]]
    print("New X: ", X)
    m = len(X)
    p = 3
    n = m + p
    T = arange(-3, m + 4, 1) # m + 4 bcs its going from [-3 to m+4[ so m+3]
    B = zeros((n-p,len(t)))
    for i in range(0,n-p):
      B[i,:] = b(t,T,i,p)
    x = X @ B
    y = Y @ B
    return x,y

   
#
# FONCTIONS A MODIFIER [end]
# ============================================================

def main() :

#
# -1- Approximation d'un rectangle :-)     
#

  X = [0,(1/4)*(10+2*5**(1/2))**(1/2),(1/4)*(10-2*5**(1/2))**(1/2),-(1/4)*(10-2*5**(1/2))**(1/2)]#,-(1/4)*(10+2*5**(1/2))**(1/2)]
  Y = [1,(1/4)*(5**(1/2)-1),-(1/4)*(5**(1/2)+1),-(1/4)*(5**(1/2)+1)]#,(1/4)*(5**(1/2)-1)]
  t = linspace(0,len(X),len(X)*100 + 1)

      
  x,y = bspline(X,Y,t)

#
# -2- Un joli dessin :-)
#

  import matplotlib.pyplot as plt
  import matplotlib 
  matplotlib.rcParams['toolbar'] = 'None'
  plt.rcParams['figure.facecolor'] = 'white'

  fig = plt.figure("Approximation avec des B-splines")
  plt.plot(X,Y,'.r',markersize=10)
  plt.plot([*X,X[0]],[*Y,Y[0]],'--r')
  plt.plot(x,y,'-b')
  plt.axis("equal"); plt.axis("off")
  #plt.show()

if __name__ == '__main__': 
  main()
