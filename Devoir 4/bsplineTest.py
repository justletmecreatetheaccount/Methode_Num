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
def b(t,T,i,p):
  if p == 0:
    return (T[i] <= t)*(t < T[i+1])
  else:
    u  = 0.0 if T[i+p ]  == T[i]   else (t-T[i])/(T[i+p]- T[i]) * b(t,T,i,p-1)
    u += 0.0 if T[i+p+1] == T[i+1] else (T[i+p+1]-t)/(T[i+p+1]-T[i+1]) * b(t,T,i+1,p-1)
    return u

def bspline(X,Y,t): 
  m = len(X)
  X = append(X, [X[0], X[1], X[2]])
  Y = append(Y, [Y[0], Y[1], Y[2]])
  T = arange(-3, len(X)-1, 1)
  a = b(t, T, 3, 3)
  x = []
  y = []
  quart_l = (len(a)-1)//m


  x_temp = zeros(quart_l)
  y_temp = zeros_like(x_temp)
  for i in range(m):
    for j in range(quart_l): #si quart l = 4
      x_temp[j] += a[j] * X[i] + a[j+quart_l] * X[i+1] + a[j+quart_l*2] * X[i+2] + a[j+quart_l*3] * X[i+3]
      y_temp[j] += a[j] * Y[i] + a[j+quart_l] * Y[i+1] + a[j+quart_l*2] * Y[i+2] + a[j+quart_l*3] * Y[i+3]

    x = insert(x,0, x_temp)
    y = insert(y,0, y_temp)
    x_temp = zeros(quart_l)
    y_temp = zeros_like(x_temp)
  
  x = append(x, x[0])
  y = append(y, y[0])
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
