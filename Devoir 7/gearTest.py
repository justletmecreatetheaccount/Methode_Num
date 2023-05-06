# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 22-23
# Problème 7
#
# Calcul de la zone de stabilité d'une méthode de Gear
#  Vincent Legat
#
# -------------------------------------------------------------------------
#

import numpy as np
from numpy import *

# ============================================================
# FONCTIONS A MODIFIER [begin]
#
# -1- La zone de stabilité de Gear d'ordre un (qui est la méthode
#     d'Euler implicite :-) vous est offerte gracieusement
#     par les descendants de la famille Gear !
#
#     Evidemment, c'est plus compliqué pour les autres ordres !
#  

def stabilityGear(x,y,order):

  coeff = {}
  coeff[1] = array([ 1,  1,                    ])   
  coeff[2] = array([ 2,  4,-  1,               ])/3
  coeff[3] = array([ 6, 18,-  9,  2,           ])/11
  coeff[4] = array([12, 48,- 36, 16,-  3,      ])/25
  coeff[5] = array([60,300,-300,200,- 75,12,   ])/137
  coeff[6] = array([60,360,-450,400,-225,72,-10])/147
  
  z = x + 1j*y
  g = zeros(shape(z))
  for i in range(len(x)) :
    for j in range(len(y)) :
       g[i,j] = max(abs(roots([(coeff[order][0]*z[i,j] - 1),*coeff[order][1:]])))

  if (order > 6) :
    return zeros(order + 1),zeros(shape(x))

  return g,coeff[order]

#
# FONCTIONS A MODIFIER [end]
# ============================================================
#
# ------------------------------------------------------------------------------------ 
#
# Script de test
#
#
# -1- Choix de l'ordre et de la précision de la carte de stabilité
#
# ------------------------------------------------------------------------------------ 

def main() : 
  order = 3
  n = 100
  x,y = np.meshgrid(np.linspace(-8,5,n),np.linspace(-8,8,n))

  
  gain,coeff = stabilityGear(x,y,order)

#
# -2- Impression des coefficients de la méthode de Gear
#     Observer l'utilisation du module "fractions" qui permet
#     d'écrire des nombres en virgule flottante sous la forme de fractions !
#     
#     Retirer la ligne np.set_printoptions
#     Noter aussi les astuces dans les paramètres :-)
#     Merci à : https://stackoverflow.com/questions/42209365/numpy-convert-decimals-to-fractions
#
#

  import fractions
  np.set_printoptions(formatter={'all':lambda x: str(fractions.Fraction(x).limit_denominator())})
  print("==== Coefficients of Gear's method for order = %d ========== " % order)
  print("     ",end='')
  print(coeff)


#
# -3- Faire le joli plot
#

  import matplotlib.pyplot as plt
  import matplotlib
  matplotlib.rcParams['toolbar'] = 'None'
  plt.rcParams['figure.facecolor'] = 'lavender'
  plt.rcParams['axes.facecolor'] = 'lavender'
  plt.figure("Stability of Gear's method : order = %d" % order)
  plt.contourf(x,y,gain,np.arange(0,1.1,0.1),cmap=plt.cm.jet_r)
  plt.contour(x,y,gain,np.arange(0,1.1,0.1),colors='black',linewidths=0.5)
  ax = plt.gca()
  ax.axhline(y=0,color='r')
  ax.axvline(x=0,color='r')
  ax.yaxis.grid(color='gray',linestyle='dashed')
  ax.xaxis.grid(color='gray',linestyle='dashed')
  ax.set_aspect('equal', 'box')
  plt.show()
  
  
main() 