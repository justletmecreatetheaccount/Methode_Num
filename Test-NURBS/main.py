import matplotlib.pyplot as plt
import numpy as np
from BSplines import *
from BSplineCurve import *

precision = 1000
degree = 3


t_axis = np.linspace(0, 2, precision)
nodes = np.array([0, 0, 0, 0, 1, 2, 2, 2, 2])
crtlPnts = np.array([(0, 0), (1, 3), (2, 0), (3, 3), (4, 0)])
rCtrlPnts = arrangeControlPoints(crtlPnts)


# plotting
fig, ax = plt.subplots()

for i in range(nodes.size-degree-1) : 
    ax.plot(t_axis, calculateBSpline(t_axis, nodes, degree, i), linewidth=2.0)
plt.show()


fig2, ax2 = plt.subplots()

x, y = calculateBSplineCurve(calculateAllBSplines(t_axis, nodes, degree), crtlPnts)

ax2.plot(x, y, linewidth=2.0)
print(x)
ax2.plot(rCtrlPnts[0], rCtrlPnts[1])


plt.show()
