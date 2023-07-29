import numpy as np
import matplotlib.pyplot as plt

def calculateBSplineCurve (BSplines : np.array, controlPoints : np.array) :
    xCurve = np.zeros_like(BSplines[0])
    yCurve = np.zeros_like(BSplines[0])
    for i in range(BSplines.shape[0]) :
        xCurve += BSplines[i] * controlPoints[i][0]
        yCurve += BSplines[i] * controlPoints[i][1]
    return (xCurve, yCurve)

def arrangeControlPoints (crtlPnts : np.array) :
    rows = np.array([np.zeros(crtlPnts.shape[0]), np.zeros(crtlPnts.shape[0])])
    for i in range(crtlPnts.shape[0]) :
        rows[0][i] = crtlPnts[i][0]
        rows[1][i] = crtlPnts[i][1]
    return rows