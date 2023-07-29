import numpy as np


calculatedBSplines = {}


def calculateBSpline (t : np.array, nodes : np.array, degree : int, pos : int) :
    if (degree, pos) not in calculatedBSplines:
        finalValues = np.zeros(t.size)

        if degree == 0 :
            for i in range(t.size):
                if t[i] >= nodes[pos] and t[i] <= nodes[pos+1] :
                    finalValues[i] = 1

                else :
                    pass

            return finalValues

        else :
            leftMonomial = t.copy()
            rightMonomial = np.negative(t).copy()

            leftMonomial -= nodes[pos]
            rightMonomial += nodes[pos+1+degree]

            leftMonomial = np.divide(leftMonomial, (nodes[pos+degree]-nodes[pos]), where=(nodes[pos+degree]-nodes[pos])!=0)
            rightMonomial = np.divide(rightMonomial, (nodes[pos+degree+1]-nodes[pos+1]), where=(nodes[pos+degree+1]-nodes[pos+1])!=0)

            leftMonomial = np.multiply(leftMonomial, calculateBSpline(t, nodes, degree-1, pos))
            rightMonomial = np.multiply(rightMonomial, calculateBSpline(t, nodes, degree-1, pos+1))

            finalValues = leftMonomial + rightMonomial
            calculatedBSplines[(degree, pos)] = finalValues

            return finalValues
    else :

        return calculatedBSplines[(degree, pos)]

    return 0

def calculateAllBSplines (t : np.array, nodes : np.array, degree : int) : 
    allBSplines = np.array([])
    for i in range(nodes.size - degree - 1) :
        allBSplines = np.vstack((allBSplines, calculateBSpline(t, nodes, degree, i))) if allBSplines.size else calculateBSpline(t, nodes, degree, i)
    return allBSplines
