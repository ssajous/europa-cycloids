import numpy as np

# evaluates cubic bezier at t, return point
def q(ctrlPoly, t):
    return (1.0-t)**3 * ctrlPoly[0] + 3*(1.0-t)**2 * t * ctrlPoly[1] + \
        3*(1.0-t)* t**2 * ctrlPoly[2] + t**3 * ctrlPoly[3]


# evaluates cubic bezier first derivative at t, return point
def qprime(ctrlPoly, t):
    return 3*(1.0-t)**2 * (ctrlPoly[1]-ctrlPoly[0]) + 6*(1.0-t) * t * \
        (ctrlPoly[2]-ctrlPoly[1]) + 3*t**2 * (ctrlPoly[3]-ctrlPoly[2])


# evaluates cubic bezier second derivative at t, return point
def qprimeprime(ctrlPoly, t):
    return 6*(1.0-t) * (ctrlPoly[2]-2*ctrlPoly[1]+ctrlPoly[0]) + 6*(t) * (ctrlPoly[3]-2*ctrlPoly[2]+ctrlPoly[1])

# evaluates linear bezier at t, return point
def bezierLinearInterp(points, time):
    return (1 - time)*points[0] + time * points[1]

# evaluates a cubic bezier at t to find the r line, return start and end points
def findCubicRPoints(points, time):
    q0 = bezierLinearInterp(np.array([points[0], points[1]]), time)
    q1 = bezierLinearInterp(np.array([points[1], points[2]]), time)
    q2 = bezierLinearInterp(np.array([points[2], points[3]]), time)
    
    r1 = bezierLinearInterp(np.array([q0, q1]), time)
    r2 = bezierLinearInterp(np.array([q1, q2]), time)
    
    return np.array([r1, r2])