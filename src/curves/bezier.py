import numpy as np


# evaluates cubic bezier at t, return point
def q(ctrl_poly, t):
    return (1.0 - t) ** 3 * ctrl_poly[0] + 3 * (1.0 - t) ** 2 * t * ctrl_poly[1] + \
           3 * (1.0 - t) * t ** 2 * ctrl_poly[2] + t ** 3 * ctrl_poly[3]


# evaluates cubic bezier first derivative at t, return point
def q_prime(ctrl_poly, t):
    return 3 * (1.0 - t) ** 2 * (ctrl_poly[1] - ctrl_poly[0]) + 6 * (1.0 - t) * t * \
           (ctrl_poly[2] - ctrl_poly[1]) + 3 * t ** 2 * (ctrl_poly[3] - ctrl_poly[2])


# evaluates cubic bezier second derivative at t, return point
def q_prime_prime(ctrl_poly, t):
    return 6 * (1.0 - t) * (ctrl_poly[2] - 2 * ctrl_poly[1] + ctrl_poly[0]) + 6 * t * (
        ctrl_poly[3] - 2 * ctrl_poly[2] + ctrl_poly[1])


# evaluates linear bezier at t, return point
def bezier_linear_interp(points, time):
    return (1 - time) * points[0] + time * points[1]


# evaluates a cubic bezier at t to find the r line, return start and end points
def find_cubic_r_points(points, time):
    q0 = bezier_linear_interp(np.array([points[0], points[1]]), time)
    q1 = bezier_linear_interp(np.array([points[1], points[2]]), time)
    q2 = bezier_linear_interp(np.array([points[2], points[3]]), time)

    r1 = bezier_linear_interp(np.array([q0, q1]), time)
    r2 = bezier_linear_interp(np.array([q1, q2]), time)

    return np.array([r1, r2])
