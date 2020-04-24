import pandas as pd
import numpy as np
import curves.fitCurves as fit
import curves.bezier as bezier
import utils

def findHeading(points, reverse = False):
    if reverse:
        origin = 1
        destination = 0
    else:
        origin = 0
        destination = 1

    rise = points[destination][1] - points[origin][1] # lats
    run = points[destination][0] - points[origin][0]  # lons

    degrees = np.degrees(np.arctan2(run, rise))
    if degrees < 0:
        degrees += 360
    return degrees if degrees > 180 else degrees + 180


def fit_arc(arc, max_error=0.05, reverse=False, startingPoint=1, tolerance=1):
    points = np.array(arc)
    controls = fit.fitCurve(points, max_error)

    if reverse:
        controls = controls[::-1]

    rows = []
    pointNumber = startingPoint
    for control in controls:
        if reverse:
            control.reverse()
        for time in np.arange(0, 1, 0.01):
            point = bezier.q(control, time)

            heading = findHeading(bezier.findCubicRPoints(control, time))
            rows.append({
                'pointNumber': pointNumber,
                'lon': point[0],
                'lat': point[1],
                'heading': heading,
                'headingCategory': utils.round_heading(heading, tolerance)
            })
            pointNumber += 1

    return pd.DataFrame(rows)
