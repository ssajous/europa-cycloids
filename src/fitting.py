import pandas as pd
import numpy as np
import StressTools as tools
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


def fit_arc(arc, max_error=0.05, reverse=False, startingPoint=1, tolerance=1, output_points=100):
    points = np.array(arc)
    controls = fit.fitCurve(points, max_error)

    if reverse:
        controls = controls[::-1]

    rows = []
    pointNumber = startingPoint
    for control in controls:
        if reverse:
            control.reverse()
        for time in np.arange(0, 1, 1 / output_points):
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


def match_orientations(curve, stresses, positive_only=True):
    merged = curve.merge(
        stresses,
        how='left',
        on=['lon', 'lat', 'headingCategory']
    )

    if positive_only:
        merged = merged.loc[merged['stress'] > 0]

    merged['maxStress'] = merged.groupby('pointNumber')['stress'].transform('max')
    merged_unique = merged.loc[merged['stress'] == merged['maxStress']]

    return merged_unique

def calculate_loss(frame, startingPoint=1, pointCount=100):
    if len(frame) == 0:
        return pointCount

    diffs = np.diff(frame.sort_values('pointNumber')['pointNumber'])

    edges = frame.iloc[[0, -1]]
    startDiff = edges.iloc[0].pointNumber - 1
    endDiff = 100 - edges.iloc[1].pointNumber

    diffs = np.append(diffs, [startDiff, endDiff])
    return max(diffs)

def test_arc(
    arc,
    phase_increment,
    interior,
    eccentricity,
    obliquity,
    nsr,
    steps=360, positive_only=True):
    results = []
    startingPoint = arc['pointNumber'].min()
    points = len(arc)

    for phase in range(0, 361, phase_increment):
        field = tools.get_simon_stress_field(
            interior,
            arc,
            phase,
            eccentricity,
            obliquity,
            nsr,
            steps=steps)
        matches = match_orientations(arc, field, positive_only)
        loss = calculate_loss(matches, startingPoint, points)

        results.append({
            'phase': phase,
            'error': loss
        })
        print(f'Calculated for phase {phase}')

    return pd.DataFrame(results)

def find_heading_error(curve, stresses, positive_only=True):
    data = stresses.loc[stresses['stress'] > 0] if positive_only else stresses
    merged = curve.merge(
        data,
        how='left',
        on=['lat', 'lon']
    )
    merged['deltaHeading'] = np.abs(merged['heading_x'] - merged['heading_y'])
    merged['minDeltaHeading'] = merged.groupby('pointNumber')['deltaHeading'].transform('min')
    merged['maxStress'] = merged.groupby('pointNumber')['stress'].transform('max')

    merged_unique = merged.loc[(merged['deltaHeading'] == merged['minDeltaHeading'])].copy()
    merged_unique['maxStress'] = merged_unique.groupby('pointNumber')['stress'].transform('max')

    merged_unique = merged_unique.loc[merged_unique['stress'] == merged_unique['maxStress']]

    return merged_unique[['pointNumber', 'lon', 'lat', 'time', 'heading_x',
                'heading_y', 'stress', 'deltaHeading',
                'deltaStress']]


def test_stress_parameters(batch, params, interior):
    min_vals = np.array([0, 0.01, 0])
    max_vals = np.array([360, 1, 360])

    if len(params) == 3:
        variables = params * (max_vals - min_vals) + min_vals # denormalize
        phase, obliquity, longitude = variables
    else:
        variables = params * (max_vals[0:2:] - min_vals[0:2:]) + min_vals[0:2:]
        phase, obliquity = variables
        longitude = 0

    test_data = batch.copy()
    test_data['lon'] = test_data['lon'] + longitude

    field = tools.get_simon_stress_field(
        interior,
        test_data,
        phase=phase,
        eccentricity=0.01,
        obliquity=np.radians(obliquity),
        nsr=0,
        is_async=True,
        steps=360)
    error = find_heading_error(test_data, field)

    result = error['deltaHeading']

    root_mean_squared_error = np.sqrt(np.sum(np.power(result, 2))) / result.shape[0]

     # calculate jacobian & gradient
    diffs = np.insert(np.diff(result), 0, result.iloc[0], axis=0)
    jac = np.array([[(-1*error)/param for param in params] for error in diffs])
    loss_vector = np.array([root_mean_squared_error/error for error in result])
    gradient = loss_vector @ jac

    return root_mean_squared_error, gradient


def match_stresses(batch, params, interior):
    if len(params) == 3:
        phase, obliquity, longitude = params
    else:
        phase, obliquity = params
        longitude = 0

    test_data = batch.copy()
    test_data['lon'] = test_data['lon'] + longitude

    field = tools.get_simon_stress_field(
        interior,
        test_data,
        phase=phase,
        eccentricity=0.01,
        obliquity=np.radians(obliquity),
        nsr=0,
        is_async=True,
        steps=360)
    error = find_heading_error(test_data, field)

    return error

class Adam:

    def __init__(self, alpha=1e-3, beta1=0.9, beta2=0.999, epsilon=1e-8):
        self.alpha = alpha
        self.beta1 = beta1
        self.beta2 = beta2
        self.epsilon = epsilon

    def minimize(
        self,
        objective_function,
        dataset,
        starting_params,
        interior,
        batch_size=32,
        threshold=0.01,
        max_iterations=1000,
        verbose=False):

        params = starting_params

        best_case = dict(loss=100, parameters=params)
        worst_case = dict(loss=0, parameters=params)
        losses = []

        moment = [np.zeros_like(params)]
        raw_moment = [np.zeros_like(params)]

        loss = 100
        time = 1
        while loss > threshold and time < max_iterations:
            batch = dataset.sample(batch_size)

            loss, gradient = objective_function(batch, params, interior)

            losses.append(loss)
            if loss < best_case['loss']:
                best_case['loss'] = loss
                best_case['parameters'] = params

            moment.append(self.beta1 * moment[time - 1] + (1. - self.beta1) * gradient)
            raw_moment.append(self.beta2 * raw_moment[time - 1]  + (1. - self.beta2) * gradient**2)

            learning_rate = self.alpha * (np.sqrt(1. - self.beta2**time)/(1. - self.beta1**time))
            params = params - learning_rate * moment[time]/(np.sqrt(raw_moment[time]) + self.epsilon)

            params[params >=  1] = 1
            params[params <= 0] = self.epsilon

            time += 1

            if verbose:
                print(f'Loss Output: {loss}')


        return np.array(losses), best_case, dict(loss=loss, parameters=params)
