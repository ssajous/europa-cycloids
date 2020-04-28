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
    find_stress_level_of_total(stresses)
    data = stresses.loc[stresses['stress'] > 0] if positive_only else stresses
    # data = data.loc[data['stress'] > 35]
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
                'deltaStress', 'overallMaxStress', 'stressPctOfMax']]

def find_stress_level_of_total(field):
    field['overallMaxStress'] = field['stress'].max()
    field['stressPctOfMax'] = field['stress'] / field['overallMaxStress']

def test_stress_parameters(batch, params, paramDiff, previousError, interior):
    min_vals = np.array([0, 0.1, 0])
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

    field = tools.build_simon_stress_field(
        interior,
        test_data,
        phase=phase,
        eccentricity=0.01,
        obliquity=np.radians(obliquity),
        nsr=0,
        is_async=True,
        steps=360)
    error = find_heading_error(test_data, field)

    result = error['deltaHeading'] * (1 + (1 - error['stressPctOfMax']))
    errorArray = np.array(result)

    root_mean_squared_error = np.sqrt(np.sum(np.power(result, 2))) / result.shape[0]

     # calculate jacobian & gradient
    # diffs = np.insert(np.diff(result), 0, result.iloc[0], axis=0)
    diffs = errorArray - previousError
    jac = np.array([[error/param if param != 0 else 0 for param in paramDiff] for error in diffs])
    loss_vector = np.array([root_mean_squared_error/error for error in result])
    gradient = loss_vector @ jac

    return root_mean_squared_error, gradient, errorArray


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

    @staticmethod
    def adjustParameters(params, contraints):
        for index in range(len(params)):
            constraint = contraints[index]
            minValue = constraint.get('minValue')
            maxValue = constraint.get('maxValue')
            wrapValue = constraint.get('wrapValue')

            if wrapValue is not None and minValue is not None and maxValue is not None:
                if wrapValue:
                    if params[index] > maxValue:
                        params[index] = minValue + (params[index] - maxValue)
                    elif params[index] < minValue:
                        params[index] = maxValue - (minValue - params[index])
            elif maxValue is not None:
                if params[index] > maxValue:
                    params[index] = maxValue
            elif minValue is not None:
                if params[index] < minValue:
                    params[index] = minValue

        return params

    def minimize(
        self,
        objective_function,
        dataset,
        starting_params,
        interior,
        constraints=[{},{},{}],
        batch_size=32,
        threshold=0.01,
        max_iterations=10000,
        verbose=False):

        params = starting_params
        oldParams = np.zeros_like(params)
        previousErrorVector = np.zeros(batch_size)

        best_case = dict(loss=100, parameters=params)
        losses = []

        moment = [np.zeros_like(params)]
        raw_moment = [np.zeros_like(params)]

        loss = 100
        time = 1
        while loss > threshold and time < max_iterations:
            batch = dataset.sample(batch_size)

            deltaParams = params - oldParams
            loss, gradient, previousErrorVector = objective_function(batch, params, deltaParams, previousErrorVector, interior)

            losses.append(loss)
            if loss < best_case['loss']:
                best_case['loss'] = loss
                best_case['parameters'] = params

            moment.append(self.beta1 * moment[time - 1] + (1. - self.beta1) * gradient)
            raw_moment.append(self.beta2 * raw_moment[time - 1]  + (1. - self.beta2) * gradient**2)

            learning_rate = self.alpha * (np.sqrt(1. - self.beta2**time)/(1. - self.beta1**time))

            oldParams = params
            params = params - learning_rate * moment[time]/(np.sqrt(raw_moment[time]) + self.epsilon)

            params = self.adjustParameters(params, constraints)

            if verbose:
                window_size = 25
                avg_loss = np.average(losses) if len(losses) < window_size else np.average(losses[-1*window_size:])
                print(f'Iteration {time}/{max_iterations} -- Loss Output: {loss} -- Moving Avg Loss: {avg_loss}')

            time += 1

        return np.array(losses), best_case, dict(loss=loss, parameters=params)

class Nesterov:
    def __init__(self, alpha=1e-5, gamma=0.9):
        self.alpha = alpha
        self.gamma = gamma

    def minimize(
        self,
        objective_function,
        dataset,
        starting_params,
        interior,
        batch_size=32,
        threshold=0.01,
        max_iterations=10000,
        verbose=False):

        params = starting_params
        oldParams = np.zeros_like(params)
        previousErrorVector = np.zeros(batch_size)

        best_case = dict(loss=100, parameters=params)
        losses = []

        velocity = np.zeros_like(params)

        loss = 100
        time = 1
        while loss > threshold and time < max_iterations:
            batch = dataset.sample(batch_size)

            deltaParams = params - oldParams
            loss, gradient, previousErrorVector = objective_function(batch, params, deltaParams, previousErrorVector, interior)

            losses.append(loss)
            if loss < best_case['loss']:
                best_case['loss'] = loss
                best_case['parameters'] = params

            velocity = self.gamma * velocity - self.alpha * gradient

            oldParams = params
            params = params + velocity

            params[params >=  1] = 1
            params[params <= 0] = 1e-8

            if verbose:
                window_size = 25
                avg_loss = np.average(losses) if len(losses) < window_size else np.average(losses[-1*window_size:])
                print(f'Iteration {time}/{max_iterations} -- Loss Output: {loss} -- Moving Avg Loss: {avg_loss}')

            time += 1

        return np.array(losses), best_case, dict(loss=loss, parameters=params)

