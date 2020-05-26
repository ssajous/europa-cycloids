import pandas as pd
import numpy as np
import StressTools as tools
import curves.fitCurves as fit
import curves.bezier as bezier
import utils
from scipy import stats
import matplotlib.pyplot as plt


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
                'headingCategory': utils.round_heading(heading, tolerance),
                'isCusp': False
            })
            pointNumber += 1

    if len(rows) > 0:
        rows[0]['isCusp'] = True
        rows[-1]['isCusp'] = True

    return pd.DataFrame(rows)


def createCycloidBezier(arcs, pointsPerCurve=100, maxError=0.05):
    all_curves = None

    for index, arc in enumerate(arcs):
        startingPoint = all_curves.shape[0] + 1 if all_curves is not None else 1
        curve = fit_arc(arc,
                        startingPoint=startingPoint,
                        output_points=pointsPerCurve,
                        max_error=maxError)
        curve['arcNumber'] = index + 1
        all_curves = curve if all_curves is None else pd.concat([all_curves, curve], ignore_index=True)

    return all_curves


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
                'heading_y', 'stress', 'deltaHeading', 'isCusp', 'arcNumber',
                'deltaStress', 'overallMaxStress', 'stressPctOfMax']]

def find_stress_level_of_total(field):
    groups = field.groupby(['lon', 'lat'])['stress']
    field['overallMaxStress'] = groups.transform('max')
    field['stressPctOfMax'] = field['stress'] / field['overallMaxStress']

def normalize_parameters(params):
    min_vals = np.array([0, 0.1, 0])
    max_vals = np.array([360, 1, 360])

    if len(params) == 3:
        variables = (params - min_vals) / (max_vals - min_vals)
    else:
        variables = (params - min_vals[0:2:]) / (max_vals[0:2:] - min_vals[0:2:])

    return variables

def calc_non_monotonic_error_rate(series):
    if series.shape[0] == 0:
        return

    forwardError = 0
    reverseError = 0
    prevValue = None

    for row in series:
        if prevValue is None:
            pass
        elif row > prevValue:
            reverseError += 1
        elif row < prevValue:
            forwardError += 1
        else: # equal values
            pass

        prevValue = row

    error = min([forwardError, reverseError])
    errorRate = error / series.shape[0]

    return errorRate

def get_time_error_coefficient(data):
    timedf = data.copy().sort_values('pointNumber')

    mask = timedf.time < 180
    timedf.loc[mask, 'time'] = timedf.loc[mask, 'time'] + 360

    error_rate = calc_non_monotonic_error_rate(timedf.time)
    return 1 + error_rate

def get_stress_error_coefficient(data):
    cusps = data.copy()

    # forward
    starts = cusps.sort_values('pointNumber', ascending=True).drop_duplicates(['arcNumber'])
    forward_error_rate = calc_non_monotonic_error_rate(starts.stress)

    # reverse
    starts = cusps.sort_values('pointNumber', ascending=False).drop_duplicates(['arcNumber'])
    reverse_error_rate = calc_non_monotonic_error_rate(starts.stress)

    return 1 + min([forward_error_rate, reverse_error_rate])

def calc_monotonic_errors(data, dataset):
    error = get_stress_error_coefficient(dataset)
    data['stressError'] = error

    error = data.groupby('arcNumber').apply(get_time_error_coefficient)
    data['timeError'] = data['arcNumber'].map(error)

    return data


def test_stress_parameters(batch, dataset, params, paramDiff, previousError, interior):
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

    cusps = dataset.loc[dataset['isCusp']]
    cuspField = tools.build_simon_stress_field(
        interior,
        cusps,
        phase=phase,
        eccentricity=0.01,
        obliquity=np.radians(obliquity),
        nsr=0,
        is_async=True,
        steps=360)

    cuspError = find_heading_error(cusps.copy(), cuspField, positive_only=False)
    error = find_heading_error(test_data, field)
    error = calc_monotonic_errors(error, cuspError)

    result = error['deltaHeading'] * (1 + (1 - error['stressPctOfMax'])) * error['stressError'] * error['timeError']
    errorArray = np.array(result)

    root_mean_squared_error = np.sqrt(np.sum(np.power(result, 2))) / result.shape[0]

     # calculate jacobian & gradient
    # diffs = np.insert(np.diff(result), 0, result.iloc[0], axis=0)
    diffs = errorArray - previousError
    jac = np.array([[error/param if param != 0 else 1 for param in paramDiff] for error in diffs])
    loss_vector = np.array([root_mean_squared_error/error for error in result])
    gradient = loss_vector @ jac

    return root_mean_squared_error, gradient, errorArray

def getArcDistance(lonSeries, latSeries, radius, reverse=False):
    if reverse:
        lons = np.radians(lonSeries[::-1])
        lats = np.radians(latSeries[::-1])
    else:
        lons = np.radians(lonSeries)
        lats = np.radians(latSeries)
    lonsPrev = np.roll(lons, 1)
    lonsPrev[0] = lonsPrev[1]

    latsPrev = np.roll(lats, 1)
    latsPrev[0] = latsPrev[1]

    c1 = (np.sin((lats - latsPrev) / 2)) ** 2
    c2 = (np.sin((lons - lonsPrev) / 2)) ** 2

    distance = 2 * radius * np.arcsin(np.sqrt(c1 + np.cos(lats) * np.cos(latsPrev) * c2))

    if reverse:
        return distance[::-1]
    else:
        return distance

def calcKStress(stresses, lengths):
    return (1.12 * stresses * 1000 * \
        np.sqrt(np.pi * lengths * 1000)) / 1000

def calcChangeRate(x, y):
    dx = np.diff(x)
    dy = np.diff(y)
    changes = np.insert(np.abs(dy/dx), 0, 0)

    return changes

def addMetrics(stress, interior):
    radiusKm = interior.radius / 1000
    stress['segLengthKm'] = getArcDistance(stress['lon'], stress['lat'], radiusKm)
    stress['segLengthKmReverse'] = getArcDistance(stress['lon'], stress['lat'], radiusKm, reverse=True)

    stress['cumulativeCycloidLength'] = np.cumsum(stress['segLengthKm'])
    stress['cumulativeCycloidLengthReverse'] = (np.cumsum(stress[::-1]['segLengthKmReverse']))[::-1]

    stress['cumulativeArcLength'] = stress.groupby('arcNumber')['segLengthKm'].cumsum()
    stress['cumulativeArcLengthReverse'] = stress[::-1].groupby('arcNumber')['segLengthKmReverse'].cumsum()[::-1]

    stress['KStress'] = calcKStress(stress['stress'], stress['cumulativeArcLength'])
    stress['KStressReverse'] = calcKStress(stress['stress'], stress['cumulativeArcLengthReverse'])
    stress['headingChangeRate'] = calcChangeRate(stress['lon'], stress['heading_x'])
    stress['headingAcceleration'] = calcChangeRate(stress['lon'], stress['headingChangeRate'])

    return stress

def match_stresses(batch, params, interior, saveStressField=False, path='./output/stressfield.csv.gz'):
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
    error = addMetrics(error, interior)

    if saveStressField:
        field.to_csv(path, index=False, compression='gzip')

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
            unstick = constraint.get('unstick')

            if wrapValue and minValue is not None and maxValue is not None:
                if params[index] > maxValue:
                    params[index] = minValue + (params[index] - maxValue)
                elif params[index] < minValue:
                    params[index] = maxValue - (minValue - params[index])
            elif unstick and minValue is not None and maxValue is not None:
                if params[index] > maxValue or params[index] < minValue:
                    params[index] = np.random.rand() * (maxValue - minValue) + minValue
            elif maxValue is not None and params[index] > maxValue:
                params[index] = maxValue
            elif minValue is not None and params[index] < minValue:
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
        verbose=False,
        infoFrequency=150):

        params = starting_params
        oldParams = np.zeros_like(params)
        previousErrorVector = np.random.randn(batch_size)

        best_case = dict(loss=100, parameters=params)
        losses = []
        history = []

        moment = [np.zeros_like(params)]
        raw_moment = [np.zeros_like(params)]

        loss = 100
        time = 1
        while loss > threshold and time < max_iterations + 1:
            batch = dataset.sample(batch_size)

            deltaParams = params - oldParams
            loss, gradient, previousErrorVector = objective_function(batch, dataset, params, deltaParams, previousErrorVector, interior)
            history.append(np.array([loss, *params]))

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
                print(f'\tParameters used: {oldParams}')
            elif time % infoFrequency == 0:
                window_size = 25
                avg_loss = np.average(losses) if len(losses) < window_size else np.average(losses[-1*window_size:])
                print(f'Iteration {time}/{max_iterations} -- Loss Output: {loss} -- Moving Avg Loss: {avg_loss}')
                print(f'\tParameters used: {oldParams}')

            time += 1

        return np.array(losses), best_case, dict(loss=loss, parameters=params), np.array(history)


def optimizationHeatMap(optimize_output, cycloid_name):
    phase = optimize_output[3].T[1]
    obliquity = optimize_output[3].T[2]
    loss = optimize_output[3].T[0]

    phase_min = phase.min()
    obl_min = obliquity.min()
    loss_min = loss.min()
    phase_max = phase.max()
    obl_max = obliquity.max()
    loss_max = loss.max()

    # First Obliquity vs loss
    dataset = np.vstack([obliquity, loss])
    X, Y = np.mgrid[obl_min:obl_max:100j, loss_min:loss_max:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    kernel = stats.gaussian_kde(dataset)
    Z = np.reshape(kernel(positions).T, X.shape)

    plt.figure(figsize=(10, 5))
    plt.imshow(Z.transpose(),
            origin="lower",
            extent=[obl_min, obl_max, loss_min, loss_max],
            aspect=(obl_max-obl_min)/(loss_max-loss_min))
    plt.ylim(ymin = loss_min, ymax = min([5, loss_max]))
    plt.title(f'{cycloid_name} - Obliquity Concentration')
    plt.xlabel('Obliquity')
    plt.ylabel('Loss')

    # Then do Phase vs loss
    dataset = np.vstack([phase, loss])
    kernel = stats.gaussian_kde(dataset)
    X, Y = np.mgrid[phase_min:phase_max:100j, loss_min:loss_max:100j]
    Z = np.reshape(kernel(positions).T, X.shape)

    plt.figure(figsize=(10, 5))
    plt.imshow(Z.transpose(),
            origin="lower",
            extent=[phase_min, phase_max, loss_min, loss_max],
            aspect=(phase_max-phase_min)/(loss_max-loss_min))
    plt.ylim(ymin = loss_min, ymax = min([5, loss_max]))
    plt.title(f'{cycloid_name} - Phase Concentration')
    plt.xlabel('Phase')
    plt.ylabel('Loss')

def findBestParametersFromFrame(optimizeHistory):
    df = optimizeHistory.loc[optimizeHistory['loss'] < 1]

    # Find loss threshold
    lossHist = np.histogram(df['loss'], bins=100, density=True)
    index = np.argmax(lossHist[0])
    threshold = lossHist[1][index + 1]

    # find param values
    bestFits = df.loc[df['loss'] <= threshold]
    paramHist = np.histogram2d(bestFits['phase'], bestFits['obliquity'], bins=50, density=True)
    index = np.unravel_index(np.argmax(paramHist[0]), paramHist[0].shape)
    phase = np.average(paramHist[1][index[0]:index[0] + 2])
    obliquity = np.average(paramHist[2][index[1]:index[1] + 2])

    return phase, obliquity

def findBestParameters(optimize_output):
    df = pd.DataFrame(optimize_output[3], columns=['loss', 'phase', 'obliquity'])

    return findBestParametersFromFrame(df)

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


