import pandas as pd
import numpy as np
import StressTools as tools
import curves.fitCurves as fit
import curves.bezier as bezier
import utils
from scipy import stats
import matplotlib.pyplot as plt


def find_heading(points, reverse=False):
    if reverse:
        origin = 1
        destination = 0
    else:
        origin = 0
        destination = 1

    rise = points[destination][1] - points[origin][1]  # lats
    run = points[destination][0] - points[origin][0]  # longitudes

    degrees = np.degrees(np.arctan2(run, rise))
    if degrees < 0:
        degrees += 360
    return degrees if degrees > 180 else degrees + 180


def fit_arc(arc, max_error=0.05, reverse=False, starting_point=1, tolerance=1, output_points=100):
    points = np.array(arc)
    controls = fit.fitCurve(points, max_error)

    if reverse:
        controls = controls[::-1]

    rows = []
    point_number = starting_point
    for control in controls:
        if reverse:
            control.reverse()
        for time in np.arange(0, 1, 1 / output_points):
            point = bezier.q(control, time)

            heading = find_heading(bezier.find_cubic_r_points(control, time))
            rows.append({
                'pointNumber': point_number,
                'lon': point[0],
                'lat': point[1],
                'heading': heading,
                'headingCategory': utils.round_heading(heading, tolerance),
                'isCusp': False
            })
            point_number += 1

    if len(rows) > 0:
        rows[0]['isCusp'] = True
        rows[-1]['isCusp'] = True

    return pd.DataFrame(rows)


def create_cycloid_bezier(arcs, points_per_curve=100, max_error=0.05):
    all_curves = None

    for index, arc in enumerate(arcs):
        starting_point = all_curves.shape[0] + 1 if all_curves is not None else 1
        curve = fit_arc(arc,
                        starting_point=starting_point,
                        output_points=points_per_curve,
                        max_error=max_error)
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


def calculate_loss(frame, point_count=100):
    if len(frame) == 0:
        return point_count

    diffs = np.diff(frame.sort_values('pointNumber')['pointNumber'])

    edges = frame.iloc[[0, -1]]
    start_diff = edges.iloc[0].pointNumber - 1
    end_diff = 100 - edges.iloc[1].pointNumber

    diffs = np.append(diffs, [start_diff, end_diff])
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
        loss = calculate_loss(matches, points)

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

    forward_error = 0
    reverse_error = 0
    prev_value = None

    for row in series:
        if prev_value is None:
            pass
        elif row > prev_value:
            reverse_error += 1
        elif row < prev_value:
            forward_error += 1
        else:  # equal values
            pass

        prev_value = row

    error = min([forward_error, reverse_error])
    error_rate = error / series.shape[0]

    return error_rate


def get_time_error_coefficient(data):
    time_df = data.copy().sort_values('pointNumber')

    mask = time_df.time < 180
    time_df.loc[mask, 'time'] = time_df.loc[mask, 'time'] + 360

    error_rate = calc_non_monotonic_error_rate(time_df.time)
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


# noinspection DuplicatedCode
def test_stress_parameters(batch, dataset, params, param_diff, previous_error, interior):
    min_vals = np.array([0, 0.1, 0])
    max_vals = np.array([360, 1, 360])

    if len(params) == 3:
        variables = params * (max_vals - min_vals) + min_vals  # denormalize
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
    cusp_field = tools.build_simon_stress_field(
        interior,
        cusps,
        phase=phase,
        eccentricity=0.01,
        obliquity=np.radians(obliquity),
        nsr=0,
        is_async=True,
        steps=360)

    cusp_error = find_heading_error(cusps.copy(), cusp_field, positive_only=False)
    error = find_heading_error(test_data, field)
    error = calc_monotonic_errors(error, cusp_error)

    result = error['deltaHeading'] * (1 + (1 - error['stressPctOfMax'])) * error['stressError'] * error['timeError']
    error_array = np.array(result)

    root_mean_squared_error = np.sqrt(np.sum(np.power(result, 2))) / result.shape[0]

    # calculate jacobian & gradient
    # diffs = np.insert(np.diff(result), 0, result.iloc[0], axis=0)
    diffs = error_array - previous_error
    jac = np.array([[error / param if param != 0 else 1 for param in param_diff] for error in diffs])
    loss_vector = np.array([root_mean_squared_error / error for error in result])
    gradient = loss_vector @ jac

    return root_mean_squared_error, gradient, error_array


def get_arc_distance(lon_series, lat_series, radius, reverse=False):
    if reverse:
        lons = np.radians(lon_series[::-1])
        lats = np.radians(lat_series[::-1])
    else:
        lons = np.radians(lon_series)
        lats = np.radians(lat_series)
    lons_prev = np.roll(lons, 1)
    lons_prev[0] = lons_prev[1]

    lats_prev = np.roll(lats, 1)
    lats_prev[0] = lats_prev[1]

    c1 = (np.sin((lats - lats_prev) / 2)) ** 2
    c2 = (np.sin((lons - lons_prev) / 2)) ** 2

    distance = 2 * radius * np.arcsin(np.sqrt(c1 + np.cos(lats) * np.cos(lats_prev) * c2))

    if reverse:
        return distance[::-1]
    else:
        return distance


def calc_k_stress(stresses, lengths):
    return (1.12 * stresses * 1000 *
            np.sqrt(np.pi * lengths * 1000)) / 1000


def calc_change_rate(x, y):
    dx = np.diff(x)
    dy = np.diff(y)
    changes = np.insert(np.abs(dy / dx), 0, 0)

    return changes


def add_metrics(stress, interior):
    radius_km = interior.radius / 1000
    stress['segLengthKm'] = get_arc_distance(stress['lon'], stress['lat'], radius_km)
    stress['segLengthKmReverse'] = get_arc_distance(stress['lon'], stress['lat'], radius_km, reverse=True)

    stress['cumulativeCycloidLength'] = np.cumsum(stress['segLengthKm'])
    stress['cumulativeCycloidLengthReverse'] = (np.cumsum(stress[::-1]['segLengthKmReverse']))[::-1]

    stress['cumulativeArcLength'] = stress.groupby('arcNumber')['segLengthKm'].cumsum()
    stress['cumulativeArcLengthReverse'] = stress[::-1].groupby('arcNumber')['segLengthKmReverse'].cumsum()[::-1]

    stress['KStress'] = calc_k_stress(stress['stress'], stress['cumulativeArcLength'])
    stress['KStressReverse'] = calc_k_stress(stress['stress'], stress['cumulativeArcLengthReverse'])
    stress['headingChangeRate'] = calc_change_rate(stress['lon'], stress['heading_x'])
    stress['headingAcceleration'] = calc_change_rate(stress['lon'], stress['headingChangeRate'])

    return stress


# noinspection DuplicatedCode
def match_stresses(batch, params, interior, save_stress_field=False, path='./output/stressfield.csv.gz'):
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

    cusps = test_data.loc[test_data['isCusp']]
    cusp_field = tools.build_simon_stress_field(
        interior,
        cusps,
        phase=phase,
        eccentricity=0.01,
        obliquity=np.radians(obliquity),
        nsr=0,
        is_async=True,
        steps=360)

    cusp_error = find_heading_error(cusps.copy(), cusp_field, positive_only=False)
    error = find_heading_error(test_data, field)
    error = calc_monotonic_errors(error, cusp_error)
    error = add_metrics(error, interior)

    result = error['deltaHeading'] * (1 + (1 - error['stressPctOfMax'])) * error['stressError'] * error['timeError']

    root_mean_squared_error = np.sqrt(np.sum(np.power(result, 2))) / result.shape[0]

    if save_stress_field:
        field.to_csv(path, index=False, compression='gzip')

    return error, root_mean_squared_error


# noinspection PyUnresolvedReferences
class Adam:

    def __init__(self, alpha=1e-3, beta1=0.9, beta2=0.999, epsilon=1e-8):
        self.alpha = alpha
        self.beta1 = beta1
        self.beta2 = beta2
        self.epsilon = epsilon

    @staticmethod
    def adjust_parameters(params, contraints):
        for index in range(len(params)):
            constraint = contraints[index]
            min_value = constraint.get('minValue')
            max_value = constraint.get('maxValue')
            wrap_value = constraint.get('wrapValue')
            unstick = constraint.get('unstick')

            if wrap_value and min_value is not None and max_value is not None:
                if params[index] > max_value:
                    params[index] = min_value + (params[index] - max_value)
                elif params[index] < min_value:
                    params[index] = max_value - (min_value - params[index])
            elif unstick and min_value is not None and max_value is not None:
                if params[index] > max_value or params[index] < min_value:
                    params[index] = np.random.rand() * (max_value - min_value) + min_value
            elif max_value is not None and params[index] > max_value:
                params[index] = max_value
            elif min_value is not None and params[index] < min_value:
                params[index] = min_value

        return params

    def minimize(
                self,
                objective_function,
                dataset,
                starting_params,
                interior,
                constraints=[{}, {}, {}],
                batch_size=32,
                threshold=0.01,
                max_iterations=10000,
                verbose=False,
                info_frequency=150):

        params = starting_params
        old_params = np.zeros_like(params)
        previous_error_vector = np.random.randn(batch_size)

        best_case = dict(loss=100, parameters=params)
        losses = []
        history = []

        moment = [np.zeros_like(params)]
        raw_moment = [np.zeros_like(params)]

        loss = 100
        time = 1
        while loss > threshold and time < max_iterations + 1:
            batch = dataset.sample(batch_size)

            delta_params = params - old_params
            loss, gradient, previous_error_vector = objective_function(batch, dataset, params, delta_params,
                                                                       previous_error_vector, interior)
            history.append(np.array([loss, *params]))

            losses.append(loss)
            if loss < best_case['loss']:
                best_case['loss'] = loss
                best_case['parameters'] = params

            moment.append(self.beta1 * moment[time - 1] + (1. - self.beta1) * gradient)
            raw_moment.append(self.beta2 * raw_moment[time - 1] + (1. - self.beta2) * gradient ** 2)

            learning_rate = self.alpha * (np.sqrt(1. - self.beta2 ** time) / (1. - self.beta1 ** time))

            old_params = params
            params = params - learning_rate * moment[time] / (np.sqrt(raw_moment[time]) + self.epsilon)

            params = self.adjust_parameters(params, constraints)

            if verbose:
                window_size = 25
                avg_loss = np.average(losses) if len(losses) < window_size else np.average(losses[-1 * window_size:])
                print(f'Iteration {time}/{max_iterations} -- Loss Output: {loss} -- Moving Avg Loss: {avg_loss}')
                print(f'\tParameters used: {old_params}')
            elif time % info_frequency == 0:
                window_size = 25
                avg_loss = np.average(losses) if len(losses) < window_size else np.average(losses[-1 * window_size:])
                print(f'Iteration {time}/{max_iterations} -- Loss Output: {loss} -- Moving Avg Loss: {avg_loss}')
                print(f'\tParameters used: {old_params}')

            time += 1

        return np.array(losses), best_case, dict(loss=loss, parameters=params), np.array(history)


# noinspection DuplicatedCode
def optimization_heat_map(optimize_output, cycloid_name):
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
    x, y = np.mgrid[obl_min:obl_max:100j, loss_min:loss_max:100j]
    positions = np.vstack([x.ravel(), y.ravel()])
    kernel = stats.gaussian_kde(dataset)
    z = np.reshape(kernel(positions).T, x.shape)

    plt.figure(figsize=(10, 5))
    plt.imshow(z.transpose(),
               origin="lower",
               extent=[obl_min, obl_max, loss_min, loss_max],
               aspect=(obl_max - obl_min) / (loss_max - loss_min))
    plt.ylim(ymin=loss_min, ymax=min([5, loss_max]))
    plt.title(f'{cycloid_name} - Obliquity Concentration')
    plt.xlabel('Obliquity')
    plt.ylabel('Loss')

    # Then do Phase vs loss
    dataset = np.vstack([phase, loss])
    kernel = stats.gaussian_kde(dataset)
    x, y = np.mgrid[phase_min:phase_max:100j, loss_min:loss_max:100j]
    z = np.reshape(kernel(positions).T, x.shape)

    plt.figure(figsize=(10, 5))
    plt.imshow(z.transpose(),
               origin="lower",
               extent=[phase_min, phase_max, loss_min, loss_max],
               aspect=(phase_max - phase_min) / (loss_max - loss_min))
    plt.ylim(ymin=loss_min, ymax=min([5, loss_max]))
    plt.title(f'{cycloid_name} - Phase Concentration')
    plt.xlabel('Phase')
    plt.ylabel('Loss')


def find_best_parameters_from_frame(optimize_history):
    df = optimize_history.loc[optimize_history['loss'] < 1]

    # Find loss threshold
    loss_hist = np.histogram(df['loss'], bins=100, density=True)
    index = np.argmax(loss_hist[0])
    threshold = loss_hist[1][index + 1]

    # find param values
    best_fits = df.loc[df['loss'] <= threshold]
    param_hist = np.histogram2d(best_fits['phase'], best_fits['obliquity'], bins=50, density=True)
    index = np.unravel_index(np.argmax(param_hist[0]), param_hist[0].shape)
    phase = np.average(param_hist[1][index[0]:index[0] + 2])
    obliquity = np.average(param_hist[2][index[1]:index[1] + 2])

    return phase, obliquity


def find_best_parameters(optimize_output):
    df = pd.DataFrame(optimize_output[3], columns=['loss', 'phase', 'obliquity'])

    return find_best_parameters_from_frame(df)
