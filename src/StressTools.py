import pandas as pd
import numpy as np
import sympy as sym
import math
import multiprocessing
from sympy.printing.theanocode import theano_function
import StressEquations as simon
import utils
from joblib import Parallel, delayed, Memory

TIME_STEPS = 360
MIN_LAT = -75
MAX_LAT = 90
MIN_LON = 0
MAX_LON = 360
LAT_STEP_SIZE = 15
LON_STEP_SIZE = 10
RAD_MULTIPLIER = np.pi / 180
DEG_MULTIPLIER = 180 / np.pi
CPUS = multiprocessing.cpu_count()
CACHE_DIR = './cache'

mem = Memory(CACHE_DIR, verbose=0)

r, θ, φ, t = sym.symbols('r θ φ t', real=True)

get_principal1 = None
get_principal2 = None
get_principal_orientation = None
get_principal_orientation2 = None


# noinspection PyCallingNonCallable
def get_stress_for_latitude(step, time, lat):
    """
    Do NOT use this function directly! Use build_stress_field instead.

    Generates a data frame with stresses for each longitude within the specified
    longitude and time step.

    Parameters
    ----------
    step: int
        The current time step to use for stress generation.  Value should be
        between 1 and 360 inclusive.
    time: float
        The number of seconds to reach the specified step in the orbit
    lat: float
        The latitude in degrees to use for stress generation

    Returns
    -------
    DataFrame
        Stress data which includes time step, latitude (degrees), longitude (degrees),
        principal1 stress, principal2 stress, principal1 orientation (degrees),
        principal2 orientation (degrees), max stress value, max stress orientation.
    """
    results = []

    lat_radians = (90 - lat) * RAD_MULTIPLIER
    time_value = time

    for lon in range(MIN_LON, MAX_LON + 1, 10):
        if lat == 90 or lon == 0:
            continue

        lon_radians = lon * RAD_MULTIPLIER

        principal1 = get_principal1(time_value, lat_radians, lon_radians)
        principal2 = get_principal2(time_value, lat_radians, lon_radians)
        principal_phi = get_principal_orientation(time_value, lat_radians, lon_radians)
        principal_phi2 = get_principal_orientation2(time_value, lat_radians, lon_radians)

        max_stress = max(principal1, principal2)
        max_stress_orientation = principal_phi if max_stress == principal1 else principal_phi2

        results.append({
            'time_step': step,
            'latitude': lat,
            'longitude': lon,
            'principal1': principal1,
            'principal2': principal2,
            'principal_orientation': principal_phi * DEG_MULTIPLIER,
            'principal_orientation2': principal_phi2 * DEG_MULTIPLIER,
            'max_stress': max_stress,
            'max_stress_orientation': max_stress_orientation * DEG_MULTIPLIER
        })

    return results


# noinspection PyUnboundLocalVariable
def build_stress_field(satellite, orbit_time_seconds, rotations=1, is_async=True):
    """
    Creates a data frame with the results of stress calculations for a range of
    latitudes and longitudes across 360 time steps

    Parameters
    ----------
    satellite : MEWtools.satellite
        The body on which stress values will be calculated
    orbit_time_seconds: int,
        The total number of seconds it takes for the satellite to orbit
        its primary
    rotations: float, optional
        The number of orbital rotations to loop over for time step
        calculations. Default value is 1
    is_async: boolean, optional
        When True (default) the calculation will be spread across multiple
        processes for performance.  When False all calculation is done in the
        current process, this most likely slower, but better for debugging.

    Returns
    -------
    DataFrame
        Stress data which includes time step, latitude (degrees), longitude (degrees),
        principal1 stress, principal2 stress, principal1 orientation (degrees),
        principal2 orientation (degrees), max stress value, max stress orientation.
    """

    # Create "lamdified" versions of each of the stress equations.  This turns
    # the symbolic python formulas into standard python functions.  Which improves their
    # performance by about 2 orders of magnitude.  The use of global variables for the
    # functions is necessary to allow the multiprocess functionality to work.  Locally
    # defined functions cannot be pickled.
    global get_principal1
    global get_principal2
    global get_principal_orientation
    global get_principal_orientation2
    get_principal1 = sym.lambdify([t, φ, θ], satellite.PC1, modules=["math", {"cot": math.atan}])
    get_principal2 = sym.lambdify([t, φ, θ], satellite.PC2, modules=["math", {"cot": math.atan}])
    get_principal_orientation = sym.lambdify([t, φ, θ], satellite.PCΨ, modules=["math", {"cot": math.atan}])
    get_principal_orientation2 = sym.lambdify([t, φ, θ], satellite.PCΨ2, modules=["math", {"cot": math.atan}])

    data = []

    def callback(items):
        data.extend(items)

    def error_callback(err):
        print(err)

    if is_async:
        pool = multiprocessing.Pool()

    # mean_motion = (2 * np.pi) / orbit_time_seconds
    total_steps = int(TIME_STEPS * rotations)
    for step in range(1, total_steps):
        # time = np.radians(step/TIME_STEPS)
        time = (step / TIME_STEPS) * orbit_time_seconds
        for lat in range(MIN_LAT, MAX_LAT + 1, LAT_STEP_SIZE):
            if is_async:
                # pool.apply_async will schedule the processing within separate python processes
                # which allows the work to be distributed to multiple CPU cores
                pool.apply_async(get_stress_for_latitude,
                                 args=(step, time, lat,),
                                 callback=callback,
                                 error_callback=error_callback)
            else:
                data.extend(get_stress_for_latitude(step, time, lat))

    if is_async:
        pool.close()
        pool.join()

    df = pd.DataFrame(data)
    return df.sort_values(['latitude', 'longitude', 'time_step'])


def build_mew_stress_field(satellite, orbit_time_seconds, point_frame, rotations=1):
    """
    Generates a dataset for tidal stress values at each lat/lon location input over the specified number of orbits.

    :type rotations: float
    :type point_frame: DataFrame
    :type orbit_time_seconds: int/float
    :type satellite: MEWTools.Satellite
    :rtype: DataFrame

    :param satellite: Definition of the interior structure of the satellite
    :param orbit_time_seconds: Total number of seconds it takes the satellite to orbit it's primary
    :param point_frame: DataFrame of lat/lon locations in radians. Lat is colat and lon is east lon
    :param rotations: The number of orbits to calculate stress
    :return: DataFrame with the stresses at each time step for all latitude/longitude points
    """
    stress_calc = theano_function([t, φ, θ], (satellite.PC1, satellite.PC2, satellite.PCΨ, satellite.PCΨ2),
                                  dims={t: 1, φ: 1, θ: 1},
                                  dtypes={t: 'float64', φ: 'float64', θ: 'float64'})

    # mean_motion = (2 * np.pi) / orbit_time_seconds
    max_time = orbit_time_seconds * rotations
    times = np.linspace(1, max_time, TIME_STEPS * rotations)
    parameters = np.array([[time, point.lon, point.lat] for point in point_frame.itertuples() for time in times])
    stress_output = stress_calc(parameters[:, 0:1].flatten(),
                                parameters[:, 1:2].flatten(),
                                parameters[:, 2:3].flatten())
    output_stacked = np.column_stack(stress_output)

    results = pd.DataFrame(np.hstack([parameters, output_stacked]),
                           columns=['time', 'lat', 'lon', 'principal1', 'principal2', 'orientation1', 'orientation2'])

    return results.sort_values(['lat', 'lon', 'time'])


def get_stresses_for_point(
                           interior,
                           lon,
                           lat,
                           phase,
                           tolerance,
                           steps,
                           eccentricity,
                           obliquity,
                           nsr):
    results = []
    previous = 0
    for step in range(steps):
        current = simon.getStress(
            interior_value=interior,
            e_in=eccentricity,
            colat=np.radians(90 - lat),
            lon=np.radians(360 - lon),
            steps=steps,
            this_step=step,
            oblq=obliquity,
            phase=np.radians(phase),
            NSRdelta=nsr)
        heading_degrees = np.degrees(current[1])
        results.append({
            'lon': lon,
            'lat': lat,
            'stress': current[0],
            'heading': heading_degrees,
            'headingCategory': utils.round_heading(heading_degrees, tolerance),
            'deltaStress': current[0] - previous,
            'time': step
        })
        previous = current[0]

    results[0]['deltaStress'] = results[0]['stress'] - results[-1]['stress']
    return results


def build_simon_stress_field(
                            interior,
                            point_frame,
                            phase_degrees,
                            eccentricity,
                            obliquity_radians,
                            nsr_radians,
                            is_async=True,
                            tolerance=1,
                            steps=360):
    """
    Generates a dataset for tidal stress values at each lat/lon location input over the duration of a satellite's full
    orbit.

    :rtype: DataFrame
    :type steps: int
    :type tolerance: float
    :type is_async: bool
    :type nsr_radians: float
    :type obliquity_radians: float
    :type eccentricity: float
    :type phase_degrees: float
    :type point_frame: DataFrame
    :type interior: Interior

    :param interior: Interior structure for the satellite
    :param point_frame: DataFrame containing the points where stress will be evaluated. Must contain columns lat and lon
    which represent degrees latitude and degrees west longitude respectively.

    :param phase_degrees: Spin pole direction in degrees
    :param eccentricity: satellite eccentricity value
    :param obliquity_radians: satellite obliquity value in radians
    :param nsr_radians: Non-synchronous rotation amount in radians
    :param is_async: When True multiple processors will be used to calculate the stress field. Using False can sometimes
    aid in debugging.
    :param tolerance: Heading tolerance in degrees used to group headings into categories.
    :param steps: The number of time steps across an orbit used to calculate stress at different points. A full orbit is
    always used, steps only impacts the granularity.
    :return: DataFrame with the stresses at each time step for all latitude/longitude points
    """
    stresses = []

    if is_async:
        point_stresses = Parallel(n_jobs=CPUS)(delayed(get_stresses_for_point)
                                               (interior, point.lon, point.lat, phase_degrees, tolerance, steps,
                                                eccentricity, obliquity_radians,
                                                nsr_radians) for point in point_frame.itertuples())
    else:
        point_stresses = [
            get_stresses_for_point(interior, point.lon, point.lat, phase_degrees, tolerance, steps, eccentricity,
                                   obliquity_radians, nsr_radians) for point in point_frame.itertuples()]

    for stress in point_stresses:
        stresses.extend(stress)

    df = pd.DataFrame(stresses).sort_values(['lat', 'lon', 'time'])

    return df


get_simon_stress_field = mem.cache(build_simon_stress_field)
get_stress_field = mem.cache(build_stress_field)
get_mew_stress_field = mem.cache(build_mew_stress_field)
