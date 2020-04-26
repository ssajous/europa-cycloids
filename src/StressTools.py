import pandas as pd
import numpy as np
import sympy as sym
import math
import multiprocessing
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

r, θ, φ, t = sym.symbols('r θ φ t', real = True)

get_principal1 = None
get_principal2 = None
get_principal_orientation = None
get_principal_orientation2 = None

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

    lat_radians = (90-lat) * RAD_MULTIPLIER #np.radians(lat)
    time_value = time

    for lon in range(MIN_LON, MAX_LON + 1, 10):
        if (lat == 90 or lon == 0):
            continue

        lon_radians = lon * RAD_MULTIPLIER # np.radians(lon)

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

def build_stress_field(satellite, orbit_time_seconds, rotations = 1, is_async = True):
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
        processes for peformance.  When False all calculation is done in the
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
    get_principal1 = sym.lambdify([t, φ, θ], satellite.PC1, modules = ["math", {"cot": math.atan}])
    get_principal2 = sym.lambdify([t, φ, θ], satellite.PC2, modules = ["math", {"cot": math.atan}])
    get_principal_orientation = sym.lambdify([t, φ, θ], satellite.PCΨ, modules = ["math", {"cot": math.atan}])
    get_principal_orientation2 = sym.lambdify([t, φ, θ], satellite.PCΨ2, modules = ["math", {"cot": math.atan}])

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
        time = (step/TIME_STEPS) * orbit_time_seconds
        for lat in range(MIN_LAT, MAX_LAT + 1, LAT_STEP_SIZE):
            if is_async:
                # pool.apply_async will schedule the processing within separate python processes
                # which allows the work to be distributed to multiple CPU cores
                pool.apply_async(get_stress_for_latitude,
                    args = (step, time, lat, ),
                    callback=callback,
                    error_callback=error_callback)
            else:
                data.extend(get_stress_for_latitude(step, time, lat))

    if is_async:
        pool.close()
        pool.join()

    df = pd.DataFrame(data)
    return df.sort_values(['latitude', 'longitude', 'time_step'])


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
                colat=np.radians(90-lat),
                lon=np.radians(360-lon),
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

get_stresses_for_point_cached = mem.cache(get_stresses_for_point)


def build_simon_stress_field(
    interior,
    pointFrame,
    phase,
    eccentricity,
    obliquity,
    nsr,
    is_async=True,
    tolerance=1,
    steps=360):
    stresses = []

    if is_async:
        pointStresses = Parallel(n_jobs=CPUS)(delayed(get_stresses_for_point_cached)\
            (interior, point.lon, point.lat, phase, tolerance, steps, eccentricity, obliquity, nsr) for point in pointFrame.itertuples())
    else:
        pointStresses = [get_stresses_for_point_cached(interior, point.lon, point.lat, phase, tolerance, steps, eccentricity, obliquity, nsr) for point in pointFrame.itertuples()]

    for stress in pointStresses:
        stresses.extend(stress)

    df = pd.DataFrame(stresses)

    return df

get_simon_stress_field = mem.cache(build_simon_stress_field)
get_stress_field = mem.cache(build_stress_field)
