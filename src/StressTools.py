import pandas as pd
import numpy as np
import sympy as sym
import math
import multiprocessing

TIME_STEPS = 360
MIN_LAT = -75
MAX_LAT = 90
MIN_LON = 0
MAX_LON = 360
LAT_STEP_SIZE = 15
LON_STEP_SIZE = 10
RAD_MULTIPLIER = np.pi / 180
DEG_MULTIPLIER = 180 / np.pi

r, θ, φ, t = sym.symbols('r θ φ t', real = True)

get_principal1 = None
get_principal2 = None
get_principal_orientation = None
get_principal_orientation2 = None

def get_stress_for_latitude(step, lat):
    results = []

    lat_radians = lat * RAD_MULTIPLIER #np.radians(lat)
    step_value = step / TIME_STEPS
    
    for lon in range(MIN_LON, MAX_LON + 1, 10):
        if (lat == 90 or lon == 0):
            continue
            
        lon_radians = lon * RAD_MULTIPLIER # np.radians(lon)

        principal1 = get_principal1(step_value, lat_radians, lon_radians)
        principal2 = get_principal2(step, lat_radians, lon_radians)
        principal_phi = get_principal_orientation(step, lat_radians, lon_radians)
        principal_phi2 = get_principal_orientation2(step, lat_radians, lon_radians)        

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

def build_stress_field(satellite, is_async = True):
    """ 
    Creates a data frame with the results of stress calculations for a range of 
    latitudes and longitudes across 360 time steps

    Parameters
    ----------
    satelite : MEWtools.satellite
        The body on which stress values will be calculated
    is_async: boolean, optional
        When True (default) the calculation will be spread across multiple
        processes for peformance.  When False all calculation is done in the 
        current process, this most likely slower, but better for debugging.
    """

    # Create "lamdified" versions of each of the stress equations.  This turns
    # the symbolic python formulas into standar python functions.  Which improves their 
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

    for step in range(TIME_STEPS):
        for lat in range(MIN_LAT, MAX_LAT + 1, LAT_STEP_SIZE):
            if is_async:
                # pool.apply_async will schedule the processing within separate python processes
                # which allows the work to be distributed to multiple CPU cores
                pool.apply_async(get_stress_for_latitude, args = (step, lat, ), callback=callback, error_callback=error_callback)
            else:
                data.extend(get_stress_for_latitude(step, lat))
    
    if is_async:
        pool.close()
        pool.join()

    df = pd.DataFrame(data)
    return df.sort_values(['latitude', 'longitude', 'time_step'])