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

r, θ, φ, t = sym.symbols('r θ φ t', real = True)


def build_stress_field(satellite, is_async = False):
    get_principal1 = sym.lambdify([t, φ, θ], satellite.PC1, modules = ["math", {"cot": math.atan}])
    get_principal2 = sym.lambdify([t, φ, θ], satellite.PC2, modules = ["math", {"cot": math.atan}])
    get_principal_orientation = sym.lambdify([t, φ, θ], satellite.PCΨ, modules = ["math", {"cot": math.atan}])
    get_principal_orientation2 = sym.lambdify([t, φ, θ], satellite.PCΨ2, modules = ["math", {"cot": math.atan}])

    def get_stress_for_latitude(step, lat):
        results = []
        lat_radians = np.radians(lat)
        step_value = step / TIME_STEPS
        
        for lon in range(MIN_LON, MAX_LON + 1, 10):
            if (lat == 90 or lon == 0):
                continue
                
            lon_radians = np.radians(lon)

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
                    'principal_orientation': np.rad2deg(principal_phi),
                    'principal_orientation2': np.rad2deg(principal_phi2),
                    'max_stress': max_stress,
                    'max_stress_orientation': np.rad2deg(max_stress_orientation)
                })
            
        return results

    data = []

    def callback(items):
        data.extend(items)

    if is_async:
        pool = multiprocessing.Pool()

    for step in range(TIME_STEPS):
        for lat in range(MIN_LAT, MAX_LAT + 1, LAT_STEP_SIZE):
            if is_async:
                pool.apply_async(get_stress_for_latitude, args = (step, lat, ), callback=callback)
            else:
                data.extend(get_stress_for_latitude(step, lat))
    
    if is_async:
        pool.close()
        pool.join()

    df = pd.DataFrame(data)
    return df.sort_values(['latitude', 'longitude', 'time_step'])