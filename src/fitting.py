import pandas as pd
import numpy as np
import StressTools as tools

def match_orientations(curve, stresses):
    merged = curve.merge(
        stresses,
        how='left',
        on=['lon', 'lat', 'headingCategory']
    )
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

def test_arc(arc, phase_increment, interior, steps=360):
    results = []
    startingPoint = arc['pointNumber'].min()
    points = len(arc)

    for phase in range(0, 361, phase_increment):
        field = tools.get_simon_stress_field(interior, arc, phase, steps=steps)
        matches = match_orientations(arc, field)
        loss = calculate_loss(matches, startingPoint, points)

        results.append({
            'phase': phase,
            'error': loss
        })
        print(f'Calculated for phase {phase}')

    return pd.DataFrame(results) 