import fitting
import pandas as pd
import numpy as np
import collections
from joblib import Memory

CACHE_DIR = './cache'
mem = Memory(CACHE_DIR, verbose=0)

Cycloid = collections.namedtuple('Cycloid', 'points curve arcs')


def convert_lon_raw(lon):
    new_lon = lon - 180
    if new_lon < 0:
        new_lon = 360 + new_lon

    return new_lon


convert_lon = np.vectorize(convert_lon_raw)


def load_delphi(points):
    delphi = pd.read_csv("./obsData/DelphiLonLatAT.txt", header=None, sep=' ', names=['lon', 'lat'])

    delphi_arcs = [
        delphi[0:8],
        delphi[7:15],
        delphi[14:24],
        delphi[23:33],
        delphi[32:]
    ]

    # delphiCurve = fitting.createCycloidBezier(delphi_arcs, maxError=0.008)
    delphi_curve = fitting.create_cycloid_bezier(delphi_arcs, max_error=0.1, points_per_curve=points)

    return delphi, delphi_curve, delphi_arcs


def load_tyrrel(points):
    tyrrel = pd.read_csv("./obsData/TyrrelLonLat2.txt", header=None, sep=' ', names=['lon', 'lat'])

    tyrrel_arcs = [
        tyrrel[0:13],
        tyrrel[12:24],
        tyrrel[23:34],
        tyrrel[33:44],
        tyrrel[43:54],
        tyrrel[53:67],
        tyrrel[66:74],
        tyrrel[73:83],
        tyrrel[82:94],
        tyrrel[93:]
    ]

    tyrrel_curve = fitting.create_cycloid_bezier(tyrrel_arcs, max_error=0.09, points_per_curve=points)

    return tyrrel, tyrrel_curve, tyrrel_arcs


def load_alex(points):
    # alex = pd.read_csv("./obsDataCorrected/Alex.csv")
    alex = pd.read_csv("./obsData/AlexLonLatCut.txt", header=None, sep=' ', names=['lon', 'lat'])

    alex_arcs = [
        alex[0:26],
        alex[25:50],
        alex[49:]
    ]

    # alexCurve = fitting.createCycloidBezier(alex_arcs, maxError=0.01135)
    alex_curve = fitting.create_cycloid_bezier(alex_arcs, max_error=0.1, points_per_curve=points)

    return alex, alex_curve, alex_arcs


def load_sidon(points):
    sidon = pd.read_csv("./obsData/SidonLonLatAT.txt", header=None, sep=' ', names=['lon', 'lat'])

    sidon_arcs = [
        sidon[0:7],
        sidon[6:19],
        sidon[18:29],
        sidon[28:39],
        sidon[38:47],
        sidon[46:54],
        sidon[53:60],
        sidon[59:66],
        sidon[65:]
    ]

    # sidonCurve = fitting.createCycloidBezier(sidon_arcs, maxError=0.012) # Smoothing the last arc wiggle
    sidon_curve = fitting.create_cycloid_bezier(sidon_arcs, max_error=0.065, points_per_curve=points)

    return sidon, sidon_curve, sidon_arcs


def load_carly(points):
    carly = pd.read_csv("./obsData/CarlyLonLatCut.txt", header=None, sep='  ', names=['lon', 'lat'], engine='python')

    carly_arcs = [
        carly[0:23],
        carly[22:38],
        carly[37:55],
        carly[54:82],
        carly[81:]
    ]

    # carlyCurve = fitting.createCycloidBezier(carly_arcs, maxError=0.0085)
    carly_curve = fitting.create_cycloid_bezier(carly_arcs, max_error=0.035, points_per_curve=points)

    return carly, carly_curve, carly_arcs


def load_dirk(points):
    dirk = pd.read_csv("./obsData/DirkLonLat.txt", header=None, sep='\t', names=['lon', 'lat'])

    dirk_arcs = [
        dirk[0:18],
        dirk[17:28],
        dirk[27:34],
        dirk[33:41],
        dirk[40:51],
        dirk[50:]
    ]

    dirk_curve = fitting.create_cycloid_bezier(dirk_arcs, points_per_curve=points)

    return dirk, dirk_curve, dirk_arcs


def load_yaphet(points):
    yaphet = pd.read_csv("./obsData/YaphetLonLat.txt", header=None, sep=' ', names=['lon', 'lat'])

    yaphet_arcs = [
        yaphet[0:12],
        yaphet[11:23],
        yaphet[22:32],
        yaphet[31:]
    ]

    yaphet_curve = fitting.create_cycloid_bezier(yaphet_arcs, max_error=0.12, points_per_curve=points)

    return yaphet, yaphet_curve, yaphet_arcs


def load_odessa(points):
    odessa = pd.read_csv("./obsData/OdessaLonLatP180.txt", header=None, sep='\t', names=['lon', 'lat'])

    odessa_arcs = [
        odessa[0:10],
        odessa[9:16],
        odessa[15:23],
        odessa[22:]
    ]

    odessa_curve = fitting.create_cycloid_bezier(odessa_arcs, max_error=0.055, points_per_curve=points)
    odessa_curve['lon'] = convert_lon(odessa_curve['lon'])
    odessa['lon'] = convert_lon(odessa['lon'])

    return odessa, odessa_curve, odessa_arcs


def load_mira(points):
    mira = pd.read_csv("./obsData/MiraLonLatP180.txt", header=None, sep='\t', names=['lon', 'lat'])

    mira_arcs = [
        mira[0:11],
        mira[10:18],
        mira[17:21],
        mira[20:31],
        mira[30:]
    ]

    mira_curve = fitting.create_cycloid_bezier(mira_arcs, max_error=0.107, points_per_curve=points)
    mira_curve['lon'] = convert_lon(mira_curve['lon'])
    mira['lon'] = convert_lon(mira['lon'])

    return mira, mira_curve, mira_arcs


def load_cilicia(points):
    cilicia = pd.read_csv("./obsData/CiliciaLonLatAT.txt", header=None, sep=' ', names=['lon', 'lat'])

    cilicia_arcs = [
        cilicia[0:15],
        cilicia[14:29],
        cilicia[28:43],
        cilicia[42:55],
        cilicia[54:65],
        cilicia[64:]
    ]

    # ciliciaCurve = fitting.createCycloidBezier(cilicia_arcs, maxError=0.015)
    cilicia_curve = fitting.create_cycloid_bezier(cilicia_arcs, max_error=0.06, points_per_curve=points)

    return cilicia, cilicia_curve, cilicia_arcs


def load_all_cycloids_raw(points_per_curve=100):
    cycloids = {
        'alex': Cycloid(*load_alex(points_per_curve)),
        'carly': Cycloid(*load_carly(points_per_curve)),
        'cilicia': Cycloid(*load_cilicia(points_per_curve)),
        'delphi': Cycloid(*load_delphi(points_per_curve)),
        'dirk': Cycloid(*load_dirk(points_per_curve)),
        'mira': Cycloid(*load_mira(points_per_curve)),
        'odessa': Cycloid(*load_odessa(points_per_curve)),
        'sidon': Cycloid(*load_sidon(points_per_curve)),
        'tyrrel': Cycloid(*load_tyrrel(points_per_curve)),
        'yaphet': Cycloid(*load_yaphet(points_per_curve))
    }

    return cycloids


load_all_cycloids = mem.cache(load_all_cycloids_raw)
