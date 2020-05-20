import fitting
import pandas as pd
import numpy as np
import collections
from joblib import Memory

CACHE_DIR = './cache'
mem = Memory(CACHE_DIR, verbose=0)

Cycloid = collections.namedtuple('Cycloid', 'points curve arcs')


def convertLonRaw(lon):
    newLon = lon - 180
    if newLon < 0:
        newLon = 360 + newLon

    return newLon

convertLon = np.vectorize(convertLonRaw)


def loadDelphi(points):
    delphi = pd.read_csv("./obsData/DelphiLonLatAT.txt", header=None, sep=' ', names=['lon', 'lat'])

    delphi_arcs = [
        delphi[0:8],
        delphi[7:14],
        delphi[13:24],
        delphi[23:33],
        delphi[32:]
    ]

    # delphiCurve = fitting.createCycloidBezier(delphi_arcs, maxError=0.008)
    delphiCurve = fitting.createCycloidBezier(delphi_arcs, maxError=0.1, pointsPerCurve=points)

    return delphi, delphiCurve, delphi_arcs

def loadTyrrel(points):
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

    tyrrelCurve = fitting.createCycloidBezier(tyrrel_arcs, maxError=0.09, pointsPerCurve=points)

    return tyrrel, tyrrelCurve, tyrrel_arcs


def loadAlex(points):
    # alex = pd.read_csv("./obsDataCorrected/Alex.csv")
    alex = pd.read_csv("./obsData/AlexLonLatCut.txt", header=None, sep=' ', names=['lon', 'lat'])

    alex_arcs = [
        alex[0:26],
        alex[25:50],
        alex[49:]
    ]

    # alexCurve = fitting.createCycloidBezier(alex_arcs, maxError=0.01135)
    alexCurve = fitting.createCycloidBezier(alex_arcs, maxError=0.1, pointsPerCurve=points)

    return alex, alexCurve, alex_arcs


def loadSidon(points):
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
    sidonCurve = fitting.createCycloidBezier(sidon_arcs, maxError=0.065, pointsPerCurve=points)

    return sidon, sidonCurve, sidon_arcs


def loadCarly(points):
    carly = pd.read_csv("./obsData/CarlyLonLatCut.txt", header=None, sep='  ', names=['lon', 'lat'], engine='python')

    carly_arcs = [
        carly[0:23],
        carly[22:38],
        carly[37:55],
        carly[54:82],
        carly[81:]
    ]

    # carlyCurve = fitting.createCycloidBezier(carly_arcs, maxError=0.0085)
    carlyCurve = fitting.createCycloidBezier(carly_arcs, maxError=0.035, pointsPerCurve=points)

    return carly, carlyCurve, carly_arcs


def loadDirk(points):
    dirk = pd.read_csv("./obsData/DirkLonLat.txt", header=None, sep='\t', names=['lon', 'lat'])

    dirk_arcs = [
        dirk[0:18],
        dirk[17:28],
        dirk[27:34],
        dirk[33:41],
        dirk[40:51],
        dirk[50:]
    ]

    dirkCurve = fitting.createCycloidBezier(dirk_arcs, pointsPerCurve=points)

    return dirk, dirkCurve, dirk_arcs


def loadYaphet(points):
    yaphet = pd.read_csv("./obsData/YaphetLonLat.txt", header=None, sep=' ', names=['lon', 'lat'])

    yaphet_arcs = [
        yaphet[0:12],
        yaphet[11:23],
        yaphet[22:32],
        yaphet[31:]
    ]

    yaphetCurve = fitting.createCycloidBezier(yaphet_arcs, maxError=0.12, pointsPerCurve=points)

    return yaphet, yaphetCurve, yaphet_arcs


def loadOdessa(points):
    odessa = pd.read_csv("./obsData/OdessaLonLatP180.txt", header=None, sep='\t', names=['lon', 'lat'])

    odessa_arcs = [
        odessa[0:10],
        odessa[9:16],
        odessa[15:23],
        odessa[22:]
    ]

    odessaCurve = fitting.createCycloidBezier(odessa_arcs, maxError=0.055, pointsPerCurve=points)
    odessaCurve['lon'] = convertLon(odessaCurve['lon'])
    odessa['lon'] = convertLon(odessa['lon'])

    return odessa, odessaCurve, odessa_arcs


def loadMira(points):
    mira = pd.read_csv("./obsData/MiraLonLatP180.txt", header=None, sep='\t', names=['lon', 'lat'])

    mira_arcs = [
        mira[0:11],
        mira[10:18],
        mira[17:21],
        mira[20:31],
        mira[30:]
    ]

    miraCurve = fitting.createCycloidBezier(mira_arcs, maxError=0.107, pointsPerCurve=points)
    miraCurve['lon'] = convertLon(miraCurve['lon'])
    mira['lon'] = convertLon(mira['lon'])

    return mira, miraCurve, mira_arcs


def loadCilicia(points):
    cilicia = pd.read_csv("./obsData/ciliciaLonLatAT.txt", header=None, sep=' ', names=['lon', 'lat'])

    cilicia_arcs = [
        cilicia[0:15],
        cilicia[14:29],
        cilicia[28:43],
        cilicia[42:55],
        cilicia[54:65],
        cilicia[64:]
    ]

    # ciliciaCurve = fitting.createCycloidBezier(cilicia_arcs, maxError=0.015)
    ciliciaCurve = fitting.createCycloidBezier(cilicia_arcs, maxError=0.06, pointsPerCurve=points)

    return cilicia, ciliciaCurve, cilicia_arcs

def loadAllCycloidsRaw(pointsPerCurve=100):
    cycloids = {
        'alex': Cycloid(*loadAlex(pointsPerCurve)),
        'carly': Cycloid(*loadCarly(pointsPerCurve)),
        'cilicia': Cycloid(*loadCilicia(pointsPerCurve)),
        'delphi': Cycloid(*loadDelphi(pointsPerCurve)),
        'dirk': Cycloid(*loadDirk(pointsPerCurve)),
        'mira': Cycloid(*loadMira(pointsPerCurve)),
        'odessa': Cycloid(*loadOdessa(pointsPerCurve)),
        'sidon': Cycloid(*loadSidon(pointsPerCurve)),
        'tyrrel': Cycloid(*loadTyrrel(pointsPerCurve)),
        'yaphet': Cycloid(*loadYaphet(pointsPerCurve))
    }

    return cycloids


loadAllCycloids = mem.cache(loadAllCycloidsRaw)

