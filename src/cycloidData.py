import fitting
import pandas as pd
import numpy as np
import collections

Cycloid = collections.namedtuple('Cycloid', 'points curve arcs')


def convertLonRaw(lon):
    newLon = lon - 180
    if newLon < 0:
        newLon = 360 + newLon

    return newLon

convertLon = np.vectorize(convertLonRaw)


def loadDelphi():
    delphi = pd.read_csv("./obsData/DelphiLonLatAT.txt", header=None, sep=' ', names=['lon', 'lat'])

    delphi_arcs = [
        delphi[0:8],
        delphi[8:14],
        delphi[14:24],
        delphi[24:33],
        delphi[33:]
    ]

    # delphiCurve = fitting.createCycloidBezier(delphi_arcs, maxError=0.008)
    delphiCurve = fitting.createCycloidBezier(delphi_arcs, maxError=0.1)

    return delphi, delphiCurve, delphi_arcs

def loadTyrrel():
    tyrrel = pd.read_csv("./obsData/TyrrelLonLat2.txt", header=None, sep=' ', names=['lon', 'lat'])

    tyrrel_arcs = [
        tyrrel[0:13],
        tyrrel[13:25],
        tyrrel[25:40],
        tyrrel[40:53],
        tyrrel[53:59],
        tyrrel[59:63],
        tyrrel[63:74],
        tyrrel[74:84],
        tyrrel[84:95],
        tyrrel[95:]
    ]

    tyrrelCurve = fitting.createCycloidBezier(tyrrel_arcs, maxError=0.01135)

    return tyrrel, tyrrelCurve, tyrrel_arcs


def loadAlex():
    # alex = pd.read_csv("./obsDataCorrected/Alex.csv")
    alex = pd.read_csv("./obsData/AlexLonLatCut.txt", header=None, sep=' ', names=['lon', 'lat'])

    alex_arcs = [
        alex[0:26],
        alex[26:50],
        alex[50:]
    ]

    # alexCurve = fitting.createCycloidBezier(alex_arcs, maxError=0.01135)
    alexCurve = fitting.createCycloidBezier(alex_arcs, maxError=0.1)

    return alex, alexCurve, alex_arcs


def loadSidon():
    sidon = pd.read_csv("./obsData/SidonLonLatAT.txt", header=None, sep=' ', names=['lon', 'lat'])

    sidon_arcs = [
        sidon[0:7],
        sidon[7:13],
        sidon[13:19],
        sidon[19:26],
        sidon[26:34],
        sidon[34:44],
        sidon[44:55],
        sidon[55:66],
        sidon[66:]
    ]

    sidonCurve = fitting.createCycloidBezier(sidon_arcs, maxError=0.012) # Smoothing the last arc wiggle

    return sidon, sidonCurve, sidon_arcs


def loadCarly():
    carly = pd.read_csv("./obsData/CarlyLonLatCut.txt", header=None, sep='  ', names=['lon', 'lat'], engine='python')

    carly_arcs = [
        carly[0:23],
        carly[23:38],
        carly[38:55],
        carly[55:82],
        carly[82:]
    ]

    carlyCurve = fitting.createCycloidBezier(carly_arcs, maxError=0.0085)

    return carly, carlyCurve, carly_arcs


def loadDirk():
    dirk = pd.read_csv("./obsData/DirkLonLat.txt", header=None, sep='\t', names=['lon', 'lat'])

    dirk_arcs = [
        dirk[0:18],
        dirk[18:28],
        dirk[28:34],
        dirk[34:41],
        dirk[41:51],
        dirk[51:]
    ]

    dirkCurve = fitting.createCycloidBezier(dirk_arcs)

    return dirk, dirkCurve, dirk_arcs


def loadYaphet():
    yaphet = pd.read_csv("./obsData/YaphetLonLat.txt", header=None, sep=' ', names=['lon', 'lat'])

    yaphet_arcs = [
        yaphet[0:9],
        yaphet[9:18],
        yaphet[18:29],
        yaphet[29:]
    ]

    yaphetCurve = fitting.createCycloidBezier(yaphet_arcs, maxError=0.12)

    return yaphet, yaphetCurve, yaphet_arcs


def loadOdessa():
    odessa = pd.read_csv("./obsData/OdessaLonLatP180.txt", header=None, sep='\t', names=['lon', 'lat'])

    odessa_arcs = [
        odessa[0:8],
        odessa[8:15],
        odessa[15:21],
        odessa[21:]
    ]

    odessaCurve = fitting.createCycloidBezier(odessa_arcs)
    odessaCurve['lon'] = convertLon(odessaCurve['lon'])
    odessa['lon'] = convertLon(odessa['lon'])

    return odessa, odessaCurve, odessa_arcs


def loadMira():
    mira = pd.read_csv("./obsData/MiraLonLatP180.txt", header=None, sep='\t', names=['lon', 'lat'])

    mira_arcs = [
        mira[0:8],
        mira[8:18],
        mira[18:21],
        mira[21:27],
        mira[27:]
    ]

    miraCurve = fitting.createCycloidBezier(mira_arcs)
    miraCurve['lon'] = convertLon(miraCurve['lon'])
    mira['lon'] = convertLon(mira['lon'])

    return mira, miraCurve, mira_arcs


def loadCilicia():
    cilicia = pd.read_csv("./obsData/ciliciaLonLatAT.txt", header=None, sep=' ', names=['lon', 'lat'])

    cilicia_arcs = [
        cilicia[0:15],
        cilicia[15:29],
        cilicia[29:43],
        cilicia[43:55],
        cilicia[55:65],
        cilicia[65:]
    ]

    ciliciaCurve = fitting.createCycloidBezier(cilicia_arcs, maxError=0.015)

    return cilicia, ciliciaCurve, cilicia_arcs

def loadAllCycloids():
    cycloids = {
        'alex': Cycloid(*loadAlex()),
        'carly': Cycloid(*loadCarly()),
        'cilicia': Cycloid(*loadCilicia()),
        'delphi': Cycloid(*loadDelphi()),
        'dirk': Cycloid(*loadDirk()),
        'mira': Cycloid(*loadMira()),
        'odessa': Cycloid(*loadOdessa()),
        'sidon': Cycloid(*loadSidon()),
        'tyrrel': Cycloid(*loadTyrrel()),
        'yaphet': Cycloid(*loadYaphet())
    }

    return cycloids


