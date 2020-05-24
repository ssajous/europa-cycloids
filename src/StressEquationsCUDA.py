import numpy as np
import math
import utils
from numba import cuda, vectorize, guvectorize

PERIOD_SEC = 306000.0
N = 2.*np.pi/PERIOD_SEC
ZETA = np.pi/2
ZETA_COS_SQUARE = math.cos(ZETA)**math.cos(ZETA)
ZETA_SIN_SQUARE = math.sin(ZETA)**math.sin(ZETA)


@cuda.jit(device=True)
def stressththe(constA, constB, e, oblq, phase, colat, lon, he, le):

    beta20ththe = 0.75*(3.*he - 10.*le)*math.cos(2.*colat) + 0.75*(he - 2.*le)
    beta21ththe = 1.5*(3.*he - 10.*le)*math.sin(2.*colat)
    beta22ththe = -1.5*(3.*he - 10.*le)*math.cos(2.*colat) + 4.5*(he - 2.*le)

    ththe = constB*(-6.*e*beta20ththe*math.cos(constA) + e*beta22ththe*(4.*math.sin(2.*lon)*math.sin(constA)+3. *
                                                                        math.cos(2.*lon)*math.cos(constA)) + 4.*math.cos(oblq)*math.sin(oblq)*beta21ththe*math.cos(lon)*math.sin(phase + constA))

    return ththe


@cuda.jit(device=True)
def stressphphe(constA, constB, e, oblq, phase, colat, lon, he, le):

    beta20phphe = 0.75*(3.*he - 8.*le)*math.cos(2.*colat) + 0.75*(he - 4.*le)
    beta21phphe = 1.5*(3.*he - 8.*le)*math.sin(2.*colat)
    beta22phphe = -1.5*(3.*he - 8.*le)*math.cos(2.*colat) + 4.5*(he - 4.*le)

    phphe = constB*(-6.*e*beta20phphe*math.cos(constA) + e*beta22phphe*(4.*math.sin(2.*lon)*math.sin(constA)+3. *
                                                                        math.cos(2.*lon)*math.cos(constA)) + 4.*math.cos(oblq)*math.sin(oblq)*beta21phphe*math.cos(lon)*math.sin(phase + constA))

    return phphe


@cuda.jit(device=True)
def stressthphe(constA, constB, e, oblq, phase, colat, lon, he, le):

    beta21thphe = 3.*le*math.sin(colat)
    beta22thphe = 3.*le*math.cos(colat)

    thphe = constB*(2.*e*beta22thphe*(4.*math.cos(2.*lon)*math.sin(constA)-3.*math.sin(2.*lon)*math.cos(
        constA)) + 4.*math.cos(oblq)*math.sin(oblq)*beta21thphe*math.sin(lon)*math.sin(phase + constA))

    return thphe


@cuda.jit(device=True)
def modeththv(n, t, _lambda, sj, e, oblq, phase, colat, lon, hv, lv):

    beta20ththv = 0.75*(3.*hv - 10.*lv)*math.cos(2.*colat) + 0.75*(hv - 2.*lv)
    beta21ththv = 1.5*(3.*hv - 10.*lv)*math.sin(2.*colat)
    beta22ththv = -1.5*(3.*hv - 10.*lv)*math.cos(2.*colat) + 4.5*(hv - 2.*lv)

    ththv = (1./math.sqrt(1. + (-n/sj)*(-n/sj)))*(-6.*e*beta20ththv*math.cos(n*t - math.atan(-n/sj) + math.atan(_lambda)) + e*beta22ththv*(4.*math.sin(2.*lon)*math.sin(n*t - math.atan(-n/sj) + math.atan(_lambda)
                                                                                                                                                                        )+3.*math.cos(2.*lon)*math.cos(n*t - math.atan(-n/sj) + math.atan(_lambda))) + 4.*math.cos(oblq)*math.sin(oblq)*beta21ththv*math.cos(lon)*math.sin(phase + n*t - math.atan(-n/sj) + math.atan(_lambda)))

    return ththv


@cuda.jit(device=True)
def modephphv(n, t, _lambda, sj, e, oblq, phase, colat, lon, hv, lv):

    beta20phphv = 0.75*(3.*hv - 8.*lv)*math.cos(2.*colat) + 0.75*(hv - 4.*lv)
    beta21phphv = 1.5*(3.*hv - 8.*lv)*math.sin(2.*colat)
    beta22phphv = -1.5*(3.*hv - 8.*lv)*math.cos(2.*colat) + 4.5*(hv - 4.*lv)

    phphv = (1./math.sqrt(1. + (-n/sj)*(-n/sj)))*(-6.*e*beta20phphv*math.cos(n*t - math.atan(-n/sj) + math.atan(_lambda)) + e*beta22phphv*(4.*math.sin(2.*lon)*math.sin(n*t - math.atan(-n/sj) + math.atan(_lambda)
                                                                                                                                                                        )+3.*math.cos(2.*lon)*math.cos(n*t - math.atan(-n/sj) + math.atan(_lambda))) + 4.*math.cos(oblq)*math.sin(oblq)*beta21phphv*math.cos(lon)*math.sin(phase + n*t - math.atan(-n/sj) + math.atan(_lambda)))

    return phphv


@cuda.jit(device=True)
def modethphv(n, t, _lambda, sj, e, oblq, phase, colat, lon, hv, lv):

    beta21thphv = 3.*lv*math.sin(colat)
    beta22thphv = 3.*lv*math.cos(colat)

    thphv = (1./math.sqrt(1. + (-n/sj)*(-n/sj)))*(8.*e*beta22thphv*math.cos(2.*lon)*math.sin(n*t - math.atan(-n/sj) + math.atan(_lambda)) - 6.*e*beta22thphv*math.sin(2.*lon) *
                                                  math.cos(n*t - math.atan(-n/sj) + math.atan(_lambda)) + 4.*math.cos(oblq)*math.sin(oblq)*beta21thphv*math.sin(lon)*math.sin(phase + n*t - math.atan(-n/sj) + math.atan(_lambda)))

    return thphv


@cuda.jit(device=True)
def stressththeNSR(constA_NSR, constB_NSR, colat, lon, he_NSR, le_NSR):
    alpha22ththe = -1.5*(3.*he_NSR - 10.*le_NSR) * \
        math.cos(2.*colat) + 4.5*(he_NSR - 2.*le_NSR)
    ththeNSR = constB_NSR*alpha22ththe*math.cos(2.*lon + constA_NSR)

    return ththeNSR


@cuda.jit(device=True)
def stressphpheNSR(constA_NSR, constB_NSR, colat, lon, he_NSR, le_NSR):

    alpha22phphe = -1.5*(3.*he_NSR - 8.*le_NSR) * \
        math.cos(2.*colat) + 4.5*(he_NSR - 4.*le_NSR)
    phpheNSR = constB_NSR*alpha22phphe*math.cos(2.*lon + constA_NSR)

    return phpheNSR


@cuda.jit(device=True)
def stressthpheNSR(constA_NSR, constB_NSR, colat, lon, he_NSR, le_NSR):

    alpha22thphe = 3.*le_NSR*math.cos(colat)
    thpheNSR = -2.*constB_NSR*alpha22thphe*math.sin(2.*lon + constA_NSR)

    return thpheNSR


@cuda.jit(device=True)
def modeththvNSR(constA_NSR, NSRrate, sj_NSR, colat, lon, hv_NSR, lv_NSR):

    alpha22ththv = -1.5*(3.*hv_NSR - 10.*lv_NSR) * \
        math.cos(2.*colat) + 4.5*(hv_NSR - 2.*lv_NSR)
    ththvNSR = (1./math.sqrt(1. + (-2.*NSRrate/sj_NSR)*(-2.*NSRrate/sj_NSR))) * \
        alpha22ththv*math.cos(2.*lon + constA_NSR -
                              math.atan(-2.*NSRrate/sj_NSR))

    return ththvNSR


@cuda.jit(device=True)
def modephphvNSR(constA_NSR, NSRrate, sj_NSR, colat, lon, hv_NSR, lv_NSR):

    alpha22phphv = -1.5*(3.*hv_NSR - 8.*lv_NSR) * \
        math.cos(2.*colat) + 4.5*(hv_NSR - 4.*lv_NSR)
    phphvNSR = (1./math.sqrt(1. + (-2.*NSRrate/sj_NSR)*(-2.*NSRrate/sj_NSR))) * \
        alpha22phphv*math.cos(2.*lon + constA_NSR -
                              math.atan(-2.*NSRrate/sj_NSR))

    return phphvNSR


@cuda.jit(device=True)
def modethphvNSR(constA_NSR, NSRrate, sj_NSR, colat, lon, hv_NSR, lv_NSR):

    alpha22thphv = 3.*lv_NSR*math.cos(colat)
    thphvNSR = -2.*(1./math.sqrt(1. + (-2.*NSRrate/sj_NSR)*(-2.*NSRrate/sj_NSR))) * \
        alpha22thphv*math.sin(2.*lon + constA_NSR -
                              math.atan(-2.*NSRrate/sj_NSR))

    return thphvNSR


@vectorize(['float64(float64, float64,float64,float64,float64,float64,float64,float64,float64,float64,float64)'], target='cuda')
def tempGetRawStressThTheVectorized(_lambda, constB, he, le, eccentricity, colat, lon, steps, this_step, oblq, phase):
    t = (this_step/steps)*PERIOD_SEC
    constA = N * t + math.atan(_lambda)

    return stressththe(constA, constB, eccentricity, oblq, phase,
                       colat, lon, he, le)


@vectorize(['float64(float64, float64,float64,float64,float64,float64,float64,float64,float64,float64,float64)'], target='cuda')
def tempGetRawStressPhPheVectorized(_lambda, constB, he, le, eccentricity, colat, lon, steps, this_step, oblq, phase):
    t = (this_step/steps)*PERIOD_SEC
    constA = N * t + math.atan(_lambda)

    return stressphphe(constA, constB, eccentricity, oblq, phase,
                       colat, lon, he, le)


@vectorize(['float64(float64, float64,float64,float64,float64,float64,float64,float64,float64,float64,float64)'], target='cuda')
def tempGetRawStressThPheVectorized(_lambda, constB, he, le, eccentricity, colat, lon, steps, this_step, oblq, phase):
    t = (this_step/steps)*PERIOD_SEC
    constA = N * t + math.atan(_lambda)

    return stressthphe(constA, constB, eccentricity, oblq, phase,
                       colat, lon, he, le)


def tempGetRawStress(interior, curve, oblq, phase):
    # Prep items which will not change with each loop
    _lambda = interior.rigidity / (interior.viscosity * N)
    constB = 0.5*((N * N * interior.radius * interior.rigidity) /
                  interior.surface_gravity) * (1. / math.sqrt(1. + _lambda * _lambda))

    coordinates = np.array(curve[['lon', 'lat']])
    times = np.arange(0, 360)
    combinations = np.concatenate(
        np.array([[(c[0], c[1], time) for time in times] for c in coordinates]))
    lons = cuda.to_device(np.ascontiguousarray(combinations[:, 0:1]))
    lats = cuda.to_device(np.ascontiguousarray(combinations[:, 1:2]))
    fullTimes = cuda.to_device(np.ascontiguousarray(combinations[:, 2:3]))

    ththe = tempGetRawStressThTheVectorized(_lambda, constB, interior.he, interior.le, interior.eccentricity,
                                            lats, lons, 360, fullTimes, oblq, phase)
    thphe = tempGetRawStressThPheVectorized(_lambda, constB, interior.he, interior.le, interior.eccentricity,
                                            lats, lons, 360, fullTimes, oblq, phase)
    phphe = tempGetRawStressPhPheVectorized(_lambda, constB, interior.he, interior.le, interior.eccentricity,
                                            lats, lons, 360, fullTimes, oblq, phase)

    return ththe, thphe, phphe


def tempGetStressVectorized(_lambda, constB, he, le, modal_strengths, hv, lv, eccentricity, colat, lon, steps, this_step, oblq, phase):
    t = (this_step/steps)*PERIOD_SEC
    constA = N * t + math.atan(_lambda)

    ththe = stressththe(constA, constB, eccentricity, oblq, phase,
                        colat, lon, he, le)
    phphe = stressphphe(constA, constB, eccentricity, oblq, phase,
                        colat, lon, he, le)
    thphe = stressthphe(constA, constB, eccentricity, oblq, phase,
                        colat, lon, he, le)

    ththv_total = 0
    phphv_total = 0
    thphv_total = 0

    for index in range(len(modal_strengths)):
        ththv_total += modeththv(N, t, _lambda, modal_strengths[index], eccentricity, oblq,
                                 phase, colat, lon, hv[index], lv[index])
        phphv_total += modephphv(N, t, _lambda, modal_strengths[index], eccentricity, oblq,
                                 phase, colat, lon, hv[index], lv[index])
        thphv_total += modethphv(N, t, _lambda, modal_strengths[index], eccentricity, oblq,
                                 phase, colat, lon, hv[index], lv[index])

    myStressThTh = ththe + constB * thphv_total
    myStressPhPh = phphe + constB * phphv_total
    myStressThPh = thphe + constB * thphv_total

    sigTheta = myStressThTh * (ZETA_COS_SQUARE) + myStressPhPh * (ZETA_SIN_SQUARE) + \
        myStressThPh * math.sin(2.*ZETA)  # Corresponds to Terry's sigma 1
    sigPhi = myStressThTh * ZETA_SIN_SQUARE + myStressPhPh * ZETA_COS_SQUARE - \
        myStressThPh * math.sin(2.*ZETA)      # Corresponds to Terry's sigma 2

    sigThetaKPA = sigTheta*1E-3
    sigPhiKPA = sigPhi*1E-3

    zetaCW = (2.*np.pi)-ZETA

    if (sigThetaKPA < sigPhiKPA):
        stress = sigPhiKPA         # In this case, sigPhi is the max stress, so the heading should be perpendicular to sigPhi. SigTheta is -| to sigPhi by definition, so heading = zetaCW
        heading = zetaCW
    else:
        stress = sigThetaKPA
        # In this case, sigTheta is the max stress, so the heading should be perpendicular to sigTheta. SigTheta is oriented along zeta, so heading = zetaCW + 90 deg
        heading = zetaCW + (np.pi/2.)

    # Making sure azimuth values fall between 0 and 360. WOULD LOVE TO CHANGE THIS TO MOD AND INCORPORATE ABOVE.
    if (heading >= (2.*np.pi)):
        # Also finds the two, complementary heading directions (e.g. 45 and 135)
        heading = heading - (2.*np.pi)
        heading2 = heading + np.pi
    else:
        heading2 = heading - np.pi

    # Determines which of the two heading directions (e.g. 45 or 135) is largest and selects that as the output heading.
    if (heading > heading2):
        bigHeading = heading
    else:
        bigHeading = heading2

    return (stress, bigHeading)


def getStress(interior_value, e_in, colat, lon, steps, this_step, oblq, phase, NSRdelta):
    if isinstance(interior_value, str):
        interior = utils.import_interior(interior_value)
    else:
        interior = interior_value

    periodInSec = 306000.0  # sec;
    n = 2.*np.pi/periodInSec  # [rad/sec]

    t = (this_step/steps)*periodInSec

    _lambda = interior.rigidity / (interior.viscosity * n)

    if (e_in is None):
        ecc = interior.eccentricity
    else:
        ecc = e_in

    constA = n*t + math.atan(_lambda)
    constB = 0.5*((n * n * interior.radius * interior.rigidity) /
                  interior.surface_gravity) * (1. / math.sqrt(1. + _lambda * _lambda))

    ththe = stressththe(constA, constB, ecc, oblq, phase,
                        colat, lon, interior.he, interior.le)
    phphe = stressphphe(constA, constB, ecc, oblq, phase,
                        colat, lon, interior.he, interior.le)
    thphe = stressthphe(constA, constB, ecc, oblq, phase,
                        colat, lon, interior.he, interior.le)

    if len(interior.modal_strengths) < 2:
        print("Number of modes is disallowed")
        return None
    ththv_vals = []
    phphv_vals = []
    thphv_vals = []

    for index in range(len(interior.modal_strengths)):
        ththv = modeththv(n, t, _lambda, interior.modal_strengths[index], ecc, oblq,
                          phase, colat, lon, interior.hv[index], interior.lv[index])
        phphv = modephphv(n, t, _lambda, interior.modal_strengths[index], ecc, oblq,
                          phase, colat, lon, interior.hv[index], interior.lv[index])
        thphv = modethphv(n, t, _lambda, interior.modal_strengths[index], ecc, oblq,
                          phase, colat, lon, interior.hv[index], interior.lv[index])

        ththv_vals.append(ththv)
        phphv_vals.append(phphv)
        thphv_vals.append(thphv)

    myStressThTh = ththe + constB * \
        np.sum(ththv_vals)
    myStressPhPh = phphe + constB * \
        np.sum(phphv_vals)
    myStressThPh = thphe + constB * \
        np.sum(thphv_vals)  # negative sign comes later

    # sum1 = np.sum(ththv_vals)
    # sum2 = np.sum(phphv_vals)
    # sum3 = np.sum(thphv_v)

    # COMPUTING NSR STRESSES FOR EUROPA

    if NSRdelta != 0:

        # RIGIDITY AND VISCOSITY were swapped in original version, now corrected and consistent with diurnal eqs above

        # delta of 0.1 for NSR period of 11419 years; 43 roughly equals 6 Myr period // 6/23: removed *year2sec
        NSRrate = (interior.rigidity / (2. * interior.viscosity * NSRdelta))

        constA_NSR = 2. * NSRrate * t + math.atan(NSRdelta)
        #print("at constB NSR")
        constB_NSR = 0.5 * ((n * n * interior.radius * interior.rigidity) /
                            interior.surface_gravity) * (1. / math.sqrt(1. + NSRdelta * NSRdelta))

        # Currently sj values are from Hermes' paper for a reference model that is very similar to the visco model he sent us (M1)

        ththeNSR = stressththeNSR(
            constA_NSR, constB_NSR, colat, lon, interior.he_NSR, interior.le_NSR)
        phpheNSR = stressphpheNSR(
            constA_NSR, constB_NSR, colat, lon, interior.he_NSR, interior.le_NSR)
        thpheNSR = stressthpheNSR(
            constA_NSR, constB_NSR, colat, lon, interior.he_NSR, interior.le_NSR)

        if (len(interior.modal_strengths_nsr) < 2):
            print("Number of NSR modes is disallowed")
            return None

        ththvNSR_vals = []
        phphvNSR_vals = []
        thphvNSR_vals = []

        for index in range(len(interior.modal_strengths_nsr)):
            ththvNSR = modeththvNSR(
                constA_NSR, NSRrate, interior.modal_strengths_nsr[index], colat, lon,
                interior.hv_NSR[index], interior.lv_NSR[index])
            phphvNSR = modephphvNSR(
                constA_NSR, NSRrate, interior.modal_strengths_nsr[index], colat, lon,
                interior.hv_NSR[index], interior.lv_NSR[index])
            thphvNSR = modethphvNSR(
                constA_NSR, NSRrate, interior.modal_strengths_nsr[index], colat, lon,
                interior.hv_NSR[index], interior.lv_NSR[index])

            ththvNSR_vals.append(ththvNSR)
            phphvNSR_vals.append(phphvNSR)
            thphvNSR_vals.append(thphvNSR)

        # Combining stresses

        myStressThThNSR = ththeNSR + constB_NSR * \
            np.sum(ththvNSR_vals)
        myStressPhPhNSR = phpheNSR + constB_NSR * \
            np.sum(phphvNSR_vals)
        myStressThPhNSR = thpheNSR + constB_NSR * \
            np.sum(thphvNSR_vals)

        myStressThThTot = myStressThTh + myStressThThNSR
        myStressPhPhTot = myStressPhPh + myStressPhPhNSR
        myStressThPhTot = myStressThPh + myStressThPhNSR

    else:
        myStressThThTot = myStressThTh
        myStressPhPhTot = myStressPhPh
        myStressThPhTot = myStressThPh

    # Converting COMBINED STRESSES to principle stresses.

    if (myStressThThTot == myStressPhPhTot):
        zeta = np.pi/2
    else:
        # Amount by which coordinates are being rotated to get principal stresses; changes CCW, orientation of sigTheta
        zeta = 0.5*math.atan((2.*myStressThPhTot) /
                             (myStressThThTot-myStressPhPhTot))

    sigTheta = myStressThThTot*(np.square(math.cos(zeta)))+myStressPhPhTot*(np.square(
        math.sin(zeta)))+myStressThPhTot*math.sin(2.*zeta)  # Corresponds to Terry's sigma 1
    sigPhi = myStressThThTot*(np.square(math.sin(zeta)))+myStressPhPhTot*(np.square(
        math.cos(zeta)))-myStressThPhTot*math.sin(2.*zeta)      # Corresponds to Terry's sigma 2

    sigThetaKPA = sigTheta*1E-3
    sigPhiKPA = sigPhi*1E-3

    zetaCW = (2.*np.pi)-zeta   # Changes to CW, still oriented along sigTheta

    # stress, The largest (i.e. most tensile) principal stress
    # heading, Direction of crack propagation/formation. Perpendicular to direction of largest principal stress (i.e. max tension).

    # Determines which stress is most tensile and largest heading.

    if (sigThetaKPA < sigPhiKPA):
        stress = sigPhiKPA         # In this case, sigPhi is the max stress, so the heading should be perpendicular to sigPhi. SigTheta is -| to sigPhi by definition, so heading = zetaCW
        heading = zetaCW
    else:
        stress = sigThetaKPA
        # In this case, sigTheta is the max stress, so the heading should be perpendicular to sigTheta. SigTheta is oriented along zeta, so heading = zetaCW + 90 deg
        heading = zetaCW + (np.pi/2.)

    # Making sure azimuth values fall between 0 and 360. WOULD LOVE TO CHANGE THIS TO MOD AND INCORPORATE ABOVE.
    if (heading >= (2.*np.pi)):
        # Also finds the two, complementary heading directions (e.g. 45 and 135)
        heading = heading - (2.*np.pi)
        heading2 = heading + np.pi
    else:
        heading2 = heading - np.pi

    # Determines which of the two heading directions (e.g. 45 or 135) is largest and selects that as the output heading.
    if (heading > heading2):
        bigHeading = heading
    else:
        bigHeading = heading2

    return (stress, bigHeading)
