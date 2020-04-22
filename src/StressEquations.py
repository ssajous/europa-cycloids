import numpy as np
import utils


def stressththe(constA, constB, e, oblq, phase, colat, lon, he, le):

    beta20ththe = 0.75*(3.*he - 10.*le)*np.cos(2.*colat) + 0.75*(he - 2.*le)
    beta21ththe = 1.5*(3.*he - 10.*le)*np.sin(2.*colat)
    beta22ththe = -1.5*(3.*he - 10.*le)*np.cos(2.*colat) + 4.5*(he - 2.*le)

    ththe = constB*(-6.*e*beta20ththe*np.cos(constA) + e*beta22ththe*(4.*np.sin(2.*lon)*np.sin(constA)+3. *
                                                                      np.cos(2.*lon)*np.cos(constA)) + 4.*np.cos(oblq)*np.sin(oblq)*beta21ththe*np.cos(lon)*np.sin(phase + constA))

    return ththe


def stressphphe(constA, constB, e, oblq, phase, colat, lon, he, le):

    beta20phphe = 0.75*(3.*he - 8.*le)*np.cos(2.*colat) + 0.75*(he - 4.*le)
    beta21phphe = 1.5*(3.*he - 8.*le)*np.sin(2.*colat)
    beta22phphe = -1.5*(3.*he - 8.*le)*np.cos(2.*colat) + 4.5*(he - 4.*le)

    phphe = constB*(-6.*e*beta20phphe*np.cos(constA) + e*beta22phphe*(4.*np.sin(2.*lon)*np.sin(constA)+3. *
                                                                      np.cos(2.*lon)*np.cos(constA)) + 4.*np.cos(oblq)*np.sin(oblq)*beta21phphe*np.cos(lon)*np.sin(phase + constA))

    return phphe


def stressthphe(constA, constB, e, oblq, phase, colat, lon, he, le):

    beta21thphe = 3.*le*np.sin(colat)
    beta22thphe = 3.*le*np.cos(colat)

    thphe = constB*(2.*e*beta22thphe*(4.*np.cos(2.*lon)*np.sin(constA)-3.*np.sin(2.*lon)*np.cos(
        constA)) + 4.*np.cos(oblq)*np.sin(oblq)*beta21thphe*np.sin(lon)*np.sin(phase + constA))

    return thphe


def modeththv(n, t, _lambda, sj, e, oblq, phase, colat, lon, hv, lv):

    beta20ththv = 0.75*(3.*hv - 10.*lv)*np.cos(2.*colat) + 0.75*(hv - 2.*lv)
    beta21ththv = 1.5*(3.*hv - 10.*lv)*np.sin(2.*colat)
    beta22ththv = -1.5*(3.*hv - 10.*lv)*np.cos(2.*colat) + 4.5*(hv - 2.*lv)

    ththv = (1./np.sqrt(1. + (-n/sj)*(-n/sj)))*(-6.*e*beta20ththv*np.cos(n*t - np.arctan(-n/sj) + np.arctan(_lambda)) + e*beta22ththv*(4.*np.sin(2.*lon)*np.sin(n*t - np.arctan(-n/sj) + np.arctan(_lambda)
                                                                                                                                                                )+3.*np.cos(2.*lon)*np.cos(n*t - np.arctan(-n/sj) + np.arctan(_lambda))) + 4.*np.cos(oblq)*np.sin(oblq)*beta21ththv*np.cos(lon)*np.sin(phase + n*t - np.arctan(-n/sj) + np.arctan(_lambda)))

    return ththv


def modephphv(n, t, _lambda, sj, e, oblq, phase, colat, lon, hv, lv):

    beta20phphv = 0.75*(3.*hv - 8.*lv)*np.cos(2.*colat) + 0.75*(hv - 4.*lv)
    beta21phphv = 1.5*(3.*hv - 8.*lv)*np.sin(2.*colat)
    beta22phphv = -1.5*(3.*hv - 8.*lv)*np.cos(2.*colat) + 4.5*(hv - 4.*lv)

    phphv = (1./np.sqrt(1. + (-n/sj)*(-n/sj)))*(-6.*e*beta20phphv*np.cos(n*t - np.arctan(-n/sj) + np.arctan(_lambda)) + e*beta22phphv*(4.*np.sin(2.*lon)*np.sin(n*t - np.arctan(-n/sj) + np.arctan(_lambda)
                                                                                                                                                                )+3.*np.cos(2.*lon)*np.cos(n*t - np.arctan(-n/sj) + np.arctan(_lambda))) + 4.*np.cos(oblq)*np.sin(oblq)*beta21phphv*np.cos(lon)*np.sin(phase + n*t - np.arctan(-n/sj) + np.arctan(_lambda)))

    return phphv


def modethphv(n, t, _lambda, sj, e, oblq, phase, colat, lon, hv, lv):

    beta21thphv = 3.*lv*np.sin(colat)
    beta22thphv = 3.*lv*np.cos(colat)

    thphv = (1./np.sqrt(1. + (-n/sj)*(-n/sj)))*(8.*e*beta22thphv*np.cos(2.*lon)*np.sin(n*t - np.arctan(-n/sj) + np.arctan(_lambda)) - 6.*e*beta22thphv*np.sin(2.*lon) *
                                                np.cos(n*t - np.arctan(-n/sj) + np.arctan(_lambda)) + 4.*np.cos(oblq)*np.sin(oblq)*beta21thphv*np.sin(lon)*np.sin(phase + n*t - np.arctan(-n/sj) + np.arctan(_lambda)))

    return thphv


def stressththeNSR(constA_NSR, constB_NSR, colat, lon, he_NSR, le_NSR):
    alpha22ththe = -1.5*(3.*he_NSR - 10.*le_NSR) * \
        np.cos(2.*colat) + 4.5*(he_NSR - 2.*le_NSR)
    ththeNSR = constB_NSR*alpha22ththe*np.cos(2.*lon + constA_NSR)

    return ththeNSR


def stressphpheNSR(constA_NSR, constB_NSR, colat, lon, he_NSR, le_NSR):

    alpha22phphe = -1.5*(3.*he_NSR - 8.*le_NSR) * \
        np.cos(2.*colat) + 4.5*(he_NSR - 4.*le_NSR)
    phpheNSR = constB_NSR*alpha22phphe*np.cos(2.*lon + constA_NSR)

    return phpheNSR


def stressthpheNSR(constA_NSR, constB_NSR, colat, lon, he_NSR, le_NSR):

    alpha22thphe = 3.*le_NSR*np.cos(colat)
    thpheNSR = -2.*constB_NSR*alpha22thphe*np.sin(2.*lon + constA_NSR)

    return thpheNSR


def modeththvNSR(constA_NSR, NSRrate, sj_NSR, colat, lon, hv_NSR, lv_NSR):

    alpha22ththv = -1.5*(3.*hv_NSR - 10.*lv_NSR) * \
        np.cos(2.*colat) + 4.5*(hv_NSR - 2.*lv_NSR)
    ththvNSR = (1./np.sqrt(1. + (-2.*NSRrate/sj_NSR)*(-2.*NSRrate/sj_NSR))) * \
        alpha22ththv*np.cos(2.*lon + constA_NSR -
                            np.arctan(-2.*NSRrate/sj_NSR))

    return ththvNSR


def modephphvNSR(constA_NSR, NSRrate, sj_NSR, colat, lon, hv_NSR, lv_NSR):

    alpha22phphv = -1.5*(3.*hv_NSR - 8.*lv_NSR) * \
        np.cos(2.*colat) + 4.5*(hv_NSR - 4.*lv_NSR)
    phphvNSR = (1./np.sqrt(1. + (-2.*NSRrate/sj_NSR)*(-2.*NSRrate/sj_NSR))) * \
        alpha22phphv*np.cos(2.*lon + constA_NSR -
                            np.arctan(-2.*NSRrate/sj_NSR))

    return phphvNSR


def modethphvNSR(constA_NSR, NSRrate, sj_NSR, colat, lon, hv_NSR, lv_NSR):

    alpha22thphv = 3.*lv_NSR*np.cos(colat)
    thphvNSR = -2.*(1./np.sqrt(1. + (-2.*NSRrate/sj_NSR)*(-2.*NSRrate/sj_NSR))) * \
        alpha22thphv*np.sin(2.*lon + constA_NSR -
                            np.arctan(-2.*NSRrate/sj_NSR))

    return thphvNSR


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

    constA = n*t + np.arctan(_lambda)
    constB = 0.5*((n * n * interior.radius * interior.rigidity) /
                  interior.surface_gravity) * (1. / np.sqrt(1. + _lambda * _lambda))

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

        constA_NSR = 2. * NSRrate * t + np.arctan(NSRdelta)
        #print("at constB NSR")
        constB_NSR = 0.5 * ((n * n * interior.radius * interior.rigidity) /
                            interior.surface_gravity) * (1. / np.sqrt(1. + NSRdelta * NSRdelta))

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
        zeta = 0.5*np.arctan((2.*myStressThPhTot) /
                             (myStressThThTot-myStressPhPhTot))

    sigTheta = myStressThThTot*(np.square(np.cos(zeta)))+myStressPhPhTot*(np.square(
        np.sin(zeta)))+myStressThPhTot*np.sin(2.*zeta)  # Corresponds to Terry's sigma 1
    sigPhi = myStressThThTot*(np.square(np.sin(zeta)))+myStressPhPhTot*(np.square(
        np.cos(zeta)))-myStressThPhTot*np.sin(2.*zeta)      # Corresponds to Terry's sigma 2

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
