////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                            //
// Viscoelastic cycloid modeling (adapted from Europa cycloid, with Burkhard Militzer)                        //
// For use with FracFits optimization code                                                                    //
// Alyssa Rhoden                                NASA Goddard 3/2013                                           //
//                                                                                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function that calculates angular distance from the bulge location
// to the current lon/lat, measured in radians west of the bulge.

#include "Standard.h"

inline double stressththe(const double constA, const double constB, const double e, const double oblq, const double phase, const double colat, const double lon, const double he, const double le){
    
    double beta20ththe = 0.75*(3.*he - 10.*le)*cos(2.*colat) + 0.75*(he - 2.*le);
    double beta21ththe = 1.5*(3.*he - 10.*le)*sin(2.*colat);
    double beta22ththe = -1.5*(3.*he - 10.*le)*cos(2.*colat) + 4.5*(he - 2.*le);
    
    double ththe = constB*(-6.*e*beta20ththe*cos(constA) + e*beta22ththe*(4.*sin(2.*lon)*sin(constA)+3.*cos(2.*lon)*cos(constA)) + 4.*cos(oblq)*sin(oblq)*beta21ththe*cos(lon)*sin(phase + constA));
    
    return ththe;
}

inline double stressphphe(const double constA, const double constB, const double e, const double oblq, const double phase, const double colat, const double lon, const double he, const double le){
    
    double beta20phphe = 0.75*(3.*he - 8.*le)*cos(2.*colat) + 0.75*(he - 4.*le);
    double beta21phphe = 1.5*(3.*he - 8.*le)*sin(2.*colat);
    double beta22phphe = -1.5*(3.*he - 8.*le)*cos(2.*colat) + 4.5*(he - 4.*le);
    
    double phphe = constB*(-6.*e*beta20phphe*cos(constA) + e*beta22phphe*(4.*sin(2.*lon)*sin(constA)+3.*cos(2.*lon)*cos(constA)) + 4.*cos(oblq)*sin(oblq)*beta21phphe*cos(lon)*sin(phase + constA));
    
    return phphe;
}

inline double stressthphe(const double constA, const double constB, const double e, const double oblq, const double phase, const double colat, const double lon, const double he, const double le){
    
    double beta21thphe = 3.*le*sin(colat);
    double beta22thphe = 3.*le*cos(colat);
    
    double thphe = constB*(2.*e*beta22thphe*(4.*cos(2.*lon)*sin(constA)-3.*sin(2.*lon)*cos(constA)) + 4.*cos(oblq)*sin(oblq)*beta21thphe*sin(lon)*sin(phase + constA));
    
    return thphe;
}

inline double modeththv(const double n, const double t, const double lambda, const double sj, const double e, const double oblq, const double phase, const double colat, const double lon, const double hv, const double lv){
    
    double beta20ththv = 0.75*(3.*hv - 10.*lv)*cos(2.*colat) + 0.75*(hv - 2.*lv);
    double beta21ththv = 1.5*(3.*hv - 10.*lv)*sin(2.*colat);
    double beta22ththv = -1.5*(3.*hv - 10.*lv)*cos(2.*colat) + 4.5*(hv - 2.*lv);
    
    double ththv = (1./sqrt(1. + (-n/sj)*(-n/sj)))*(-6.*e*beta20ththv*cos(n*t - atan(-n/sj) + atan(lambda)) + e*beta22ththv*(4.*sin(2.*lon)*sin(n*t - atan(-n/sj) + atan(lambda))+3.*cos(2.*lon)*cos(n*t - atan(-n/sj) + atan(lambda))) + 4.*cos(oblq)*sin(oblq)*beta21ththv*cos(lon)*sin(phase + n*t - atan(-n/sj) + atan(lambda)));
    
    return ththv;
}

inline double modephphv(const double n,const double t,const double lambda,const double sj,const double e,const double oblq,const double phase,const double colat, const double lon, const double hv,const double lv){
    
    double beta20phphv = 0.75*(3.*hv - 8.*lv)*cos(2.*colat) + 0.75*(hv - 4.*lv);
    double beta21phphv = 1.5*(3.*hv - 8.*lv)*sin(2.*colat);
    double beta22phphv = -1.5*(3.*hv - 8.*lv)*cos(2.*colat) + 4.5*(hv - 4.*lv);
    
    double phphv = (1./sqrt(1. + (-n/sj)*(-n/sj)))*(-6.*e*beta20phphv*cos(n*t - atan(-n/sj) + atan(lambda)) + e*beta22phphv*(4.*sin(2.*lon)*sin(n*t - atan(-n/sj) + atan(lambda))+3.*cos(2.*lon)*cos(n*t - atan(-n/sj) + atan(lambda))) + 4.*cos(oblq)*sin(oblq)*beta21phphv*cos(lon)*sin(phase + n*t - atan(-n/sj) + atan(lambda)));
    
    return phphv;
}

inline double modethphv(const double n,const double t,const double lambda,const double sj,const double e,const double oblq,const double phase,const double colat, const double lon,const double hv,const double lv){
    
    double beta21thphv = 3.*lv*sin(colat);
    double beta22thphv = 3.*lv*cos(colat);
    
    double thphv = (1./sqrt(1. + (-n/sj)*(-n/sj)))*(8.*e*beta22thphv*cos(2.*lon)*sin(n*t - atan(-n/sj) + atan(lambda)) - 6.*e*beta22thphv*sin(2.*lon)*cos(n*t - atan(-n/sj) + atan(lambda)) + 4.*cos(oblq)*sin(oblq)*beta21thphv*sin(lon)*sin(phase + n*t - atan(-n/sj) + atan(lambda)));
    
    return thphv;
}

inline double stressththeNSR(const double constA_NSR, const double constB_NSR, const double colat, const double lon, const double he_NSR, const double le_NSR){
    
    double alpha22ththe = -1.5*(3.*he_NSR - 10.*le_NSR)*cos(2.*colat) + 4.5*(he_NSR - 2.*le_NSR);
    double ththeNSR = constB_NSR*alpha22ththe*cos(2.*lon + constA_NSR);
    
    return ththeNSR;
}

inline double stressphpheNSR(const double constA_NSR, const double constB_NSR, const double colat, const double lon, const double he_NSR, const double le_NSR){
    
    double alpha22phphe = -1.5*(3.*he_NSR - 8.*le_NSR)*cos(2.*colat) + 4.5*(he_NSR - 4.*le_NSR);
    double phpheNSR = constB_NSR*alpha22phphe*cos(2.*lon + constA_NSR);
    
    return phpheNSR;
}

inline double stressthpheNSR(const double constA_NSR, const double constB_NSR, const double colat, const double lon, const double he_NSR, const double le_NSR){
    
    double alpha22thphe = 3.*le_NSR*cos(colat);
    double thpheNSR = -2.*constB_NSR*alpha22thphe*sin(2.*lon + constA_NSR); // changed cos to sin on 6/23/17
    
    return thpheNSR;
}

inline double modeththvNSR(const double constA_NSR, const double NSRrate, const double sj_NSR, const double colat, const double lon, const double hv_NSR, const double lv_NSR){
    
    double alpha22ththv = -1.5*(3.*hv_NSR - 10.*lv_NSR)*cos(2.*colat) + 4.5*(hv_NSR - 2.*lv_NSR);
    double ththvNSR = (1./sqrt(1. + (-2.*NSRrate/sj_NSR)*(-2.*NSRrate/sj_NSR)))*alpha22ththv*cos(2.*lon + constA_NSR - atan(-2.*NSRrate/sj_NSR));
    
    return ththvNSR;
}

inline double modephphvNSR(const double constA_NSR, const double NSRrate, const double sj_NSR, const double colat, const double lon, const double hv_NSR, const double lv_NSR){
    
    double alpha22phphv = -1.5*(3.*hv_NSR - 8.*lv_NSR)*cos(2.*colat) + 4.5*(hv_NSR - 4.*lv_NSR);
    double phphvNSR = (1./sqrt(1. + (-2.*NSRrate/sj_NSR)*(-2.*NSRrate/sj_NSR)))*alpha22phphv*cos(2.*lon + constA_NSR - atan(-2.*NSRrate/sj_NSR));
    
    return phphvNSR;
}

inline double modethphvNSR(const double constA_NSR, const double NSRrate, const double sj_NSR, const double colat, const double lon, const double hv_NSR, const double lv_NSR){
    
    double alpha22thphv = 3.*lv_NSR*cos(colat);
    double thphvNSR = -2.*(1./sqrt(1. + (-2.*NSRrate/sj_NSR)*(-2.*NSRrate/sj_NSR)))*alpha22thphv*sin(2.*lon + constA_NSR - atan(-2.*NSRrate/sj_NSR));
    
    return thphvNSR;
}

double GetStress(const double colat, const double lon, const int steps, const double thisStep, const double oblq, const double phase, const double NSRdelta, const double libAmp, const double libPhase, double & bigHeading) {

    int interior = 1; // THIS IS TERRIBLE. MUST FIX.
    
    double kyear = 1./3.156E10; // kyear in sec (inverted to ease convertion of s_j from inverse kyears to sec
    double year2sec = 3600.*24.*365.;  // conversion for NSR sj values
    
    //  // COMPUTING DIURNAL STRESSES
    
    double period, periodInSec, n, t, ecc, sgrav, radius;
    double lonEast;
    double visc, ridg, lambda;
    double he, le;
    double sj1, sj2, sj3, sj4, sj5, sj6;
    double hv1, hv2, hv3, hv4, hv5, hv6;
    double lv1, lv2, lv3, lv4, lv5, lv6;
    double he_NSR, le_NSR;
    double sj1_NSR, sj2_NSR, sj3_NSR, sj4_NSR;
    double hv1_NSR, hv2_NSR, hv3_NSR, hv4_NSR;
    double lv1_NSR, lv2_NSR, lv3_NSR, lv4_NSR;
    
    int numModes, numModes_NSR;
    
    double sum1,sum2,sum3;
    
    period = 0.; //days;
    periodInSec = 0;
    n = 0.; // [rad/sec]
    t = 0.;
    visc = 0.;          // surface, Pa*s
    ridg = 0.;           // surface, Pa
    lambda = ridg/(visc*n);
    sgrav = 0.0;           // m/s^2
    radius = 0.*1000.;   // m
    ecc = 0.0;
    numModes = 0;
    numModes_NSR = 0;
    
    he = 0.;
    le = 0.;
    
    sj1 = kyear*0; //(*M0*)
    sj2 = kyear*0; //(*M2*)
    sj3 = kyear*0; //(*C0*)
    sj4 = kyear*0;   //(*M3*)
    sj5 = kyear*0;   //(*T1*)
    sj6 = kyear*0;   //(*T2*)
    
    hv1 = 0;
    hv2 = 0;
    hv3 = 0;
    hv4 = 0;
    hv5 = 0;
    hv6 = 0;
    
    lv1 = 0;
    lv2 = 0;
    lv3 = 0;
    lv4 = 0;
    lv5 = 0;
    lv6 = 0;
    
    he_NSR = 0.;
    le_NSR = 0.;
    
    sj1_NSR = kyear*0; //(*M0*)
    sj2_NSR = kyear*0; //(*M2*)
    sj3_NSR = kyear*0; //(*C0*)
    sj4_NSR = kyear*0;   //(*M3*)
    
    hv1_NSR = 0;
    hv2_NSR = 0;
    hv3_NSR = 0;
    hv4_NSR = 0;
    
    lv1_NSR = 0;
    lv2_NSR = 0;
    lv3_NSR = 0;
    lv4_NSR = 0;

//  // COMPUTING DIURNAL STRESSES
    
    periodInSec = 306000.0; // sec;
    n = 2.*pi/periodInSec; // [rad/sec]
    t = (thisStep/steps)*periodInSec;
    visc = 1.E21;          // surface, Pa*s
    ridg = 3.487E9;           // surface, Pa    Hermes: 3.487E9     Thin shell: 3.52E9
    lambda = ridg/(visc*n);
    lonEast = (2.*pi)-lon;  // switching to east longitude
    sgrav = 1.313;           // m/s^2
    radius = 1562.0*1000.;   // m               Hermes: 1562.0      Thin shell: 1561.5
    ecc = 0.01;    //                  Hermes: 0.0094      Thin shell: 0.01

    
    if (interior == 1) {
        
        numModes = 6;
        numModes_NSR = 4;
        
        he = 1.15106;
        le = 0.308014;
        
        // Model 1 values, from Hermes =VISCO=
        sj1 = kyear*-1.08100E-1; //(*M0*)
        sj2 = kyear*-3.23708E-1; //(*M2*)
        sj3 = kyear*-4.49843E-1; //(*C0*)
        sj4 = kyear*-3.44856;   //(*M3*)
        sj5 = kyear*-1.82043E5;   //(*T1*)
        sj6 = kyear*-1.02489E6;   //(*T2*)
        
        hv1 = 3.69312E-2;
        hv2 = 1.38303E-3;
        hv3 = 5.17520E-2;
        hv4 = 7.16927E-1;
        hv5 = 7.19007E-5;
        hv6 = 8.49126E-2;
        
        lv1 = 1.00431E-2;
        lv2 = 2.35498E-4;
        lv3 = 1.40299E-2;
        lv4 = 1.94722E-1;
        lv5 = 6.00150E-3;
        lv6 = 2.26341E-2;
        
        // These s_j's almost match the ones in Hermes' paper, but they migt be specific to the *slightly* different values he gave me later via email.
        
//        he_NSR = 1.85155;
//        le_NSR = 4.95366E-1;
//         
//        sj1_NSR = year2sec*-3.09117E3; //(*M0*)
//        sj2_NSR = year2sec*-9.24992E3; //(*M2*)
//        sj3_NSR = year2sec*-5.49324E-3;   //(*T1*)
//        sj4_NSR = year2sec*-9.84029E-4;   //(*T2*)
//         
//        hv1_NSR = 1.80231E-3;
//        hv2_NSR = 3.60537E-2;
//        hv3_NSR = 1.16432E-6;
//        hv4_NSR = 1.53522E-1;
//         
//        lv1_NSR = 3.07295E-3;
//        lv2_NSR = 9.80448E-3;
//        lv3_NSR = 9.28203E-3;
//        lv4_NSR = 4.09211E-2;
        
        
        // These are the s_j values that Wade computed for this interior model:
        // NSR HERMES EUROPA Visc ductile 1e14, Shell 30km = 25km ductile, 5km brittle L1=600 km, L2=832 km, L3=100 km, L4=25 km, L5=5 km
        
        he_NSR = 1.90474;
        le_NSR = 0.509593;
        
        sj1_NSR = -1./(3.37991E6);
        sj2_NSR = -1./(1.1289E6);
        sj3_NSR = -1./(2.00708);
        sj4_NSR = -1./(0.359672);
        
        hv1_NSR = 0.0382668;
        hv2_NSR = 0.00204669;
        hv3_NSR = 0.000668372;
        hv4_NSR = 0.159271;
        
        lv1_NSR = 0.0104063;
        lv2_NSR = 0.000385665;
        lv3_NSR = 0.00968342;
        lv4_NSR = 0.0424533;
        
        //        } else if (interior == "elastic") {
    } else if (interior == 2) {
        
        numModes = 6;
        numModes_NSR = 0;
        
        he = 1.151065;
        le = 0.3080142;
        
        // Model 2 values, from Hermes =ELASTIC=
        sj1 = kyear*-4.44842E-4; //(*M0*)
        sj2 = kyear*-1.08230E-1; //(*M2*)
        sj3 = kyear*-4.49843E-1; //(*C0*)
        sj4 = kyear*-3.44846;   //(*M3*)
        sj5 = kyear*-1.82134E2;   //(*T1*)
        sj6 = kyear*-1.02492E3;   //(*T2*)
        
        hv1 = 4.31352E-3;
        hv2 = 3.43781E-2;
        hv3 = 5.15023E-2;
        hv4 = 7.17001E-1;
        hv5 = 7.10977E-5;
        hv6 = 8.47114E-2;
        
        lv1 = 1.60251E-1;
        lv2 = 9.34250E-3;
        lv3 = 1.39915E-2;
        lv4 = 1.94817E-1;
        lv5 = 5.93813E-3;
        lv6 = 2.25805E-2;
        
        //Write3(sj1,hv1,lv1);
        
        
    } else if (interior == 3) {
        
        numModes = 6;
        numModes_NSR = 4;
        
        he = 1.0796;
        le = 0.28983;
        
        // Europa - Case 3 (Visc ductile 1e13, Shell 25km = 20km ductile, 5km brittle; L1=600 km, L2=832 km, L3=105 km, L4=20 km, L5=5 km)
        sj1 = -1./(2.92024E11);
        sj2 = -1./(6.9823E10);
        sj3 = -1./(1.91402E10);
        sj4 = -1./(8.73645E9);
        sj5 = -1./(14441.4);
        sj6 = -1./(3035.89);
        
        hv1 = 0.0402921;
        hv2 = 0.058432;
        hv3 = 0.00227749;
        hv4 = 0.881648;
        hv5 = 0.00079317;
        hv6 = 0.0644347;
        
        lv1 = 0.0109549;
        lv2 = 0.01589;
        lv3 = 0.000535483;
        lv4 = 0.2394;
        lv5 = 0.00459703;
        lv6 = 0.0172331;
        
        he_NSR = 1.94845;
        le_NSR = 0.52298;
        
        sj1_NSR = -1./(3.37991E6);
        sj2_NSR = -1./(221530.);
        sj3_NSR = -1./(0.167146);
        sj4_NSR = -1./(0.0354193);
        
        hv1_NSR = 0.0399374;
        hv2_NSR = 0.0029286;
        hv3_NSR = 0.00144285;
        hv4_NSR = 0.131424;
        
        lv1_NSR = 0.0108584;
        lv2_NSR = 0.000591633;
        lv3_NSR = 0.00796846;
        lv4_NSR = 0.0351483;
        
    } else if (interior == 4) {
        
        numModes = 6;
        numModes_NSR = 4;
        
        he = 1.0796;
        le = 0.28983;
        
        // Europa - Case 4 (Visc ductile 1e14, Shell 25km = 20km ductile, 5km brittle; L1=600 km, L2=832 km, L3=105 km, L4=20 km, L5=5 km)
        sj1 = -1./(2.92608E11);
        sj2 = -1./(1.91094E11);
        sj3 = -1./(6.9823E10);
        sj4 = -1./(8.73645E9);
        sj5 = -1./(144325.);
        sj6 = -1./(30370.6);
        
        hv1 = 0.0470372;
        hv2 = 0.000363229;
        hv3 = 0.0587308;
        hv4 = 0.880277;
        hv5 = 0.0000815638;
        hv6 = 0.0657785;
        
        lv1 = 0.0127975;
        lv2 = 0.000163256;
        lv3 = 0.0159528;
        lv4 = 0.239157;
        lv5 = 0.00440269;
        lv6 = 0.0175925;
        
        he_NSR = 1.94845;
        le_NSR = 0.52298;
        
        sj1_NSR = -1./(3.38329E6);
        sj2_NSR = -1./(2.21173E6);
        sj3_NSR = -1./(1.67043);
        sj4_NSR = -1./(0.354329);
        
        hv1_NSR = 0.0399628;
        hv2_NSR = 0.000517585;
        hv3_NSR = 0.000144936;
        hv4_NSR = 0.133819;
        
        lv1_NSR = 0.0108727;
        lv2_NSR = 0.000195862;
        lv3_NSR = 0.00761448;
        lv4_NSR = 0.0357888;
        
    } else if (interior == 5) {
        
        numModes = 6;
        numModes_NSR = 4;
        
        he = 1.0796;
        le = 0.28983;
        
        // Europa - Case 5 (Visc ductile 1e15, Shell 25km = 20km ductile, 5km brittle; L1=600 km, L2=832 km, L3=105 km, L4=20 km, L5=5 km)
        sj1 = -1./(1.90596E12);
        sj2 = -1./(2.91733E11);
        sj3 = -1./(6.9823E10);
        sj4 = -1./(8.73645E9);
        sj5 = -1./(1.44381E6);
        sj6 = -1./(303519.);
        
        hv1 = 0.00475959;
        hv2 = 0.0380883;
        hv3 = 0.0586485;
        hv4 = 0.880203;
        hv5 = 0.000530642;
        hv6 = 0.0636394;
        
        lv1 = 0.00401839;
        lv2 = 0.0103509;
        lv3 = 0.0159355;
        lv4 = 0.239145;
        lv5 = 0.0045248;
        lv6 = 0.0170204;
        
        he_NSR = 1.94845;
        le_NSR = 0.52298;
        
        sj1_NSR = -1./(2.20597E7);
        sj2_NSR = -1./(3.37654E6);
        sj3_NSR = -1./(16.7107);
        sj4_NSR = -1./(3.54465);
        
        hv1_NSR = 0.00476059;
        hv2_NSR = 0.0377858;
        hv3_NSR = 0.000964149;
        hv4_NSR = 0.136254;
        
        lv1_NSR = 0.00402118;
        lv2_NSR = 0.0102686;
        lv3_NSR = 0.00783787;
        lv4_NSR = 0.03644;
        
    } else if (interior == 6) {
        
        numModes = 6;
        numModes_NSR = 4;
        
        he = 1.12665;
        le = 0.306313;
        
        // Europa - Case 6 (Visc ductile 1e13, Shell 5km = 4km ductile, 1km brittle; L1=600 km, L2=832 km, L3=125 km, L4=4 km, L5=1 km)
        sj1 = -1./(2.441096E12);
        sj2 = -1./(2.87973E11);
        sj3 = -1./(6.98927E10);
        sj4 = -1./(8.96607E9);
        sj5 = -1./(14369.5);
        sj6 = -1./(2902.59);
        
        hv1 = 0.00123524;
        hv2 = 0.0117981;
        hv3 = 0.0647433;
        hv4 = 0.959494;
        hv5 = 0.000758626;
        hv6 = 0.014568;
        
        lv1 = 0.00106205;
        lv2 = 0.00321539;
        lv3 = 0.0176442;
        lv4 = 0.261483;
        lv5 = 0.00106861;
        lv6 = 0.00395794;
        
        he_NSR = 2.13243;
        le_NSR = 0.579738;
        
        sj1_NSR = -1./(2.82533E7);
        sj2_NSR = -1./(3.33302E6);
        sj3_NSR = -1./(0.166314);
        sj4_NSR = -1./(0.0336619);
        
        hv1_NSR = 0.00123468;
        hv2_NSR = 0.011725;
        hv3_NSR = 0.00143917;
        hv4_NSR = 0.0323838;
        
        lv1_NSR = 0.00106246;
        lv2_NSR = 0.00319547;
        lv3_NSR = 0.00193662;
        lv4_NSR = 0.00879825;
        
    } else if (interior == 7) {
        
        numModes = 6;
        numModes_NSR = 4;
        
        he = 1.12665;
        le = 0.306313;
        
        // Europa - Case 7 (Visc ductile 1e14, Shell 5km = 4km ductile, 1km brittle; L1=600 km, L2=832 km, L3=125 km, L4=4 km, L5=1 km)
        sj1 = -1./(2.41296E13);
        sj2 = -1./(2.87973E11);
        sj3 = -1./(6.98927E10);
        sj4 = -1./(8.96607E9);
        sj5 = -1./(143607.);
        sj6 = -1./(29037.);
        
        hv1 = 0.00121084;
        hv2 = 0.0116411;
        hv3 = 0.0647423;
        hv4 = 0.959491;
        hv5 = 0.0000567636;
        hv6 = 0.0160851;
        
        lv1 = 0.00792541;
        lv2 = 0.00317263;
        lv3 = 0.0176439;
        lv4 = 0.261482;
        lv5 = 0.000877126;
        lv6 = 0.00437013;
        
        he_NSR = 2.13243;
        le_NSR = 0.579738;
        
        sj1_NSR = -1./(2.79278E8);
        sj2_NSR = -1./(3.33302E6);
        sj3_NSR = -1./(1.66212);
        sj4_NSR = -1./(0.336748);
        
        hv1_NSR = 0.00121079;
        hv2_NSR = 0.0115661;
        hv3_NSR = 0.000107574;
        hv4_NSR = 0.035255;
        
        lv1_NSR = 0.00792613;
        lv2_NSR = 0.00315219;
        lv3_NSR = 0.00157342;
        lv4_NSR = 0.00957833;
        
    } else if (interior == 8) {
        
        numModes = 6;
        numModes_NSR = 4;
        
        he = 1.12665;
        le = 0.306313;
        
        // Europa - Case 8 (Visc ductile 1e15, Shell 5km = 4km ductile, 1km brittle; L1=600 km, L2=832 km, L3=125 km, L4=4 km, L5=1 km)
        sj1 = -1./(2.16512E14);
        sj2 = -1./(2.87973E11);
        sj3 = -1./(6.98927E10);
        sj4 = -1./(8.96607E9);
        sj5 = -1./(1.43662E6);
        sj6 = -1./(290192.);
        
        hv1 = 0.00248642;
        hv2 = 0.0116274;
        hv3 = 0.0647422;
        hv4 = 0.959491;
        hv5 = 0.000499664;
        hv6 = 0.0137214;
        
        lv1 = 0.0692122;
        lv2 = 0.0031689;
        lv3 = 0.0176439;
        lv4 = 0.261482;
        lv5 = 0.000997851;
        lv6 = 0.00372794;
        
        he_NSR = 2.13243;
        le_NSR = 0.579738;
        
        sj1_NSR = -1./(2.50592E9);
        sj2_NSR = -1./(3.33302E6);
        sj3_NSR = -1./(16.6276);
        sj4_NSR = -1./(3.36542);
        
        hv1_NSR = 0.00248637;
        hv2_NSR = 0.0115523;
        hv3_NSR = 0.000947997;
        hv4_NSR = 0.0307656;
        
        lv1_NSR = 0.069213;
        lv2_NSR = 0.00314843;
        lv3_NSR = 0.00180264;
        lv4_NSR = 0.00835862;
        
    } else if (interior == 9) {
        
        numModes = 6;
        numModes_NSR = 4;
        
        he = 1.09287;
        le = 0.294329;
        
        // Europa - Case 9 (DIURNAL ThinLowViscLowerLayer  5km 1e13 Ductile, 15km 1e21 Brittle; L1=600 km, L2=832 km, L3=110 km, L4=5 km, L5=15 km)
        
        sj1 = -1./(1.2177E12);
        sj2 = -1./(3.01499E11);
        sj3 = -1./(6.98927E10);
        sj4 = -1./(8.90366E9);
        sj5 = -1./(3834.33);
        sj6 = -1./(2908.39);
        
        hv1 = 0.00582762;
        hv2 = 0.110008;
        hv3 = 0.0585987;
        hv4 = 0.877646;
        hv5 = 0.000782287;
        hv6 = 0.0153222;
        
        lv1 = 0.0020604;
        lv2 = 0.0297098;
        lv3 = 0.0158261;
        lv4 = 0.237032;
        lv5 = 0.00123089;
        lv6 = 0.00408376;
        
        he_NSR = 2.01125;
        le_NSR = 0.541582;
        
        sj1_NSR = -1./(1.40937E7);
        sj2_NSR = -1./(3.48609E6);
        sj3_NSR = -1./(0.0443788);
        sj4_NSR = -1./(0.0337628);
        
        hv1_NSR = 0.0058192;
        hv2_NSR = 0.102652;
        hv3_NSR = 0.00145207;
        hv4_NSR = 0.0384067;
        
        lv1_NSR = 0.00205788;
        lv2_NSR = 0.0277231;
        lv3_NSR = 0.0022563;
        lv4_NSR = 0.0102352;
        
    } else {
        cout << "Undefined interior";
        Quit();
    }

    double myStressThTh, myStressPhPh, myStressThPh;
    double myStressThThTot, myStressPhPhTot, myStressThPhTot;
    double ththe, ththv1, ththv2, ththv3, ththv4, ththv5, ththv6;
    double phphe, phphv1, phphv2, phphv3, phphv4, phphv5, phphv6;
    double thphe, thphv1, thphv2, thphv3, thphv4, thphv5, thphv6;
    
    ththv3 = 0;
    ththv4 = 0;
    ththv5 = 0;
    ththv6 = 0;
    phphv3 = 0;
    phphv4 = 0;
    phphv5 = 0;
    phphv6 = 0;
    thphv3 = 0;
    thphv4 = 0;
    thphv5 = 0;
    thphv6 = 0;
    
    double constA = n*t + atan(lambda);
    double constB = 0.5*((n*n*radius*ridg)/sgrav)*(1./sqrt(1. + lambda*lambda));
    
    ththe = stressththe(constA,constB,ecc,oblq,phase,colat,lonEast,he,le);
    phphe = stressphphe(constA,constB,ecc,oblq,phase,colat,lonEast,he,le);
    thphe = stressthphe(constA,constB,ecc,oblq,phase,colat,lonEast,he,le);
    
    if (numModes < 2 || numModes > 6){
        cout << "Number of modes is disallowed";
        Quit();
    }
    
    ththv1 = modeththv(n,t,lambda,sj1,ecc,oblq,phase,colat,lonEast,hv1,lv1);
    phphv1 = modephphv(n,t,lambda,sj1,ecc,oblq,phase,colat,lonEast,hv1,lv1);
    thphv1 = modethphv(n,t,lambda,sj1,ecc,oblq,phase,colat,lonEast,hv1,lv1);
        
    ththv2 = modeththv(n,t,lambda,sj2,ecc,oblq,phase,colat,lonEast,hv2,lv2);
    phphv2 = modephphv(n,t,lambda,sj2,ecc,oblq,phase,colat,lonEast,hv2,lv2);
    thphv2 = modethphv(n,t,lambda,sj2,ecc,oblq,phase,colat,lonEast,hv2,lv2);
    
    if (numModes > 2){
        
        ththv3 = modeththv(n,t,lambda,sj3,ecc,oblq,phase,colat,lonEast,hv3,lv3);
        phphv3 = modephphv(n,t,lambda,sj3,ecc,oblq,phase,colat,lonEast,hv3,lv3);
        thphv3 = modethphv(n,t,lambda,sj3,ecc,oblq,phase,colat,lonEast,hv3,lv3);
        
    }
    
    if (numModes > 3){
        
        ththv4 = modeththv(n,t,lambda,sj4,ecc,oblq,phase,colat,lonEast,hv4,lv4);
        phphv4 = modephphv(n,t,lambda,sj4,ecc,oblq,phase,colat,lonEast,hv4,lv4);
        thphv4 = modethphv(n,t,lambda,sj4,ecc,oblq,phase,colat,lonEast,hv4,lv4);
        
    }
    
    if (numModes > 4){
        
        ththv5 = modeththv(n,t,lambda,sj5,ecc,oblq,phase,colat,lonEast,hv5,lv5);
        phphv5 = modephphv(n,t,lambda,sj5,ecc,oblq,phase,colat,lonEast,hv5,lv5);
        thphv5 = modethphv(n,t,lambda,sj5,ecc,oblq,phase,colat,lonEast,hv5,lv5);
        
    }
    
    if (numModes > 5){
        
        ththv6 = modeththv(n,t,lambda,sj6,ecc,oblq,phase,colat,lonEast,hv6,lv6);
        phphv6 = modephphv(n,t,lambda,sj6,ecc,oblq,phase,colat,lonEast,hv6,lv6);
        thphv6 = modethphv(n,t,lambda,sj6,ecc,oblq,phase,colat,lonEast,hv6,lv6);
        
    }
    
    myStressThTh = ththe + constB*(ththv1 + ththv2 + ththv3 + ththv4 + ththv5 + ththv6);
    myStressPhPh = phphe + constB*(phphv1 + phphv2 + phphv3 + phphv4 + phphv5 + phphv6);
    myStressThPh = thphe + constB*(thphv1 + thphv2 + thphv3 + thphv4 + thphv5 + thphv6); // negative sign comes later
    
    
    // COMPUTING NSR STRESSES ===== NOT UPDATED FROM EUROPA!!
    
    double myStressThThNSR, myStressPhPhNSR, myStressThPhNSR;
    double NSRrate, NSRperiod, constA_NSR, constB_NSR;
    
    if (NSRdelta != 0) {
        
        double ththv1NSR, ththv2NSR, ththv3NSR, ththv4NSR;
        double phphv1NSR, phphv2NSR, phphv3NSR, phphv4NSR;
        double thphv1NSR, thphv2NSR, thphv3NSR, thphv4NSR;
        
        ththv3NSR= 0;
        ththv4NSR= 0;
        phphv3NSR = 0;
        phphv4NSR = 0;
        thphv3NSR = 0;
        thphv4NSR = 0;
        
        // RIGIDITY AND VISCOSITY were swapped in original version, now corrected and consistent with diurnal eqs above
        
        NSRrate = (ridg/(2.*visc*NSRdelta)); // delta of 0.1 for NSR period of 11419 years; 43 roughly equals 6 Myr period // 6/23: removed *year2sec
        NSRperiod = 2.*pi/(NSRrate); // testing
        
        constA_NSR = 2.*NSRrate*t + atan(NSRdelta);
        constB_NSR = 0.5*((n*n*radius*ridg)/sgrav)*(1./sqrt(1. + NSRdelta*NSRdelta));
        
        //Write6(ridg, visc, NSRdelta, year2sec, NSRrate, NSRperiod); // testing
        
        // Currently sj values are from Hermes' paper for a reference model that is very similar to the visco model he sent us (M1)
        
        double ththeNSR = stressththeNSR(constA_NSR,constB_NSR,colat,lonEast,he_NSR,le_NSR);
        double phpheNSR = stressphpheNSR(constA_NSR,constB_NSR,colat,lonEast,he_NSR,le_NSR);
        double thpheNSR = stressthpheNSR(constA_NSR,constB_NSR,colat,lonEast,he_NSR,le_NSR);
        
        if (numModes_NSR < 2 || numModes_NSR > 4){
            cout << "Number of NSR modes is disallowed";
            Quit();
        }
        
        ththv1NSR = modeththvNSR(constA_NSR, NSRrate,sj1_NSR,colat,lonEast,hv1_NSR,lv1_NSR);
        phphv1NSR = modephphvNSR(constA_NSR, NSRrate,sj1_NSR,colat,lonEast,hv1_NSR,lv1_NSR);
        thphv1NSR = modethphvNSR(constA_NSR, NSRrate,sj1_NSR,colat,lonEast,hv1_NSR,lv1_NSR);
        
        ththv2NSR = modeththvNSR(constA_NSR, NSRrate,sj2_NSR,colat,lonEast,hv2_NSR,lv2_NSR);
        phphv2NSR = modephphvNSR(constA_NSR, NSRrate,sj2_NSR,colat,lonEast,hv2_NSR,lv2_NSR);
        thphv2NSR = modethphvNSR(constA_NSR, NSRrate,sj2_NSR,colat,lonEast,hv2_NSR,lv2_NSR);
        
        
        if (numModes_NSR > 2){
            
            ththv3NSR = modeththvNSR(constA_NSR, NSRrate,sj3_NSR,colat,lonEast,hv3_NSR,lv3_NSR);
            phphv3NSR = modephphvNSR(constA_NSR, NSRrate,sj3_NSR,colat,lonEast,hv3_NSR,lv3_NSR);
            thphv3NSR = modethphvNSR(constA_NSR, NSRrate,sj3_NSR,colat,lonEast,hv3_NSR,lv3_NSR);
            
        }
        
        if (numModes_NSR > 3){
            
            ththv4NSR = modeththvNSR(constA_NSR, NSRrate,sj4_NSR,colat,lonEast,hv4_NSR,lv4_NSR);
            phphv4NSR = modephphvNSR(constA_NSR, NSRrate,sj4_NSR,colat,lonEast,hv4_NSR,lv4_NSR);
            thphv4NSR = modethphvNSR(constA_NSR, NSRrate,sj4_NSR,colat,lonEast,hv4_NSR,lv4_NSR);
            
        }
        // Combining stresses
        
        myStressThThNSR = ththeNSR + constB_NSR*(ththv1NSR + ththv2NSR + ththv3NSR + ththv4NSR);
        myStressPhPhNSR = phpheNSR + constB_NSR*(phphv1NSR + phphv2NSR + phphv3NSR + phphv4NSR);
        myStressThPhNSR = thpheNSR + constB_NSR*(thphv1NSR + thphv2NSR + thphv3NSR + thphv4NSR);
        
        myStressThThTot = myStressThTh + myStressThThNSR;
        myStressPhPhTot = myStressPhPh + myStressPhPhNSR;
        myStressThPhTot = myStressThPh + myStressThPhNSR;  // 
        
        
    } else {
        myStressThThTot = myStressThTh;
        myStressPhPhTot = myStressPhPh;
        myStressThPhTot = myStressThPh;  //
    }

    
  // Converting COMBINED STRESSES to principle stresses.

    double zeta;
    if (myStressThThTot == myStressPhPhTot) {
        zeta = pi/2;
    } else {
        zeta = 0.5*atan(2.*myStressThPhTot/(myStressThThTot-myStressPhPhTot)); // Amount by which coordinates are being rotated to get principal stresses; changes CCW, orientation of sigTheta
    }
    
    double sigTheta   = myStressThThTot*(sqr(cos(zeta)))+myStressPhPhTot*(sqr(sin(zeta)))+myStressThPhTot*sin(2.*zeta);  // Corresponds to Terry's sigma 1
    double sigPhi = myStressThThTot*(sqr(sin(zeta)))+myStressPhPhTot*(sqr(cos(zeta)))-myStressThPhTot*sin(2.*zeta);      // Corresponds to Terry's sigma 2
    
    double sigThetaBars = sigTheta*1E-5;
    double sigPhiBars   = sigPhi*1E-5;
    
    //double sigThetaKPA = sigTheta*1E-3;
    //double sigPhiKPA   = sigPhi*1E-3;
    
    double zetaCW = (2.*pi)-zeta;   // Changes to CW, still oriented along sigTheta
    
    double stress;   // The largest (i.e. most tensile) principal stress
    double heading;  // Direction of crack propagation/formation. Perpendicular to direction of largest principal stress (i.e. max tension).
    double heading2;
    
    //Write3(sigThetaBars,sigPhiBars,zeta);
  
// Determines which stress is most tensile and largest heading.

    if (sigThetaBars < sigPhiBars) {
        stress = sigPhiBars;         // In this case, sigPhi is the max stress, so the heading should be perpendicular to sigPhi. SigTheta is -| to sigPhi by definition, so heading = zetaCW
        heading = zetaCW;
    } else {
        stress = sigThetaBars;
        heading = zetaCW +(pi/2.);  // In this case, sigTheta is the max stress, so the heading should be perpendicular to sigTheta. SigTheta is oriented along zeta, so heading = zetaCW + 90 deg
    }
    
    if (heading >= (2.*pi)) {       // Making sure azimuth values fall between 0 and 360. WOULD LOVE TO CHANGE THIS TO MOD AND INCORPORATE ABOVE.
        heading = heading -(2.*pi); // Also finds the two, complementary heading directions (e.g. 45 and 135)
        heading2 = heading +pi;
    } else {
        heading2 = heading -pi;
    }
    
    if (heading > heading2) {       // Determines which of the two heading directions (e.g. 45 or 135) is largest and selects that as the output heading.
        bigHeading = heading;
    } else {
        bigHeading = heading2;
    }
    
    return stress;
    
}
