////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                            //
// Europa cycloid                                                                                             //
//                                                                                                            //
// Alyssa Sarid, Burkhard Militzer                                UC Berkeley, 07-08-08                       //
//                                                                                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _CYCLOIDPARAMETERS_
#define _CYCLOIDPARAMETERS_

#include "Array.h"

class CycloidParameters {
 public:
  static const int nPar = 10; // originally 10
  int nOpt;

  //static const bool eccAmpDegFlag   = true;
  //static const bool eccPhaseDegFlag  = true;
  //static const bool eccFlag  = true;
  static const bool oblqDegFlag   = true;
  static const bool phaseDegFlag  = true;
  static const bool thresholdFlag = true;
  static const bool speedFlag     = true;
  static const bool firstStepFlag = true;
  static const bool NSRdeltaFlag    = true;
  static const bool lonDegFlag    = true;
  static const bool libAmpDegFlag    = false;
  static const bool libPhaseDegFlag  = false;  
  static const bool deltaPhaseDegFlag = false;

  //double eccAmpDeg;        // ecc lib amp
  //double eccPhaseDeg;   // ecc phase lag
  //double eccVal;           // eccentricity
  double oblqDeg;       // obliquity
  double phaseDeg;      // direction of the spin pole
  double threshold;       // not fitted directly - minimum tensil stress for crack to keep propagating
  double speed;         // crack propagation speed = distance/step
  static const double speedMax;    // change in CycloidParameters.C
  double firstStep;        // point in orbit where we are starting, orbit has been divided into 850 steps, 'firststep' indicate step that we are going to start (integer)
  double NSRdelta;        // nonsyncronous rotation
  double lonDeg;         // Starting longitude
  double libAmpDeg;		// amplitude of libration
  double libPhaseDeg;	// phase of libration
  double deltaPhaseDeg; // precession rate

  Array1 < double* > v;

  double latDeg;         // Starting latitude (not a fitting parameter)
  int steps;            // number of steps in an orbit
  int nArcs;

  bool dir;             // true if east-to-west, false if west-to-east
  double pixDiff;       // accuracy of data ?
  double res;           // resolution in m/pixel
  double radius;        // Europa's radius (in km)

  void SetPointers() {
    nOpt=0;
//    if (eccAmpDegFlag) {
//      v.PushBack(&eccAmpDeg);
//      nOpt++;
//    }
//    if (eccPhaseDegFlag) {
//      v.PushBack(&eccPhaseDeg);
//      nOpt++;
//    }
//    if (eccFlag) {
//      v.PushBack(&eccVal);
//      nOpt++;
//    }
    if (oblqDegFlag) {
      v.PushBack(&oblqDeg);
      nOpt++;
    }
    if (phaseDegFlag) {
      v.PushBack(&phaseDeg);
      nOpt++;
    }
    if (thresholdFlag) {
      v.PushBack(&threshold);
      nOpt++;
    }
    if (speedFlag) {
      v.PushBack(&speed);
      nOpt++;
    }
    if (firstStepFlag) {
      v.PushBack(&firstStep);
      nOpt++;
    }
    if (NSRdeltaFlag) {
      v.PushBack(&NSRdelta);
      nOpt++;
    }
    if (lonDegFlag) {
      v.PushBack(&lonDeg);
      nOpt++;
    }
	if (libAmpDegFlag) {
      v.PushBack(&libAmpDeg);
      nOpt++;
    }
	if (libPhaseDegFlag) {
      v.PushBack(&libPhaseDeg);
      nOpt++;
    }
	if (deltaPhaseDegFlag) {
      v.PushBack(&deltaPhaseDeg);
      nOpt++;
    }
    if (nOpt==0) error("No parameter to optimize");
  }
  
  CycloidParameters() {
    SetPointers();
  }

  CycloidParameters(const CycloidParameters & cp) {
    //eccAmpDeg    = cp.eccAmpDeg;
    //eccPhaseDeg  = cp.eccPhaseDeg;
    //eccVal    = cp.eccVal;
    oblqDeg   = cp.oblqDeg;
    phaseDeg  = cp.phaseDeg;
    threshold = cp.threshold;
    speed     = cp.speed;
    firstStep = cp.firstStep;
    NSRdelta  = cp.NSRdelta;
    lonDeg    = cp.lonDeg;
	libAmpDeg    = cp.libAmpDeg;
	libPhaseDeg    = cp.libPhaseDeg;
	deltaPhaseDeg = cp.deltaPhaseDeg;

    latDeg    = cp.latDeg;
    steps     = cp.steps;
    nArcs     = cp.nArcs;
    dir       = cp.dir;
    pixDiff   = cp.pixDiff;
    res       = cp.res;
    radius    = cp.radius;
    SetPointers();
  }

  CycloidParameters & operator=(const CycloidParameters & cp) { 
    //eccAmpDeg    = cp.eccAmpDeg;
    //eccPhaseDeg  = cp.eccPhaseDeg;
    //eccVal    = cp.eccVal;
    oblqDeg   = cp.oblqDeg;
    phaseDeg  = cp.phaseDeg;
    threshold = cp.threshold;
    speed     = cp.speed;
    firstStep = cp.firstStep;
    NSRdelta    = cp.NSRdelta;
    lonDeg    = cp.lonDeg;
	libAmpDeg    = cp.libAmpDeg;
	libPhaseDeg    = cp.libPhaseDeg;
	deltaPhaseDeg = cp.deltaPhaseDeg;

    latDeg    = cp.latDeg;
    steps     = cp.steps;
    nArcs     = cp.nArcs;
    dir       = cp.dir;
    pixDiff   = cp.pixDiff;
    res       = cp.res;
    radius    = cp.radius;

    return *this;
  }

  string GetName(int i) {
    if (i<0 || i>nPar) error("Parameter invalid",i);
    //if (eccAmpDegFlag   && i--==0) return "eccAmpDeg";
    //if (eccPhaseDegFlag   && i--==0) return "eccPhaseDeg";
    //if (eccFlag   && i--==0) return "eccVal";
    if (oblqDegFlag   && i--==0) return "oblq";
    if (phaseDegFlag  && i--==0) return "phase";
    if (thresholdFlag && i--==0) return "threshold";
    if (speedFlag     && i--==0) return "speed";
    if (firstStepFlag && i--==0) return "first";
    if (NSRdeltaFlag    && i--==0) return "NSR";
    if (lonDegFlag    && i--==0) return "lon";
	if (libAmpDegFlag    && i--==0) return "libAmp";
	if (libPhaseDegFlag    && i--==0) return "libPhase";
	if (deltaPhaseDegFlag    && i--==0) return "deltaPhase";
    error("wrong i in GetName()",i);
    return " ";
  }

  // parameter range for scanning
  static void ParameterRange(const int i, double & pMin, double & pMax) {
    int ii=i;
    if (i<0 || i>nPar) error("Parameter invalid",i);
    //if (eccAmpDegFlag   && ii--==0) { pMin = 0.;  pMax =   2.0; } // eccAmpDeg
    //if (eccPhaseDegFlag  && ii--==0) { pMin = 0.0;  pMax = 360.0; } // eccPhaseDeg
    //if (eccFlag   && ii--==0) { pMin = 0.;  pMax =   0.1; } // ecc
    if (oblqDegFlag   && ii--==0) { pMin = 0.0;  pMax =   3.0; } // oblqDeg
    if (phaseDegFlag  && ii--==0) { pMin = 0.0;  pMax = 360.0; } // phaseDeg
    if (thresholdFlag && ii--==0) { pMin = 0.0; pMax =   10.0; } // threshold WARNING!!!
    if (speedFlag     && ii--==0) { pMin = 0.05; pMax = speedMax; } // speed, before was [0.1 ... 2]
    if (firstStepFlag && ii--==0) { pMin = 0.0;  pMax = 1700.0; } // first step in orbit
    if (NSRdeltaFlag    && ii--==0) { pMin = 10.0;  pMax = 1000.0; } // NSR delta
    if (lonDegFlag    && ii--==0) { pMin = 0.0;  pMax = 360.0; } // lonDeg
	if (libAmpDegFlag    && ii--==0) { pMin = 0.0;  pMax = 2.0; } // libration amplitude
	if (libPhaseDegFlag    && ii--==0) { pMin = 0.0;  pMax = 360.0; } // libration phase
	if (deltaPhaseDegFlag    && ii--==0) { pMin = 0.0;  pMax = 2.5; } // rate of change of the phase
  }
  
  // parameter range for scanning
  static double ParameterRange(const int i) { 
    double pMin, pMax;
    ParameterRange(i,pMin,pMax);
    return abs(pMax-pMin);
  }

  // is parameter valid? Careful: No restriction on 'lon' and 'phase' 
  static bool ParameterValid(const double p, const int i) {
    int ii=i;
    if (i<0 || i>nPar) error("Parameter invalid",i);
    //if (eccAmpDegFlag   && ii--==0) { return p>=0. && p<=2.0;  } // eccAmpDeg
    //if (eccPhaseDegFlag  && ii--==0) { return true; } // eccPhaseDeg
    //if (eccFlag   && ii--==0) { return p>=0. && p<=0.1;  } // ecc
    if (oblqDegFlag   && ii--==0) { return p>=0.0 && p<=3.0;  } // oblqDeg
    if (phaseDegFlag  && ii--==0) { return true;              } // phaseDeg, any value ok
    if (thresholdFlag && ii--==0) { return p>=0.0 && p<=10.0; } // threshold (it easily gets stuck at boundary 0, avoid it
    if (speedFlag     && ii--==0) { return p>=0.05 && p<=speedMax; } // speed, before was [0.1 ... 2]
    if (firstStepFlag && ii--==0) { return true;              } // first step can be anything
    if (NSRdeltaFlag    && ii--==0) { return p>=10.0 && p<=1000.0;  } // NSR delta
    if (lonDegFlag    && ii--==0) { return true;              } // lonDeg, any value ok
    if (libAmpDegFlag    && ii--==0) { return p>=0.0 && p<=2.0;  } // libration amplitude	
    if (libPhaseDegFlag    && ii--==0) { return true;              } // libration phase
	if (deltaPhaseDegFlag    && ii--==0) {return p>=0.0 && p<=2.5; } // rate of change of the phase
    error("Parameter index invalid",i);
    return true;
  }

  friend ostream& operator<<(ostream & os, const CycloidParameters & cp) {
    //os << ((cp.eccAmpDegFlag) ? " eccAmpDeg= "     : " <eccAmpDeg>= "    ) << cp.eccAmpDeg;
    //os << ((cp.eccPhaseDegFlag) ? " eccPhaseDeg= "     : " <eccPhaseDeg>= "    ) << cp.eccPhaseDeg;
    //os << ((cp.eccFlag) ? " ecc= "     : " <ecc>= "    ) << cp.eccVal;
    os << ((cp.firstStepFlag) ? " first= "     : " <first>= "    ) << cp.firstStep;
    os << ((cp.lonDegFlag)    ? " lon= "       : " <lon>= "      ) << cp.lonDeg;
    os << ((cp.oblqDegFlag)   ? " oblq= "      : " <oblq>= "     ) << cp.oblqDeg;
    os << ((cp.phaseDegFlag)  ? " phase= "     : " <phase>= "    ) << cp.phaseDeg;
    os << ((cp.NSRdeltaFlag)    ? " NSR= "       : " <NSR>= "      ) << cp.NSRdelta;
    os << ((cp.libAmpDegFlag)   ? " libAmp= "      : " <libAmp>= "     ) << cp.libAmpDeg;
    os << ((cp.libPhaseDegFlag)  ? " libPhase= "     : " <libPhase>= "    ) << cp.libPhaseDeg;
    os << ((cp.thresholdFlag) ? " threshold= " : " <threshold>= ") << cp.threshold;
    os << ((cp.speedFlag)     ? " speed= "     : " <speed>= "    ) << cp.speed;
	os << ((cp.deltaPhaseDegFlag) ? " deltaPhase= " : " <deltaPhase>= " ) << cp.deltaPhaseDeg;
    return os;
  }
    

};

#endif
