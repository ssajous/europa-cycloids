////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                            //
// Europa cycloid                                                                                             //
//                                                                                                            //
// Alyssa Sarid, Burkhard Militzer                                UC Berkeley, 07-08-08                       //
//                                                                                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Standard.h"
#include "CycloidParameters.h"
#include "Array.h"
#include "StressCalcModuleVisco.h"

// This module generates a cycloid using the passed parameters.

void MoveCrack(const double colat, const double lon, const double deltaT, 
	       const double stress, const double heading, const double speed,
	       double & newColat, double & newLon) {  

  // deltaT = 1 for a full step
  // deltaT < 1 is a fraction of a step

  double distance = deltaT*speed/1562.0; // Hermes' radius value
  newColat = acos( (cos(colat)*cos(distance)) + (sin(colat)*sin(distance)*cos(heading)) );
  double lonNumer = cos(distance)-cos(newColat)*cos(colat);
  double lonDenom = sin(newColat)*sin(colat);

  double tempLon = lonNumer / lonDenom;
  if (tempLon > +1.0) tempLon = +1.0; // no change in lon below
  if (tempLon < -1.0) tempLon = -1.0; 

  /*
  double tempLon;
  if (abs(abs(lonNumer)-abs(lonDenom)) < 1E-15) {

    if (lonNumer/lonDenom < 0) {
      tempLon = -1.0;
    } else {
      tempLon = +1.0;
    }

  } else {
    tempLon = (lonNumer)/(lonDenom);
  }
  */

  double temp2 = acos(tempLon);
    
  if (heading <= pi) {
    if (lon-temp2 < 0) {
      newLon = lon-temp2 + 2.0*pi;
    } else {
      newLon = lon-temp2;
    }
  } else {
    if (lon+temp2 >= 2.0*pi) {
      newLon = lon+temp2 - 2.0*pi;
    } else {
      newLon = lon+temp2;
    }
  }

  //  Write5(speed,colat,lon,newColat,newLon);

}

// rounding disabled
inline double MyRound(const double x, const int n) {
  return x;
  //  double f = pow(10.0,-n);
  //  return f*round(x/f);
}

void MakeCrack(const CycloidParameters & cp, const double lengthLimit, 
	       Array1 <double> & lonDegs, Array1 <double> & latDegs,  
	       bool & startStressTooLow, bool & crackTooShort, const bool print) {  
  //double eccAmp = Radians(cp.eccAmpDeg);
  //double eccPhase = Radians(cp.eccPhaseDeg);
  double colat = Radians(90.-cp.latDeg);
  double lon = Radians(cp.lonDeg);
  double oblq = Radians(cp.oblqDeg);
  double phase = Radians(cp.phaseDeg);
  double NSRdelta = cp.NSRdelta;
  double libAmp = Radians(cp.libAmpDeg);
  double libPhase = Radians(cp.libPhaseDeg);
  double deltaPhase = Radians(cp.deltaPhaseDeg);

  lonDegs.PushBack(cp.lonDeg);
  latDegs.PushBack(cp.latDeg);

  crackTooShort = false;
  startStressTooLow = false;

  double totLength = 0.0;
  int iArcs = 0; // ***** added by Alyssa -- will be compared to real number of arcs inputted somewhere else
  bool firstTimeMoving = true;
  int notMoving = 0; // counts the steps where has not propagated recently

  double thisStep = cp.firstStep;
    
  // Calculates the largest tensile stress and the largest
  // associated principle stress direction (mainheading is always 
  // 180-360).

  //  cout << endl << endl;
  //  Write3(thisStep,iArcs,cp.nArcs);

  double mainHeading;
  double stress = GetStress(colat, lon, cp.steps, thisStep, oblq, phase, NSRdelta, libAmp, libPhase, mainHeading);
  //  Write(mainHeading);
  double startStress = stress;
  //  double threshold = startStress*cp.threshPercent;
  //  double threshold = 0.5;
  //  warning("DD");

  // To enter the loop. the crack must be able to begin propagating; 
  // the stress has to be higher than the cracking threshold.

  if (startStress < cp.threshold ) {
    startStressTooLow=true; 
    //    warning("Stress too low");
  } else {

    // The first restriction keeps the crack from propagating
    // farther than the real crack; the second checks whether
    // the crack will ever reinitiate after stopping at a cusp.

    double oldHeading;
    // ***** while (totLength <= lengthLimit && notMoving <= cp.steps) { // loops over the whole time
    while (iArcs < cp.nArcs && totLength <= 5*lengthLimit && notMoving <= cp.steps) { // loops over the whole time  

      //      Write2(iArcs,cp.nArcs);

     // Arc formation loop. 
      double deltaT = 1.0;
      // ***** while (totLength <= lengthLimit && stress >= cp.threshold && deltaT==1.0) { // only move if stress > threshold!!!
      while (iArcs < cp.nArcs && totLength <= 5.*lengthLimit && stress >= cp.threshold && deltaT==1.0) { // only move if stress > threshold!!!

	// Forces the crack to begin or continue propagating
	// in the specified direction at various times.
	// 1 means E2W propagation, 2 means W2E.
            
	double newColat, newLon;

	//	if (print) write3(thisStep,mainHeading,oldHeading);

	if (firstTimeMoving) { // start very first time only

	  // When the crack in the FIRST arc begins propagating...
	  if (cp.dir) { // formerly == 1 --> east-to-west
	    oldHeading = mainHeading;
	  } else {
	    oldHeading = mainHeading-pi;
	  }
	  firstTimeMoving = false;

	} else if (notMoving > 0) { // restart after cusp

	  // When crack reinitiates after a cusp...
	  if (oldHeading < pi) {
	    oldHeading = mainHeading;
	  } else {
	    oldHeading = mainHeading-pi;
	  }

	} else { // you are in an arc
	  
	  // As the crack propagates...
	  if (abs(oldHeading-mainHeading) > pi/2.0 && abs(oldHeading-mainHeading)< 3.0*pi/2.0) {
	    oldHeading = mainHeading-pi;
	  } else {
	    oldHeading = mainHeading;
	  }

	}

	//	if (print) Write(oldHeading);

	// Calculates lat/lon of next point.
	MoveCrack(colat, lon, deltaT, stress, oldHeading, cp.speed, newColat, newLon);
        
	double newStress = GetStress(newColat, newLon, cp.steps, thisStep+deltaT, oblq, phase+(deltaT*deltaPhase/cp.steps), NSRdelta, libAmp, libPhase, mainHeading);
	if (newStress < cp.threshold) { // partial step
	  deltaT = (stress-cp.threshold) / (stress - newStress); // how much of partial step?
	  if (deltaT<0.0 || deltaT>1.0) error("deltaT problem",deltaT);
	  MoveCrack(colat, lon, deltaT, stress, oldHeading, cp.speed, newColat, newLon);
	}

	//	double latDeg   = 90.0 - Degrees(colat);
	double colatDeg = Degrees(colat);
	double lonDeg   = Degrees(lon);
	double newLatDeg   = 90. - Degrees(newColat);
	double newColatDeg = Degrees(newColat);
	double newLonDeg   = Degrees(newLon);

	// Adds to list of lat/lons for generated cycloid.
	lonDegs.PushBack(newLonDeg);
	latDegs.PushBack(newLatDeg);
	
	//	Write3(thisStep,newLonDeg,newLatDeg);
        
	if (abs(newLonDeg-lonDeg) > 180.0) {
	  if (newLonDeg > lonDeg) {
	    lonDeg += 360.0;
	  } else {
	    newLonDeg += 360.0;
	  }
	}
	
	// Calculates new arc length as new lat/lon points
	// are determined.
	double segLength = sqrt(sqr(newLonDeg-lonDeg)+sqr(newColatDeg-colatDeg));
	totLength += segLength;
        
	colat = newColat;
	lon   = newLon;
	
	// Increments time step in orbit.
	thisStep += deltaT;
	phase += deltaT*deltaPhase/cp.steps;
	if (thisStep >= double(cp.steps)) {
	  thisStep -= double(cp.steps);
	  }
	if (phase >= 360.0) {
	  phase -= 360.0;
	  }
	
	// Gets stress and heading for new lat/lon and time step.
	stress = GetStress(colat, lon, cp.steps, thisStep, oblq, phase, NSRdelta, libAmp, libPhase, mainHeading);
	notMoving = 0;
	
      } // ... (totLength <= lengthLimit && stress >= threshold) 

      deltaT = 1.0;
      // Cusp formation loop -- sit and crack is NOT moving
	  
      iArcs++; // *****
      //      Write3(thisStep,iArcs,cp.nArcs);
	  
      while (iArcs < cp.nArcs && totLength <= 5*lengthLimit && stress < startStress && notMoving <= cp.steps && deltaT==1.0) { // ***** 
	// while (totLength <= lengthLimit && stress < startStress && notMoving <= cp.steps && deltaT==1.0) { // (notMoving <= cp.steps) have we waited for an entire orbit?
	// Increments time and checks stress.  If stress >
	// threshold, the crack reinitiates and enters the
	// arc formation loop until it reaches the length
	// limit.  If the time step loops through an entire
	// orbit, it is assumed that the crack never 
	// reinitiates.
	
	//	if (print) write3(thisStep,mainHeading,oldHeading);

	double oldStress = stress;
	stress = GetStress(colat, lon, cp.steps, thisStep+deltaT, oblq, phase+(deltaT*deltaPhase/cp.steps), NSRdelta, libAmp, libPhase, mainHeading);
	if (stress >= startStress) { // partial wait step
	  //	  double oldStress = GetStress(colat, lon, cp.steps, thisStep, oblq, phase, NSRdelta, mainHeading);
	  deltaT = (oldStress-startStress) / (oldStress - stress); // how much of partial step?
	  if (deltaT<0.0 || deltaT>1.0) error("deltaT problem",deltaT);
	}

	if (abs(oldHeading-mainHeading) > pi/2.0 && abs(oldHeading-mainHeading) < (3.0/2.0*pi)) {
	  oldHeading = mainHeading-pi;
	} else {
	  oldHeading = mainHeading;
	}

	thisStep += deltaT;
	phase += deltaT*deltaPhase/cp.steps;
	if (thisStep >= double(cp.steps)) {
	  thisStep -= double(cp.steps);
	  }
	if (phase >= 360.0) {
	  phase -= 360.0;
	  }
	

	//	if (print) Write2(oldHeading,"NOT_MOVING_PART");

	notMoving += 1;
      } // ... while (totLength <= lengthLimit && stress < startStress && notMoving <= steps)


    } // ... while (totLength <= lengthLimit && stress >= threshold) {
  } // else ... (startStress < threshold )

  if (notMoving > cp.steps) {
    crackTooShort = true;
    //    warning("too short");
    //    cout << cp << endl;
  }

//    Write(lonDegs);
//    Write(latDegs);
//	Quit("In MakeCrack");

} // BEFORE    return(points, start_stress, start_stress_too_low, crackTooShort)
