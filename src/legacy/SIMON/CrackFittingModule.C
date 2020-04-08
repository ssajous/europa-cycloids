////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                            //
// Europa cycloid                                                                                             //
//                                                                                                            //
// Alyssa Sarid, Burkhard Militzer                                UC Berkeley, 07-08-08                       //
//                                                                                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include "Array.h"
#include "CycloidParameters.h"

inline double DiffernceLonDeg(const double lon1, const double lon2) { // make diff=lon1-lon2 between -180 and +180
  double diff = lon1 - lon2;
  diff -= 360.0*floor(diff/360.0+0.5);
  //  Write3(lon1,lon2,diff);
  return diff;
}

// lon(end) - lon(begin)           
double TotalLonDegDifference(const Array1 <double> & lonDegs) {
  if (lonDegs.Size()<2) error("LonDeg too short");
  double d=0.0;
  for(int i=1; i<lonDegs.Size(); i++) {
    d += DiffernceLonDeg(lonDegs[i],lonDegs[i-1]);
  }
  return d;
}

// west to east --> log is deceasing --> diff(end-begin)<0 --> dir=false
bool DetermineDirection(const Array1 <double> & lonDeg) {
  double d=TotalLonDegDifference(lonDeg);
  if (d==0.0) error("No diff in lonDeg");
  return d>0.0;
}

// Not anymore -- Careful: This is routine now returns the fractional arc length
Array1 <double> CalcArcLengths(const Array1 <double> & lonDeg, const Array1 <double> & latDeg) {
  int n = lonDeg.Size();
  Array1 <double> arcLength(n);
  
  arcLength[0]=0.0;

  double totLength = 0.0;

  for(int i=1; i<n; i++) {
        
    double lonPt   = lonDeg[i-1];
    double colatPt = 90.0-latDeg[i-1];
    double lonNextPt   = lonDeg[i];
    double colatNextPt = 90.0-latDeg[i];

    //    double lonDiff = abs(lonNextPt-lonPt);
    //    lonDiff -= floor(lonDiff / 360.0)*360.0; // makes lonDiff between 0 and 360
    double lonDiff = DiffernceLonDeg(lonNextPt,lonPt);
    //    if (i==10) Write3(lonDiff,lonNextPt,lonPt);
    //    warning("Check BB");

    double segLength = sqrt( lonDiff*lonDiff + sqr(colatNextPt-colatPt) );
    totLength += segLength;
    arcLength[i] = totLength;
  }
  //  Write3(lonDeg.Size(),latDeg.Size(),arcLength.Size());
  //  cout << arcLength;
  //  Write(arcLength.Last());

  //  Write2(arcLength,totLength);

  //  arcLength /= totLength;

  return arcLength;
}

double CalcChi2(const Array1 <double> & simLonDegs, const Array1 <double> & simLatDegs, 
		const Array1 <double> & dataLonDegs, const Array1 <double> & dataLatDegs, const Array1 <double> & dataArcLengths, 
		const double pixDiff, const double res, const bool print) {
  
  //  Write(DiffernceLonDeg(15,5));
  //  Write(DiffernceLonDeg(355,5));
  //  Write(DiffernceLonDeg(-5,365));
  //  Write(DiffernceLonDeg(-355,-345));
  //  Quit("Q");

  // Calculates estimated error in selecting the data points of 
  // the real cycloid.
    
  double radius = 1562.0;        //km
  double kmDiff = (pixDiff*res)/1000.0;    // so it's in km
  double latErrorDeg = (kmDiff/radius)*(180.0/pi);

  // Determines diffference in starting longitude between real 
  // and generated cycloids.

  //  Write3(simLonDegs[0],dataLonDegs[0],simLonDegs[0]-dataLonDegs[0]);
  double lonShift = dataLonDegs[0] - simLonDegs[0]; // add 'lonShift' to 'simLonDegs' to overlay with 'dataLonDeg'
  //  Write(lonShift);
  //  if (print) Write3(simLonDegs[0], dataLonDegs[0], lonShift);

  Array1 <double> dataLonDegsUnrolled(dataLonDegs.Size());
  dataLonDegsUnrolled[0] = dataLonDegs[0];
  //  Write3(dataLonDegsUnrolled[0],dataLonDegs[0],lonShift);

  for(int i=1; i<dataLonDegs.Size(); i++) {
    dataLonDegsUnrolled[i] = dataLonDegsUnrolled[i-1] + DiffernceLonDeg(dataLonDegs[i],dataLonDegs[i-1]);
  }
  //  Write(dataLonDegsUnrolled[1]);

  Array1 <double> simLonDegsUnrolled(simLonDegs.Size());
  simLonDegsUnrolled[0] = simLonDegs[0] + lonShift; //
  for(int i=1; i<simLonDegs.Size(); i++) {
    simLonDegsUnrolled[i] = simLonDegsUnrolled[i-1] + DiffernceLonDeg(simLonDegs[i],simLonDegs[i-1]);
    //    if (print) Write3(i,simLonDegsUnrolled[i],simLatDegs[i]);
  }

  Array1 <double> simArcLengths = CalcArcLengths(simLonDegsUnrolled,simLatDegs);
  double simTotLength = simArcLengths.Last();
  simArcLengths /= simTotLength;

  int i=0;
  double summedChi2 = 0.0;
  for(int j=1; j<dataArcLengths.Size(); j++) { // before 'each' is now  dataArcLengths[j]    do not use first point - chi2=0 because of shift

    double each = dataArcLengths[j] / dataArcLengths.Last();
    if (simArcLengths.Last() < each) error("Simulated cycloid too short",simArcLengths.Last(),each);

    //    Write(simArcLengths);

    while (simArcLengths[i] < each) {
      //      Write2(simArcLengths[i],each);

      // Finds first point along real cycloid where arc length
      // is larger than that of the generated cycloid - assures
      // that interpolated point is bounded.
      i += 1;
    }

    if (i>=simArcLengths.Size()) error("sim arc l loop error");

    double f = (each - simArcLengths[i-1]) / (simArcLengths[i]-simArcLengths[i-1]);
    double interpLon   = simLonDegsUnrolled[i-1]   + f * ( simLonDegsUnrolled[i] - simLonDegsUnrolled[i-1] );
    double interpLat   = simLatDegs[i-1]           + f * ( simLatDegs[i]         - simLatDegs[i-1]         );
    double interpColat = 90.0 - interpLat;
    
    double lonErrorDeg = kmDiff / (radius*cos(Radians(90.-dataLatDegs[j]))) * 180.0/pi; 
    
    double chi2Point = sqr((interpLon-dataLonDegsUnrolled[j])/lonErrorDeg) + sqr( (interpColat - (90.0-dataLatDegs[j])) /latErrorDeg );

    if (print) Write6(j,i,dataLonDegsUnrolled[j],interpLon,dataLatDegs[j],interpLat);

    summedChi2 += chi2Point;
  }			     
  
  double chi2 = summedChi2/double(dataArcLengths.Size());
    
  return chi2;
}

double CalcChi2Max(const CycloidParameters & cp) {
  double lMax = cp.steps * cp.speedMax * cp.nArcs;
  return 4.0/3.0*lMax*lMax*lMax;
}
