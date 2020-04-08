////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                            //
// Viscoelastic cycloid modeling (adapted from Europa cycloid, with Burkhard Militzer                         //
//                                                                                                            //
// Alyssa Rhoden                                NASA Goddard 3/2013                                           //
//                                                                                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "ReadInTable.h"
#include "CycloidParameters.h"
#include "StressCalcModuleVisco.h"
#include "MakeCrackModule.h"
#include "CrackFittingModule.h"
#include "BFGSOptimization.h"
#include "Random.h"

void ReadCycloidFile(const string & name, int & n, Array1 <double> & lonDeg, Array1 <double> & latDeg) {
  ReadInTable(name, 1, 2, n, lonDeg, latDeg);
  cout << n << " values read from file " << name << endl;
}

double GetChi2(const CycloidParameters & cpIn, 
	       const Array1<double> & dataLonDegs, const Array1<double> & dataLatDegs, const Array1<double> & dataArcLengths,
	       const double invalid, const bool print) {

  CycloidParameters cp(cpIn);
  cp.lonDeg   -= floor(cp.lonDeg  /360.0)*360.0;
  cp.phaseDeg -= floor(cp.phaseDeg/360.0)*360.0;
  cp.libPhaseDeg -= floor(cp.libPhaseDeg/360.0)*360.0;
  
  /*
  while (cp.phaseDeg<0.0) {
    warning("Phase mapping",cp.lonDeg, cp.phaseDeg);
    cp.phaseDeg += 180.0;
    cp.lonDeg   += 180.0;
  }
  */
/* **************************** THIS SECTION USED TO CHECK STRESS FIELD AGAINST TERRY'S DATA
  double notUsedMainHeading;
  
  cp.latDeg = 80.0; // overwritten below with data from file
  cp.pixDiff = 2.0;
  cp.res     = 1000;
  cp.radius  = 1561.5; 
  cp.steps     = 4.0; // *10 how many steps in orbit

  cp.lonDeg   = 0.0;  // 210.0;
  cp.oblqDeg  = 0.0; // 0.375;
  cp.phaseDeg = 0.0; // 135.0;

  cp.firstStep = 1.0; // *10 must change with steps 
  cp.threshold = 0.629684;
  cp.speed     = 0.669482; // *10 CAREFUL: distance per step , use 0.25 for 850 steps
  cp.NSRdelta    = 0.000;
  
  cp.libAmpDeg = 1.0;
  cp.libPhaseDeg = 270.0;

  Write3(cp.oblqDeg, cp.libAmpDeg, cp.libPhaseDeg);
  Write2(cp.latDeg, cp.lonDeg);
  double startStress = GetStress(Radians(90.-cp.latDeg), Radians(cp.lonDeg), cp.steps, cp.firstStep,
				 Radians(cp.oblqDeg), Radians(cp.phaseDeg), cp.NSRdelta, Radians(cp.libAmpDeg), Radians(cp.libPhaseDeg), 
				 notUsedMainHeading);
  
  for(int i=0; i<17; i++) {
		for(int j=0; j<36; j++) {
			Write2(cp.latDeg, cp.lonDeg);
			double startStress = GetStress(Radians(90.-cp.latDeg), Radians(cp.lonDeg), cp.steps, cp.firstStep,
				 Radians(cp.oblqDeg), Radians(cp.phaseDeg), cp.NSRdelta, Radians(cp.libAmpDeg), Radians(cp.libPhaseDeg), 
				 notUsedMainHeading);
			cp.lonDeg += 10.0;
	}
	cp.latDeg -= 10.0;
	cp.lonDeg = 0.0;
  }
  
  Quit("W");
****************************** */


  double notUsedMainHeading;
  double startStress = GetStress(Radians(90.-cp.latDeg), Radians(cp.lonDeg), cp.steps, cp.firstStep,
				 Radians(cp.oblqDeg), Radians(cp.phaseDeg), cp.NSRdelta, Radians(cp.libAmpDeg), Radians(cp.libPhaseDeg), 
				 notUsedMainHeading);  
  //  double lengthLimit=dataArcLengths.Last();
  double lengthLimit=2.0*dataArcLengths.Last(); // for testing purpose, allow longer cycloids

  //Write2(startStress, notUsedMainHeading);
  
  double chi2;
  if (startStress <= 0.0) {
    //    chi2 = 1E50; // starting stress is compressive
    chi2 = invalid;
  } else {
    //    double threshold = cp.threshPercent*startStress;
    //    threshold = 0.5;
    //    warning("CC",threshold);
    //    Write3(startStress,cp.threshPercent,threshold);
    //    Quit("WW");

    Array1 <double> simLonDegs, simLatDegs;
    bool startStressTooLow, crackTooShort;
    MakeCrack(cp, lengthLimit, simLonDegs, simLatDegs, startStressTooLow, crackTooShort, print); 
    //    Write2(startStressTooLow, crackTooShort);

	
    if (!startStressTooLow && !crackTooShort) {
      chi2 = CalcChi2(simLonDegs, simLatDegs, dataLonDegs, dataLatDegs, dataArcLengths, cp.pixDiff, cp.res, print);
    } else if (startStressTooLow) {
      //      chi2 = 1E50; // starting stress is below threshold
      chi2 = invalid;
    } else {
      //      chi2 = 1E20; // crack too short
      chi2 = invalid;
    }
    //    Quit("Q");
  }
 
//  cout.precision(10);
//  Write3(cp.latDeg, startStress, notUsedMainHeading);
//  Write(chi2);
//  Quit("W");


  return chi2;
}

double TunnelingScan(OptFunction & of, const int n, bool & found) {
  cout << "Before " << n << " point 1D tunneling scan: chi=" << of();

  int nF = of.GetNumber();
  Array1<double> grad(nF);
  of.GetGradient(grad);

  Array1<double> pChi2Min(nF);
  Array1<double> pOld(nF);
  of.GetParameters(pOld);

  double aMax = 0.0;
  for(int i=0; i<nF; i++) {
    double p = *(of.cp.v[i]);
    double pMin,pMax;
    of.cp.ParameterRange(i,pMin,pMax);
    double a;
    if (grad[i]==0.0) continue;
    if (grad[i]>0) { 
      a = (pMax-p)/grad[i];
    } else {
      a = (pMin-p)/grad[i];
    }
    if (aMax==0.0 || a<aMax) aMax = a;
  }
      

  double chi2Min=of();
  bool foundHere   = false;
  for(int i=1; i<=n; i++) {
    double a = aMax*double(i)/double(n);

    Array1<double> delta(grad);
    delta *= a;
    of.SetParameters(pOld);
    of.ApplyPerturbations(delta);
    double chi2i = of();
    //    Write3(i,a,chi2i);
    if (chi2i<chi2Min) {
      of.GetParameters(pChi2Min);
      chi2Min  = chi2i;
      found    = true;
      foundHere= true;
    }
  }
  if (foundHere) 
    of.SetParameters(pChi2Min);
  else 
    of.SetParameters(pOld);

  cout << ((foundHere) ? " +++" : " ") << endl;
  cout << "After  " << n << " point 1D tunneling scan: chi=" << of() 
       << ((foundHere) ? " +++" : " ")
       << endl;
  return chi2Min;

}

double OneDimensionalScan(OptFunction & of, const int ii, const int n, bool & found) {
  cout << "Before " << n << " point 1D scan of variable " << of.cp.GetName(ii) << "=" << *(of.cp.v[ii]) << " chi=" << of();

  double pMin,pMax;
  of.cp.ParameterRange(ii,pMin,pMax);

  double chi2Min=of();
  double pChi2Min= *(of.cp.v[ii]);

  Array1 <double> chi2(n+1),p(n+1);
  bool foundHere   = false;
  //  double di = Random()*0.0;
  double di = Random();
  for(int i=0; i<n; i++) {
    *(of.cp.v[ii]) = p[i] = pMin+(double(i)+di)/double(n)*(pMax-pMin);
    double chi2i = chi2[i] = of();
    if (chi2i<chi2Min) {
      pChi2Min = *(of.cp.v[ii]);
      chi2Min  = chi2i;
      found    = true;
      foundHere= true;
    }
  }
  cout << ((foundHere) ? " +++" : " ") << endl;

  /*
  //  to investigate problems with BFGS, uncomment and restart --> it will print the values from the scan
  if (chi2Min<1059.85) {
    int nn=20000;
    for(int i=0; i<nn; i++) {
      *(of.cp.v[ii]) = pMin+double(i)/double(nn)*(pMax-pMin);
      double chi2i = of();
      Write3(i,*(of.cp.v[ii]),chi2i);
    }
    Quit("Quit during 1D scan");
  } 
  */

  /*if (!foundBefore && foundHere && n>=200) 
    for(int i=0; i<n; i++) 
      Write3(i,p[i],chi2[i]);
  */

  //  Write3(ii,chi2Min,pChi2Min);
  *(of.cp.v[ii]) = pChi2Min;
  cout << "After  " << n << " point 1D scan of variable " << of.cp.GetName(ii) << "=" << *(of.cp.v[ii]) << " chi=" << of() 
       << ((foundHere) ? " +++" : " ")
       << endl;
  return chi2Min;
}

double MonteCarloSearch(OptFunction & of, const int n, bool & found) {
  //  cout << "Before " << n << " point 1D scan of variable v[" << ii+1 << "]=" << of.cp.v[ii] << " chi=" << of() << endl;
  cout << endl << "Before " << n << " point MC scan chi=" << of() << endl;

  CycloidParameters cpOrg = of.cp;
  double chi2Min=of();
  bool foundHere = false;

  for(int i=0; i<n; i++) { // try so many MC points
    for(int ii=0; ii<of.cp.nOpt; ii++) { // lopp over all parameters that will be changed during the optimization
      double range = 0.00005*of.cp.ParameterRange(ii); // 0.005%
      double delta = range*(Random()-0.5);
      *(of.cp.v[ii]) = *(cpOrg.v[ii]) + delta;
      if (!of.cp.ParameterValid(*(of.cp.v[ii]),ii)) 
	*(of.cp.v[ii]) = *(cpOrg.v[ii]);
    }
    //    Write(cpOrg);
    //    Write(of.cp);
    double chi2 = of();
    //    Write3(chi2,chi2Min,chi2-chi2Min);
    if (chi2<chi2Min) {
      cpOrg = of.cp;
      chi2Min = chi2;
      found    = true;
      foundHere= true;
    }
  }

  of.cp = cpOrg;
  cout << "After  " << n << " point MC scan chi=" << of()
       << ((foundHere) ? " +++" : " ")
       << endl;
  return chi2Min;
}

double TwoDimensionalScan(OptFunction & of, const int ii, const int jj, const int n, bool & found) {
  cout << "Before " << n*n << " point 2D scan of variable " << of.cp.GetName(ii) << "=" << *(of.cp.v[ii]) 
       << " and " << of.cp.GetName(jj) << "=" << *(of.cp.v[jj]) << " chi=" << of(); // << endl;
  
  if (ii==jj) error("ii==jj",ii,jj);

  double pMin1,pMax1;
  of.cp.ParameterRange(ii,pMin1,pMax1);
  double pMin2,pMax2;
  of.cp.ParameterRange(jj,pMin2,pMax2);

  double chi2Min=of();
  double pChi2Min1=*(of.cp.v[ii]);
  double pChi2Min2=*(of.cp.v[jj]);
  //  double di = Random()*0.0;
  //  double dj = Random()*0.0;
  double di = Random();
  double dj = Random();
  bool foundHere = false;

  for(int i=0; i<=n; i++) {
    *(of.cp.v[ii]) = pMin1+(double(i)+di)/double(n)*(pMax1-pMin1);
    for(int j=0; j<=n; j++) {
      *(of.cp.v[jj]) = pMin2+(double(j)+dj)/double(n)*(pMax2-pMin2);
      double chi2 = of();
      if (chi2<chi2Min) {
	pChi2Min1 = *(of.cp.v[ii]);
	pChi2Min2 = *(of.cp.v[jj]);
	chi2Min   = chi2;
	found     = true;
	foundHere = true;
      }
    }
  }
  cout << ((foundHere) ? " +++" : " ") << endl;
  *(of.cp.v[ii]) = pChi2Min1;
  *(of.cp.v[jj]) = pChi2Min2;
  cout << "After  " << n*n << " point 2D scan of variable " << of.cp.GetName(ii) << "=" << *(of.cp.v[ii]) 
       << " and " << of.cp.GetName(jj) << "=" << *(of.cp.v[jj]) << " chi=" << of() 
       << ((foundHere) ? " +++" : " ")
       << endl;
  return chi2Min;
}

double LocalTwoDimensionalScan(OptFunction & of, const int ii, const int jj, const int n, bool & found) {
  int p = cout.precision();
  cout.precision(p+2);
  cout << "Before " << n*n << " point local 2D scan of variable " << of.cp.GetName(ii) << "=" << *(of.cp.v[ii]) 
       << " and " << of.cp.GetName(jj) << "=" << *(of.cp.v[jj]) << " chi=" << of(); // << endl;
  
  if (ii==jj) error("ii==jj",ii,jj);
  int nn = n/2;

  CycloidParameters cpOrg = of.cp;
  double pRange1 = of.cp.ParameterRange(ii);
  double pRange2 = of.cp.ParameterRange(jj);

  double chi2Min=of();
  double pChi2Min1=*(of.cp.v[ii]);
  double pChi2Min2=*(of.cp.v[jj]);
  bool foundHere = false;
  double epsStart = 1e-2;

  for(int i=-nn; i<=nn; i++) {
    if (i!=0) {
      *(of.cp.v[ii]) = *(cpOrg.v[ii]) + pRange1*epsStart*pow(10.0,-abs(i)) * ((i<0) ? -1.0 : 1.0);
      for(int j=-nn; j<=nn; j++) {
	if (j!=0) {
	  *(of.cp.v[jj]) = *(cpOrg.v[jj]) + pRange2*epsStart*pow(10.0,-abs(j)) * ((j<0) ? -1.0 : 1.0);
	  double chi2 = of();
	  if (chi2<chi2Min) {
	    pChi2Min1 = *(of.cp.v[ii]);
	    pChi2Min2 = *(of.cp.v[jj]);
	    chi2Min   = chi2;
	    found     = true;
	    foundHere = true;
	  }
	}
      }
    }
  }
  cout << ((foundHere) ? " +++" : " ") << endl;
  *(of.cp.v[ii]) = pChi2Min1;
  *(of.cp.v[jj]) = pChi2Min2;
  cout << "After  " << n*n << " point local 2D scan of variable " << of.cp.GetName(ii) << "=" << *(of.cp.v[ii]) 
       << " and " << of.cp.GetName(jj) << "=" << *(of.cp.v[jj]) << " chi=" << of() 
       << ((foundHere) ? " +++" : " ")
       << endl;
  cout.precision(p);
  return chi2Min;
}

double ThreeDimensionalScan(OptFunction & of, const int ii, const int jj, const int kk, const int n, bool & found) {
  cout << "Before " << n*n*n << " point 3D scan of variables " << of.cp.GetName(ii) << "=" << *(of.cp.v[ii]) 
       << " and " << of.cp.GetName(jj) << "=" << *(of.cp.v[jj]) 
       << " and " << of.cp.GetName(kk) << "=" << *(of.cp.v[kk]) 
       << " chi=" << of(); // << endl;

  if (ii==jj) error("ii==jj",ii,jj);
  if (ii==kk) error("ii==kk",ii,kk);
  if (jj==kk) error("jj==kk",jj,kk);

  double pMin1,pMax1;
  of.cp.ParameterRange(ii,pMin1,pMax1);
  double pMin2,pMax2;
  of.cp.ParameterRange(jj,pMin2,pMax2);
  double pMin3,pMax3;
  of.cp.ParameterRange(kk,pMin3,pMax3);

  double chi2Min=of();
  double pChi2Min1=*(of.cp.v[ii]);
  double pChi2Min2=*(of.cp.v[jj]);
  double pChi2Min3=*(of.cp.v[kk]);
  double di = Random();
  double dj = Random();
  double dk = Random();

  bool foundHere = false;
  for(int i=0; i<n; i++) {
    *(of.cp.v[ii]) = pMin1+(double(i)+di)/double(n)*(pMax1-pMin1);
    for(int j=0; j<n; j++) {
      *(of.cp.v[jj]) = pMin2+(double(j)+dj)/double(n)*(pMax2-pMin2);
      for(int k=0; k<n; k++) {
	*(of.cp.v[kk]) = pMin3+(double(k)+dk)/double(n)*(pMax3-pMin3);
	double chi2 = of();
	if (chi2<chi2Min) {
	  pChi2Min1 = *(of.cp.v[ii]);
	  pChi2Min2 = *(of.cp.v[jj]);
	  pChi2Min3 = *(of.cp.v[kk]);
	  chi2Min  = chi2;
	  found    = true;
	  foundHere= true;
	}
      }
    }
  }
  cout << ((foundHere) ? " +++" : " ") << endl;
  //  Write3(ii,chi2Min,pChi2Min);
  *(of.cp.v[ii]) = pChi2Min1;
  *(of.cp.v[jj]) = pChi2Min2;
  *(of.cp.v[kk]) = pChi2Min3;
  cout << "After  " << n*n*n << " point 3D scan of variables " << of.cp.GetName(ii) << "=" << *(of.cp.v[ii]) 
       << " and " << of.cp.GetName(jj) << "=" << *(of.cp.v[jj]) 
       << " and " << of.cp.GetName(kk) << "=" << *(of.cp.v[kk]) 
       << " chi=" << of() 
       << ((foundHere) ? " +++" : " ")
       << endl;
  return chi2Min;
}

void Randomize(CycloidParameters & cp) {
  for(int i=0; i<cp.nOpt; i++) {
    double pMin,pMax;
    cp.ParameterRange(i,pMin,pMax);
    *(cp.v[i]) = pMin+Random()*(pMax-pMin);
  }
}

///////////////////////////////////////////////////////////////////////

void AnalyzeDerivatives(OptFunction & of) {
  CycloidParameters cpOrg=of.cp;
  double chi2 = of();
  for(int ii=0; ii<of.cp.nOpt; ii++) { // lopp over all parameters that will be changed d  
    double eps1=1e-3;
    double eps2=1e-14;
    int n=200;
    for(int i=0; i<=n; i++) {
      double eps=eps1*exp(double(i)/double(n)*log(eps2/eps1));
      double delta = eps*of.cp.ParameterRange(ii); // 0.005%
      *(of.cp.v[ii]) = *(cpOrg.v[ii]) + delta;
      double chi2i = of();
      Write7(ii,i,eps,delta,chi2i,chi2i-chi2,(chi2i-chi2)/delta);
      delta *= -1.0;
      *(of.cp.v[ii]) = *(cpOrg.v[ii]) + delta;
      chi2i = of();
      Write7(ii,i,eps,delta,chi2i,chi2i-chi2,(chi2i-chi2)/delta);
    }
  }
  Quit("WW");
}

int main(int argc, char *argv[]) {
  string filename = "TyrrelLonLat2.txt";
  bool reorder=false;

  CycloidParameters cp;

  if (argc!=3 && argc!=4) {
    cout << "Europa cycloid C++ fitting program" << endl << endl;
    cout << argv[0] << " filename  number_of_arcs<int> " << endl;
    cout << argv[0] << " filename  number_of_arcs<int> r " << endl << endl;
    return 0;
  } else {
    filename = argv[1];
    cp.nArcs = atoi(argv[2]);
    if (cp.nArcs<=0 || cp.nArcs>=100) error("Wrong number of arcs specified",cp.nArcs);
    if (argc==4) {
      string s = argv[3];
      if (s!="r") error("wrong second argument");
      reorder = true;
    }
  }

  // Print input parameters
  cout << "Command line: ";
  for(int i=0; i<argc; i++) {
    cout << argv[i] << " ";
  }
  cout << endl;
  //  Write3(filename,cp.nArcs,reorder);

  cp.latDeg = 15.0; // overwritten below with data from file
  cp.pixDiff = 0.2;
  cp.res     = 1000;
  cp.radius  = 1562.0; // Hermes: 1562000;  Thin shell: 1561500
  cp.steps     = 850.0; // *10 how many steps in orbit

  //cp.eccAmpDeg = 1.46344291;  // 1.42000947
  //cp.eccPhaseDeg = 17.71933217; // 173.5926446
  
  InitRandom(true); // use different random numbers each time
  Write(Random());
  Randomize(cp);
  
  //cp.NSRdelta = 0.00000000;
  //cp.lonDeg   = 32.115;
  
// Alex Visco 1r
//  cp.firstStep = 92.12022206;
//  cp.lonDeg   = 32.115; //309.7813303
//  cp.oblqDeg  = 0.1; //0.7540096499
//  cp.phaseDeg =  48.25395781;
//  cp.NSRdelta = 0.000000001;
//  cp.libAmpDeg = 0.0;
//  cp.libPhaseDeg = 0.0;
//  cp.threshold = 0.8383560519;
//  cp.speed     = 0.9418831271;
//  cp.deltaPhaseDeg = 0.0;
  
    //  // Alex OBLQ
//  cp.firstStep = 90.53813847;
//  cp.lonDeg   =  307.7666949;
//  cp.oblqDeg  = 0.8318765055;
//  cp.phaseDeg =  45.22795869;
//  cp.NSRdelta = 0.0;
//  cp.libAmpDeg = 0.0;
//  cp.libPhaseDeg = 0.0;
//  cp.threshold = 0.8306583025;
//  cp.speed     = 0.9110724955;
//  cp.deltaPhaseDeg = 0.0;
  
//  // Alex new oblq case with interior 1
//  cp.firstStep = 694.8163379;
//  cp.lonDeg   = 184.4627102; //309.7813303
//  cp.oblqDeg  = 0.6074544763; //0.7540096499
//  cp.phaseDeg =  31.89838603;
//  cp.NSRdelta = 0.000000;
//  cp.libAmpDeg = 0.0;
//  cp.libPhaseDeg = 0.0;
//  cp.threshold = 0.531878199;
//  cp.speed     = 1.026173492;
//  cp.deltaPhaseDeg = 0.0;

  // Alex w/our ecc and libration
//  cp.firstStep = 824.6824977;
//  cp.lonDeg   = 35.69133413; 
//  cp.oblqDeg  = 0.872043615;
//  cp.phaseDeg =  317.045978;
//  cp.NSRdelta = 0.0;
//  cp.libAmpDeg = 0.0;
//  cp.libPhaseDeg = 0.0;
//  cp.threshold = 0.5032441716;
//  cp.speed     = 1.236430096;
//  cp.deltaPhaseDeg = 0.0;

  // Alex w/Hermes' ecc and libration
//  cp.firstStep = 824.4186003;
//  cp.lonDeg   = 35.82122934; 
//  cp.oblqDeg  = 0.875664976;
//  cp.phaseDeg =  316.9347476;
//  cp.NSRdelta = 0.0;
//  cp.libAmpDeg = 0.0;
//  cp.libPhaseDeg = 0.0;
//  cp.threshold = 0.4986486012;
//  cp.speed     = 1.232501119;
//  cp.deltaPhaseDeg = 0.0;
  //  
//  // Cil OBLQ
//  cp.firstStep = 505.0993888;
//  cp.lonDeg   = 52.76128398;
//  cp.oblqDeg  = 0.6853921333;
//  cp.phaseDeg =  8.54002165;
//  cp.NSRdelta    = 0.0;
//  cp.libAmpDeg = 0.0;
//  cp.libPhaseDeg = 0.0;
//  cp.threshold = 0.3476305393;
//  cp.speed     = 0.560803445;
//  cp.deltaPhaseDeg = 0.0;

//    // Del OBLQ
//    cp.firstStep = 435.4357187;
//    cp.lonDeg   = 208.7878697;
//    cp.oblqDeg  = 0.4429288147;
//    cp.phaseDeg =  187.7444267;
//    cp.NSRdelta    = 0.0;
//    cp.libAmpDeg = 0.0;
//    cp.libPhaseDeg = 0.0;
//    cp.threshold = 0.6259688114;
//    cp.speed     = 0.6617663166;
//    cp.deltaPhaseDeg = 0.0;
//
//      // Sid OBLQ
//  cp.firstStep = 451.6097507;
//  cp.lonDeg   = 245.6485646;
//  cp.oblqDeg  = 0.5792415798;
//  cp.phaseDeg =  122.4561786;
//  cp.NSRdelta    = 0.0;
//  cp.libAmpDeg = 0.0;
//  cp.libPhaseDeg = 0.0;
//  cp.threshold = 0.4192731443;
//  cp.speed     = 0.4502617851;
//  cp.deltaPhaseDeg = 0.0;
//
//    // Tyrrel (EQ1) OBLQ
//  cp.firstStep = 1589.446335;
//  cp.lonDeg   = 220.759044;
//  cp.oblqDeg  = 0.1753590155;
//  cp.phaseDeg =  11.83844637;
//  cp.NSRdelta    = 0.0;
//  cp.libAmpDeg = 0.0;
//  cp.libPhaseDeg = 0.0;
//  cp.threshold = 0.4508149621;
//  cp.speed     = 0.8257468642;
//  cp.deltaPhaseDeg = 0.0;
  
//      // Carly (EQ3) OBLQ
//  cp.firstStep = 662.0860227;
//  cp.lonDeg   = 7.780883413;
//  cp.oblqDeg  = 0.5664911897;
//  cp.phaseDeg =  193.5736321;
//  cp.NSRdelta    = 0.0;
//  cp.libAmpDeg = 0.0;
//  cp.libPhaseDeg = 0.0;
//  cp.threshold = 0.6018510168;
//  cp.speed     = 0.9350019446;
//  cp.deltaPhaseDeg = 0.0;
//  

  int nCycloidPoints;
  Array1 <double> cycloidLonDeg, cycloidLatDeg;
  //  ReadCycloidFile("SynthDataNoNSR.txt", nCycloidPoints, cycloidLonDeg, cycloidLatDeg);
  //  ReadCycloidFile("TyrrelLonLat2.txt", nCycloidPoints, cycloidLonDeg, cycloidLatDeg);
  ReadCycloidFile(filename, nCycloidPoints, cycloidLonDeg, cycloidLatDeg);
  cout << endl;

  if (reorder) {
    cycloidLonDeg.Invert();
    cycloidLatDeg.Invert();
  }

  cp.dir     = DetermineDirection(cycloidLonDeg);

  Write(nCycloidPoints);
  cp.latDeg = cycloidLatDeg[0];
  //cp.lonDeg = cycloidLonDeg[0];

  Array1 <double> cycloidArcLength = CalcArcLengths(cycloidLonDeg, cycloidLatDeg);
  
  Write2(cp,GetChi2(cp, cycloidLonDeg, cycloidLatDeg, cycloidArcLength));

  //Quit();

  OptFunction of(cp, cycloidLonDeg, cycloidLatDeg, cycloidArcLength);
  double chi2 = of();

  cout << endl;
  cout << "Starting " << of.cp.nOpt << " dimensional search from: " << endl;
  Write2(chi2,of.cp);
  cout << endl;
  cout << endl;

  int nScan    = 50;
  //  int nScanMax = 12800; 
  int nScanMax = 1000; // ***** Changed from 400 by Alyssa
  //  int nScanMax = 4000; 


  double chi2Target=10.0;
  bool flag = true; // yes=run BFGS, either true at beginning, or when the scans have found a new starting point

  // If old results cannto be reproduced because if insufficient digits then run MC first. Uncomment the following line
  //  chi2 = MonteCarloSearch(of, 1000, flag);

  do {
  
    cout << endl << "Start 1D scans ";
    write1(chi2); cout << " <-- " << of.cp << endl;
    for(int i=0; i<of.cp.nOpt; i++) {
      chi2 = OneDimensionalScan(of, i, nScan, flag);
    }
    if (flag) {
      chi2 = BFGSOptimization(of);
      of.PrintFit();
      cout << "BFGS +++"; write1(chi2); cout << " <-- " << of.cp << endl;
      if (chi2<chi2Target) break;
      flag = false;
    }

    chi2 = MonteCarloSearch(of, nScan, flag);
    cout.precision(8);
	write1(chi2);
	
    if (flag) {
      chi2 = BFGSOptimization(of);
      of.PrintFit();
      cout << "BFGS +++"; write1(chi2); cout << " <-- " << of.cp << endl;
      if (chi2<chi2Target) break;
      flag = false;
    }

    //    chi2 = TunnelingScan(of, nScan, flag);

    cout << endl << "Start 2D scans ";
    write1(chi2); cout << " <-- " << of.cp << endl;
    for(int i=0; i<of.cp.nOpt-1; i++) {
      for(int j=i+1; j<of.cp.nOpt; j++) {
	chi2 = TwoDimensionalScan(of, i, j, int(sqrt(nScan)), flag);
      }
    }
    if (flag) {
      chi2 = BFGSOptimization(of);
      of.PrintFit();
      cout << "BFGS +++"; write1(chi2); cout << " <-- " << of.cp << endl;
      flag = false;
    }

    cout << endl << "Start local 2D scans ";
    write1(chi2); cout << " <-- " << of.cp << endl;
    for(int i=0; i<of.cp.nOpt-1; i++) {
      for(int j=i+1; j<of.cp.nOpt; j++) {
	chi2 = LocalTwoDimensionalScan(of, i, j, int(sqrt(nScan)), flag);
      }
    }
    if (flag) {
      chi2 = BFGSOptimization(of);
      of.PrintFit();
      cout << "BFGS +++"; write1(chi2); cout << " <-- " << of.cp << endl;
      flag = false;
    }

    cout << endl << "Start 3D scans ";
    write1(chi2); cout << " <-- " << of.cp << endl;
    for(int i=0; i<of.cp.nOpt-2; i++) {
      for(int j=i+1; j<of.cp.nOpt-1; j++) {
	for(int k=j+1; k<of.cp.nOpt; k++) {
	  chi2 = ThreeDimensionalScan(of, i, j, k, int(pow(nScan,1.0/3.0)), flag);
	}
      }
    }
    if (flag) {
      chi2 = BFGSOptimization(of);
      of.PrintFit();
      cout << "BFGS +++"; write1(chi2); cout << " <-- " << of.cp << endl;
      flag = false;
    }


    nScan *= 2; // improve scan quality assuming we are stuck somewhere
  } while (chi2>chi2Target && nScan<nScanMax);

  cout << endl;
  if (chi2<=chi2Target) {
    cout << "Hurray! Found a cycloid with a chi2 below the target chi2 of " << chi2Target << endl;
  } else {
    cout << "Sorry. Did not find a cycloid with a chi2 below the target chi2 of " << chi2Target << endl;
  }
  cout << endl;
  cout.precision(10);
  cout << "FINAL ANSWER: Best fit obtained for: " << endl;
  write1(chi2); cout << " <-- " << of.cp << endl;
  cout << endl;

}


