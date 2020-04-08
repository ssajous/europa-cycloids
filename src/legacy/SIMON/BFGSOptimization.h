////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                            //
// BFSG Opimization algorithm                                                                                 //
//                                                                                                            //
// adopted to C++ by Burkhard Militzer                                                    Washington, DC 2007 //
//                                                                                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _BFGSOPTIMZATION_
#define _BFGSOPTIMZATION_

#include "Array.h"
#include "EuropaCycloids.h"
#include "CrackFittingModule.h"

class OptFunction {
public:
  CycloidParameters cp;
  const Array1<double> & dataLonDegs; 
  const Array1<double> & dataLatDegs; 
  const Array1<double> & dataArcLengths;
  int n;
  const double chi2Max;

  OptFunction(const CycloidParameters & cp_, 
	      const Array1<double> & dataLonDegs_, const Array1<double> & dataLatDegs_, const Array1<double> & dataArcLengths_):
    cp(cp_), dataLonDegs(dataLonDegs_), dataLatDegs(dataLatDegs_), dataArcLengths(dataArcLengths_),n(cp.nOpt),
    chi2Max(CalcChi2Max(cp)) 
    {}
    
  int GetNumber() const { // how many fit parameters
    return n;
  }

  void ApplyPerturbations(const Array1 <double> d) {
    if (d.Size()!=n) error("Size problem",d.Size(),n);
    for(int i=0; i<n; i++) {
      FitParameter(i) += d[i];
    }
  }

  // symmetric formule gives better fit at low P
  //  void GetGradient(Array1 <double> & grad) {
  void GetGradient(Array1 <double> & grad, const bool modify=false) {
    if (grad.Size()!=n) error("Size problem",grad.Size(),n);
    double f = FitQuality();
    if (f==chi2Max) error("Cannot compute gradient at an invalid point",cp);

    for(int i=0; i<n; i++) {
      double p = FitParameter(i);
      //      double eps=1e-3;
      //      double eps=max(1e-8,abs(delta[i]*a));
      //      double eps=1e-6*cp.ParameterRange(i); // tested it: 1e-6 ot 1e-07 seemed good
      double eps=1e-8*cp.ParameterRange(i); // tested it: 1e-6 ot 1e-07 seemed good
      double d   = eps; // absolute
      FitParameter(i) += d;
      //      cout << cp << endl;
      double fp = FitQuality();
      FitParameter(i) -= 2.0*d;
      //      cout << cp << endl;
      double fm = FitQuality();
      FitParameter(i) = p;
      //      cout << cp << endl;

      double gp = (fp-f)/d;
      double gm = (f-fm)/d;
      if (fp==chi2Max || abs(gp)>10.0*abs(gm)) {
	if (fm==chi2Max && fp==chi2Max) error("Both points a invalid",cp.GetName(i),p);
	grad[i] = gm;
	//	if (modify && gm<0.0 && gp>0.0) grad[i]=0.0; // runs BFGS into a wall
      } else if  (fm==chi2Max || abs(gm)>10.0*abs(gp)) {
	grad[i] = gp;
	//	if (modify && gp<0.0 && gm>0.0) grad[i]=0.0; // runs BFGS into a wall
      } else {
	grad[i] = (gp+gm)/2.0; // original line
	//	if (modify && gm<0.0 && gp>0.0) grad[i]=0.0; // runs BFGS into a wall
	//	if (modify && gp<0.0 && gm>0.0) grad[i]=0.0; // runs BFGS into a wall
      }
      if (modify && gm<0.0 && gp>0.0) grad[i]=0.0; // runs BFGS into a wall
      //      if (gp*gm<0.0) grad[i]=0.0; // if you cannot agree of sign of gradient then set to zero
      if (modify) Write10(i,p,eps,fm,f,fp,gm,gp,chi2Max,grad[i]);
    }

    /*
    double feps=1e-03;
    for(int j=0; j<8; j++) {
      for(int i=0; i<n; i++) {
	double p = FitParameter(i);
	//      double eps=1e-3;
	//      double eps=max(1e-8,abs(delta[i]*a));
	double eps=feps*cp.ParameterRange(i);
	double d   = eps; // absolute
	FitParameter(i) += d;
	//      cout << cp << endl;
	double fp = FitQuality();
	FitParameter(i) -= 2.0*d;
	//      cout << cp << endl;
	double fm = FitQuality();
	FitParameter(i) = p;
	//      cout << cp << endl;
	grad[i] = (fp-fm)/(2.0*d);
	//      if (i==4) Write5(p,eps,fp,fm,grad[i]);
	Write7(i,p,eps,fm,f,fp,grad[i]);
      }
      feps*=0.1;
    }
    */

    //    Array1 <double> da = delta;
    //    da *= a;
    //    Write2(grad,da);
    //    Quit("Q Gradient");
  }

  /*
  void GetGradient(Array1 <double> & grad) {
    if (grad.Size()!=n) error("Size problem",grad.Size(),n);
    double f = FitQuality();
    for(int i=0; i<n; i++) {
      double p = FitParameter(i);
      double eps = 1e-6;
      double d   = p*eps;
      FitParameter(i) += d;
      double fp = FitQuality();
      grad[i] = (fp-f)/d;
      FitParameter(i) = p;
     }
  }
  */

  bool AllParametersValid() {
    for(int i=0; i<cp.nOpt; i++) {
      if (!cp.ParameterValid(*(cp.v[i]),i)) return false;
    }
    return true;
  }

  double FitQuality() {
    //    cout << "### " << cp << endl;
    if (!AllParametersValid()) return chi2Max;
    double chi2 = GetChi2(cp, dataLonDegs, dataLatDegs, dataArcLengths, chi2Max);
    if (chi2>chi2Max) error("chi2>chi2Max",chi2,chi2Max);
    //    write1(chi2); cout << " <-- " << cp << endl;
    return chi2;
  }

  void PrintFit() {
    GetChi2(cp, dataLonDegs, dataLatDegs, dataArcLengths, chi2Max, true);
  }

  double operator()() {
    return FitQuality();
  }

  void GetParameters(Array1 <double> & p) const {
    if (p.Size()!=n) error("Size problem",p.Size(),n);
    for(int i=0; i<n; i++) 
      p[i] = FitParameter(i);
  }

  Array1 <double> GetParameters() const {
    Array1 <double> p(n);
    GetParameters(p);
    return p;
  }

  void SetParameters(const Array1 <double> & p) {
    if (p.Size()!=n) error("Size problem",p.Size(),n);
    for(int i=0; i<n; i++) 
      FitParameter(i) = p[i];
    //    Write(cp);
    //    Write(p);
    //    Write(cp);
    //    Quit("QQ");
  }

  void Print() {
    double chi2 = FitQuality();
    write1(chi2); cout << " <-- " << cp << endl;

    /*
    CycloidParameters cp2(cp);
    cp.lonDeg   += 180.0;
    cp.phaseDeg += 180.0;
    chi2 = FitQuality();
    write1(chi2); cout << " <-- " << cp << endl;
    cp = cp2;
    */
  }

private:
  const double & FitParameter(const int i) const {
    return *(cp.v[i]);
  }

  double & FitParameter(const int i) {
    return *(cp.v[i]);
  }
};

double BFGSOptimization(OptFunction & f, double a = 1e-3);

#endif
