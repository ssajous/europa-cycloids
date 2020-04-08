////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                            //
// BFSG Opimization algorithm                                                                                 //
//                                                                                                            //
// adopted to C++ by Burkhard Militzer                                                    Washington, DC 2007 //
//                                                                                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Array.h"
#include "BFGSOptimization.h"
#include "CycloidParameters.h"
#include "EuropaCycloids.h"
#include "MatrixAlgebra.h"

inline double Norm(const Array1 <double> & a) {
  double norm = 0.0;
  for(int i=0; i<a.Size(); i++) {
    norm += a[i]*a[i];
  }
  return sqrt(norm);
}

inline int MaxNormIndex(const Array1 <double> & a) {
  double n2 = 0.0;
  int iMax = 0;
  for(int i=0; i<a.Size(); i++) {
    double n2i = a[i]*a[i];
    if (n2i>n2) {
      n2 = n2i;
      iMax = i;
    }
  }
  return iMax;
}

void PrintLineSearchFunction(OptFunction & f, const Array1 <double> & delta, double & a) {
  Write(delta);
  cout.precision(12);
  int nn=100;
  Array1 <double> pOrg = f.GetParameters();
  //  f.PrintFit();
  for(int ii=-nn; ii<=nn; ii++) {
    //    double aa = 100/Norm(delta)*a*double(ii)/double(nn);
    //    double aa = 1.0/Norm(delta)*double(ii)/double(nn);
    double aa = 0.01/Norm(delta)*double(ii)/double(nn);

    Array1 <double> d = delta;
    d *= aa;

    f.SetParameters(pOrg);
    f.ApplyPerturbations(d);
    //    double e = f();
    write2(ii,aa);
    f.Print();
    // XCXC run without any of three lines
    // XCXC 2 switch on one of the three statements to see if there is a big jump,
    //    if (ii==-1) f.PrintFit();
    //    if (ii==0) f.PrintFit(); // compare data and sim cycloid
    //    if (ii==1) f.PrintFit();
  }
  Quit("PrintLineSearchFunction");
}

double LineSearch(OptFunction & f, const Array1 <double> & delta, double & a) {
  Array1 <double> pOrg = f.GetParameters();
  double e = f();

  Array1 <double> d = delta;
  d *= -a;

  Array1 <double> p = f.GetParameters();
  f.ApplyPerturbations(d);
  double e1 = f();

  double aa;
  if (e1<e) {
    f.ApplyPerturbations(d); // apply perturbations again
    double e2 = f();
    cout << " LS++ "; Write5(e,e1,e2,e1-e,e2-e);
    if (e2<=e1) {
      aa=2.0*a;
      a*=2.0; // larger steps next time
    } else {
      double d2 = e2-e;
      double d1 = e1-e;
      aa = a * (d2-4.0*d1) / (2.0*d2-4.0*d1); // estimate optimal future step size
    }
  } else { // did not get better, go in opposite direction
    f.SetParameters(p);
    d *= -1.0;
    f.ApplyPerturbations(d);
    double em1 = f();
    cout << " LS-- "; Write5(em1,e,e1,em1-e,e1-e);
    if (em1<=e) {
      aa =-a; 
    } else {
      double d2 = e1-em1;
      double d1 = e -em1;
      aa = a * (-d2) / (2.0*d2-4.0*d1);
    }
    a /= 2.0; // smaller steps next time
  }
  //  Write3(aa,aa/aOrg,a/aOrg);

  f.SetParameters(pOrg); // restore original parameter set
  return aa;
}

double LineSearchLog(OptFunction & f, const Array1 <double> & delta, double & a) {
  Array1 <double> pOrg = f.GetParameters();
  double e = f();
  //  Array1<double> ee;

  double eMin,aMin;
  int iMin;
  int nMin = -2;
  int nMax = +2;
  for(int i=nMin; i<=nMax; i++) {
    Array1 <double> d = delta;
    double ai = +a*pow(2.0,i);
    d *= -ai;
    f.SetParameters(pOrg);
    f.ApplyPerturbations(d);
    double ei = f();
    //    ee.PushBack(ei);
    if (i==nMin || ei<eMin) { eMin=ei; aMin=ai; iMin=i; } 
  }
  //  double dd = ee[0];
  //  ee -= dd;
  //  Write(ee);

  //  Write(delta);

  if (eMin<e) {
    cout << " LS++ "; Write4(e,eMin,eMin-e,iMin);
    a = aMin;
  } else {
    cout << " LS+B "; Write4(e,eMin,eMin-e,iMin);
    Array1 <double> d = delta;

    //    Array1<double> ee;
    for(int i=nMin; i<=nMax; i++) {
      Array1 <double> d = delta;
      double ai = -a*pow(2.0,i);
      d *= -ai;
      f.SetParameters(pOrg);
      f.ApplyPerturbations(d);
      double ei = f();
      //      ee.PushBack(ei);
      if (i==nMin || ei<eMin) { eMin=ei; aMin=ai; iMin=i; } 
    }
    //    double dd = ee[0];
    //    ee -= dd;
    //    ee -= ee[0]; // does not work
    //    Write(ee);

    if (eMin<e) {
      cout << " LS-G "; Write5(e,eMin,eMin-e,a,iMin);
      a = abs(aMin); // 'a' should never change sign
    } else {
      cout << " LS-- "; Write5(e,eMin,eMin-e,a,iMin);
      a /= 2.0; // 'a' should never change sign
      //      warning("Line search did not improve anything",a);
    }

  } 

  f.SetParameters(pOrg); // restore orginal parameters
  return aMin;
}

// relax one degree of freedom after the other
double BFGSOptimization(OptFunction & f, double a) {
  cout << endl << "/----------------BFGS optimization started------------------\\" << endl << endl;

  int nF = f.GetNumber();
  //  double a=1.0e-3;

  Array1<double> grad(nF);
  Array2<double> H = IdentityMatrix(nF);
  Array1<double> delta(nF);
  Array1<double> pOld(nF);

  double e = f();
  int n=0;
  int nEnd=0;
  int nReset=0;
  int nResetMax=5;
  bool getNewDeltaFlag = true;
  do {
   
    f.GetGradient(grad);
    if (getNewDeltaFlag) MatrixTimesVector(H,grad,delta,nF);
    //    Write(delta);
    //    Write(grad);

    //    double k = LineSearch(f,delta,a); // return the current step size 'k' (can be neg) and adjust the future step size 'a' (always pos)
    double k = LineSearchLog(f,delta,a); // return the current step size 'k' and adjust the future step size 'a'

    // Careful, this could yield k=0.0 and then a rotation routine gives "Singular matrix in routine ludcmp" upon matrxi inversion
    if (k==0.0) {
      warning("Problem: LineSearch yield k=0. Check whether the BFGS relaxation has converged");
      break;
    }
    
    delta *= -k;

    f.GetParameters(pOld);
    f.ApplyPerturbations(delta);
    double eNew = f();
    //    if (abs(eNew-e)<1e-2 || k<1e-4) { warning("Lousy lousy BFGS relaxation"); nEnd++; }
    //    if (abs(eNew-e)<1e-4 || k<1e-6) { warning("Lousy BFGS relaxation"); nEnd++; }
    if (abs(eNew-e)<1e-14 || abs(k)<1e-16) {
      nEnd++;

      if (nReset++<nResetMax) {
	IdentityMatrix(H); // reset H matrix, (delta gets overwritten)
	a *= 10;
	f.SetParameters(pOld);

	f.GetGradient(grad,true);
	int ii = MaxNormIndex(grad);
	double normDelta = Norm(delta);
	delta = 0.0;
	if (grad[ii]==0.0) { 
	  //	  error("grad is all zero",ii,grad); 
	  cout << "BFGS cannot reset direction anymore. Quit BFGS now." << endl;
	  break; 
	}
	delta[ii] = grad[ii]; // normDelta;
	a = abs(1e-6*f.cp.ParameterRange(ii)/delta[ii]);
	getNewDeltaFlag = false;
	//      Write(delta);
 	cout << " BFGS:: Reset direction "; write2(nReset,nEnd); 
	cout << endl;
	cout << endl;
	continue;
      }
    } else { 
      nEnd=0;
      if (eNew<e-1e-14) nReset=0; // new
    }
    //    Write6(++n,a,k,nEnd,eNew,eNew-e);
    Write8(++n,a,k,nEnd,e,eNew,eNew-e,grad);
    //    Write(delta);
    /*
    if (abs(delta.Max())+abs(delta.Min())<1e-9) {
      PrintLineSearchFunction(f,delta,a);
      Quit("XX");
    }
    */

    if (eNew<e) {
      e = eNew;
      Array1<double> grad1(nF);
      f.GetGradient(grad1);
      
      Array1 <double> gamma(nF);
      double dbg=0.0;
      for(int i=0; i<nF; i++) {
	gamma[i] = grad1[i]-grad[i];
	dbg     += gamma[i]*delta[i];
	grad[i]  = grad1[i];
      }
      
      Array1 <double> hg(nF);
      MatrixTimesVector(H,gamma,hg,nF);
      double gphg=0.0;
      for(int i=0; i<nF; i++) {
	for(int j=0; j<nF; j++) {
	  gphg += H(i,j)*gamma[i]*gamma[j];
	}
      }
      for(int i=0; i<nF; i++) {
	for(int j=0; j<nF; j++) {
	  H(i,j) += 
	    -(hg(i)*delta(j)/dbg)
	    -(delta(i)*hg(j)/dbg)
	    +(1.0+gphg/dbg)*(delta[i]*delta[j]/dbg);
	}
      }
    } else {
      f.SetParameters(pOld);
      //      a /= 10.0;
      a /= 4.0;
      Write2("Rejected",a);
      // XCXC 1 switch on next line
      //      if (a<=1e-10 && e<800.0) PrintLineSearchFunction(f,delta,a);
    }

    f.Print();
    cout << endl;
    getNewDeltaFlag = true;
  } while (nEnd<5); 
  //  } while (nEnd<5 && n<10); <<<-- quit early

  //  PrintLineSearchFunction(f,delta,a);
  //  delta = 0.0;
  //  delta[2] = 1.0;
  //  PrintLineSearchFunction(f,delta,a);
  //  f.GetGradient(grad);
  //  PrintLineSearchFunction(f,grad,a);
  f.GetGradient(grad,true);

  cout << "\\_____________ End of BFGS optimization.___________/" << endl << endl;
  //  Quit("Exit at end of BFGS Relaxation");

  return f();
}

