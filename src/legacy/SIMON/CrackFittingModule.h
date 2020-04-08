////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                            //
// Europa cycloid                                                                                             //
//                                                                                                            //
// Alyssa Sarid, Burkhard Militzer                                UC Berkeley, 07-08-08                       //
//                                                                                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _CRACKFITTINGMODULE_
#define _CRACKFITTINGMODULE_

Array1 <double> CalcArcLengths(const Array1 <double> & lonDeg, const Array1 <double> & latDeg);
double CalcChi2(const Array1 <double> & simLonDegs, const Array1 <double> & simLatDegs, 
		const Array1 <double> & dataLonDegs, const Array1 <double> & dataLatDegs, const Array1 <double> & dataArcLengths, 
		const double pixDiff, const double res, const bool print=false);
double CalcChi2Max(const CycloidParameters & cp);
//double CalcChi2Max(const Array1 <double> & dataArcLengths, 
//		   const double pixDiff, const double res);
bool DetermineDirection(const Array1 <double> & lonDeg);

#endif
