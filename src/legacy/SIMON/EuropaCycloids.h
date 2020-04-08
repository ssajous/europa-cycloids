////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                            //
// Europa cycloid                                                                                             //
//                                                                                                            //
// Alyssa Sarid, Burkhard Militzer                                UC Berkeley, 07-08-08                       //
//                                                                                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _EUROPACYCLOIDS_
#define _EUROPACYCLOIDS_

#include "CycloidParameters.h"
#include "Array.h"

double GetChi2(const CycloidParameters & cp, 
	       const Array1<double> & dataLonDegs, const Array1<double> & dataLatDegs, const Array1<double> & dataArcLengths, 
	       const double invalid = 1e50, const bool print=false);

#endif _EUROPACYCLOIDS_
