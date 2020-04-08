////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                            //
// Europa cycloid                                                                                             //
//                                                                                                            //
// Alyssa Sarid, Burkhard Militzer                                UC Berkeley, 07-08-08                       //
//                                                                                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _MAKECRACKMODULE_
#define _MAKECRACKMODULE_

void MakeCrack(const CycloidParameters & cp, const double lengthLimit, 
	       Array1 <double> & lonDegs, Array1 <double> & latDegs,  // formerly 'points'
	       bool & startStressTooLow, bool & crackTooShort, const bool print=false);

#endif // _MAKECRACKMODULE_
