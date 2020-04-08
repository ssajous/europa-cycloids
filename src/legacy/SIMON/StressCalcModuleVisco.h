////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                            //
// Viscoelastic cycloid modeling (adapted from Europa cycloid, with Burkhard Militzer                         //
//                                                                                                            //
// Alyssa Rhoden                                NASA Goddard 3/2013                                           //
//                                                                                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _STRESSCALCMODULEVISCO_
#define _STRESSCALCMODULEVISCO_

double GetStress(const double colat, const double lon, const int steps, const double thisStep, const double oblq, const double phase, const double NSRdelta, const double libAmp, const double libPhase, double & bigHeading);

#endif
