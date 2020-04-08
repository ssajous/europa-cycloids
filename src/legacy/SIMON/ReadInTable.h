/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Read text files                                                         // 
//                                                                         //
// Burkhard Militzer                                  Livermore 4-27-01    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _READINTABLE_
#define _READINTABLE_

#include "Parser.h"
#include "Array.h"


void ReadInTable(const string & fileName, 
		 const Array1<int> & nPos,
		 const int nVar,
		 int & n,
		 Array1 <Array1 <double> *> & x);

void ReadInTable(const string & fileName,
		 const int nx,
		 const int ny,
		 int & n,
		 Array1 <double> & x,
		 Array1 <double> & y);

void ReadInTable(const string & fileName,
		 const int nx,
		 const int ny,
		 const int nz,
		 int & n,
		 Array1 <double> & x,
		 Array1 <double> & y,
		 Array1 <double> & z);

void ReadInTable(const string & fileName,
		 const int n1,
		 const int n2,
		 const int n3,
		 const int n4,
		 int & n,
		 Array1 <double> & x1,
		 Array1 <double> & x2,
		 Array1 <double> & x3, 
		 Array1 <double> & x4);

void ReadInTable(const string & fileName,
		 const int n1,
		 const int n2,
		 const int n3,
		 const int n4,
		 const int n5,
		 int & n,
		 Array1 <double> & x1,
		 Array1 <double> & x2,
		 Array1 <double> & x3, 
		 Array1 <double> & x4,
		 Array1 <double> & x5);

void ReadInTable(const string & fileName,
		 const int n1,
		 const int n2,
		 const int n3,
		 const int n4,
		 const int n5,
		 const int n6,
		 int & n,
		 Array1 <double> & x1,
		 Array1 <double> & x2,
		 Array1 <double> & x3, 
		 Array1 <double> & x4,
		 Array1 <double> & x5,
		 Array1 <double> & x6);

void ReadInTable(const string & fileName,
		 const int n1,
		 const int n2,
		 const int n3,
		 const int n4,
		 const int n5,
		 const int n6,
		 const int n7,
		 int & n,
		 Array1 <double> & x1,
		 Array1 <double> & x2,
		 Array1 <double> & x3, 
		 Array1 <double> & x4,
		 Array1 <double> & x5,
		 Array1 <double> & x6,
		 Array1 <double> & x7);

int FindInTable(const double x, const Array1 <double> & t);
int FindInReverseTable(const double x, const Array1 <double> & t);

#endif // _READINTABLE_
