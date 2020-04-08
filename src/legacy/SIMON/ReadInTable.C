/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Read text files                                                         // 
//                                                                         //
// Burkhard Militzer                                  Livermore 4-27-01    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "ReadInTable.h"

#define NO_WARNINGS

void ReadInTable(const string & fileName, 
		 const Array1<int> & nPos,
		 const int nVar,
		 int & n,
		 Array1 <Array1 <double> *> & x) {
  n = 0;
  int nMax=0;
  for(int i=0; i<nVar; i++) {
    const int & nP(nPos[i]);
    if (nP<0) 
      error("Incorrect column index",i,nP);
    if (nP>nMax) nMax=nP;
  }

  cout << "Opening file \"" << fileName << "\" to read " << nVar << " columns of floating point values: Columns= "; 
  for(int i=0; i<nVar; i++) {
    cout << nPos[i] << " ";
  }
  cout << endl;

  ifstream is(fileName.c_str());
  if (!is) 
    error("Cannot open file",fileName);
  Parser p(is);

  while(p.ReadLine()) {

    const int nW=p.GetNWords();
    if (nW<nMax) {
      warning("Input line has too few words",p.GetLineNumber(),nW,"\""+p.GetLineString()+"\"");
    } else {
      n++;
      //      cout << "Read table: " << n << ":  ";
      for(int i=0; i<nVar; i++) {
	double q=p.GetDouble(nPos[i]-1);
	(x[i])->PushBack(q);
	//	cout << "  " << q ;
      }
      //      cout << endl;
    }
  }
  cout << "Read " << n << " lines from file \"" << fileName << "\"." << endl;
  is.close();
}

void ReadInTable(const string & fileName,
		 const int nx,
		 const int ny,
		 int & n,
		 Array1 <double> & x,
		 Array1 <double> & y) {
  const int nn=2;
  Array1 <Array1 <double> *> xx(nn);

  xx[0] = & x;
  xx[1] = & y;

  Array1 <int> nPos(nn);

  nPos[0] = nx;
  nPos[1] = ny;

  ReadInTable(fileName,nPos,nn,n,xx);
}

void ReadInTable(const string & fileName,
		 const int nx,
		 const int ny,
		 const int nz,
		 int & n,
		 Array1 <double> & x,
		 Array1 <double> & y,
		 Array1 <double> & z) {
  Array1 <Array1 <double> *> xx(3);

  xx[0] = & x;
  xx[1] = & y;
  xx[2] = & z;

  Array1 <int> nPos(3);

  nPos[0] = nx;
  nPos[1] = ny;
  nPos[2] = nz;

  ReadInTable(fileName,nPos,3,n,xx);
}

void ReadInTable(const string & fileName,
		 const int n1,
		 const int n2,
		 const int n3,
		 const int n4,
		 int & n,
		 Array1 <double> & x1,
		 Array1 <double> & x2,
		 Array1 <double> & x3, 
		 Array1 <double> & x4) {
  const int ii=4;
  Array1 <Array1 <double> *> xx(ii);

  xx[0] = & x1;
  xx[1] = & x2;
  xx[2] = & x3;
  xx[3] = & x4;

  Array1 <int> nPos(ii);

  nPos[0] = n1;
  nPos[1] = n2;
  nPos[2] = n3;
  nPos[3] = n4;

  ReadInTable(fileName,nPos,ii,n,xx);
}

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
		 Array1 <double> & x5) {
  const int ii=5;
  Array1 <Array1 <double> *> xx(ii);

  xx[0] = & x1;
  xx[1] = & x2;
  xx[2] = & x3;
  xx[3] = & x4;
  xx[4] = & x5;

  Array1 <int> nPos(ii);

  nPos[0] = n1;
  nPos[1] = n2;
  nPos[2] = n3;
  nPos[3] = n4;
  nPos[4] = n5;

  ReadInTable(fileName,nPos,ii,n,xx);
}

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
		 Array1 <double> & x6) {
  const int ii=6;
  Array1 <Array1 <double> *> xx(ii);

  xx[0] = & x1;
  xx[1] = & x2;
  xx[2] = & x3;
  xx[3] = & x4;
  xx[4] = & x5;
  xx[5] = & x6;

  Array1 <int> nPos(ii);

  nPos[0] = n1;
  nPos[1] = n2;
  nPos[2] = n3;
  nPos[3] = n4;
  nPos[4] = n5;
  nPos[5] = n6;

  ReadInTable(fileName,nPos,ii,n,xx);
}

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
		 Array1 <double> & x7) {
  const int ii=7;
  Array1 <Array1 <double> *> xx(ii);

  xx[0] = & x1;
  xx[1] = & x2;
  xx[2] = & x3;
  xx[3] = & x4;
  xx[4] = & x5;
  xx[5] = & x6;
  xx[6] = & x7;

  Array1 <int> nPos(ii);

  nPos[0] = n1;
  nPos[1] = n2;
  nPos[2] = n3;
  nPos[3] = n4;
  nPos[4] = n5;
  nPos[5] = n6;
  nPos[6] = n7;

  ReadInTable(fileName,nPos,ii,n,xx);
}

// assumes t[i]>=t[i+1]
// find i for t[i] <= x < t[i+1] with i <= n-2
int FindInTable(const double x, const Array1 <double> & t) {
  const int n=t.Size();
  if (n==0)
    error("Empty table");
  
  if (x<t[0]) {
#ifndef NO_WARNINGS
    warning("Look up beyond lower table boundaries",x,t[0],t[n-1]);
#endif
    return 0;
  }

  int ix;
  for(ix=1; ix<n; ix++) {
    if (x<t[ix]) {
      return ix-1;
    }
  }
  
  if (ix==n) {
#ifndef NO_WARNINGS
    warning("Look up beyond upper table boundaries",x,t[0],t[n-1]);
#endif
  }

  ix=n-2;
  return ix;
}

// assumes t[i]<=t[i+1]
// find i for t[i+1] > x > t[i] with i <= n-2
int FindInReverseTable(const double x, const Array1 <double> & t) {
  const int n=t.Size();
  if (n==0)
    error("Empty table");
  
  if (x>t[0]) {
    warning("Look up beyond lower table boundaries",x,t[0],t[n-1]);
    return 0;
  }

  int ix;
  for(ix=1; ix<n; ix++) {
    if (x>t[ix]) {
      return ix-1;
    }
  }
  
  if (ix==n) {
    warning("Look up beyond upper table boundaries",x,t[0],t[n-1]);
  }

  ix=n-2;
  return ix;
}


