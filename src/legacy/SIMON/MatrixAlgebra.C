// Matrix operations
//
// Burkhard Militzer                                    Urbana 4-1-99
//
#include <math.h>
#include "Standard.h"
#include "Array.h"
const double TINY=1.0e-20;

Array2 <double> tmp2;
Array1 <double> tmp11,tmp12;
Array1 <int> itmp1;

// LU decomposition of matrix a, which is overwritten
void ludcmp(Array2 <double> & a, const int n, Array1 <int> & indx, double & d) {
  Array1 <double>& vv(tmp11);
  vv.Resize(n);

  d=1.0;
  for (int i=0;i<n;++i) {
    double big=0.0;
    for (int j=0;j<n;++j) {
      double temp;
      if ((temp=abs(a(i,j))) > big) big=temp;
    }
    if (big == 0.0) error("Singular matrix in routine ludcmp");
    vv[i]=1.0/big;
  }
  
  for (int j=0;j<n;++j) {
    int imax;
    for (int i=0;i<j;++i) {
      double sum=a(i,j);
      for (int k=0;k<i;++k) {
	sum -= a(i,k)*a(k,j);
      }
      a(i,j)=sum;
    }
    double big(0.0);
    for (int i=j;i<n;++i) {
      double sum=a(i,j);
      for (int k=0;k<j;++k) {
	sum -= a(i,k)*a(k,j);
      }
      a(i,j)=sum;
      double dum;
      if ( (dum=vv[i]*abs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (int k=0;k<n;++k) {
	double dum(a(imax,k));
	a(imax,k)=a(j,k);
	a(j,k)=dum;
      }
      d = -d;
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a(j,j) == 0.0) {
      a(j,j)=TINY;
      /*
      Write(n);
      for (int q=0;q<n;++q)
	for (int qq=0;qq<n;++qq)
	  Write3(q,qq,a(q,qq));
      */
      error("Singular matrix in routine ludcmp II.");
    }
    if (j != n-1) {
      double dum=1.0/(a(j,j));
      for (int i=j+1;i<n;++i) a(i,j) *= dum;
    }
  }
}

void lubksb(const Array2 <double> & a, const int n, const Array1 <int> & indx, Array1 <double> & b) {
  int ii(-1);
  for (int i=0;i<n;++i) {
    const int ip=indx[i];
    double sum(b[ip]);
    b[ip]=b[i];
    if (ii>=0)
      for (int j=ii;j<i;++j) sum -= a(i,j)*b[j];
    else 
      if (sum!=0.0) ii=i;
    b[i]=sum;
  }
  
  for (int i=n-1;i>=0;--i) {
    double sum(b[i]);
    for (int j=i+1;j<n;++j) sum -= a(i,j)*b[j];
    b[i]=sum/a(i,i);
  }
}

// Invert matrix a, but do not change a
// Inverse in a1
void InvertMatrix(const Array2 <double> & a, Array2 <double> & a1, const int n) {
  Array2 <double>& temp(tmp2);
  temp.Resize(n,n);
  Array1 <int>& indx(itmp1);
  indx.Resize(n);
  double d;

  // a(i,j) first index i is row index (convention)
  // elements of column vectors are stored contiguous in memory in C style arrays
  // a(i) refers to a column vector

  // calculate the inverse of the transposed matrix because this
  // allows to pass a column vector to lubksb() instead of a row

  // put the transposed matrix in temp
  for(int i=0;i<n;++i) {
    for(int j=0;j<n;++j) {
      temp(i,j)=a(j,i);
      a1(i,j)=0.0;
    }
    a1(i,i)=1.0;
  }

  ludcmp(temp,n,indx,d);
  
  for(int j=0;j<n;++j) {
    // get column vector
    Array1 <double> yy(a1,j);
    lubksb(temp,n,indx,yy);
  }
}

// Solve matrix equation A x = b for x
// Destroy contents of A
// Return the solution for x in the b argument
void SolveMatrixEquation(Array2 <double> & a, Array1 <double> & b, const int n) {
  Array1 <int>& indx(itmp1);
  indx.Resize(n);
  double d;
  ludcmp(a,n,indx,d);
  lubksb(a,n,indx,b);

//   cout << "ND= ";
//   for(int j=0;j<n;++j) {
//     cout <<  a(j,j) << "  " ;
//   }
//   cout << endl;
}

// Solve matrix equation A x = b for x
// Destroy contents of A
// Return the solution in x
void SolveMatrixEquation(Array2 <double> & a, const Array1 <double> & b, Array1 <double> & x, const int n) {
  for(int i=0; i<n; i++) 
    x[i]=b[i];
  SolveMatrixEquation(a,x,n);
}

// Assumes A is a diagponal matrix 
void SolveDiagonalMatrixEquation(Array2 <double> & a, const Array1 <double> & b, Array1 <double> & x, const int n) {
  for(int i=0; i<n; i++) 
    x[i]=b[i]/a(i,i);
}

// Assumes A is a diagponal matrix 
void SolveDiagonalMatrixEquation(Array2 <double> & a, Array1 <double> & x, const int n) {
  for(int i=0; i<n; i++) 
    x[i] /=a(i,i);
}

// Solve matrix equation A x = b for x
// Keep contents of A
// Return the solution in x
void SolveMatrixEquationPreserveA(const Array2 <double> & a, const Array1 <double> & b, 
				  Array1 <double> & x, const int n) {
  Array2 <double>& temp(tmp2);
  temp.Resize(n,n);
  for(int i=0;i<n;++i) {
    for(int j=0;j<n;++j) {
      temp(i,j)=a(i,j);
    }
  }
  
  for(int i=0; i<n; i++) 
    x[i]=b[i];
  SolveMatrixEquation(temp,x,n);
}

// Calculate the transpose inverse of matrix a
// and return the determinante
double TransposeInverseMatrix(const Array2 <double> & a, Array2 <double> & a1, const int n) {
  Array2 <double>& temp(tmp2);
  temp.Resize(n,n);
  Array1 <int>& indx(itmp1);
  indx.Resize(n);
  double d;

  // a(i,j) first index i is row index (convention)
  // elements of column vectors are stored contiguous in memory in C style arrays
  // a(i) refers to a column vector

  // calculate the inverse of the transposed matrix because this
  // allows to pass a column vector to lubksb() instead of a row

  // put the transposed matrix in temp
  for(int i=0;i<n;++i) {
    for(int j=0;j<n;++j) {
      temp(i,j)=a(i,j);
      a1(i,j)=0.0;
    }
    a1(i,i)=1.0;
  }

  ludcmp(temp,n,indx,d);
  
  for(int j=0;j<n;++j) {
    // get column vector
    Array1 <double> yy(a1,j);
    lubksb(temp,n,indx,yy);
  }

  // return the determinante as well
  for(int j=0;j<n;++j) {
    d *= temp(j,j);
  }
  return d;
}

// Return det|a| and leave matrix a constant
double Determinant(const Array2 <double> & a, const int n) {
  Array2 <double>& temp(tmp2);
  temp.Resize(n,n);
  Array1 <int>& indx(itmp1);
  indx.Resize(n);
  double d;

  // a(i,j) first index i is row index (convention)
  // elements of column vectors are stored contiguous in memory in C style arrays
  // a(i) refers to a column vector

  // calculate the inverse of the transposed matrix because this
  // allows to pass a column vector to lubksb() instead of a row

  // put the transposed matrix in temp
  for(int i=0;i<n;++i) {
    for(int j=0;j<n;++j) {
      temp(i,j)=a(i,j);
    }
  }

  ludcmp(temp,n,indx,d);
  
  for(int j=0;j<n;++j) {
    d *= temp(j,j);
  }
  return d;
}

// Return log(|det(A)|), sign=sgn(det(A))
// sign<0 indicates that det(A)<0, which normally would have cuased an error 
// when calculating log(det(A))
// matrix A gets destroyed
void LogDeterminant(Array2 <double> & a, const int n, double & logD, double & sign) {
  Array1 <int>& indx(itmp1);
  indx.Resize(n);
  double d;

  ludcmp(a,n,indx,d);
  if (d<0.0) {
    sign = -1.0;
    d    *=-1.0;
  } else {
    sign = 1.0;
  }

  logD = log(d);
  
  for(int j=0;j<n;++j) {
    d = a(j,j);
    if (d<0) {
      sign *= -1.0;
      d    *= -1.0;
    }
    logD += log(d);
  }
}

// transpose matrix a
void TransposeMatrix(Array2 <double> & a, const int n) {
  for(int i=0;i<n;++i) {
    for(int j=i+1;j<n;++j) {
      double & r(a(i,j));
      double & rt(a(j,i));
      double x=r; 
      r=rt;
      rt=x;
    }
  }
}

// multiply c=a*b
void MultiplyMatrices(const Array2 <double> & a, const Array2 <double> & b,
		      Array2 <double> & c, const int n) {
  for(int i=0; i<n; ++i) {
    for(int j=0; j<n; ++j) {
      double & cc(c(i,j));
      cc=0.0;
      for(int k=0; k<n; ++k) {
	//	if(i==0 && j==1 && fabs(a(i,k)*b(k,j))>1e-10) {
	//	  Write5("01 contrib=",k,a(i,k),b(k,j),a(i,k)*b(k,j));
	//	}
	cc +=a(i,k)*b(k,j);
      }
    }
  }
}

// multiply b=A*x
void MatrixTimesVector(const Array2 <double> & a, const Array1 <double> & x,
		       Array1 <double> & b, int n) {
  for(int i=0;i<n;++i) {
    const Array1 <double> ai (a,i);
    b[i]=0.0;
    for(int j=0;j<n;++j) {
      b[i] += ai[j]*x[j];
      //      b[i] += a(i,j)*x[j];
    }
  }
}

// multiply b=Atrans*x
void TransposeMatrixTimesVector(const Array2 <double> & a, const Array1 <double> & x,
				Array1 <double> & b, int n) {
  for(int i=0; i<n; ++i) {
    b[i] =0.0;
  }
  for(int j=0; j<n; ++j) {
    const Array1 <double> aj (a,j);
    for(int i=0; i<n; ++i) {
      //      b[i] += aj[i]*x[i];
      b[i] += a(j,i)*x[j];
    }
  }
}

//////////////////////// Update Inverse /////////////////////


// Update inverse a1 after row in matrix a has changed
// get new row  out of matrix a
//
// a1 = old inverse
// a  = new matrix with new row lRow
// returns Det(a_old)/Det(a_new)
double InverseUpdateRow(Array2 <double> & a1, const Array2 <double> & a, 
			const int lRow, const int n) {
  Array1 <double> & tmpColL(tmp11);
  tmpColL.Resize(n);
  Array1 <double> & prod(tmp12);
  prod.Resize(n);

  double f=0.0;
  for(int i=0;i<n;++i) {
    f += a(lRow,i)*a1(i,lRow);
  }
  f =-1.0/f;

  for(int j=0;j<n;++j) {
    prod[j]   =0.0;
    tmpColL[j]=a1(j,lRow);
    for(int i=0;i<n;++i) {
      prod[j] += a(lRow,i)*a1(i,j);
    }
    prod[j] *= f;
  }

  for(int ii=0;ii<n;++ii) {
    double & t(tmpColL[ii]);
    for(int j=0;j<n;++j) {
      a1(ii,j) += t*prod[j];
    }
  }

  f = -f;
  for(int i=0;i<n;++i) {
    a1(i,lRow) = f*tmpColL[i];
  }
  return f; 
}
  
// Update inverse a1 after row in matrix a has changed
// get new row out of array1 newRow
//
// a1 = old inverse
// newRow  = new row lRow in new matrix a_new
// returns Det(a_old)/Det(a_new)
double InverseUpdateRow(Array2 <double> & a1, const Array1 <double> & newRow, 
			const int lRow, const int n) {
  Array1 <double> & tmpColL(tmp11);
  tmpColL.Resize(n);
  Array1 <double> & prod(tmp12);
  prod.Resize(n);

  double f=0.0;
  for(int i=0;i<n;++i) {
    f += newRow[i]*a1(i,lRow);
  }
  f =-1.0/f;

  for(int j=0;j<n;++j) {
    prod[j]   =0.0;
    tmpColL[j]=a1(j,lRow);
    for(int i=0;i<n;++i) {
      prod[j] += newRow[i]*a1(i,j);
    }
    prod[j] *= f;
  }

  for(int i=0;i<n;++i) {
    double & t(tmpColL[i]);
    for(int j=0;j<n;++j) {
      a1(i,j) += t*prod[j];
    }
  }

  f = -f;
  for(int i=0;i<n;++i) {
    a1(i,lRow) = f*tmpColL[i];
  }
  return f; 
}
  

// Update inverse a1 after column in matrix a has changed
// get new column out of matrix a
//
// a1= old inverse
// a = new matrix with new column lCol
// returns Det(a_old)/Det(a_new)
double InverseUpdateColumn(Array2 <double> & a1, const Array2 <double> & a, 
			   const int lCol, const int n) {
  Array1 <double> & tmpColL(tmp11);
  tmpColL.Resize(n);
  Array1 <double> & prod(tmp12);
  prod.Resize(n);

  double f=0.0;
  for(int i=0;i<n;++i) {
    f += a1(lCol,i)*a(i,lCol);
  }
  f =-1.0/f;

  for(int j=0;j<n;++j) {
    tmpColL[j]=a1(lCol,j);
    prod[j]   =0.0;
    for(int i=0;i<n;++i) {
      prod[j] += a1(j,i)*a(i,lCol);
    }
    prod[j] *= f;
  }

  for(int i=0;i<n;++i) {
    double & p(prod[i]);
    for(int j=0;j<n;++j) {
      a1(i,j) += tmpColL[j]*p;
    }
  }

  f = -f;
  for(int j=0;j<n;++j) {
    a1(lCol,j) = f*tmpColL[j];
  }
  return f;
}
  
// Update inverse a1 after column in matrix a has changed
// get new column out of array1 newCol
//
// a1= old inverse
// newCol = new column lCol in the new matrix a_new
// returns Det(a_old)/Det(a_new)
double InverseUpdateColumn(Array2 <double> & a1, const Array1 <double> & newCol, 
			   const int lCol, const int n) {
  Array1 <double> & tmpColL(tmp11);
  tmpColL.Resize(n);
  Array1 <double> & prod(tmp12);
  prod.Resize(n);

  double f=0.0;
  for(int i=0;i<n;++i) {
    f += a1(lCol,i)*newCol[i];
  }
  f =-1.0/f;

  for(int j=0;j<n;++j) {
    tmpColL[j]=a1(lCol,j);
    prod[j]   =0.0;
    for(int i=0;i<n;++i) {
      prod[j] += a1(j,i)*newCol[i];
    }
    prod[j] *= f;
  }

  for(int i=0;i<n;++i) {
    double & p(prod[i]);
    for(int j=0;j<n;++j) {
      a1(i,j) += tmpColL[j]*p;
    }
  }

  f = -f;
  for(int j=0;j<n;++j) {
    a1(lCol,j) = f*tmpColL[j];
  }
  return f;
}
  
///////////////// Update Transpose Inverse /////////////////////


// Update transpose inverse a1 after row in matrix a has changed
// get new row out of matrix a
// (This is actually the modified routine InverseUpdateColumn)
//
// a1= old inverse
// a = new matrix with new column lCol
// returns Det(a_old)/Det(a_new)
double TransposeInverseUpdateRow(Array2 <double> & a1, const Array2 <double> & a, 
				 const int lCol, const int n) {
  Array1 <double> & tmpColL(tmp11);
  tmpColL.Resize(n);
  Array1 <double> & prod(tmp12);
  prod.Resize(n);

  double f=0.0;
  for(int i=0;i<n;++i) {
    f += a1(lCol,i)*a(lCol,i);
  }
  f =-1.0/f;

  for(int j=0;j<n;++j) {
    tmpColL[j]=a1(lCol,j);
    prod[j]   =0.0;
    for(int i=0;i<n;++i) {
      prod[j] += a1(j,i)*a(lCol,i);
    }
    prod[j] *= f;
  }

  for(int i=0;i<n;++i) {
    double & p(prod[i]);
    for(int j=0;j<n;++j) {
      a1(i,j) += tmpColL[j]*p;
    }
  }

  f = -f;
  for(int j=0;j<n;++j) {
    a1(lCol,j) = f*tmpColL[j];
  }
  return f;
}

// Update inverse a1 after column in matrix a has changed
// get new column out of matrix a
// (This is actually the modified routine InverseUpdateRow)
//
// a1 = old inverse
// a  = new matrix with new row lRow
// returns Det(a_old)/Det(a_new)
double TransposeInverseUpdateColumn(Array2 <double> & a1, const Array2 <double> & a, 
				    const int lRow, const int n) {
  Array1 <double> & tmpColL(tmp11);
  tmpColL.Resize(n);
  Array1 <double> & prod(tmp12);
  prod.Resize(n);

  double f=0.0;
  for(int i=0;i<n;++i) {
    f += a(i,lRow)*a1(i,lRow);
  }
  f =-1.0/f;

  for(int j=0;j<n;++j) {
    prod[j]   =0.0;
    tmpColL[j]=a1(j,lRow);
    for(int i=0;i<n;++i) {
      prod[j] += a(i,lRow)*a1(i,j);
    }
    prod[j] *= f;
  }

  for(int i=0;i<n;++i) {
    double & t(tmpColL[i]);
    for(int j=0;j<n;++j) {
      a1(i,j) += t*prod[j];
    }
  }

  f = -f;
  for(int i=0;i<n;++i) {
    a1(i,lRow) = f*tmpColL[i];
  }
  return f; 
}

// #define PRINT_NON_ZERO
#ifdef PRINT_NON_ZERO
double limit=1e-10;
void WriteMatrix(const Array2 <double> & a, int n1=-1, int n2=-1) {
  if (n1<0 || n1>a.dim[0]) n1=a.dim[0];
  if (n2<0 || n2>a.dim[1]) n2=a.dim[1];
  for(int i=0; i<n1; i++) {
    for(int j=0; j<n2; j++) {
      if (abs(a(i,j))>limit) Write3(i,j,a(i,j));
    }
  }
}

void WriteVector(const Array1 <double> & a, int n=-1) {
  if (n<0 || n>a.size) n=a.size;
  for(int i=0; i<n; i++) {
    if (abs(a[i])>limit) Write2(i,a[i]);
  }
}
#else 
void WriteDouble(const double a) {
  const ios::fmtflags current = cout.flags();
  cout.setf(ios::scientific);
  cout.setf(ios::showpos);
  cout << a << " ";
  cout.unsetf(ios::scientific);
  cout.unsetf(ios::showpos);
  cout.setf(current);
}

  
void WriteMatrix(const Array2 <double> & a, int n1 =-1, int n2 =-1) {
  if (n1<0 || n1>a.dim[0]) n1=a.dim[0];
  if (n2<0 || n2>a.dim[1]) n2=a.dim[1];
  for(int i=0; i<n1; i++) {
    cout << i << " ( ";
    for(int j=0; j<n2; j++) {
      WriteDouble(a(i,j));
    }
    cout << " )" << endl;
  }
  cout << endl;
}

void WriteVector(const Array1 <double> & a, int n=-1) {
  if (n<0 || n>a.size) n=a.size;
  cout << " [ ";
  for(int i=0; i<n; i++) {
    WriteDouble(a[i]);
  }
  cout << " ]" << endl;
}
#endif

void SaveVector(const string & filename, const Array1 <double> & a) {
  FILE * file = fopen(filename.c_str(),"w");
  if (file==0) error("Could not open file for writing",filename);

  int nWrite,n;
  n=a.Size();
  nWrite = fwrite(&n,sizeof(int),1,file);
  if (nWrite!=1) error("Could not write one integer");

  nWrite = fwrite(a.v,sizeof(double),n,file);
  if (nWrite!=n) error("Could not write double vector",nWrite,n);
  fclose(file);
}

void SaveMatrix(const string & filename, const Array2 <double> & a) {
  FILE * file = fopen(filename.c_str(),"w");
  if (file==0) error("Could not open file for writing",filename);

  int nWrite,n1,n2;
  n1=a.dim[0];
  n2=a.dim[1];
  nWrite = fwrite(&n1,sizeof(int),1,file);
  if (nWrite!=1) error("Could not write one integer");
  nWrite = fwrite(&n2,sizeof(int),1,file);
  if (nWrite!=1) error("Could not write one integer");

  nWrite = fwrite(a.v,sizeof(double),n1*n2,file);
  if (nWrite!=n1*n2) error("Could not write double array",nWrite,n1*n2);
  fclose(file);
}

void LoadVector(const string & filename, Array1 <double> & a) {
  FILE * file = fopen(filename.c_str(),"r");
  if (file==0) error("Could not open file to read",filename);
  
  int nRead,n;
  nRead = fread(&n,sizeof(int),1,file);
  if (nRead!=1) error("Could not read one integer (1)");

  a.Resize(n);
  nRead = fread(a.v,sizeof(double),n,file);
  if (nRead!=n) error("Could not read double vector",nRead,n);

  fclose(file);
}

void LoadMatrix(const string & filename, Array2 <double> & a) {
  FILE * file = fopen(filename.c_str(),"r");
  if (file==0) error("Could not open file to read",filename);
  
  int nRead,n1,n2;
  nRead = fread(&n1,sizeof(int),1,file);
  if (nRead!=1) error("Could not read one integer (1)");
  nRead = fread(&n2,sizeof(int),1,file);
  if (nRead!=1) error("Could not read one integer (2)");

  a.Resize(n1,n2);
  nRead = fread(a.v,sizeof(double),n1*n2,file);
  if (nRead!=n1*n2) error("Could not read double array",nRead,n1*n2);

  fclose(file);
}

Array2<double> IdentityMatrix(const int nF) {
  Array2 <double> a(nF,nF,0.0);
  for(int i=0; i<nF; i++) a(i,i)=1.0;
  return a;
}

void IdentityMatrix(Array2 <double> & a) {
  a = 0.0;
  for(int i=0; i<min(a.dim[0],a.dim[1]); i++) a(i,i)=1.0;
}

/*
int main() {
  const int n=4;
  Array2 <double> a(n,n),a1(n,n),one(n,n);
  for(int i=0;i<n;++i) {
    for(int j=0;j<n;++j) {
      a(i,j)=sin((i*2+1)*(j+1)*0.314);
    }
  }
  WriteMatrix(a);
  InvertMatrix(a,a1,n);
  WriteMatrix(a1);
  MultiplyMatrices(a,a1,one,n);
  WriteMatrix(one);

  cout << "---- Solve matrix equation ----\n\n";
  
  Array1 <double> b(n),x(n);
  for(int j=0;j<n;++j) {
    x[j]=0.0;
    b[j]=cos((j+1)*2.45);
  }
  WriteVector(b);
  SolveMatrixEquationPreserveA(a,b,x,n);
  WriteVector(x);
  MatrixTimesVector(a,x,b,n);
  WriteVector(b);

  exit(1);

  cout << " Changing row 1\n";
  int i=1;
  for(int j=0;j<n;++j) {
    a(i,j)+=sin((i*2+1)*(j+1)*0.111);
  }

  InverseUpdateRow(a1,a,i,n);
  MultiplyMatrices(a,a1,one,n);
  WriteMatrix(one);

  cout << "Changing column 0\n";
  i=0;
  for(int j=0;j<n;++j) {
    a(j,i)+=sin((i*2+1)*(j+1)*0.234);
  }
  InverseUpdateColumn(a1,a,i,n);
  MultiplyMatrices(a,a1,one,n);
  WriteMatrix(one);

  cout << "-----Transposed inverse -----\n\n";
  TransposeInverseMatrix(a,a1,n);
  TransposeMatrix(a,n);
  MultiplyMatrices(a,a1,one,n);
  WriteMatrix(one);
  TransposeMatrix(a,n);

  cout << " Changing row 1\n";
  i=1;
  for(int j=0;j<n;++j) {
    a(i,j)+=sin((i*2+1)*(j+1)*0.696);
  }
  TransposeInverseUpdateRow(a1,a,i,n);
  TransposeMatrix(a,n);
  MultiplyMatrices(a,a1,one,n);
  WriteMatrix(one);
  TransposeMatrix(a,n);

  cout << " Changing column 0\n";
  i=0;
  for(int j=0;j<n;++j) {
    a(j,i)+=sin((i*2+1)*(j+1)*0.5757);
  }
  TransposeInverseUpdateColumn(a1,a,i,n);
  TransposeMatrix(a,n);
  MultiplyMatrices(a,a1,one,n);
  WriteMatrix(one);
  TransposeMatrix(a,n);

  return 0;
}
*/
