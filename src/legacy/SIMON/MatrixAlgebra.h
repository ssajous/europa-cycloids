// Matrix operations
//
// Burkhard Militzer                                    Urbana 4-1-99
//
#ifndef _MATRIXALGEBRA_
#define _MATRIXALGEBRA_

#include "Array.h"

void ludcmp(Array2 <double> & a, const int n, Array1 <int> & indx, double & d);
void lubksb(const Array2 <double> & a, const int n, const Array1 <int> & indx, Array1 <double> & b);

double Determinant(const Array2 <double> & a, const int n);
void LogDeterminant(Array2 <double> & a, const int n, double & logD, double & sign);

void MultiplyMatrices(const Array2 <double> & a, const Array2 <double> & b,
		      Array2 <double> & c, int n);

void InvertMatrix(const Array2 <double> & a, Array2 <double> & a1, const int n);
double InverseUpdateRow(Array2 <double> & a1, const Array2 <double> & a, 
			const int lRow, const int n);
double InverseUpdateRow(Array2 <double> & a1, const Array1 <double> & newRow, 
			const int lRow, const int n);
double InverseUpdateColumn(Array2 <double> & a1, const Array2 <double> & a, 
			   const int lCol, const int n);
double InverseUpdateColumn(Array2 <double> & a1, const Array1 <double> & newCol, 
			   const int lCol, const int n);
 
void TransposeMatrix(Array2 <double> & a, const int n);
double TransposeInverseMatrix(const Array2 <double> & a, Array2 <double> & a1, const int n);
double TransposeInverseUpdateRow(Array2 <double> & a1, const Array2 <double> & a, 
				 const int lCol, const int n);
double TransposeInverseUpdateColumn(Array2 <double> & a1, const Array2 <double> & a, 
				    const int lRow, const int n);
void SolveMatrixEquation(Array2 <double> & a, Array1 <double> & b, const int n);
void SolveMatrixEquation(Array2 <double> & a, const Array1 <double> & b, Array1 <double> & x, const int n);
void SolveDiagonalMatrixEquation(Array2 <double> & a, const Array1 <double> & b, Array1 <double> & x, const int n);
void SolveDiagonalMatrixEquation(Array2 <double> & a, Array1 <double> & x, const int n);
void SolveMatrixEquationPreserveA(const Array2 <double> & a, const Array1 <double> & b, 
				  Array1 <double> & x, const int n);
void MatrixTimesVector(const Array2 <double> & a, const Array1 <double> & x,
		       Array1 <double> & b, int n);
void TransposeMatrixTimesVector(const Array2 <double> & a, const Array1 <double> & x,
				Array1 <double> & b, int n);

void WriteMatrix(const Array2 <double> & a, int n1= -1, int n2= -1);
void WriteVector(const Array1 <double> & a, int n= -1);
void SaveVector(const string & filename, const Array1 <double> & a);
void LoadVector(const string & filename, Array1 <double> & a);
void SaveMatrix(const string & filename, const Array2 <double> & a);
void LoadMatrix(const string & filename, Array2 <double> & a);
Array2<double> IdentityMatrix(const int nF);
void IdentityMatrix(Array2 <double> & a);

#endif // _MATRIXALGEBRA_
