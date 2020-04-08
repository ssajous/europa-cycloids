/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Contains a short vector of fixed length                                 //
//                                                                         //
// Burkhard Militzer                               Urbana 4-1-99           //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _SHORTVECTOR_
#define _SHORTVECTOR_

#include "Standard.h"

template <class T, int n> 
class ShortVector {
 public:
  T v[n];
  
  // Create without initialization
  ShortVector(){};

  // construction & initialization with one agrument is below
  inline ShortVector(const T & t1, const T & t2) {
    // this is not nice but I did not get the partical specialization <T,2> to work.
    if (n!=2) error("Template constructor problem");
    v[0]=t1;
    v[1]=t2;
  }

  inline ShortVector(const T & t1, const T & t2, const T & t3) {
    // this is not nice but I did not get the partical specialization <T,3> to work.
    if (n!=3) error("Template constructor problem");
    v[0]=t1;
    v[1]=t2;
    v[2]=t3;
  }

  // Create and make a copy 
  ShortVector(const ShortVector& x) {
    for(int i=0; i<n; i++) {
      v[i]=x[i];
    }
  }

  // Create and make a copy 
  template<class X> 
  ShortVector(const ShortVector<X,n> & x) {
    for(int i=0; i<n; i++) {
      v[i]=T(x[i]);
    }
  }

  // Creat and fill with one element
  //  template<class X> 
  //  ShortVector(const X & vv) {
  //    for(int i=0; i<n; i++) {
  //      v[i]=T(vv);
  //    }
  //  }
  // Creat and fill with one element
  ShortVector(const T & vv) {
    for(int i=0; i<n; i++) {
      v[i]=vv;
    }
  }
  
  T Norm2() const {
    T x2=T(0);
    for(int i=0; i<n; i++) {
      x2 += v[i]*v[i];
    }
    return x2;
  }

  double Norm() const {
    return sqrt(Norm2());
  }

  void Normalize(const double nn=1.0) {
    *this *= nn/Norm();
  }

  T Sum() const {
    T x=T(0);
    for(int i=0; i<n; i++) {
      x += v[i];
    }
    return x;
  }

  T AbsSum() const {
    T x=T(0);
    for(int i=0; i<n; i++) {
      x += abs(v[i]);
    }
    return x;
  }

  const T& operator() (const int i) const {     // Read element
    Limits(i,n);
    return v[i];
  };

  const T& operator[] (const int i) const {     // Read element
    Limits(i,n);
    return v[i];
  };

  T& operator() (const int i) {                 // Access element
    Limits(i,n);    return v[i];
  };

  T& operator[] (const int i) {                 // Access element
    Limits(i,n);
    return v[i];
  };

  ShortVector & operator=(const ShortVector & x) {      // Copy into predefined array
    for(int i=0; i<n; ++i) v[i]=x.v[i];
    return *this;
  }

  inline friend ShortVector operator-(const ShortVector & x) {     
    ShortVector c;
    for(int i=0; i<n; ++i) c.v[i]=-x.v[i];
    return c;
  }

  ShortVector & operator=(const T & x) {           // Fill with one element
    for(int i=0; i<n; ++i) v[i]=x;
    return *this;
  }

  ShortVector & operator+=(const T & x) {           // Add x to all elements
    for(int i=0; i<n; ++i) v[i]+=x;
    return *this;
  }

  ShortVector & operator+=(const ShortVector & x) {           // Add x to all elements
    for(int i=0; i<n; ++i) v[i]+=x.v[i];
    return *this;
  }

  ShortVector & operator-=(const T & x) {           // Subtract x from all elements
    for(int i=0; i<n; ++i) v[i]-=x;
    return *this;
  }

  ShortVector & operator-=(const ShortVector & x) {           // Subtract x from all elements
    for(int i=0; i<n; ++i) v[i]-=x.v[i];
    return *this;
  }

  ShortVector & operator*=(const T & x) {           // Multiply all elements by x
    for(int i=0; i<n; ++i) v[i]*=x;
    return *this;
  }

  /*
  ShortVector & operator*=(const ShortVector & x) {           // Multiply all elements by x
    for(int i=0; i<n; ++i) v[i]*=x.v[i];
    return *this;
  }
  */

  ShortVector & operator/=(const T & x) {           // Divide all elements by x
    for(int i=0; i<n; ++i) v[i]/=x;
    return *this;
  }

  /*
  ShortVector & operator/=(const ShortVector & x) {           // Divide all elements by x
    for(int i=0; i<n; ++i) v[i]/=x.v[i];
    return *this;
  }
  */

  bool operator==(const ShortVector & x) const {      // Are Array identical?
    for(int i=0; i<n; ++i) {
      if (v[i]!=x.v[i]) return false;
    }
    return true;
  }

  bool operator!=(const ShortVector & a) const {      // Are Array identical?
    return !(*this==a);
  }

  int Size() const { return n;}

  // Return the point but do not give up the ownership
  T* GetPointer() const {
    return v;
  }

/*
  template<class X> 
  inline friend ShortVector operator+(const ShortVector & a, const ShortVector<X,n> & b) { 
    ShortVector <T,n> c;
    for(int i=0; i<n; ++i) {
      c[i] = a[i]+T(b[i]);
    }
    return c;
  }
*/

  inline friend ShortVector operator+(const ShortVector & a, const ShortVector & b) { 
    ShortVector c;
    for(int i=0; i<n; ++i) {
      c[i] = a[i]+b[i];
    }
    return c;
  }

  /*
  template<class X> 
  inline friend ShortVector operator-(const ShortVector & a, const ShortVector<X,n> & b) { 
    ShortVector <T,n> c;
    for(int i=0; i<n; ++i) {
      c[i] = a[i]-T(b[i]);
    }
    return c;
  }
  */

  inline friend ShortVector operator-(const ShortVector & a, const ShortVector & b) { 
    ShortVector c;
    for(int i=0; i<n; ++i) {
      c[i] = a[i]-b[i];
    }
    return c;
  }

  // scalar product
  inline friend T operator*(const ShortVector & a, const ShortVector & b) { 
    T s=T(0);
    for(int i=0; i<n; ++i) {
      s += a[i]*b[i];
    }
    return s;
  }

  inline friend ShortVector operator*(const ShortVector & a, const T & x) { 
    ShortVector <T,n> c;
    for(int i=0; i<n; ++i) {
      c[i] = a[i]*x;
    }
    return c;
  }

  inline friend ShortVector operator*(const T & x, const ShortVector & a) { 
    ShortVector <T,n> c;
    for(int i=0; i<n; ++i) {
      c[i] = a[i]*x;
    }
    return c;
  }

  inline friend ShortVector operator/(const ShortVector & a, const T & x) { 
    ShortVector <T,n> c;
    for(int i=0; i<n; ++i) {
      c[i] = a[i]/x;
    }
    return c;
  }
  friend ShortVector operator/(const T & x, const ShortVector & a) { // disabled
    error("scalar/vector calculation is undefined.");
  }

  friend ostream& operator<<(ostream & os, const ShortVector & x) {
    os << "(";
    for(int i=0; i<n; ++i) {
      if (i>0) 
	os << ",";
      os << x.v[i];
    }
    os << ")";
    return os;
  }

  void CycleUp() {
    T x=v[n-1];
    for(int i=n-1; i>0; i--) {
      v[i]=v[i-1];
    }
    v[0]=x;
  }

  void CycleUp(const int m) {
    for(int i=0; i<m; i++) {
      CycleUp();
    }
  }

  // move elements down, [i+1] goes to [i] and [0] goes to the end
  void CycleDown() {
    T x=v[0];
    for(int i=0; i<n-1; ++i) {
      v[i]=v[i+1];
    }
    v[n-1]=x;
  }

  void CycleDown(const int m) {
    for(int i=0; i<m; i++) {
      CycleDown();
    }
  }

};

template <>
inline double ShortVector<double,1>::Norm2() const {
  return v[0]*v[0];
}

template <>
inline double ShortVector<double,2>::Norm2() const {
  return v[0]*v[0]+v[1]*v[1];
}

template <>
inline double ShortVector<double,3>::Norm2() const {
  return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
}

#endif // _SHORTVECTOR_
