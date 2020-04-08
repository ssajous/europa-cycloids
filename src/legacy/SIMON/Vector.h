/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Vector with range checking                                              //
//                                                                         //
// Burkhard Militzer                                    Urbana 4-1-99      //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _RANGE_CHECKING_VECTOR_
#define _RANGE_CHECKING_VECTOR_

#include <vector>
#include "Standard.h"

template <class T>
class Vector : public vector <T> {
  void Limits(int n) const {
    if (n<0 || n>=int(this->size()) )
      error("Array Error:  0<= ",n," < ",this->size());
  }
 public:
  Vector(){}
  Vector(int i):vector <T>(i){}  
  Vector(const Vector <T> & v ):vector <T>(v){}
  const T& operator[](int i) const {
    Limits(i); 
    const vector<T>& tmp=*this; 
    return tmp[i]; 
  }
  T& operator[](int i) {
    Limits(i); 
    vector<T>& tmp=*this; 
    return tmp[i]; 
  }
};

#endif // _RANGE_CHECKING_VECTOR_
