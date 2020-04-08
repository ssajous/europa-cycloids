/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Array class                                                             //
//                                                                         //
// Burkhard Militzer                                    Urbana 4-1-99      //
// Updated later                                   Livermore 01-11-03      //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _ARRAY_
#define _ARRAY_

#include "Standard.h"

// Define SLOW_WARNING to warn if arrays are copied or [][] are used.
// #define SLOW_WARNING

// Array1, Array2, and Array3 are template classes for 1, 2, and 3
// dimensional arrays storing elements contiguously in memory.  Array
// boundary checking as well as resizing functions are implemented.

// Array1 can accessed by () or []. [] make it somewhat
// compatible with STL's <vector> class but does not show its
// performance problems.

// If you access Array2 a[i][j] then a temporary Array1 must be
// generated. It works but is slower than a(i,j)

// We could define a virtual base class Array providing a universal access to rank
// and GetDims(). However, we have not had much use for it.
// RArray get be read from a stream e.g. cin >> mat but Array cannot because
// many template elements do not provide stream handling.

// Constraints of this class:
// (1) Reading must be as fast as possible: 
//     a(i)     <- a.v+i                and not "a.v+i+shift"
//     a(i,j)   <- a.v+i*step+j         and not "a.v+i*step+j*step2"         nor "a[i][j]"
//     a(i,j,k) <- a.v+i*step+j*step2+k and not "a.v+i*step+j*step2+k*step3" nor "a[i][j][k]"
// (2) Can have arrays that own and some that do not own the memory, see "bool b".

// Functions that work for all arrays:
//   Copy constructor
//   Constructor with given dimensions
//   Constructor with given dimensions and initialization
//   Constructor with given dimensions and external pointer (do not deallocate memory later)
//   Read and write elements (...)
//   CopyElements(const Array & a)
//   Assignment operator=(const Array & a) A=B
//   Assignment operator=(const T & x)     A=x
//   A += x
//   A += B
//   A -= x
//   A -= B
//   A *= x
//   A /= x
//   A == B
//   A != B
//   Max(), Min(), Sum()
//   Resize(...)
//   Resize(...,mn)
//   CopyResize(...)
//   CopyResize(...,mn)
//   Free()
//   SetPointer(v,...,bool b=true);
//   T* GetPointer()
//   Check() 
template <class T> class Array1;
template <class T> class Array2;
template <class T> class Array3;

const int stInc=10; // Default increment for PushBack
const int notFound = -1;

template <class T> class Array1 {
 public:
  int size;      // logical boundary of array: 0<= index < size
  int mSize;     // physical boundary of allocated memory 
  int inc;       // current increment for pushback
  T* v;          // pointer to the elements
  bool b;        // true if the objects owns the memory and mst deallocated memory at the end

  //default constructor
  Array1():
    size(0),
    mSize(0),
    inc(stInc),
    v(0),
    b(true)
    {};

  // What should never be done: Make a copy of the pointer and prevent 
  // allocation and copyingL fast but NOT SAFE because copy and original share the 
  // same memory indicate by v. Do never use this!!!
  //  Array1(const Array1& a):size(a.size),mSize(a.mSize),v(a.v),b(false) {
  //  }
  
  // Make a real copy in the standard way, slow and safe. Copy also the extra memory
  // so that Array keeps 'all' its properties when copied.
  Array1(const Array1 & a):
    size(a.size),
    mSize(a.mSize),
    inc(a.inc),
    b(a.b) {
    //    warning("do not copy Array1",a.Size());
    if (b) {
      v = new T[mSize];
      if (mSize>0 && v==0) error("Could not allocate memory",mSize);
      for(int i=0; i<size; i++)
	v[i]=a.v[i];
    } else v=a.v; // memory is neither owned by old array and nor by new array.
  }

  // Create with certain size but without initialization
  Array1(const int n):
    size(n),
    mSize(n),
    inc(stInc),
    v(new T[n]),
    b(true) {          
    if (mSize>0 && v==0) error("Could not allocate memory",n);
  }

  // Create with certain size and initialize
  Array1(const int n, const T & x):size(0),mSize(0),inc(stInc),v(0),b(true) {             
    Resize(n);                        
    for(int i=0;i<size;++i) v[i]=x;  
  }

  // Create with certain size, add some extra memory and initialize
  Array1(const int n, const int mn, const T & x):size(0),mSize(0),inc(stInc),v(0),b(true) { 
    Resize(n,mn);
    for(int i=0;i<size;++i) v[i]=x;
  }

  // Create if memory is already allocated, do not deallocate memory later
  Array1(const int n, T* vv):size(n),mSize(n),v(vv),b(false) {}

  // Create from existing bigger arrays
  // These constructors ignore const of the input array because
  // there is no way to enforce that the created array is const too.
  inline Array1(const Array1 <T> & a1, const int n1, const int n2);
  inline Array1(const Array2 <T> & a2, const int n1);
  inline Array1(const Array2 <T> & a2, const int n1, const int n21, const int n22);
  inline Array1(const Array3 <T> & a3, const int n1, const int n2);

  ~Array1() {
    if (v && b) delete[] v;
  };

  const T& operator() (int n) const {     // Read element
    Limits(n,size);
    return v[n];
  };

  const T& operator[] (int n) const {     // Read element
    Limits(n,size);
    return v[n];
  };

  T& operator() (int n) {                 // Access element
    Limits(n,size);
    return v[n];
  };

  T& operator[] (int n) {                 // Access element
    Limits(n,size);
    return v[n];
  };

  const T& First() const {                // Read first element
    return (*this)[0];
  }
  
  T& First() {                            // Access first element
    return (*this)[0];
  }
  
  const T& Second() const {                // Read 2nd element
    return (*this)[1];
  }
  
  T& Second() {                            // Access 2nd element
    return (*this)[1];
  }
  
  const T& BeforeLast() const {                 
    return (*this)[size-2];
  }
  
  T& BeforeLast() {                             
    return (*this)[size-2];
  }
  
  const T& Last() const {                 // Read last element
    return (*this)[size-1];
  }
  
  T& Last() {                             // Access last element
    return (*this)[size-1];
  }
  
  // Copy all elements and do not check size unless RANGE_CHECKING is defined
  void CopyElements(const Array1 & x) {     
    EqualLimits(size,x.size);
    for(int i=0;i<size;++i) v[i]=x.v[i];
  }

  // Make a real copy of the array like in the copy constructor, safe but slow
  Array1 & operator=(const Array1 & a) {   
    inc=a.inc;
    b  =a.b;
    if (b) {
#ifdef SLOW_WARNING
      warning("Avoid copying Array1",a.Size());
#endif
      Resize(a.size,a.mSize);
      for(int i=0;i<size;++i) v[i]=a.v[i];
    } else {
      size = a.size;
      mSize = a.mSize;
    }
    return *this;
  }

  Array1 & operator=(const T & x) {            // Fill with one value
    for(int i=0;i<size;++i) v[i]=x;
    return *this;
  }

  Array1 & operator+=(const T & x) {           // Add x to all elements
    for(int i=0;i<size;++i) v[i]+=x;
    return *this;
  }

  Array1 & operator+=(const Array1 <T> & a) {  // Add array of equal dimensions
    for(int i=0;i<size;++i) v[i]+=a[i];
    return *this;
  }

  Array1 & operator-=(const T & x) {           // Subtract x from all elements
    for(int i=0;i<size;++i) v[i]-=x;
    return *this;
  }

  Array1 & operator-=(const Array1 <T> & a) {  // Subtract array of equal dimensions   
    for(int i=0;i<size;++i) v[i]-=a[i];
    return *this;
  }

  Array1 & operator*=(const T & x) {           // Multiply all elements by x
    for(int i=0;i<size;++i) v[i]*=x;
    return *this;
  }

  Array1 & operator/=(const T & x) {           // Divide all elements by x
    for(int i=0;i<size;++i) v[i]/=x;
    return *this;
  }

  bool operator==(const Array1 <T> & a) const {      // Are Array identical?
    if (size != a.Size()) 
      return false;
    for(int i=0;i<size;++i) {
      if (v[i]!=a[i]) return false;
    }
    return true;
  }

  bool operator!=(const Array1 <T> & a) const {      // Are Array identical?
    return !(*this==a);
  }

  T Max() const {
    Limits(0,size);
    T x=v[0];
    for(int i=1; i<size; i++) {
      if (v[i]>x) x=v[i];
    }
    return x;
  }

  T Min() const {
    Limits(0,size);
    T x=v[0];
    for(int i=1; i<size; i++) {
      if (v[i]<x) x=v[i];
    }
    return x;
  }

  T Sum() const {
    T sum=0;
    for(int i=0; i<size; i++) {
      sum += v[i];
    }
    return sum;
  }

  int MaxIndex() const {
    Limits(0,size);
    int ii=0;
    for(int i=1; i<size; i++) {
      if (v[i]>v[ii]) ii=i;
    }
    return ii;
  }
  
  int MinIndex() const {
    Limits(0,size);
    int ii=0;
    for(int i=1; i<size; i++) {
      if (v[i]<v[ii]) ii=i;
    }
    return ii;
  }

  // Adjust logical size of array, only adjust physical size (and allocate memory) if necessary
  void Resize(const int n) { 
    Resize(n,max(n,mSize));
  }

  // Adjust physical and logical size, allocate memory if nm!=mSize
  // Do not copy existing elements over!!!
  void Resize(const int n, const int nm) { 
    if (!b) error ("Cannot resize fixed-size array.");
    if (nm!=mSize) {
      mSize=nm;
      if (v) delete[] v;
      v = new T[mSize];
      if (mSize>0 && v==0) error("Could not allocate memory",mSize);
    }
    size = n;
  }

  void Set(const T & x) {                   // Set all elements, equivalent to 
    for(int i=0; i<size; ++i) v[i]=x;       // Array1=x; but nothing is returned
  }

  
  // Resize and copy elements, keep extra memory
  void CopyResize(const int n) {
    CopyResize(n,max(n,mSize));
  }

  // Resize and copy elements and reserve memory for exactly nm elements (for later CopyResize)
  void CopyResize(const int n, const int nm) {
    if (!b) error ("Cannot copy and resize fixed-size array.");
    if (nm<n) error("Copy resize limit error",n,nm);
    if (nm!=mSize) { //if we use nm>mSize we would lose some memory but no memory leak.
      mSize=nm;
      T* vv = new T[mSize];
      if (mSize>0 && vv==0) error("Could not allocate memory",mSize);
      for(int i=0;i<size;i++) vv[i] = v[i];
      if (v) delete[] v;
      v = vv;
    }
    size = n;
  }

  // Ask for more memory without changing the logical size of the array yet
  void Allocate(const int nm) {
    CopyResize(size,nm);
  }
 
  // make an extra element at the end available, increase size by 'inc' if necessary
  int PushBack() {
    if (!b) error ("Cannot push back in a fixed array");
    if (size==mSize) {
      CopyResize(size+1,size+inc);
      if (inc< (1<<10)*stInc) inc*=2;
    } else size++;
    return size-1;
  }      

  // put x end the end of current array, increase size by 'inc' if necessary
  int PushBack(const T& x) {
    int i = PushBack();
    v[i]=x;
    return i;
  }      

  // Pushback every element of a
  void PushBack(const Array1 <T> & a) {
    int s=size;
    CopyResize(size+a.size);
    for(int i=0; i<a.Size(); i++) {
      v[s+i]=a[i];
    }
  }      

  // make an extra element at the end available, increase size by 'inc' if necessary
  T PopBack() {
    if (!b) error("Cannot pop back in a fixed array");
    Limits(0,size);
    //    if (size<=0) error("Cannot pop back if array is empty");
    //careful: the reduces variable size first and then still accesses the last element
    //Only if there is a size check in the overloaded v[i] there would be a porblem.
    return v[--size]; 
  }      

  void SetIncrement(const int i) {
    if (i<=0) error("Increment invalid",i);
    inc = i;
  }

  void Free() {
    if (!b) error ("Cannot free memory of a fixed array");
    if (mSize>0 && v) {
      delete[] v;
      // Must set to zero because otherwise deallocation is tried again (for egcs at least)
      v = 0;
      mSize = size = 0;
    }
  }

  int Size() const { return size;}
  int Rank() const { return 1; }

  // Assume the memory has be allocated by some other routine
  // if (b1==true) take over the ownership and deallocated the memory at the end
  void SetPointer(T* vv, const int n, const bool b1=true) {
    if (v && b) delete[] v;
    size=n;
    mSize=n;
    v=vv;
    b=b1;
  }
  // Return the point but do not give up the ownership
  T* GetPointer() const {
    return v;
  }

  int Find(const T & x) const {
    for(int i=0; i<size; i++) {
      if (v[i]==x) return i;
    }
    return notFound;
  }

  void Delete(const int i, const int n=1) {
    for(int j=i; j<size-n; j++) 
      v[j] = v[j+n];
    size -= n;
  }

  // do not allocate memory, just overlay a 2D array
  inline Array2 <T> MakeArray2(const int n1);
  inline const Array2 <T> MakeConstArray2(const int n1);
    
  void Check() const {
    if (size<0)     error("in Array dimension",size);
    if (mSize<size) error("mSize is wrong",mSize,size);
    if (inc<0)      error("inc is wrong",inc);
  }
  
  friend ostream& operator<<(ostream & os, const Array1 <T> & a) {
    os << "[" << a.size << "] (";
    for(int i=0; i<a.size; ++i) {
      if (i>0) 
	os << ",";
      os << a[i];
    }
    os << ")";
    return os;
  }

  void Insert(const int ii) {
    InsertN(ii,1);
  }

  // ii==0      means insert at the beginning. If x[0] exists it is moved up.
  // ii==size-1 means x[size-1] is moved up.
  // in==size   would mean to append nn elements BUT THIS IS NOT ALLOWED
  void InsertN(const int ii, const int nn) {
    LimitsInclusive(ii,size);
    if (!b) error ("Cannot copy and resize fixed-size array.");
    if (nn<0) error("Insert limit error",ii,nn);
    const int newSize = size+nn; 
    if (newSize>mSize) { // yes, we must reallocate 
      mSize=newSize;
      T* vv = new T[mSize];
      if (mSize>0 && vv==0) error("Could not allocate memory",mSize);
      for(int i=0;  i<ii;   i++) vv[i]    = v[i]; // this assumes ii<=size
      for(int i=ii; i<size; i++) vv[i+nn] = v[i]; // "wrong" dir in "i" works because vv and v are different
      if (v) delete[] v;
      v = vv;
    } else { // No, we have enough memory
      for(int i=size-1; i>=ii; i--) v[i+nn] = v[i];
    }      
    size = newSize;
  }

  void Insert(const T & x, const int ii) {
    InsertN(x,ii,1);
  }
  void InsertN(const T & x, const int ii, const int nn) {
    InsertN(ii,nn);
    for(int i=ii; i<ii+nn; i++) 
      v[i]=x;
  }

  void Invert() { // reorder elements
    for(int i=0; i<size/2; i++) {
      int j = size-1-i;
      T x = v[i];
      v[i] = v[j];
      v[j] = x;
    }
   }

};

template <class T> class Array2 {
 public:
  int size;       // == dim[0]*dim[1], logical memory
  int mSize;      // >=size, physical memory
  T* v;
  bool b;
  int dim[2];     // array dimensions: Array2(i,j) with i<dim[0] and j<dim[1] 
                  // only one obvious reason why ordering could not be reversed (dim[0]<->dim[1])
                  // because this definition associates "Limits(n1,dim[0]); Limits(n2,dim[1]);"
  Array2():
    size(0),
    mSize(0),
    v(0),
    b(true)
    {
      dim[0]=dim[1]=0;
    };

  // Make a real copy; slow and save in the standard way
  Array2(const Array2 & a):
    size(a.size),
    mSize(a.mSize),
    b(a.b) {
    if (b) { 
#ifdef SLOW_WARNING
      warning("Avoid copying Array2",size);
#endif
      v = new T[mSize];
      if (mSize>0 && v==0) error("Could not allocate memory",mSize);
      for(int i=0; i<size; i++)
	v[i]=a.v[i];
    } else v=a.v; // memory is neither owned by the old array nor by the new array.
    dim[0]=a.dim[0];
    dim[1]=a.dim[1];
  }

  // Create without initialization
  Array2(const int n1, const int n2):
    size(n1*n2),
    mSize(size),
    v(new T[mSize]),
    b(true) {                      
    if (mSize>0 && v==0) error("Could not allocate memory",mSize);
    dim[0] = n1; 
    dim[1] = n2; 
  }

  // Create with initialization
  Array2(const int n1, const int n2, const T & x):size(0),mSize(0),v(0),b(true) {      
    Resize(n1,n2);
    for(int i=0;i<size;++i) v[i]=x;
  }

  // Create if memory is already allocated, do not get ownership
  // do not deallocate memory later
  Array2(int n1, int n2, T* vv):size(n1*n2),mSize(size),v(vv),b(false) {
    dim[0] = n1; 
    dim[1] = n2; 
  }

  // Create Array2 from Array1, do not get ownership, do not deallocate memory later
  Array2(Array1 <T> & a1):size(a1.size),mSize(a1.size),v(a1.v),b(false) {
    dim[0] = a1.size; 
    dim[1] = 1;
  }

  // Create from existing 3d array
  inline Array2(const Array3 <T> & a3, const int n);

  ~Array2() {
    if (v && b) delete[] v;
  }

  int Size() const { return size; }
  int Rank() const { return 2; }
  int GetDim(const int i) const {
    Limits(i,Rank());
    return dim[i];
  }

  const T& operator() (const int n1, const int n2) const {  // Read element     
    Limits(n1,dim[0]); 
    Limits(n2,dim[1]);
    return v[n1*dim[1]+n2];
  };

  T& operator() (const int n1, const int n2) {              // Access element
    Limits(n1,dim[0]); 
    Limits(n2,dim[1]);
    return v[n1*dim[1]+n2];
  };

  const Array1 <T> operator[](const int n1) const {
#ifdef SLOW_WARNING
      warning("Avoid C style access Array2 a[i][j]");
#endif
    Limits(n1,dim[0]);
    Array1 <T> a(*this,n1);
    return a;
  }

  Array1 <T> operator[](const int n1) {
#ifdef SLOW_WARNING
      warning("Avoid C style access of Array2 a[i][j]");
#endif
    Limits(n1,dim[0]);
    Array1 <T> a(*this,n1);
    return a;
  }

  // Copy array like 1D in predefined array
  void CopyElements(const Array2 & a) {
    EqualLimits(a.dim[0],dim[0]); 
    EqualLimits(a.dim[1],dim[1]);
    for(int i=0;i<size;++i) v[i]=a.v[i];
  }

  // Make a real copy of the array like in the copy constructor, safe but slow
  Array2 & operator=(const Array2 & a) {
    b = a.b;
    if (b) {  // 'a' owns the memory
#ifdef SLOW_WARNING
      warning("Avoid assigning Array2",a.Size());
#endif
      Resize(a.dim[0],a.dim[1],a.mSize);
      for(int i=0; i<size; ++i) v[i]=a.v[i];
    } else { // 'a' does not own the memory
      size   = a.size;
      mSize  = a.mSize;
      dim[0] = a.dim[0];
      dim[1] = a.dim[1];
      v      = a.v;
    }
    return *this;
  }

  Array2 & operator=(const T & x) {            // Fill with one value
    for(int i=0;i<size;++i) v[i]=x;
    return *this;
  }

  Array2 & operator+=(const T & x) {           // Add x to all elements
    for(int i=0;i<size;++i) v[i]+=x;
    return *this;
  }

  Array2 & operator+=(const Array2 <T> & a) {  // Add array of equal dimensions
    for(int i=0;i<dim[0];++i) 
      for(int j=0;j<dim[1];++j) 
	(*this)(i,j) += a(i,j);
    return *this;
  }

  Array2 & operator-=(const T & x) {           // Subtract x from all elements
    for(int i=0;i<size;++i) v[i]-=x;
    return *this;
  }

  Array2 & operator-=(const Array2 <T> & a) {  // Add array of equal dimensions
    for(int i=0;i<dim[0];++i) 
      for(int j=0;j<dim[1];++j) 
	(*this)(i,j) -= a(i,j);
    return *this;
  }

  Array2 & operator*=(const T & x) {           // Multiply all elements by x
    for(int i=0;i<size;++i) v[i]*=x;
    return *this;
  }

  Array2 & operator/=(const T & x) {           // Divide all elements by x
    for(int i=0;i<size;++i) v[i]/=x;
    return *this;
  }

  bool operator==(const Array2 <T> & a) const {      // Are Array identical?
    if (dim[0]!=a.dim[0] || dim[1]!=a.dim[1])
      return false;
    for(int i=0;i<size;++i) {
      if (v[i]!=a.v[i]) return false;
    }
    return true;
  }

  bool operator!=(const Array2 <T> & a) const {      // Are Array identical?
    return !(*this==a);
  }

  T Max() const {
    Limits(0,size);
    T x=v[0];
    for(int i=1; i<size; i++) {
      if (v[i]>x) x=v[i];
    }
    return x;
  }

  T Min() const {
    Limits(0,size);
    T x=v[0];
    for(int i=1; i<size; i++) {
      if (v[i]<x) x=v[i];
    }
    return x;
  }

  T Sum() const {
    T sum=0;
    for(int i=0; i<size; i++) {
      sum += v[i];
    }
    return sum;
  }

  // Adjust logical size of array, only adjust physical size (and allocate memory) if necessary.
  // Careful, elements are NOT moved around and some of the functions above like +=Array2 assume 
  // equal dimensions and same dim[1] of the two arrays.
  void Resize(const int n1, const int n2) {
    Resize(n1,n2,max(n1*n2,mSize));
  }

  // Adjust physical and logical size, allocate memory if nm!=mSize
  void Resize(const int n1, const int n2, const int nm) {
    if (!b) error ("Cannot resize fixed-size array");
    if (nm!=mSize) {
      mSize=nm;
      if (v) delete[] v;
      v = new T[mSize];
      if (mSize>0 && v==0) error("Could not allocate memory",mSize);
    }
    size   = n1*n2;
    dim[0] = n1;
    dim[1] = n2;
  }

  // CopyResize with extra memory has not been implemented
  void CopyResize(const int n1, const int n2) {
    CopyResize(n1,n2,max(n1*n2,mSize));
  }

  // Resize and copy elements and reserve memory for exactly nm elements (for later CopyResize)
  void CopyResize(const int n1, const int n2, const int nm) {
    if (!b) error ("Cannot resize fixed-size array");
    int nn= n1*n2;
    if (nn>nm) error("Copy resize limit error",n1,n2,nm);
    if (n2!=dim[1] || nm!=mSize) {
      mSize  = nm;
      T* vv=new T[mSize];
      if (mSize>0 && vv==0) error("Could not allocate memory",mSize);
      for(int i1=0;i1<min(n1,dim[0]);++i1) {
	for(int i2=0;i2<min(n2,dim[1]);++i2) {
	  vv[i1*n2+i2] = v[i1*dim[1]+i2];
	}
      }
      if (v) delete[] v;
      v      = vv;
      dim[1] = n2;
    }
    size   = nn;
    dim[0] = n1;
  }

  void Free() {
    if (!b) error ("Cannot free memory of a fixed array");
    if (mSize>0 && v) {
      delete[] v;
      v = 0;
      mSize = size = dim[0] = dim[1] = 0;
    }
  }
  
  // Assume the memory has be allocated by some other routine
  // if (b1==true) take over the ownership and deallocated the memory at the end
  void SetPointer(T* vv, const int n1, const int n2, const bool b1=true) {
    if (v && b) delete[] v;
    int nn=n1*n2;
    size  = mSize = nn;
    v     = vv;
    b     = b1;
    dim[0]= n1;
    dim[1]= n2;
  }
  // Return the point but do not give up the ownership
  T* GetPointer() {
    return v;
  }
    
  void Check() const {
    int nn=1;
    for(int i=0; i<Rank(); i++) {
      if (dim[i]<0) error("in Array dimension",dim[i]);
      nn *= dim[i];
    }
    if (nn!=size) error("size is wrong");
    if (mSize<size) error("mSize is wrong");
  }

  friend ostream & operator<<(ostream & os, const Array2 <T> & a) {
    os << "[" << a.dim[0] << "," << a.dim[1] << "] (" << endl;
    for(int i=0; i<a.dim[0]; ++i) {
      cout << "   (";
      for(int j=0; j<a.dim[1]; ++j) {
	if (j>0) 
	  os << ",";
	os << a(i,j);
      }
      os << ")" << endl;
    }
    os << ")" << endl;
    return os;
  }

  void Insert(const int d, const int ii, const int nn=1) {
    if (d==0) { // insert a row
      CopyResize(dim[0]+nn,dim[1]); // for simplicity, copy elements twice 
      for(int i=dim[0]-nn-1; i>=ii; i--) {
	for(int j=0; j<dim[1]; j++) {
	  (*this).operator()(i+nn,j) = (*this).operator()(i,j);
	}
      }
    } else { // insert a column
      CopyResize(dim[0],dim[1]+nn); // for simplicity, copy elements twice 
      for(int i=0; i<dim[0]; i++) {
	for(int j=dim[1]-nn-1; j>=ii; j--) {
	  (*this).operator()(i,j+nn) = (*this).operator()(i,j);
	}
      }
    }		
  }

  void Insert(const T & x, const int d, const int ii, const int nn=1) {
    Insert(d,ii,nn);
    if (d==0) { 
      for(int i=ii; i<ii+nn; i++) {
	for(int j=0; j<dim[1]; j++) {
	  (*this).operator()(i,j)=x;
	}
      }
    } else {
      for(int i=0; i<dim[0]; i++) {
	for(int j=ii; j<ii+nn; j++) {
	  (*this).operator()(i,j)=x;
	}
      }
    }
  }

};

template <class T> class Array3 {
 public:
  int size;       // == dim[0]*dim[1]*dim[2], logical memory
  int mSize;      // >= size, physical memory
  T* v;
  bool b;
  int dim[3];     // array dimensions: Array3(i,j,k) with i<dim[0], j<dim[1], k<dim[2]

  Array3():
    size(0),
    mSize(0),
    v(0),
    b(true) { 
    dim[0]=dim[1]=dim[2]=0;
  }

  // Make a real copy; slow and save in the standard way
  Array3(const Array3 & a):
    size(a.size),
    mSize(a.mSize),
    b(a.b) {          
    if (b) {
#ifdef SLOW_WARNING
      warning("Avoid copying Array3",a.Size());
#endif
      v = new T[mSize];
      if (mSize>0 && v==0) error("Could not allocate memory",mSize);
      for(int i=0; i<size; i++)
	v[i]=a.v[i];
    } else v=a.v; // memory is neither owned by the old array nor by the new array.
    dim[0]=a.dim[0];
    dim[1]=a.dim[1];
    dim[2]=a.dim[2];
  }

  // Create without initialization
  Array3(const int n1, const int n2, const int n3):
    size(n1*n2*n3),
    mSize(size),
    v(new T[mSize]), 
    b(true) {
    if (mSize>0 && v==0) error("Could not allocate memory",mSize);
    dim[0] = n1; 
    dim[1] = n2; 
    dim[2] = n3;
  }

  // Create with initialization
  Array3(const int n1, const int n2, const int n3, const T & vv):size(0),mSize(size),v(0),b(true) {    
    Resize(n1,n2,n3);
    for(int i=0;i<size;++i) v[i]=vv;
  }

  // Create with initialization
  Array3(const int n1, const int n2, const int n3, T* vv):
    size(n1*n2*n3),
    mSize(size),
    v(vv),
    b(false) {    
    dim[0] = n1; 
    dim[1] = n2; 
    dim[2] = n3;
  }

  ~Array3() {
    if (b && v) delete[] v;
  } 
  
  int Size() const { return size; }
  int Rank() const { return 3; }

  const T& operator() (const int n1, const int n2, const int n3) const { // Read element
    Limits(n1,dim[0]); 
    Limits(n2,dim[1]); 
    Limits(n3,dim[2]);
    return v[(n1*dim[1]+n2)*dim[2]+n3];
  }

  T& operator() (const int n1, const int n2, const int n3) {       // Access element
    Limits(n1,dim[0]); 
    Limits(n2,dim[1]); 
    Limits(n3,dim[2]);
    return v[(n1*dim[1]+n2)*dim[2]+n3];
  }

  const Array2 <T> operator[](const int n1) const {
#ifdef SLOW_WARNING
      warning("Avoid C style access Array3 a[i][j][k]");
#endif
    Limits(n1,dim[0]);
    Array2 <T> a(*this,n1);
    return a;
  }

  Array2 <T> operator[](const int n1) {
#ifdef SLOW_WARNING
      warning("Avoid C style access Array3 a[i][j][k]");
#endif
    Limits(n1,dim[0]);
    Array2 <T> a(*this,n1);
    return a;
  }

  void CopyElements(const Array3 & a) {     // Copy array like 1D
    EqualLimits(a.dim[0],dim[0]); 
    EqualLimits(a.dim[1],dim[1]);
    EqualLimits(a.dim[2],dim[2]);
    for(int i=0;i<size;++i) v[i]=a.v[i];
  }

  Array3 operator=(const Array3 & a) {     // Copy array like 1D
    b = a.b;
    if (b) { // 'a' owns the memory
      Resize(a.dim[0],a.dim[1],a.dim[2],a.mSize);
      for(int i=0;i<size;++i) v[i]=a.v[i];
#ifdef SLOW_WARNING
      warning("Avoid assigning Array3",a.Size());
#endif
    } else {
      size   = a.size;
      mSize  = a.mSize;
      dim[0] = a.dim[0];
      dim[1] = a.dim[1];
      dim[2] = a.dim[2];
      v      = a.v;
    }
    return *this;
  }

  Array3 & operator=(const T & x) {        // Fill with one element
    for(int i=0;i<size;++i) v[i]=x;
    return *this;
  }

  Array3 & operator+=(const T & x) {           // Add x to all elements
    for(int i=0;i<size;++i) v[i]+=x;
    return *this;
  }

  Array3 & operator+=(const Array3 <T> & a) {  // Add array of equal dimensions
    for(int i=0;i<dim[0];++i) 
      for(int j=0;j<dim[1];++j) 
	for(int k=0;k<dim[2];++k) 
	  (*this)(i,j,k) += a(i,j,k);
    return *this;
  }

  Array3 & operator-=(const T & x) {           // Subtract x from all elements
    for(int i=0;i<size;++i) v[i]-=x;
    return *this;
  }

  Array3 & operator-=(const Array3 <T> & a) {  // Add array of equal dimensions
    for(int i=0;i<dim[0];++i) 
      for(int j=0;j<dim[1];++j) 
	for(int k=0;k<dim[2];++k) 
	  (*this)(i,j,k) -= a(i,j,k);
    return *this;
  }

  Array3 & operator*=(const T & x) {           // Multiply all elements by x
    for(int i=0;i<size;++i) v[i]*=x;
    return *this;
  }

  Array3 & operator/=(const T & x) {           // Divide all elements by x
    for(int i=0;i<size;++i) v[i]/=x;
    return *this;
  }

  bool operator==(const Array3 <T> & a) const {      // Are Array identical?
    if (dim[0]!=a.dim[0] || dim[1]!=a.dim[1] || dim[2]!=a.dim[2])
      return false;
    for(int i=0;i<size;++i) {
      if (v[i]!=a.v[i]) return false;
    }
    return true;
  }

  bool operator!=(const Array3 <T> & a) const {      // Are Array identical?
    return !(*this==a);
  }

  T Max() const {
    Limits(0,size);
    T x=v[0];
    for(int i=1; i<size; i++) {
      if (v[i]>x) x=v[i];
    }
    return x;
  }

  T Min() const {
    Limits(0,size);
    T x=v[0];
    for(int i=1; i<size; i++) {
      if (v[i]<x) x=v[i];
    }
    return x;
  }

  T Sum() const {
    T sum=0;
    for(int i=0; i<size; i++) {
      sum += v[i];
    }
    return sum;
  }

  // Adjust logical size of array, only adjust physical size (and allocate memory) if necessary
  // some of the functions above like +=Array2 do not work any more if dim[1]*dim[2] has been modified!!!
  void Resize(const int n1, const int n2, const int n3) {
    Resize(n1,n2,n3,max(n1*n2*n3,mSize));
  }

  // Adjust physical and logical size, allocate memory if nm!=mSize
  void Resize(const int n1, const int n2, const int n3, const int nm) {
    if (!b) error ("Cannot resize fixed-size array.");
    if (nm!=mSize) {
      mSize=nm;
      if (v) delete[] v;
      v = new T[mSize];
      if (mSize>0 && v==0) error("Could not allocate memory",mSize);
    }
    size   = n1*n2*n3;
    dim[0] = n1; 
    dim[1] = n2;
    dim[2] = n3;
  }

  // CopyResize with extra memory has not been implemented
  void CopyResize(const int n1, const int n2, const int n3) {
    CopyResize(n1,n2,n3,max(n1*n2*n3,mSize));
  }

  // Resize and copy elements and reserve memory for exactly nm elements (for later CopyResize)
  void CopyResize(const int n1, const int n2, const int n3, const int nm) {
    if (!b) error ("Cannot resize fixed-size array");
    int nn= n1*n2*n3;
    if (nn>nm) error("Copy resize limit error",n1,n2,n3,nm);
    if (n3!=dim[2] || n2!=dim[1] || nm!=mSize) {
      mSize = nm;
      T* vv=new T[mSize];
      if (mSize>0 && vv==0) error("Could not allocate memory",mSize);
      for(int i1=0;i1<min(n1,dim[0]);++i1) {
	for(int i2=0;i2<min(n2,dim[1]);++i2) {
	  for(int i3=0;i3<min(n3,dim[2]);++i3) {
	    vv[(i1*n2+i2)*n3+i3] = v[(i1*dim[1]+i2)*dim[2]+i3];
	  }
	}
      }
      if (v) delete[] v;
      v      = vv;
      dim[1] = n2;
      dim[2] = n3;
    }
    size   = nn;
    dim[0] = n1;
  }

  void Free() {
    if (!b) error ("Cannot free memory of a fixed array");
    if (mSize>0 && v) {
      delete[] v;
      v = 0;
      size = dim[0] = dim[1] = dim[2] = 0;
    }
  }

  // Assume the memory has be allocated by some other routine
  // if (b1==true) take over the ownership and deallocated the memory at the end
  void SetPointer(T* vv, const int n1, const int n2, const int n3, const bool b1=true) {
    if (v && b) delete[] v;
    int nn=n1*n2*n3;
    size  = mSize = nn;
    v     = vv;
    b     = b1;
    dim[0]= n1;
    dim[1]= n2;
    dim[2]= n3;
  }

  // Return the point but do not give up the ownership
  T* GetPointer() {
    return v;
  }
  
  void Check() const {
    int nn=1;
    for(int i=0; i<Rank(); i++) {
      if (dim[i]<0) error("in Array dimension",dim[i]);
      nn *= dim[i];
    }
    if (nn!=size) error("size is wrong");
    if (mSize<size) error("mSize is wrong");
  }

  void Insert(const int d, const int ii, const int nn=1) {
    if (d==0) { 
      CopyResize(dim[0]+nn,dim[1],dim[2]); // for simplicity, copy elements twice 
      for(int i=dim[0]-nn-1; i>=ii; i--) {
	for(int j=0; j<dim[1]; j++) {
	  for(int k=0; k<dim[2]; k++) {
	    (*this).operator()(i+nn,j,k) = (*this).operator()(i,j,k);
	  }
	}
      }
    } else if (d==1) {
      CopyResize(dim[0],dim[1]+nn,dim[2]); // for simplicity, copy elements twice 
      for(int i=0; i<dim[0]; i++) {
	for(int j=dim[1]-nn-1; j>=ii; j--) {
	  for(int k=0; k<dim[2]; k++) {
	    (*this).operator()(i,j+nn,k) = (*this).operator()(i,j,k);
	  }
	}
      }
    } else {
      CopyResize(dim[0],dim[1],dim[2]+nn); // for simplicity, copy elements twice 
      for(int i=0; i<dim[0]; i++) {
	for(int j=0; j<dim[1]; j++) {
	  for(int k=dim[2]-nn-1; k>=ii; k--) {
	    (*this).operator()(i,j,k+nn) = (*this).operator()(i,j,k);
	  }
	}
      }
    }		
  }		
  
  void Insert(const T & x, const int d, const int ii, const int nn=1) {
    Insert(d,ii,nn);
    if (d==0) { 
      for(int i=ii; i<ii+nn; i++) {
	for(int j=0; j<dim[1]; j++) {
	  for(int k=0; k<dim[2]; k++) {
	    (*this).operator()(i,j,k)=x;
	  }
	}
      }
    } else if (d==1) {
      for(int i=0; i<dim[0]; i++) {
	for(int j=ii; j<ii+nn; j++) {
	  for(int k=0; k<dim[2]; k++) {
	    (*this).operator()(i,j,k)=x;
	  }
	}
      }
    } else {
      for(int i=0; i<dim[0]; i++) {
	for(int j=0; j<dim[1]; j++) {
	  for(int k=ii; k<ii+nn; k++) {
	    (*this).operator()(i,j,k)=x;
	  }
	}
      }
    }
  }
  
};

////////////////////////////////////////////////////////////////////////////////////////////////////

template <class T> Array2<T> Array1<T>::MakeArray2(const int n1) {
  const int n2=size/n1;
  return Array2<T>(n1,n2,v);
}

template <class T> const Array2<T> Array1<T>::MakeConstArray2(const int n1) {
  const int n2=size/n1;
  return Array2<T>(n1,n2,v);
}

// create from section of 1d array
template <class T> inline Array1<T>::Array1(const Array1 <T> & a1, const int n1, const int n2):
  size(n2-n1),
  mSize(size),
  inc(stInc),
  v(a1.v+n1),
  b(false) {
  BiggerLimit(n1,n2);
  Limits(n1,a1.size);
  Limits(n2,a1.size);
}

// create from 1d vector in 2d array
template <class T> inline Array1<T>::Array1(const Array2 <T> & a2, const int n1):
  size(a2.dim[1]),
  mSize(size),
  inc(stInc),
  v(a2.v+n1*a2.dim[1]),
  b(false) {
  Limits(n1,a2.dim[0]);
}

// create from section of 1d vector in 2d array
template <class T> inline Array1<T>::Array1(const Array2 <T> & a2, const int n1, const int n21, const int n22):
  size(n22-n21),
  mSize(size),
  inc(stInc),
  v(a2.v+n1*a2.dim[1]+n21),
  b(false) {
  Limits(n1,a2.dim[0]);
  Limits(n21,a2.dim[1]);
  LimitsInclusive(n22,a2.dim[1]);
  BiggerLimit(n21,n22);
}

// create from 1d vector in 3d array
template <class T> inline Array1<T>::Array1(const Array3 <T> & a3, const int n1, const int n2):
  size(a3.dim[2]),
  mSize(size),
  inc(stInc),
  v(a3.v+(n1*a3.dim[1]+n2)*a3.dim[2]),
  b(false) {
  Limits(n1,a3.dim[0]);
  Limits(n2,a3.dim[1]);
}

// create from 2d array in 3d array
template <class T> inline Array2<T>::Array2(const Array3 <T> & a3, const int n1):
  size(a3.dim[1]*a3.dim[2]),
  mSize(size),
  v(a3.v+n1*a3.dim[1]*a3.dim[2]),
  b(false) {
  dim[0] = a3.dim[1];
  dim[1] = a3.dim[2];
  Limits(n1,a3.dim[0]);
}

template <class T>
void WriteArray1(const Array1 <T> & a) {
  cout << "[ ";
  for(int i=0; i<a.Size(); i++) {
    cout << a[i] << " ";
  }
  cout << "]";
}

template <class T>
void WriteArray2(const Array2 <T> & a, const string & s="") {
  if (s.length()==0) cout << "2D array of size " << a.dim[0] << " x " << a.dim[1] << endl;
  for(int i=0; i<a.dim[0]; i++) {
    cout << s << i << "   ";
    Array1 <T> a1(a,i);
    WriteArray1(a1);
    cout << endl;
  }
}

template <class T>
void WriteArray3(const Array3 <T> & a) {
  cout << "3D array of size " << a.dim[0] << " x " << a.dim[1] << " x " << a.dim[2] << endl;
  for(int i=0; i<a.dim[0]; i++) {
    string s = IntToString(i)+"   ";
    Array2 <T> a2(a,i);
    WriteArray2(a2,s);
  }
}

#ifdef _MAIN_

template <class T>
void Equal(const Array1 <T> & a, const Array1<T> & b) {
  a.Check();
  b.Check();
  if (a!=b) error("Array1's not equal");
  if (a.mSize != b.mSize)
    error("mSize problem");
  if (a.b != b.b)
    error("b problem");
  if (a.inc != b.inc)
    error("inc problem");
}

template <class T>
void Equal(const Array2 <T> & a, const Array2<T> & b) {
  a.Check();
  b.Check();
  if (a!=b) error("Array2's not equal");
  if (a.mSize != b.mSize)
    error("mSize problem");
  if (a.b != b.b)
    error("b problem");
}

template <class T>
void Equal(const Array3 <T> & a, const Array3<T> & b) {
  a.Check();
  b.Check();
  if (a!=b) error("Array3's not equal");
  if (a.mSize != b.mSize)
    error("mSize problem");
  if (a.b != b.b)
    error("b problem");
}

int main() {
  Array1 <int> a,b;
  Equal(a,b);

  int n=10;
  int x=314;
  a.Resize(n);
  a=x;
  Array1 <int> c(n);
  c=x;
  Array1 <int> d(n,x);
  Equal(a,c);
  Equal(a,d);
  
  a += 1;
  c = x+1;
  Equal(a,c);
  d = 2*c[0];
  d -= a;
  Equal(a,d);
  
  Array1 <int> e;
  for (int i=0; i<2*n; i++) {
    e.PushBack(e.Size());
  }
  if (e.Min()!=0) error();
  if (e.Max()!=e.Size()-1) error();
  Write2(e.Sum(),(e.Size()*(e.Size()-1))/2);
  if (e.Sum()!=(e.Size()*(e.Size()-1))/2) error();

  int nn=4*n;
  Array1 <int> f(nn);
  f.Last() = 5;
  Write2(f.v,f.Last());
  f = e;
  Equal(e,f);
  //  if (f.Last() != 5) error();
  f.CopyResize(e.Size());
  Equal(e,f);
  f.CopyResize(e.Size(),nn);
  if (f!=e) error();
  f.CopyResize(nn,nn);
  int *v = f.GetPointer();
  f.CopyResize(e.Size());
  if (e!=f) error();
  f.CopyResize(nn,nn);
  if (v != f.GetPointer()) error();
  f.PushBack(4);
  if (v == f.GetPointer()) error();
  v = f.GetPointer();
  f.PushBack(4);
  f.PushBack(4);
  f.PushBack(4);
  if (v != f.GetPointer()) error();
  f.CopyResize(e.Size());
  if (e!=f) error();
  
  cout << "Non-fixed Array1 tests passed" << endl;

  Array1 <int> g(e.Size(),e.GetPointer());
  Array1 <int> g2(g);
  if (e[0]!=g2[0]) error();
  if (e.GetPointer()!=g2.GetPointer()) error();
  if (e.Sum()!=g2.Sum()) error();
  if (e!=g2) error();
  cout << "Fixed Array1 tests passed" << endl;

  // Test making Array2 by sub-dividing Array1
  Array2 <int> e2 = e.MakeArray2(2);
  Array2 <int> e2c(e2);
  Array2 <int> e2cc;
  e2cc = e2c;
  if (e2cc.b) error();
  if (e.Last()!=e2cc(e2cc.dim[0]-1,e2cc.dim[1]-1)) error();
  
  // Make a new Array1 containing a section of the old one
  Array1 <int> ef(e,3,5);
  if (ef.Size()!=2) error();
  if (ef[0]!=e[3]) error();
  if (ef.b) error();

  ef[1]=100;
  if (ef[1]!=e[4]) error();
  
  cout << "--- all Array1 tests completed ---" << endl << endl;

  int n1=6;
  int n2=13;
  Array2 <int> h1(n1,n2,2);
  Array2 <int> h2(h1);
  for(int i=0; i<n1; i++) {
    for(int j=0; j<n2; j++) {
      if (h2(i,j)!=2) error();
      h2(i,j) = i*h2.dim[1]+j;
      //      Write3(i,j,h2(i,j));
      if (h1(i,j)!=2) error();
    }
  }

  Array2 <int> h3;
  h3 = h2;
  Equal(h2,h3);
  h3.CopyResize(h2.dim[0]*2,h2.dim[1]*3);
  h3.CopyResize(h2.dim[0],h2.dim[1]);
  if (h2!=h3) error();
  if (h2.Sum()!=h3.Sum()) error();
  if (h2.Min()!=h3.Min()) error();
  if (h2.Min()!=0) error();
  Write(h2.Max());
  if (h2.Max()!=h2.size-1) error();
  if (h2.Max()!=h3.Max()) error();

  Array1 <int> hh(h3,3);
  if (hh.Last()!=h3(3,h3.dim[1]-1)) error();

  Array1 <int> hhh(h3,3,3,5);
  if (hhh.Last()!=h3(3,4)) error();

  if (h2[5][4] != h2(5,4)) error();
  Array1 <int> q(h2[5]);
  q[4]=-21;
  (h2[5])[4] = -21;
  h2(5,4) = -21;
  if (h2[5][4] != h2(5,4)) error();
  if (h2[5][4] == h2(4,5)) error();

  h2[4][5] = -15;
  if (h2[4][5] != -15 || h2[4][5] !=  h2(4,5) ) error("[][]");

  cout << "--- all Array2 tests completed ---" << endl << endl;

  int n3=21;
  Array3 <int> m0(n1,n2,n3);
  m0.Check();
  m0 = 3;
  Array3 <int> m(n1,n2,n3,3);
  Equal(m0,m);
  Array3 <int> m2(m);
  Array3 <int> m3;
  m3 = m;
  for(int i=0; i<n1; i++) {
    for(int j=0; j<n2; j++) {
      for(int k=0; k<n3; k++) {
      if (m3(i,j,k)!=3) error();
      m3(i,j,k) = (i*m3.dim[0]+j)*m3.dim[1]+k;
      //      Write3(i,j,h2(i,j));
      if (m(i,j,k)!=3) error();
      }
    }
  }
  if (m3[5][2][7] != m3(5,2,7))
    error();
  m3[5][2][7] = m3(5,2,7) = -11;
  if (m3[5][2][7] != m3(5,2,7))
    error();
  
  cout << "--- all Array3 tests completed ---" << endl << endl;

  int tn=10;
  Array1 <double> t1(tn,0.0);
  //  Array1 <double> t1(tn,2*tn,0.0);
  Write(t1.mSize);
  for(int i=0; i<tn; i++) 
    t1[i] = i;
  WriteArray1(t1); cout << endl;
  t1.Insert(tn/2,3);
  WriteArray1(t1); cout << endl;
  t1[tn/2]   = -10.0;
  t1[tn/2+1] = -10.0;
  t1[tn/2+2] = -10.0;
  t1.Insert(-100.0,tn/2,2);
  t1.Insert(+1000.0,t1.Size(),2);
  WriteArray1(t1); cout << endl;
  if (abs(t1.Sum() - tn*(tn-1)/2+230.0-2000.0) > 1e-12) error("t1",t1.Sum());
  cout << endl;
  cout << endl;

  int tn2=4;
  Array2 <double> t2(tn2,tn2+2);
  for(int i=0; i<t2.dim[0]; i++) 
    for(int j=0; j<t2.dim[1]; j++) 
      t2(i,j) = i*10+j;
  double s=t2.Sum();
  WriteArray2(t2);  cout << endl;

  t2.Insert(50.0,0,0);
  WriteArray2(t2);  cout << endl;

  t2.Insert(-10.0,1,tn2/2,2);
  WriteArray2(t2);  cout << endl;

  if (abs(t2.Sum() - (s+(tn2+2)*50.0-(tn2+1)*2*10.0)) > 1e-12) error("t2");
  cout << endl;

  int tn3=3;
  Array3 <double> t3(tn3,tn3,tn3-1);
  for(int i=0; i<t3.dim[0]; i++) 
    for(int j=0; j<t3.dim[1]; j++) 
      for(int k=0; k<t3.dim[2]; k++) 
	t3(i,j,k) = i*100+j*10+k;
  double s3=t3.Sum();
  WriteArray3(t3);  cout << endl;

  t3.Insert(-50.0,0,1,2);
  WriteArray3(t3);  cout << endl;
  
  t3.Insert(-10.0,1,0);
  WriteArray3(t3);  cout << endl;

  t3.Insert(-1.0,2,t3.dim[2],2);
  WriteArray3(t3);  cout << endl;

  cout << "Program finished properly." << endl;
}
#endif

#endif // _ARRAY_
