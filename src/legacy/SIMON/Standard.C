////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Standard definitions and functions that are used almost everywhere     //
//                                                                        //
// Burkhard Militzer                                    Urbana 4-9-99     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "Standard.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

// ofstream COUT;

// Does this need "ss << ends" ???
string IntToString(const int i) {
  ostringstream ss;
  ss << i;
  return ss.str();
}

string DoubleToString(const double d) {
  ostringstream ss;
  ss << d;
  return ss.str();
}

int StringToInt(const string & s) {
  const string intChars(" +-0123456789");
  if (int j=s.find_first_not_of(intChars) != string::npos)
    error("Parser: Not an integer ",s,j+1);
 
  istringstream ss(s);
  int value;
  ss >> value;
  return value;
}

double StringToDouble(const string & s) {
  const string doubleChars(" +-0123456789.edDE");
  if (int j=s.find_first_not_of(doubleChars) != string::npos)
    error("Parser: Not a double ",s,j+1);

  istringstream ss(s);
  double value;
  ss >> value;
  return value;
}

string UpperCase(const string & s) {
  string sl;
  for(string::const_iterator p=s.begin();p!=s.end();++p) {
    sl += toupper(*p);
  }
  return sl;
}
 
string LowerCase(const string & s) {
  string sl;
  for(string::const_iterator p=s.begin();p!=s.end();++p) {
    sl += tolower(*p);
  }
  return sl;
}

void Terminate() {
  
#ifdef USE_MPI
  int errorcode=1;
  MPI_Abort(MPI_COMM_WORLD, errorcode);
#endif

   cerr << "Exiting now\n";
   exit(0);
}

