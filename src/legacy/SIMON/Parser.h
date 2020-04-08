/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Parser class                                                            //
//                                                                         //
// Burkhard Militzer                                    Urbana 4-1-99      //
// based on earlier version by John Shumway                                //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _PARSER_
#define _PARSER_

#include <string>
using namespace std;
#include "Vector.h"

class Parser {
public:
  Parser(istream& is_);
  Parser(const string & fn);
  
  //Read the next non-empty line and store in "words" (strips out comments)
  //Returns istream status (true=OK)
  bool ReadLine();
  void ReadLineSavely() {
    if (!ReadLine())
      error("Parser encounted unexpected end of file");
  }

  //Return "words" as a string
  const string& GetLineString() const;

  //Return the number of words in the current line ("words").
  int GetNWords() const;
  int GetLineNumber() const {
    return nLinesRead;
  }
  int GetNumberOfNonEmptyLines() const {
    return nNonEmptyLinesRead;
  }

  //Return a word as a string, int, or double. (no some error checking done)
  const string& GetString(int) const;
  const string& GetNameString(int) const;
  const string GetStringLowerCase(int) const;
  const string GetNameStringLowerCase(int) const;
  const string GetStringUpperCase(int) const;
  const string GetNameStringUpperCase(int) const;
  int GetInt(int) const;
  double GetDouble(int) const;

  void Get(int & i, const int n) const {
    i=GetInt(n);
  }
  void Get(double & d, const int n) const {
    d=GetDouble(n);
  }
  void Get(string & s, const int n) const {
    s=GetString(n);
  }

  void SetIgnoreComments() {
    ignoreComments=true;
  }
  void UnSetIgnoreComments() {
    ignoreComments=false;
  }
  void SetIgnoreEmptyLines() {
    ignoreEmptyLines=true;
  }
  void UnSetIgnoreEmptyLines() {
    ignoreEmptyLines=false;
  }

  void ProcessKeyword(const string & keyword, const int n);
  void ProcessKeyword(const string & keyword);
  void ProcessLine(const int n);

  void PutLineBack();
  bool Status() {
    return is;
  }

private:
  void CheckNotNumber(const string & s) const;

  ifstream ifs;
  istream& is;
  Vector<string> words;
  string lineString;
  bool ignoreComments;
  bool ignoreEmptyLines;
  int nLinesRead;
  int nNonEmptyLinesRead;
};

#endif // _PARSER_
