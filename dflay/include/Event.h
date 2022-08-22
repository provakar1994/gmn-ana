#ifndef EVENT_H
#define EVENT_H

// generic event class for ROOT files

#include <cstdlib>
#include <iostream>
#include <string>

namespace util_df { 
   template <typename T>
      class Event {
	 private:
	    int fID;                          // event index 
	    std::vector<std::string> fName;   // names of all variables associated with event
	    std::vector<T> fValue;            // values of all variables 

	 public:
	    Event();
	    ~Event();

	    void Print();
	    void ClearData();
	    void SetID(int v) { fID = v; }
	    void PushBackData(const char *tag,T v);
	    void SetVariableNames(std::vector<std::string> var);

	    int SetData_byIndex(int i,T v);
	    int SetData_byIndex(int i,const char *tag,T v);

	    T GetData_byIndex(int i)          const { return fValue[i]; }
	    T GetData_byName(const char *tag) const;
      };

} // ::util 

#endif
