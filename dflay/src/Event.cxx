#include "../include/Event.h"
//______________________________________________________________________________
namespace util_df { 
   //______________________________________________________________________________
   template <typename T>
      Event<T>::Event(){
	 fID = 0;
      }
   //______________________________________________________________________________
   template <typename T>
      Event<T>::~Event(){
	 ClearData();
      }
   //______________________________________________________________________________
   template <typename T>
      void Event<T>::ClearData(){
	 fName.clear();
	 fValue.clear();
      }
   //______________________________________________________________________________
   template <typename T>
      void Event<T>::Print(){
	 const int N = fName.size();
	 std::cout << "event " << fID << ": " << std::endl;
	 for(int i=0;i<N;i++){
	    std::cout << "   " << fName[i] << " = " << fValue[i] << std::endl;;
	 }
      }
   //______________________________________________________________________________
   template <typename T>
      void Event<T>::SetVariableNames(std::vector<std::string> varName){
	 // set all variable names and initialize values to zero
	 const int N = varName.size();
	 fName.resize(N);
	 fValue.resize(N);
	 for(int i=0;i<N;i++){
	    fName[i]  = varName[i];
	    fValue[i] = 0;
	 }
      }
   //______________________________________________________________________________
   template <typename T>
      int Event<T>::SetData_byIndex(int i,T v){
	 // set only the value of the data for index i
	 int NV = fName.size();
	 if(i<NV){
	    fValue[i] = v;
	 }else{
	    std::cout << "[Event::SetData]: Variable index out of range! Data NOT set" << std::endl;
	    return 1;
	 }
	 return 0;
      }
   //______________________________________________________________________________
   template <typename T>
      int Event<T>::SetData_byIndex(int i,const char *tag,T v){
	 // set name and value of the data for index i
	 int NV = fName.size();
	 std::string name = tag;
	 if(i<NV){
	    fName[i]  = name;
	    fValue[i] = v;
	 }else{
	    std::cout << "[Event::SetData]: Variable index out of range! Data NOT set" << std::endl;
	    return 1;
	 }
	 return 0;
      }
   //______________________________________________________________________________
   template <typename T>
      void Event<T>::PushBackData(const char *tag,T v){
	 std::string name = tag;
	 fName.push_back(name);
	 fValue.push_back(v);
      }
   //______________________________________________________________________________
   template <typename T>
      T Event<T>::GetData_byName(const char *tag) const{
	 int j=-1;
	 T val=0;
	 std::string theName = tag;
	 const int N = fName.size();
	 for(int i=0;i<N;i++){
	    if(fName[i].compare(theName)==0) j = i;
	 }
	 if(j<0){
	    std::cout << "[Event::GetData_byName]: Invalid name = " << theName << "! Returning 0." << std::endl;
	 }else{
	    val =fValue[j];
	 }
	 return val;
      }
} // ::util 
