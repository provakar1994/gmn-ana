#ifndef UTIL_DF_JSON_MANAGER_HH
#define UTIL_DF_JSON_MANAGER_HH

// A wrapper type class to streamline handling of JSON file data 
// based on the include file: 
// /u/site/12gev_phys/2.5/Linux_CentOS7.7.1908-gcc9.2.0/root/6.24.06/include/nlohmann/json.hpp
// since this is in ROOT already, don't need to explicitly include it here

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "/u/site/12gev_phys/2.5/Linux_CentOS7.7.1908-gcc9.2.0/root/6.24.06/include/nlohmann/json.hpp"

enum dataType{ 
   kINT,
   kDBL,
   kSTR
};

using json = nlohmann::json;

class JSONManager {

   private:
      json fObject;

         template <typename T>
            int CheckType(T data) const{
               int rc=0;
               if( std::is_arithmetic<T>::value ){
                  if( std::is_integral<T>::value ){
		     // it's an integer or boolean 
		     rc = kINT;
		  }else{
		     // it's a double or float 
		     rc = kDBL;
		  }
	       }else{
		  rc = kSTR; // it's a string
	       }
	       return rc;
	    }


   public:
      JSONManager(const char *filepath="NONE");
      ~JSONManager();

      int Print();
      int ReadFile(const char *filepath);

      bool DoesKeyExist(std::string keyName) const;

      int GetVectorFromKey_str(std::string key,std::vector<std::string> &data);

      std::string GetValueFromKey_str(std::string key) const;
      std::string GetValueFromSubKey_str(std::string key,std::string subKey) const;

      // templated methods 
      template <typename T> 
	 T GetValueFromKey(std::string key) const {
	    // this function is to retrieve a number, but could be stored as a string in the JSON object 
	    T val;
	    int outDataType = CheckType<T>(val);
	    bool exist = DoesKeyExist(key);
	    std::string DATA="";
	    if(exist){
	       // key exists 
	       if( fObject[key].is_string() ){
		  // data is a string type, convert to int or double 
		  DATA = fObject[key].get<std::string>();
		  if(outDataType==kINT) val = std::atoi( DATA.c_str() );
		  if(outDataType==kDBL) val = std::atof( DATA.c_str() );
	       }else{
		  // data isn't a string, just typecast and store the value  
		  val = (T)(fObject[key]);
	       }
	    }else{
	       // key doesn't exist, return 0
	       val = (T)(0);
	    }
	    return val;
	 }

      template <typename T> 
	 T GetValueFromSubKey(std::string key,std::string subKey) const {
	    // this function is to retrieve a number, but could be stored as a string in the JSON object 
	    T val;
	    int outDataType = CheckType<T>(val);
	    bool exist = DoesKeyExist(key);
	    std::string DATA="";
	    if(exist){
	       // key exists 
	       if( fObject[key][subKey].is_string() ){
		  // data is a string type, convert to int or double 
		  DATA = fObject[key][subKey].get<std::string>();
		  if(outDataType==kINT) val = std::atoi( DATA.c_str() );
		  if(outDataType==kDBL) val = std::atof( DATA.c_str() );
	       }else{
		  // data isn't a string, just typecast and store the value  
		  val = (T)(fObject[key][subKey]);
	       }
	    }else{
	       // key doesn't exist, return 0
	       val = (T)(0);
	    }
	    return val;
	 }

      template <typename T>
	 int GetVectorFromKey(std::string key,std::vector<T> &data){
	    T arg;
            int N=0;
	    int outDataType = CheckType<T>(arg);
	    std::string DATA="";
	    // first check if the key exists
	    bool exist = DoesKeyExist(key);
	    if(exist){
	       // found the key, fill the vector
	       N= fObject[key].size();
	       for(int i=0;i<N;i++){
		  if( fObject[key][i].is_string() ){
		     // data is a string type, convert to int or double 
		     DATA = fObject[key][i].get<std::string>();
		     if(outDataType==kINT) arg = std::atoi( DATA.c_str() );
		     if(outDataType==kDBL) arg = std::atof( DATA.c_str() );
		  }else{
		     // data isn't a string, just typecast and store the value  
		     arg = (T)(fObject[key][i]);
		  }
		  data.push_back(arg);
	       }
	    }else{
	       return 1;
	    }
	    return 0;
	 }

};

#endif
