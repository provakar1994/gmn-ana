#include "../include/JSONManager.h"
//______________________________________________________________________________
JSONManager::JSONManager(const char *filepath){
   std::string inpath = filepath;
   if(inpath.compare("NONE")!=0) ReadFile(filepath);
}
//______________________________________________________________________________
JSONManager::~JSONManager(){

}
//______________________________________________________________________________
bool JSONManager::DoesKeyExist(std::string keyName) const{
   auto it_key = fObject.find(keyName);  // this is an iterator 
   if (it_key!=fObject.end() ){
      return true;  // not at the end of fObject -- found the key 
   }else{
      std::cout << "[JSONManager::DoesKeyExist]: Key '" << keyName << "' does not exist!" << std::endl;
      return false;  // at the end of fObject -- didn't find the key 
   }
}
//______________________________________________________________________________
std::string JSONManager::GetValueFromKey_str(std::string key) const{
   bool exist = DoesKeyExist(key);
   if(exist){
      return fObject[key];
   }else{
      return "DOES_NOT_EXIST";
   }
}
//______________________________________________________________________________
std::string JSONManager::GetValueFromSubKey_str(std::string key,std::string subKey) const{
   bool exist = DoesKeyExist(key);
   if(exist){
      return fObject[key][subKey];
   }else{
      return "DOES_NOT_EXIST";
   }
}
//______________________________________________________________________________
int JSONManager::GetVectorFromKey_str(std::string key,std::vector<std::string> &data){
   int N=0;
   bool exist = DoesKeyExist(key);
   if(exist){
      // found the key, fill the vector
      N = fObject[key].size();
      for(int i=0;i<N;i++){
	 data.push_back(fObject[key][i]);
      }
   }else{
      return 1;
   }
   return 0;
}
//______________________________________________________________________________
int JSONManager::ReadFile(const char *filepath){
   std::ifstream infile;
   infile.open(filepath,std::ifstream::binary);
   if( infile.fail() ){
      std::cout << "[JSONManager::ReadFile]: Cannot open the file: " << filepath << std::endl;
      return 1;
   }else{
      infile >> fObject;
   }
   return 0;
}
//______________________________________________________________________________
int JSONManager::Print(){
   std::cout << fObject << std::endl;
   return 0;
}
