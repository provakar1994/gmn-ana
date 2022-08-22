#ifndef CSV_MANAGER_H
#define CSV_MANAGER_H

// a generic CSV file reader/writer 

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

class CSVManager {

   private:
      int fVerbosity;
      int fNumRow,fNumCol;
      bool fHeaderExists;
      std::string fDelimiter;
      std::vector<std::string> fHeader;
      std::vector< std::vector<std::string> > fData;

      int SplitString(const char delim,const std::string inStr,std::vector<std::string> &out);
      int SplitString_whiteSpace(const std::string myStr,std::vector<std::string> &out);

      template <typename T>
         int CheckType(T data){
            int rc=0;
            if( std::is_arithmetic<T>::value ){
               if( std::is_integral<T>::value ){
     		  // it's an integer or boolean 
     		  rc = 0;
     	  }else{
     		  // it's a double or float 
     		  rc = 1;
     	  }
            }else{
     	       rc = 2; // it's a string
            }
            return rc;
         }

   public:
      CSVManager(const char *delim="csv",int verbosity=0);
      ~CSVManager();

      int Print();
      int PrintHeader();
      int PrintColumns(std::string cols);
      int PrintMetaData(); 

      int ClearData();
      int ReadFile(const char *inpath,bool header=false,int lineSkip=0);
      int WriteFile(const char *outpath);
      int InitTable(int NROW,int NCOL);

      // setter methods
      void SetVerbosity(int v) {fVerbosity = v;}
      int SetHeader(std::string fullHeader);
      int SetHeader(std::vector<std::string> header);
      int SetElement_str(int row,int col,std::string x);
      int SetColumn_str(int col,std::vector<std::string> x){
	 const int N = x.size();
	 for(int i=0;i<N;i++) SetElement_str(i,col,x[i]);
	 return 0;
      }

      // templated setter methods 
      template <typename T>
         int SetElement(int row,int col,T x){
            char data[200];
            int type = CheckType<T>(x);
            if(type==0) sprintf(data,"%d"  ,(int)x   );
            if(type==1) sprintf(data,"%.15E",(double)x);
            std::string DATA = data;
            if(row<fNumRow && col<fNumCol){
                 fData[row][col] = DATA;
            }else{
                 std::cout << "[CSVManager::SetElement]: ERROR!  Invalid indices! row = "
                         << row << ", col = " << col << std::endl;
                 std::cout << "[CSVManager::SetElement]: Max dimensions are NRow = "
                         << fNumRow << ", NCol = " << fNumCol << std::endl;
                 return 1;
            }
            return 0;
         }

       template <typename T>
          int SetColumn(int col,std::vector<T> x){
             const int N = x.size();
             for(int i=0;i<N;i++) SetElement<T>(i,col,x[i]);
             return 0;
          }

       // getter methods
       int GetNumRows()    const { return fNumRow;    }
       int GetNumColumns() const { return fNumCol;    }
       int GetVerbosity()  const { return fVerbosity; }
       int GetColumnIndex_byName(std::string colName);

       int GetHeader(std::vector<std::string> &header);
       int GetColumn_byIndex_str(int colIndex,std::vector<std::string> &data);
       int GetColumn_byName_str(std::string colName,std::vector<std::string> &data);

       std::string GetElement_str(int rowIndex,int colIndex);

       // templated getter methods (to obtain arithmetic types)
       template <typename T>
          int GetColumn_byIndex(int colIndex,std::vector<T> &data){
             // find the data by col index 
             int NROW = fData.size();
             T elem;
             for(int i=0;i<NROW;i++){
                  elem = GetElement<T>(i,colIndex);
                  data.push_back(elem);
             }
             return 0;
          }

         template <typename T>
            int GetColumn_byName(std::string colName,std::vector<T> &data){
               // find the column index by name 
               int NCOL=0,NROW=0,k=-1;
               if(fHeaderExists){
                  NCOL = fHeader.size();
                  for(int i=0;i<NCOL;i++) if(fHeader[i].compare(colName)==0) k = i;
                  if(k>=0){
                       GetColumn_byIndex<T>(k,data);
                  }else{
                       std::cout << "[CSVManager::GetColumn_byName]: Cannot find the key '"
                   	    << colName << "' in header!" << std::endl;
                       return 1;
                  }
               }else{
                  std::cout << "[CSVManager::GetColumn_byName]: No header to search!";
                  std::cout << "  Try CSVManager::GetColumn_byIndex" << std::endl;
                  return 1;
               }
               return 0;
	    }

         template <typename T>
            T GetElement(int rowIndex,int colIndex){
               int NROW = fData.size();
               int NCOL = fData[0].size();
               std::string data="NONE";
               if(NROW>0 && NCOL>0){
                  data = fData[rowIndex][colIndex];
               }else{
                  std::cout << "[CSVManager::GetElement]: NO data for row "
                     << rowIndex << ", col " << colIndex << std::endl;
               }
	       // now parse the string 
	       T x;
	       if( std::is_arithmetic<T>::value ){
                  if( std::is_integral<T>::value ){
                          // it's an integer or boolean 
                          x = std::atoi( data.c_str() );
                  }else{
                          // it's a double or float 
                          x = std::atof( data.c_str() );
                  }
	       }
	       return x;
	    }

}; // ::CSVManager

#endif 
