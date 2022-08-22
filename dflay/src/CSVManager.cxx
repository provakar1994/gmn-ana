#include "../include/CSVManager.h"
//______________________________________________________________________________
CSVManager::CSVManager(const char *delim,int verbosity){
   fHeaderExists = false;
   fNumCol       = 0;
   fNumRow       = 0;
   fVerbosity    = verbosity;
   fDelimiter    = delim;
}
//______________________________________________________________________________
CSVManager::~CSVManager(){
   ClearData();
   fNumCol = 0;
   fNumRow = 0;
}
//______________________________________________________________________________
int CSVManager::ClearData(){
   // delete all data 
   fHeader.clear();
   int NROW = fData.size();
   for(int i=0;i<NROW;i++){
   	fData[i].clear();
   }
   fData.clear();
   return 0;
}
//______________________________________________________________________________
int CSVManager::InitTable(int NROW,int NCOL){
   fData.resize(NROW);
   for(int i=0;i<NROW;i++) fData[i].resize(NCOL);
   for(int i=0;i<NROW;i++){
      for(int j=0;j<NCOL;j++) fData[i][j] = "EMPTY"; 
   }
   fNumRow = NROW;
   fNumCol = NCOL;
   return 0;
}
//______________________________________________________________________________
int CSVManager::SetElement_str(int row,int col,std::string x){
   if(row<fNumRow && col<fNumCol){
   	fData[row][col] = x;
   }else{
   	std::cout << "[CSVManager::SetElement]: ERROR!  Invalid indices! row = "
   		<< row << ", col = " << col << std::endl;
   	std::cout << "[CSVManager::SetElement]: Max dimensions are NRow = "
   		<< fNumRow << ", NCol = " << fNumCol << std::endl;
   	return 1;
   }
   return 0;
}
//______________________________________________________________________________
int CSVManager::SplitString(const char delim,const std::string inStr,std::vector<std::string> &out){
   // split a string by a delimiter
   std::stringstream ss(inStr);
   while( ss.good() ){
      std::string substr;
      std::getline(ss,substr,delim);
      out.push_back(substr);
   }
   return 0;
}
//______________________________________________________________________________
int CSVManager::SplitString_whiteSpace(const std::string myStr,std::vector<std::string> &out){
   // split string on white space 
   std::istringstream buffer(myStr);
   std::copy(std::istream_iterator<std::string>(buffer),
   		std::istream_iterator<std::string>(),
   		std::back_inserter(out));
   return 0;
}
//______________________________________________________________________________
std::string CSVManager::GetElement_str(int rowIndex,int colIndex){
   int NROW = fData.size();
   int NCOL = fData[0].size();
   char msg[200];
   std::string data="NONE";
   if(NROW>0 && NCOL>0){
      data = fData[rowIndex][colIndex];
   }else{
      sprintf(msg,"[CSVManager::GetElement]: NO data for row %d, col %d",rowIndex,colIndex);
      std::cout << msg << std::endl;
   }
   return data;
}
//______________________________________________________________________________
int CSVManager::GetColumn_byIndex_str(int colIndex,std::vector<std::string> &data){
   // find the data by col index 
   int NROW = fData.size();
   std::string elem;
   for(int i=0;i<NROW;i++){
      elem = GetElement_str(i,colIndex);
      data.push_back(elem);
   }
   return 0;
}
//______________________________________________________________________________
int CSVManager::GetColumn_byName_str(std::string colName,std::vector<std::string> &data){
   // find the column index by name 
   int NCOL=0,NROW=0,k=-1;
   if(fHeaderExists){
      NCOL = fHeader.size();
      for(int i=0;i<NCOL;i++) if(fHeader[i].compare(colName)==0) k = i;
      if(k>=0){
         GetColumn_byIndex_str(k,data);
      }else{
         std::cout << "[CSVManager::GetColumn_byName]: Cannot find the key '"
         		<< colName << "' in header!" << std::endl;
         return 1;
      }
   }else{
      std::cout << "[CSVManager::GetColumn_byName]: No header to search!";
      std::cout << "  Try CSVManager::GetColumn_byName" << std::endl;
      return 1;
   }
   return 0;
}
//______________________________________________________________________________
int CSVManager::SetHeader(std::string fullHeader){
   // set a new header line
   // the input is a single string with commas separating entries
   std::vector<std::string> header;
   int rc = SplitString(',',fullHeader,header);
   // clear existing header
   // don't check for existing column sizes in case the user sets this first  
   fHeader.clear();
   int N = header.size();
   for(int i=0;i<N;i++) fHeader.push_back(header[i]);
   fNumCol = fHeader.size();
   fHeaderExists = true;
   return 0;
}
//______________________________________________________________________________
int CSVManager::SetHeader(std::vector<std::string> header){
   // set a new header line
   // don't check for existing column sizes in case the user sets this first  
   fHeader.clear();
   int N = header.size();
   for(int i=0;i<N;i++) fHeader.push_back(header[i]);
   fNumCol = fHeader.size();
   fHeaderExists = true;
   return 0;
}
//______________________________________________________________________________
int CSVManager::GetHeader(std::vector<std::string> &header){
   int N=0;
   if(fHeaderExists){
      N = fHeader.size();
      for(int i=0;i<N;i++) header.push_back(fHeader[i]);
   }else{
      std::cout << "[CSVManager::GetHeader]: No header!" <<  std::endl;
      return 1;
   }
   return 0;
}
//______________________________________________________________________________
int CSVManager::GetColumnIndex_byName(std::string colName){
   // find the inded of a given column
   int k=-1;
   int NH = fHeader.size();
   for(int i=0;i<NH;i++) if( fHeader[i].compare(colName)==0 ) k = i;
   if(k<0) std::cout << "[CSVManager::GetColumnIndex_byName]: No column named '" << colName << "'" << std::endl;
   return k;
}
//______________________________________________________________________________
int CSVManager::PrintMetaData(){
   std::cout << "[CSVManager::PrintMetaData]: CSV object data" << std::endl;
   std::cout << "Number of rows:    " << fNumRow << std::endl; 
   std::cout << "Number of columns: " << fNumCol << std::endl;
   if(fHeaderExists) PrintHeader();
   return 0;
}
//______________________________________________________________________________
int CSVManager::PrintHeader(){
   std::cout << "[CSVManager::PrintHeader]: Header data: " << std::endl;
   const int N = fHeader.size();
   for(int i=0;i<N;i++) std::cout << fHeader[i] << std::endl;
   return 0;
}
//______________________________________________________________________________
int CSVManager::Print(){
   // print to screen
   int NROW = fData.size();
   int NCOL = fHeader.size();
 
   char myStr[200];
   if(fHeaderExists){
      sprintf(myStr,"%s",fHeader[0].c_str() );
      for(int i=1;i<NCOL;i++) sprintf(myStr,"%s,%s",myStr,fHeader[i].c_str() );
      std::cout << myStr << std::endl;
   }
   
   for(int i=0;i<NROW;i++){
      sprintf(myStr,"%s",fData[i][0].c_str());
      NCOL = fData[i].size();
      for(int j=1;j<NCOL;j++) sprintf(myStr,"%s,%s",myStr,fData[i][j].c_str() );
      std::cout << myStr << std::endl;
   }
   return 0;
}
//______________________________________________________________________________
int CSVManager::PrintColumns(std::string cols){

   std::vector<std::string> colName;
   int rc  = SplitString(',',cols,colName);
   int NCN = colName.size();
   
   std::vector<int> colIndex;
   int NH = fHeader.size();
   for(int i=0;i<NCN;i++){
      for(int j=0;j<NH;j++){
         if( colName[i].compare(fHeader[j])==0 ){
            // std::cout << Form("[CSVManager::PrintColumns]: Found column '%s', index = %d",colName[i].c_str(),j) << std::endl;
            colIndex.push_back(j);
         }
      }
   }
   
   int NC = colIndex.size();
   if(NCN!=NC){
      std::cout << "[CSVManager::PrintColumns]: Did not find all columns!  Printing only those found" << std::endl;
   }
   
   std::cout << "[CSVManager::PrintColumns]: Column data: " << std::endl;
   for(int i=0;i<NC-1;i++) std::cout << colName[i] << ",";
   std::cout << colName[NC-1] << std::endl;
   
   int k=0;
   std::string sx;
   char msg[200];
   sprintf(msg,"");
   
   int NROW = fData.size();
   for(int i=0;i<NROW;i++){
      // build the row
      for(int j=0;j<NC-1;j++){
         k  = colIndex[j];
         // std::cout << "NOW GETTING COLUMN " << colName[j] << ", index = " << k << std::endl;
         sx = GetElement_str(i,k);
         if(j==0){
         	sprintf(msg,"%s",sx.c_str());
         }else{
         	sprintf(msg,"%s     %s",msg,sx.c_str());
         }
      }
      k  = colIndex[NC-1];
      sx = GetElement_str(i,k);
      sprintf(msg,"%s     %s",msg,sx.c_str());
      // print to screen 
      std::cout << msg << std::endl;
      // clear
      sprintf(msg,"");
   }
   
   return 0;
}
//______________________________________________________________________________
int CSVManager::WriteFile(const char *outpath){
   // write the data to a file
   int NROW = fData.size();
   int NCOL = fData[0].size();
   char line[2048];
   char delim = ',';
   if(fDelimiter.compare("tsv")==0) delim = ' ';
   std::ofstream outfile;
   outfile.open(outpath);
   if( outfile.fail() ){
      std::cout << "[CSVManager::WriteFile]: Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      if(fHeaderExists){
	 // print header if necessary  
	 sprintf(line,"%s",fHeader[0].c_str() );
	 for(int i=1;i<NCOL;i++) sprintf(line,"%s%c%s",line,delim,fHeader[i].c_str() );
	 outfile << line << std::endl;
      }
      // print data
      for(int i=0;i<NROW;i++){
	 sprintf(line,"%s",fData[i][0].c_str());
	 NCOL = fData[i].size();
	 for(int j=1;j<NCOL;j++) sprintf(line,"%s%c%s",line,delim,fData[i][j].c_str() );
	 outfile << line << std::endl;
      }
      std::cout << "[CSVManager::WriteFile]: Data written to the file: " << outpath << std::endl;
   }
   return 0;
}
//______________________________________________________________________________
int CSVManager::ReadFile(const char *inpath,bool headerExists,int lineSkip){
   // update variable 
   fHeaderExists = headerExists;
   
   std::string aLine;
   std::vector<std::string> line;
   
   std::ifstream infile;
   infile.open(inpath);
   
   int l=0;

   if( infile.fail() ){
      std::cout << "[CSVManager::ReadFile]: Cannot open the file: " << inpath << std::endl;
      return 1;
   }else{
      if(fVerbosity>0) std::cout << "[CSVManager::ReadFile]: Opened the file: " << inpath << std::endl;
      while( !infile.eof() ){
      	std::getline(infile,aLine);
      	if(fVerbosity>1) std::cout << aLine << std::endl;
      	if(lineSkip==0){
      		line.push_back(aLine);
      	}else if(lineSkip>0){
      		if(l>lineSkip) line.push_back(aLine);
      	}
      	l++;
      }
      line.pop_back();
      infile.close();
   }

   int NROW = line.size();
   std::vector<std::string> col;
   
   // now parse the data
   int rc=0,k=0,NCOL=0;
   for(int i=0;i<NROW;i++){
      // split the line into a vector.  This is a single row 
      if(fDelimiter.compare("csv")==0){
      	rc = SplitString(',',line[i],col);
      }else if(fDelimiter.compare("tsv")==0){
      	rc = SplitString_whiteSpace(line[i],col);
      }
      // get number of columns 
      NCOL = col.size();
      if(i==0 && fHeaderExists){
      	// this is the header
      	for(int j=0;j<NCOL;j++) fHeader.push_back(col[j]);
      }else{
      	// not the header, fill the data 
      	fData.push_back(col);
      }
      // clean up
      col.clear();
   }

   // remove hashtag if necessary
   int NVS=0;
   std::string firstEntry="";
   std::vector<std::string> vStr;
   if(fDelimiter.compare("tsv")==0 && fHeaderExists){
      // first entry is the #, so just delete it 
      if(fHeader[0].compare("#")==0) fHeader.erase( fHeader.begin() );
   }else if(fDelimiter.compare("csv")==0 && fHeaderExists){
      firstEntry = fHeader[0];
      SplitString('#',firstEntry,vStr);
      // replace first entry with characters aside from leading # (if found) 
      NVS = vStr.size();
      if(NVS>1) fHeader[0] = vStr[1];
   }

   NROW    = fData.size();
   fNumCol = fHeader.size();
   fNumRow = NROW;
   
   // column check
   if(NCOL!=fNumCol && fHeaderExists){
      std::cout << "[CSVManager::ReadFile]: ERROR!  Number of headers doesn't match number of columns!" << std::endl;
   }
   
   char msg[200];
   if(fHeaderExists)  sprintf(msg,"[CSVManager::ReadFile]: Found header, %d rows, %d columns"   ,NROW,NCOL);
   if(!fHeaderExists) sprintf(msg,"[CSVManager::ReadFile]: Found NO header, %d rows, %d columns",NROW,NCOL);
   if(fVerbosity>0) std::cout << msg << std::endl;
   
   return 0;
}
