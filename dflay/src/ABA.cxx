#include "../include/ABA.h"
//______________________________________________________________________________
ABA::ABA(){
   fUseTimeWeight = false;
   fVerbosity     = 0;
}
//______________________________________________________________________________
ABA::~ABA(){

}
//______________________________________________________________________________
int ABA::GetDifference(std::vector<double> A_time,std::vector<double> A,std::vector<double> A_err,
                       std::vector<double> B_time,std::vector<double> B,std::vector<double> B_err,
                       std::vector<double> &diff_aba,std::vector<double> &diff_aba_err){
   // WARNING: This assumes that the A measurement comes first!
   char msg[200];
   double w=0,w_prev=0,dt_tot=0;
   double diff=0,diff_prev=0;
   double arg=0,arg_err=0,arg_err_sq=0;
   const int N = B.size();
   for(int i=1;i<N;i++){
      // set up time weights 
      dt_tot    = A_time[i] - A_time[i-1];
      w         = fabs(A_time[i-1]-B_time[i-1])/dt_tot;
      w_prev    = fabs(A_time[i]-B_time[i-1])/dt_tot;
      // first compute the difference 
      diff_prev = A[i-1] - B[i-1];
      diff      = A[i]   - B[i-1];
      // now get the ABA difference
      if(fUseTimeWeight){
	 arg        = w*diff + w_prev*diff_prev;
	 // arg_err_sq = B_err[i-1]*B_err[i-1] + w*w*A_err[i]*A_err[i] + w_prev*w_prev*A_err[i-1]*A_err[i-1];
	 // account for correlations in measurements 
	 arg_err_sq = (1. - 2.*w*w_prev)*B_err[i-1]*B_err[i-1] + w*w*A_err[i]*A_err[i]
	    + w_prev*w_prev*A_err[i-1]*A_err[i-1];
	 arg_err    = sqrt(arg_err_sq);
      }else{
	 arg     = 0.5*(diff + diff_prev);
	 arg_err = sqrt( B_err[i-1]*B_err[i-1] + 0.25*A_err[i]*A_err[i] + 0.25*A_err[i-1]*A_err[i-1] );
      }
      if(fVerbosity>0){
	 sprintf(msg,"[ABA::GetDifference]: Trial %02d: ",i);
	 std::cout << msg << std::endl;
	 sprintf(msg,"A1 = %.3lf +/- %.3lf Hz, ",A[i-1],A_err[i-1]);
	 std::cout << msg << std::endl;
	 sprintf(msg,"B = %.3lf +/- %.3lf Hz, " ,B[i-1],B_err[i-1]);
	 std::cout << msg << std::endl;
	 sprintf(msg,"A2 = %.3lf +/- %.3lf Hz, ",A[i]  ,A_err[i]  );
	 std::cout << msg << std::endl;
	 sprintf(msg,"diff = %.3lf +/- %.3lf Hz",arg   ,arg_err   );
	 std::cout << msg << std::endl;
      }
      // store result 
      diff_aba.push_back(arg);
      diff_aba_err.push_back(arg_err);
   }
   return 0;
}
