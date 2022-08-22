#ifndef ABA_HH
#define ABA_HH

// a class that calculates the ABA difference between a series of measurements
// removes zero-point drift from data 

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

class ABA {

   private:
      int fVerbosity;
      bool fUseTimeWeight;

   public:
      ABA();
      ~ABA();

      void UseTimeWeight(bool v=true) { fUseTimeWeight = v; }
      void SetVerbosity(int v)        { fVerbosity     = v; }

      int GetVerbosity()          const { return fVerbosity;     }
      bool GetTimeWeightStatus()  const { return fUseTimeWeight; }

      int GetDifference(std::vector<double> A_time   ,std::vector<double> A,std::vector<double> A_err,
	                std::vector<double> B_time   ,std::vector<double> B,std::vector<double> B_err,
	                std::vector<double> &diff_aba,std::vector<double> &diff_aba_err);

};

#endif
