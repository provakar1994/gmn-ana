#ifndef UTIL_CONSTANTS_H
#define UTIL_CONSTANTS_H 

namespace Constants { 
   enum statsType {
      kBesselsCorrectionDisabled = 0,
      kBesselsCorrectionEnabled  = 1
   };

   static const int FTSIZE           = 1E+7;

   const static double second        = 1.;
   const static double minute        = 60.*second;
   const static double hour          = 60.*minute;
   const static double day           = 24.*hour;

   const static double joule_per_gev = 1.60218E-10;

   static const double alpha         = 1.0/137.0359895;
   static const double PI            = 3.14159265;

   static const double h_bar         = 1.054571817E-34;  // [Joule-sec]
   static const double k_Boltzmann   = 1.380649E-23;     // [Joule/Kelvin] 
   static const double N_Avogadro    = 6.022E+23;

   static const double radian        = 1.;
   static const double deg           = 0.017453293;      // 1 deg = pi/180 radians 
   static const double electron_mass = 0.511E-3;         // GeV 
   static const double proton_mass   = 0.9383;           // GeV 
   static const double pion_mass     = 0.140;            // GeV

   static const double mu_0          = 1.2566371e-06;    // 4pi*1E-7

   static const double mu_N          = 5.050783699E-27;  // nuclear magneton [Joule/Tesla] ± 31 in last two digits  
   static const double mu_B          = 9.274009994E-24;  // Bohr magneton [Joule/Tesla] ± 57 in last two digits  

   static const double Gauss         = 1E-4;             // 1 Gauss    = 1E-4 Tesla 
   static const double Oerstead      = 79.577471546;     // 1 Oerstead = 79.6 A/m 

   static const double g_p           = 5.5856946893;     // dg_p/g_p = ± 2.9E-10  

}

#endif 
