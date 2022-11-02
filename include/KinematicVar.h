/* This namespace holds functions to calculate interesting physics quantities for GMn-nTPE 
   Reaction: e + e' -> N + N' (fixed target)
   4-momentum conservation: Pe + Peprime -> PN + PNprime
   Please try to be consistent with the above naming convension.
   -----
   Created P. Datta <pdbforce@jlab.org> 05-27-2022
*/

#ifndef KINE_VAR_H
#define KINE_VAR_H

#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"

#include "Constants.h"
#include "ExpConstants.h"

namespace kine{

  // central scattered e- momentum
  double pcentral(double ebeam, double etheta, std::string Ntype);       // GeV 
  double etheta(TLorentzVector Peprime);  // Scattering angle (rad)
  double ephi(TLorentzVector Peprime);    // Angle of scattering plane (rad)

  // Constructing target nucleon 4-momentum (assuming at rest)
  void SetPN(std::string Ntype, TLorentzVector &PN);

  // expected recoil nucleon momentum
  double pN_expect(double nu,             // energy of the virtual photon (GeV)
		   std::string Ntype);    // type of nucleon in the reaction
		   
  // projected nucleon 3 vector using elastic kinematics (pNhat)
  TVector3 qVect_unit(double Ntheta,      // recoil nucleon theta (rad)
		      double Nphi);       // recoil nucleon phi (rad) 

  // Constructing HCAL co-ordinate system (CoS) in terms of Hall CoS
  void SetHCALaxes(double sbstheta_rad,           // SBS angle (rad)
		   vector<TVector3> &HCAL_axes);  // HCAL axes in order: X, Y, Z (Output)

  // Returns HCAL origin offset vector: A vector pointing from HCAL center 
  // defined by DB xpos and ypos to real HCAL origin.
  TVector3 HCALOriginOffset(vector<TVector3> HCAL_axes, // HCAL CoS axes [in Hall CoS]
			    std::string dataOrsimu);    // Flag to choose "data" or "simu"
  // Constructs HCAL origin vector from vertex in Hall CoS
  void SetHCALorigin(double sbsdist,                                // SBS distance (m)
		     vector<TVector3> HCAL_axes,                    // HCAL CoS axes [in Hall CoS]
		     std::string dataOrsimu,                        // Choose "Data" or "Simulation"
		     TVector3 &HCAL_origin);                        // Output: HCAL origin vector from vertex

  // Get the expected vertical (x) and horizontal (y) positions of the recoil 
  // nucleon at the face of HCAL.
  void GetxyHCALexpect(TVector3 vertex,                 // vertex vector [in Hall CoS]
		       TVector3 pNhat,                  // projected q vector
		       TVector3 HCAL_origin,            // HCAL origin vector [in Hall CoS]
		       vector<TVector3> HCAL_axes,      // HCAL CoS axes [in Hall CoS]
		       vector<double> &xyHCALexpect);   // expected x and y positions (Output)       

  double Q2(double ebeam, double eeprime, double etheta);                // GeV, GeV, rad
  double W2(double ebeam, double eeprime, double Q2, std::string Ntype); // GeV, GeV, GeV2
  double W(double ebeam, double eeprime, double Q2, std::string Ntype);  // GeV, GeV, GeV2

  double Luminosity(double ibeam, std::string targetType);               // A
  
}

#endif
