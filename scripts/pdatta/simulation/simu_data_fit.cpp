/* 
   This macro uses single parameter R to fit simulation dxHCAL plot to data.
   -----
   P. Datta  Created  08-15-2022
*/

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TEventList.h"
#include "TCut.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TLatex.h"
#include <set>
#include <map>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <TStopwatch.h>

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "../../../include/Constants.h"
#include "../../../src/SetROOTVar.cpp"
#include "../../../src/ExpConstants.cpp"

double PI = TMath::Pi();

double Mp = constant::Mp;
double Mn = constant::Mn;

int SBSM=30;

void simu_data_fit( const char *configfilename,
		    double R=1.,
		    const char *outputfilename="simu_data_fit_SBS50.root" )
{
  
  TStopwatch *sw = new TStopwatch();
  sw->Start();

  ifstream configfile(configfilename);

  // so as usual the main things to define are:
  // 1. List of files
  // 2. Global cuts
  // 3. Nominal beam energy (perhaps corrected for average energy loss along half target thickness
  // 4. File name for old momentum coefficients (optional?)
  // 5. BigBite central angle
  // 6. SBS central angle
  // 7. HCAL distance

  // lets read in the real data histograms
  TFile *fdata = new TFile("../pdout/sbs4_sbs50p_ld2.root");

  TChain *C = new TChain("T");
  
  TString currentline;

  cout << endl << "Chaining all the ROOT files.." << endl;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline);
    }
  }

  TCut globalcut="";
  
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
    }
  }
  
  int model=0; //model for calculation
  double ebeam=1.916;
  double bbtheta = 51.0*TMath::DegToRad();
  double sbstheta = 34.5*TMath::DegToRad();
  double sbsdist = 2.25; //m
  double hcaldist = 13.5; //m
  double Ltgt = 15.0; //cm
  double rho_tgt = 0.0723; //g/cc
  double rho_Al = 2.7; //g/cc

  double dipGap = 1.22; //m
  double sbsfield_frac = 1.;
  double sbsmaxfield = 3.1 * atan( 0.85/(11.0 - 2.25 - dipGap/2.0 ))/0.3/dipGap/0.7;
  
  double celldiameter = 1.6*2.54; //cm, right now this is a guess

  //Eventually we will grab HALLA:p from the EPICs tree for the beam energy:
  double Ztgt = 1.0;
  double Atgt = 1.0;
  double Mmol_tgt = 1.008; //g/mol
  
  //For energy-loss correction to beam energy:
  double dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy:
  double dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV:
  double uwallthick_LH2 = 0.0145; //cm
  double dwallthick_LH2 = 0.015; //cm
  double cellthick_LH2 = 0.02; //cm, this is a guess;
  int useAlshield = 0; //Use 1/8" aluminum shield on scattering chamber exit?
  double Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch 
  //Mean energy loss of the beam prior to the scattering:
  double MeanEloss = Ltgt/2.0 * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV:
  double MeanEloss_outgoing = celldiameter/2.0/sin(bbtheta) * rho_tgt * dEdx_tgt; //Approximately 1 MeV  


  int usehcalcut = 0, usepcut=0, usencut=0; //default to false:
  double dx0=0.0, dy0=0.0, dx0_p=0.0, dy0_p=0.0, dx0_n=0.0, dy0_n=0.0; //zero offsets for deltax and deltay cuts:
  double dxsigma = 0.06, dysigma = 0.1, dxsigma_p = 0.0, dysigma_p = 0.0, dxsigma_n = 0.0, dysigma_n = 0.0; //sigmas for deltax and deltay cuts:

  double GEMpitch = 10.0*TMath::DegToRad();

  double dpelmin_fit = -0.03;
  double dpelmax_fit = 0.03;
  double Wmin = 0.86;
  double Wmax = 1.02; 

  //To calculate a proper dx, dy for HCAL, we need
  double hcalheight = 0.0; //-0.2897; //m (we are guessing that this is the height of the center of HCAL above beam height:
  //The following are the positions of the "first" row and column from HCAL database (top right block as viewed from upstream)
  double xoff_hcal = -2.090; //0.92835;
  double yoff_hcal = -0.825; //0.47305; 
  double blockspace_hcal = 0.152; //0.15254; //don't need to make this configurable
  int nrows_hcal=24;
  int ncols_hcal=12;

  //By default we fix the zero-order coefficient for pth and the first-order pth coefficient for xfp
  int fix_pth0_flag = 1;
  int fix_pthx_flag = 1;
  double pth0 = 0.275;
  double pthx = 0.102;
  
  double BBdist = 1.85;
  
  double thtgtmin_fit = -0.15;
  double thtgtmax_fit = 0.15;

  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endconfig") ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");

      int ntokens = tokens->GetEntries();

      if( ntokens >= 2 ){
	TString skey = ( (TObjString*) (*tokens)[0] )->GetString();
	if( skey == "SBSM" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  SBSM = stemp.Atoi();
	}
	
	if( skey == "model" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  model = stemp.Atoi();
	}
		
	if( skey == "dEdx" ){ //assumed to be given in GeV/(g/cm^2)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dEdx_tgt = stemp.Atof();
	}

	if( skey == "dEdx_Al" ){ //assumed to be given in GeV/(g/cm^2)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dEdx_Al = stemp.Atof();
	}

	if( skey == "useAlShield" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  useAlshield = stemp.Atoi();
	}
	
	if( skey == "ebeam" ){ //primary beam energy (GeV):
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ebeam = stemp.Atof();
	}

	if( skey == "bbtheta" ){ //BigBite central angle (deg) (assumed to be on beam left)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  bbtheta = stemp.Atof() * TMath::DegToRad();
	}

	if( skey == "bbdist" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  BBdist = stemp.Atof();
	}

	if( skey == "sbstheta" ){ //SBS/HCAL central angle (deg) (assumed to be on beam right)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sbstheta = stemp.Atof() * TMath::DegToRad();
	}
	
	if( skey == "sbsmaxfield" ){ //100% 48D48 field: 1.5T.m/1.22m = 1.23T 
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sbsmaxfield = stemp.Atof();
	}

	if( skey == "sbsfield_frac" ){ //fraction of 100% 48D48 field 
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sbsfield_frac = stemp.Atof();
	}

	if( skey == "dipGap" ){ //internal gap in 48D48 (m)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dipGap = stemp.Atof();
	}
	
	if( skey == "sbsdist" ){ //48D48 magnet distance from target (m)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sbsdist = stemp.Atof();
	}

	if( skey == "hcaldist" ){ //HCAL distance from target (m)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  hcaldist = stemp.Atof();
	}

	if( skey == "Ltgt" ){ //target length (cm):
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Ltgt = stemp.Atof();
	}

	if( skey == "rho_tgt" ){ //target density (g/cm^3)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  rho_tgt = stemp.Atof();
	}

	if( skey == "celldiameter" ){ //cell diameter (cm)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  celldiameter = stemp.Atof();
	}
	
	
	if( skey == "usehcalcut" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  usehcalcut = stemp.Atoi();
	}

	// **** ---------
	if( skey == "usepcut" ){ //calculate dxHCAL, dyHCAL for proton and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  usepcut = stemp.Atoi();
	}
	if( skey == "usencut" ){ //calculate dxHCAL, dyHCAL for neutron and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  usencut = stemp.Atoi();
	}
	if( skey == "dx0_p" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dx0_p = stemp.Atof();
	}

	if( skey == "dy0_p" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dy0_p = stemp.Atof();
	}

	if( skey == "dxsigma_p" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dxsigma_p = stemp.Atof();
	}

	if( skey == "dysigma_p" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dysigma_p = stemp.Atof();
	}
	if( skey == "dx0_n" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dx0_n = stemp.Atof();
	}

	if( skey == "dy0_n" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dy0_n = stemp.Atof();
	}

	if( skey == "dxsigma_n" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dxsigma_n = stemp.Atof();
	}

	if( skey == "dysigma_n" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dysigma_n = stemp.Atof();
	}	
	// **** ---------
	
	if( skey == "GEMpitch" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  GEMpitch = stemp.Atof()*TMath::DegToRad();
	}

	if( skey == "dpel_min" ){ //p/pel - 1 minimum to include in fit
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dpelmin_fit = stemp.Atof();
	}

	if( skey == "dpel_max" ){ //p/pel - 1 maximum to include in fit
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dpelmax_fit = stemp.Atof();
	}

	if( skey == "Wmin" ){ //Wmin to include in fit
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Wmin = stemp.Atof();
	}

	if( skey == "Wmax" ){ //Wmax to include in fit
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Wmax = stemp.Atof();
	}

	if( skey == "fix_pth0" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  fix_pth0_flag = stemp.Atoi();
	}

	if( skey == "fix_pthx" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  fix_pthx_flag = stemp.Atoi();
	}

	if( skey == "pth0" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pth0 = stemp.Atof();
	}

	if( skey == "pthx" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pthx = stemp.Atof();
	}

	if( skey == "fit_thtgt_min" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  thtgtmin_fit = stemp.Atof();
	}

	if( skey == "fit_thtgt_max" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  thtgtmax_fit = stemp.Atof();
	}
	
      }
      
      tokens->Delete();
    }
  }

  // if(SBSM==30)
  //   //outputfilename = "test_dig_sbs4_30p_ld2.root";
  //   outputfilename = "test.root";
  // if(SBSM==50)
  //   outputfilename = "simu_sbs4_50p_ld2_imac.root";
  
  
  //Note that both of these calculations neglect the Aluminum end windows and cell walls:
  
  MeanEloss = Ltgt/2.0 * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //hopefully the user has provided these quantities with correct units
  
  MeanEloss_outgoing = celldiameter/2.0/sin(bbtheta) * rho_tgt * dEdx_tgt +
    cellthick_LH2/sin(bbtheta) * rho_Al * dEdx_Al;
  if( useAlshield != 0 ) MeanEloss_outgoing += Alshieldthick * rho_Al * dEdx_Al;

  cout << "-----" << endl;
  cout << " Mean eloss in aluminum shield = " << Alshieldthick * rho_Al * dEdx_Al << endl;
  cout << " Use Aluminum shield = " << useAlshield << endl;
  cout << " Mean E loss (beam electron) = " << MeanEloss << endl;
  cout << " BigBite theta = " << bbtheta * TMath::RadToDeg() << endl;
  cout << " Mean E loss (outgoing electron) = " << MeanEloss_outgoing << endl;
  cout << " Ebeam corrected (GeV) = " << ebeam - MeanEloss << endl;
  cout << "-----" << endl;

  // In order to get 

  TEventList *elist = new TEventList("elist","Event list for BigBite momentum calibration");
  
  cout << endl << "Creating a list of events that passed globalcut.." << endl;
  C->Draw(">>elist",globalcut);

  int maxNtr=1000;
  C->SetBranchStatus("*",0);
  //bbsh clus var
  double ESH, xSH, ySH;
  std::vector<std::string> shclvar = {"e","x","y"};
  std::vector<void*> shclvar_mem = {&ESH,&xSH,&ySH};
  setrootvar::setbranch(C,"bb.sh",shclvar,shclvar_mem);

  //bbps clus var
  double EPS, xPS, yPS;
  std::vector<std::string> psclvar = {"e","x","y"};
  std::vector<void*> psclvar_mem = {&EPS,&xPS,&yPS};
  setrootvar::setbranch(C,"bb.ps",psclvar,psclvar_mem);

  //hcal clus var
  double EHCAL, xHCAL, yHCAL;
  std::vector<std::string> hcalclvar = {"e","x","y"};
  std::vector<void*> hcalclvar_mem = {&EHCAL,&xHCAL,&yHCAL};
  setrootvar::setbranch(C,"sbs.hcal",hcalclvar,hcalclvar_mem);

  //track var
  double ntrack, p[maxNtr],px[maxNtr],py[maxNtr],pz[maxNtr];
  double vx[maxNtr],vy[maxNtr],vz[maxNtr];
  double xfp[maxNtr],yfp[maxNtr],thfp[maxNtr],phfp[maxNtr];
  double xtgt[maxNtr],ytgt[maxNtr],thtgt[maxNtr],phtgt[maxNtr];
  std::vector<std::string> trvar = {"n","p","px","py","pz","vx","vy","vz",
  				    "r_x","r_y","r_th","r_ph",
  				    "tg_x","tg_y","tg_th","tg_ph"};
  std::vector<void*> trvar_mem = {&ntrack,&p,&px,&py,&pz,&vx,&vy,&vz,
				  &xfp,&yfp,&thfp,&phfp,
  				  &xtgt,&ytgt,&thtgt,&phtgt};
  setrootvar::setbranch(C,"bb.tr",trvar,trvar_mem);

  //tdctrig variable
  // int tdcElemN;
  // double tdcTrig[maxNtr], tdcElem[maxNtr];
  // std::vector<std::string> tdcvar = {"tdcelemID","tdcelemID","tdc"};
  // std::vector<void*> tdcvar_mem = {&tdcElem,&tdcElemN,&tdcTrig};
  // setrootvar::setbranch(C,"bb.tdctrig",tdcvar,tdcvar_mem,1);

  //MC variables
  double mc_sigma, mc_omega, mc_fnucl;
  std::vector<std::string> mc = {"mc_sigma","mc_omega","mc_fnucl"};
  std::vector<void*> mc_mem = {&mc_sigma,&mc_omega,&mc_fnucl};
  setrootvar::setbranch(C,"MC",mc,mc_mem);

  double pcentral = ebeam/(1.+ebeam/Mp*(1.-cos(bbtheta)));

  TFile *fout = new TFile(outputfilename,"RECREATE");

  //defining all the interesting histograms  
  TH1D *h_W = new TH1D("h_W",";W (GeV);",250,0,2);
  TH1D *h_dpel = new TH1D("h_dpel",";p/p_{elastic}(#theta)-1;",250,-0.25,0.25);

  TH1D *h_dxHCAL_n = new TH1D("h_dxHCAL_n","Simulation n; x_{HCAL} - x_{exp} (m);", 90, -2.5,2.0);
  TH1D *h_dxHCAL_p = new TH1D("h_dxHCAL_p","Simulation p; x_{HCAL} - x_{exp} (m);", 90, -2.5,2.0);
  TH1D *h_n = new TH1D("h_n","Simulation n; x_{HCAL} - x_{exp} (m);", 90, -2.5,2.0);
  TH1D *h_p = new TH1D("h_p","Simulation p; x_{HCAL} - x_{exp} (m);", 90, -2.5,2.0);
  TH1D *h_comb_MC_ap = new TH1D("h_comb_MC_ap","Simulation; x_{HCAL} - x_{exp} (m);", 90, -2.5,2.0);
  TH1D *h_comb_MC = new TH1D("h_comb_MC","Simulation; x_{HCAL} - x_{exp} (m);", 90, -2.5,2.0);

  TH1D *h_dxHCAL_data;
  fdata->GetObject("h_dxHCAL_data",h_dxHCAL_data);
  
  long nevent=0;
  long nevents=elist->GetN();

  // looping through all the events
  while( C->GetEntry(elist->GetEntry(nevent++)) ){
    if( nevent % 1000 == 0 ) cout << nevent << "/" << nevents << "\r";
    cout.flush();

    // Calculating weight
    double ngen_total = 100000; //shouldn't be hard-coded
    double I_beam = 1.0E-6; //A, shouldn't be hard-coded
    double lumi = ((I_beam/constant::qe)*expconst::tarlen*expconst::ld2tarrho*(constant::N_A/constant::D2_Amass));
    double weight = mc_sigma*mc_omega*lumi/ngen_total;
      
    //The first thing we want to do is to calculate the "true" electron momentum incident on BigBite:
    double Ebeam_corrected = ebeam - MeanEloss;
      
    double etheta = acos(pz[0]/p[0]);
    double ephi = atan2( py[0], px[0] );    
      
    // Calculate the expected momentum of an elastically scattered electron at the reconstructed scattering angle and then correct it for the mean energy loss of the
    // electron on its way out of the target:
    double pelastic = Ebeam_corrected/(1.+(Ebeam_corrected/Mp)*(1.0-cos(etheta))); 
      
    double precon = p[0] + MeanEloss_outgoing; //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
      
    double nu_recon = Ebeam_corrected - precon;
      
    double Q2recon = 2.0*Ebeam_corrected*precon*(1.0-cos(etheta));
    double W2recon = pow(Mp,2) + 2.0*Mp*nu_recon - Q2recon;
    double Wrecon; //= sqrt( std::max(0.0,W2recon) );
    if(W2recon>0) Wrecon = sqrt(W2recon); 
    h_W->Fill( Wrecon );
      
    double dpel = precon/pelastic-1.0;
    h_dpel->Fill( dpel );
 
    TVector3 vertex( 0, 0, vz[0] );
    TLorentzVector Pbeam(0,0,Ebeam_corrected,Ebeam_corrected);
    TLorentzVector Kprime(px[0],py[0],pz[0],p[0]);
    TLorentzVector Ptarg(0,0,0,Mp);
    TLorentzVector q = Pbeam - Kprime;

    // Let's try using TLorentzVector
    TLorentzVector Pnprime = q + Ptarg;

    // Calculation of theta_pq for elastic, expected to be 0
    TVector3 q_3v = q.Vect();
    TVector3 Pnprime_3v = Pnprime.Vect();
    double thetapq_elastic = acos(( Pnprime_3v.Mag()*q_3v.Mag() )/Pnprime_3v.Dot(q_3v));
    // T_thetapq = theta_pq; //57.29578*theta_pq;
    // ****

    // Usually when we are running this code, the angle reconstruction is already well calibrated, but the momentum reconstruction is
    // unreliable; use pel(theta) as electron momentum for kinematic correlation:
    // Model 0 = uses reconstructed p as independent variable
    // model 1 = uses reconstructed angles as independent variable
    // model 2 = same as model 1 but uses TLorentzVector for the calculation of pp_expect
    double nu;
    if(model==0)
      nu = Ebeam_corrected - p[0];
    else if(model==1)
      nu = Ebeam_corrected - pelastic;
    else if(model==2)
      nu = q.E();
      
    double pp_expect = sqrt(pow(nu,2)+2.*Mp*nu); //sqrt(nu^2 + Q^2) = sqrt(nu^2 + 2Mnu)
    if(model==2) pp_expect = Pnprime.P();
      
    double pphi_expect = ephi + PI;
      
    double ptheta_expect;
    if(model==0)
      //will this give better resolution than methods based on electron angle only? Not entirely clear
      ptheta_expect = acos( (Ebeam_corrected-pz[0])/pp_expect ); 
    else if(model==1 || model==2)
      ptheta_expect = acos( (Ebeam_corrected-pelastic*cos(etheta))/pp_expect );
      
    TVector3 pNhat( sin(ptheta_expect)*cos(pphi_expect), sin(ptheta_expect)*sin(pphi_expect), cos(ptheta_expect) );
    // In simulation HCAL co-ordinate axes follows the same convension as Hall co-ordinate system
    TVector3 HCAL_zaxis(sin(-sbstheta),0,cos(-sbstheta));
    // TVector3 HCAL_yaxis(0,1,0);
    // TVector3 HCAL_xaxis = HCAL_yaxis.Cross(HCAL_zaxis).Unit();
    TVector3 HCAL_xaxis(0,-1,0);
    TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
      
    TVector3 HCAL_origin = hcaldist * HCAL_zaxis + hcalheight * HCAL_xaxis;

    // Expected position of the q vector at HCAL
    double sintersect = (HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot(HCAL_zaxis));
    TVector3 HCAL_intersect = vertex + sintersect * pNhat; //Straight-line proj to the surface of HCAL:
    double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
    double xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );

    // Calculation of thetapq for n using HCAL information for QE scattering
    TVector3 HCAL_pos = HCAL_origin + xHCAL*HCAL_xaxis + yHCAL*HCAL_yaxis;
    TVector3 n_dir = ( HCAL_pos - vertex ).Unit();
 
    // For p theta_pq calculation we need to estimate the deflection
    double BdL = sbsfield_frac*sbsmaxfield*dipGap;
    double proton_thetabend = 0.3*BdL/q_3v.Mag();
    double proton_deflection = tan(proton_thetabend)*( hcaldist - (sbsdist + dipGap/2.0) );
    TVector3 p_dir = ( HCAL_pos + proton_deflection*HCAL_xaxis - vertex ).Unit();


    //**** fiducial cut ****
    bool fidu_cut, HCAL_active;
    double deltax_shift;
    HCAL_active = yHCAL>-0.75 && yHCAL<0.75 && xHCAL>-2.015 && xHCAL<1.285;
    if(SBSM==30){
      deltax_shift=0.65;
      fidu_cut = (yexpect_HCAL-0.26)>-.75 && (yexpect_HCAL+0.26)<.75 && (xexpect_HCAL-deltax_shift-0.175)>-2.015 && (xexpect_HCAL+0.187)<1.285;
    }
    else if(SBSM==50){
      deltax_shift=1.10;
      fidu_cut = (yexpect_HCAL-0.26)>-.75 && (yexpect_HCAL+0.26)<.75 && (xexpect_HCAL-deltax_shift-0.164)>-2.015 && (xexpect_HCAL+0.182)<1.285;
    }
      

    if( Wrecon >= Wmin && Wrecon <= Wmax && fidu_cut && HCAL_active)
      {
	if(int(mc_fnucl)==0){ // choosing n
	  h_dxHCAL_n->Fill( xHCAL-xexpect_HCAL, weight );
	  h_comb_MC->Fill( xHCAL-xexpect_HCAL, R*weight );
	} 
	if(int(mc_fnucl)==1){ // choosing p
	  h_dxHCAL_p->Fill( xHCAL-xexpect_HCAL, weight );
	  h_comb_MC->Fill( xHCAL-xexpect_HCAL, weight );
	}
      }
   
  } // Event loop

  cout << endl;

  // Scaling simulation histo which we got using my method
  h_comb_MC->Scale(1./h_comb_MC->Integral());

  // leaving the tails for now. Need to implement radiative correction.
  int lbin = h_dxHCAL_data->FindBin(-1.7);
  int hbin = h_dxHCAL_data->FindBin(0.5);

  double norm = 1./(h_dxHCAL_p->Integral() + R*h_dxHCAL_n->Integral());
  // Looping over bins to create combine simulation histo using Andrew's method
  for(int ibin=lbin;ibin<hbin;ibin++){
    double simu = norm*(h_dxHCAL_p->GetBinContent(ibin) + R*h_dxHCAL_n->GetBinContent(ibin));
    h_comb_MC_ap->SetBinContent(ibin, simu);
    h_n->SetBinContent(ibin, norm*(R*h_dxHCAL_n->GetBinContent(ibin)));
    h_p->SetBinContent(ibin, norm*(h_dxHCAL_p->GetBinContent(ibin)));
  }

  TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
  c1->Divide(2,1);

  c1->cd(1);
  h_dxHCAL_data->Draw();
  h_comb_MC->Draw("same");
  h_comb_MC->SetLineColor(1);
  h_comb_MC->SetMarkerStyle(8);
  h_n->Draw("same");
  h_p->Draw("same");

  c1->cd(2);
  h_dxHCAL_data->Draw();
  h_comb_MC_ap->Draw("same");

  // looping over bins again to calculate chi2
  double chi2 = 0., chi2_ap = 0.;
  for(int ibin=lbin;ibin<hbin;ibin++)
    {
      double simu = h_comb_MC->GetBinContent(ibin);
      double simu_ap = h_comb_MC_ap->GetBinContent(ibin);
      double data = h_dxHCAL_data->GetBinContent(ibin);
      if(data>0)
	{ 
	  double dataErr = sqrt(data)/sqrt(h_dxHCAL_data->GetEntries());
	  chi2 += (data-simu)*(data-simu) / ((dataErr)*(dataErr));
	  chi2_ap += (data-simu_ap)*(data-simu_ap) / ((dataErr)*(dataErr));
	}
    }

  cout << endl << "---" << endl
       << "R = " << R << " chi2 = " << chi2 << " chi2_ap " 
       << chi2_ap << " lbin = " << lbin << " hbin = " << hbin
       << endl << "---" << endl;

  ofstream datafile1("chi2_vs_R_sbs50_.82T.csv", ios_base::app | ios_base::out);
  datafile1 << R << "," << chi2 << endl;
  ofstream datafile2("chi2_vs_R_sbs50_.82T_ap.csv", ios_base::app | ios_base::out);
  datafile2 << R << "," << chi2_ap << endl;

  h_dxHCAL_n->Scale(1./h_dxHCAL_n->Integral());
  h_dxHCAL_p->Scale(1./h_dxHCAL_p->Integral());
  c1->cd(2);
  h_dxHCAL_data->Draw();
  h_comb_MC_ap->Draw("same");
  h_n->Draw("same");
  h_p->Draw("same");

  TString plotsfilename = outputfilename;
  plotsfilename.ReplaceAll(".root",".pdf");
  
  c1->Print(plotsfilename.Data(),"pdf");
  plotsfilename.ReplaceAll(".pdf",".png");
  c1->Print(plotsfilename.Data(),"png");

  fout->Write(); //fout->Close();
  
  elist->Delete();
  //fout->Delete();

  sw->Stop();
  cout << "CPU time " << sw->CpuTime() << endl;

}
