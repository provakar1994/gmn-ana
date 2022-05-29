/* 
   This macro plots proton and neutron spots on HCal and
   helps determining the QE event selection cuts for p and n.
   -----
   P. Datta  Created  03-20-2022 (Based on AJR Puckett's momentum calibration script)
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

#include "../../src/SetROOTVar.cpp"

double PI = TMath::Pi();

double Mp = 0.938272;
double Mn = 0.939565;

int SBSM=30;

//gaussian background function for deltax
bool reject_bg = false;
double bg_fit(double *x, double *par){
  if (SBSM==30 && reject_bg && x[0]>-1.2 && x[0]<0.48) { 
    TF1::RejectPoint();
    return 0;
  }
  if (SBSM==50 && reject_bg && x[0]>-1.77 && x[0]<0.51){ //x[0]>-1.65 && x[0]<0.556) {  
    TF1::RejectPoint();
    return 0;
  }
  return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2.));
}

// double W_n_fit_total(double *x, double *par){
//   if (x[0]>0.07 && x[0]<0.76) {
//     return par[0] + par[1]*x[0] + par[2]*pow(x[0],2.) + par[3]*pow(x[0],3.);
//   }else if (x[0]>0.76 && x[0]<1.09) {
//     return par[4]*exp(-0.5*pow((x[0]-par[5])/par[6],2.));
//   }else if (x[0]>1.09 && x[0]<1.75) {
//     return par[7] + par[8]*x[0] + par[9]*pow(x[0],2.) + par[10]*pow(x[0],3.);
//   }else{
//     return 0;
//   }
// }

void hcal_dxdy_simu( const char *configfilename,
		    const char *outputfilename="hcal_dxdy_ld2.root" )
{
  
  TStopwatch *sw = new TStopwatch();
  sw->Start();

  ifstream configfile(configfilename);

  //so as usual the main things to define are:
  // 1. List of files
  // 2. Global cuts
  // 3. Nominal beam energy (perhaps corrected for average energy loss along half target thickness
  // 4. File name for old momentum coefficients (optional?)
  // 5. BigBite central angle
  // 6. SBS central angle
  // 7. HCAL distance

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

  if(SBSM==30)
    //outputfilename = "test_dig_sbs4_30p_ld2.root";
    outputfilename = "test.root";
  if(SBSM==50)
    outputfilename = "test_dig_sbs4_50p_ld2.root";
  
  
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

  int MAXNTRACKS=1000;
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
  double ntrack, p[MAXNTRACKS],px[MAXNTRACKS],py[MAXNTRACKS],pz[MAXNTRACKS];
  double vx[MAXNTRACKS],vy[MAXNTRACKS],vz[MAXNTRACKS];
  double xfp[MAXNTRACKS],yfp[MAXNTRACKS],thfp[MAXNTRACKS],phfp[MAXNTRACKS];
  double xtgt[MAXNTRACKS],ytgt[MAXNTRACKS],thtgt[MAXNTRACKS],phtgt[MAXNTRACKS];
  std::vector<std::string> trvar = {"n","p","px","py","pz",
  				    "vx","vy","vz",
  				    "r_x","r_y","r_th","r_ph",
  				    "tg_x","tg_y","tg_th","tg_ph"};
  std::vector<void*> trvar_mem = {&ntrack,&p,&px,&py,&pz,
  				  &vx,&vy,&vz,&xfp,&yfp,&thfp,&phfp,
  				  &xtgt,&ytgt,&thtgt,&phtgt};
  setrootvar::setbranch(C,"bb.tr",trvar,trvar_mem);

  //tdctrig variable
  // int tdcElemN;
  // double tdcTrig[MAXNTRACKS], tdcElem[MAXNTRACKS];
  // std::vector<std::string> tdcvar = {"tdcelemID","tdcelemID","tdc"};
  // std::vector<void*> tdcvar_mem = {&tdcElem,&tdcElemN,&tdcTrig};
  // setrootvar::setbranch(C,"bb.tdctrig",tdcvar,tdcvar_mem,1);
  

  double pcentral = ebeam/(1.+ebeam/Mp*(1.-cos(bbtheta)));

  TFile *fout = new TFile(outputfilename,"RECREATE");

  //defining all the interesting histograms
  TH1D *h_coin_time = new TH1D("h_coin_time",";HCAL trigTime - BBCAL trigTime (ns);",150,450,600);
  
  TH1D *h_W = new TH1D("h_W",";W (GeV);",250,0,2);
  TH1D *h_W_n_cut = new TH1D("h_W_n_cut",";W (GeV);",250,0,2);
  TH1D *h_W_p_cut = new TH1D("h_W_p_cut",";W (GeV);",250,0,2);

  TH1D *h_dpel = new TH1D("h_dpel",";p/p_{elastic}(#theta)-1;",250,-0.25,0.25);
  TH1D *h_dpel_n_cut = new TH1D("h_dpel_n_cut",";p/p_{elastic}(#theta)-1;",250,-0.25,0.25);
  TH1D *h_dpel_p_cut = new TH1D("h_dpel_p_cut",";p/p_{elastic}(#theta)-1;",250,-0.25,0.25);

  TH1D *h_dxHCAL = new TH1D("h_dxHCAL","; x_{HCAL} - x_{exp} (m);", 250, -2.5,2.5);
  TH1D *h_dyHCAL = new TH1D("h_dyHCAL","; y_{HCAL} - y_{exp} (m);", 250, -1.25,1.25);

  TH1D *h_dxHCAL_cut = new TH1D("h_dxHCAL_cut","; x_{HCAL} - x_{exp} (m);", 250, -2.5,2.5);
  // TH1D *h_dxHCAL_cut_1 = new TH1D("h_dxHCAL_cut_1","; x_{HCAL} - x_{exp} (m);", 250, -2.5,2.5);
  // TH1D *h_dxHCAL_cut_2 = new TH1D("h_dxHCAL_cut_2","; x_{HCAL} - x_{exp} (m);", 250, -2.5,2.5);
  // TH1D *h_dxHCAL_cut_3 = new TH1D("h_dxHCAL_cut_3","; x_{HCAL} - x_{exp} (m);", 250, -2.5,2.5);
  // TH1D *h_dxHCAL_cut_4 = new TH1D("h_dxHCAL_cut_4","; x_{HCAL} - x_{exp} (m);", 250, -2.5,2.5);
  // TH1D *h_dxHCAL_cut_5 = new TH1D("h_dxHCAL_cut_5","; x_{HCAL} - x_{exp} (m);", 250, -2.5,2.5);
  TH1D *h_dyHCAL_cut = new TH1D("h_dyHCAL_cut","; y_{HCAL} - y_{exp} (m);", 250, -1.25,1.25);  

  TH2D *h2_dxdyHCAL = new TH2D("h2_dxdyHCAL","; y_{HCAL} - y_{exp} (m); x_{HCAL} - x_{exp} (m)",
			       250, -1.25,1.25, 250,-2.5,2.5);
  TH2D *h2_dxdyHCAL_cut = new TH2D("h2_dxdyHCAL_cut","; y_{HCAL} - y_{exp} (m); x_{HCAL} - x_{exp} (m)",
				   250, -1.25,1.25, 250,-2.5,2.5);
  
  TH2D *h2_xyHCAL = new TH2D("h2_xyHCAL",";y_{HCAL} (m);x_{HCAL} (m)",12,-0.9,0.9,24,-2.165,1.435);
  TH2D *h2_xyHCAL_W_cut = new TH2D("h2_xyHCAL_W_cut",";y_{HCAL} (m);x_{HCAL} (m)"
				   ,12,-0.9,0.9,24,-2.165,1.435);
  TH2D *h2_xyHCAL_p_cut = new TH2D("h2_xyHCAL_p_cut","x_{HCAL} vs y_{HCAL} w/ W, p spot, & Fiducial cut;y_{HCAL} (m);x_{HCAL} (m)"
				   ,12,-0.9,0.9,24,-2.165,1.435);
  TH2D *h2_xyHCAL_n_cut = new TH2D("h2_xyHCAL_n_cut","x_{HCAL} vs y_{HCAL} w/ W, n spot, & Fiducial cut;y_{HCAL} (m);x_{HCAL} (m)"
				   ,12,-0.9,0.9,24,-2.165,1.435);

  TH2D *h2_xyBB_p_exp = new TH2D("h2_xyBB_p_exp","x_{exp}-0.65 vs y_{exp} w/ W & Fiducial cut (exp. p pos.);y_{exp} (m);x_{exp}-0.65 (m)"
				   ,12,-0.9,0.9,24,-2.165,1.435);
  TH2D *h2_xyBB_n_exp = new TH2D("h2_xyBB_n_exp","x_{exp} vs y_{exp} w/ W & Fiducial cut (exp. n pos.);y_{exp} (m);x_{exp} (m)"
				   ,12,-0.9,0.9,24,-2.165,1.435);

  TH2D *h2_W_vs_dxHCAL = new TH2D("h2_W_vs_dxHCAL","; x_{HCAL} - x_{exp} (m); W (GeV)", 250,-2.5,2.5,250,0,2);

  // TH2D *h2_dpel_xfp = new TH2D("h2_dpel_xfp", "; xfp (m); p/p_{elastic}(#theta)-1",
  // 			       250,-0.75,0.75,250,-0.125,0.125);
  // TH2D *h2_dpel_yfp = new TH2D("h2_dpel_yfp", "; yfp (m); p/p_{elastic}(#theta)-1",
  // 			       250,-0.25,0.25,250,-0.125,0.125);
  // TH2D *h2_dpel_xpfp = new TH2D("h2_dpel_xpfp", "; xpfp (m); p/p_{elastic}(#theta)-1",
  // 				250,-0.4,0.4,250,-0.125,0.125);
  // TH2D *h2_dpel_ypfp = new TH2D("h2_dpel_ypfp", "; ypfp (m); p/p_{elastic}(#theta)-1",
  // 				250,-0.15,0.15,250,-0.125,0.125);
  // TH2D *h2_dpel_xptar = new TH2D("h2_dpel_xptar", "; xptar; p/p_{elastic}(#theta)-1",
  // 				 250,-0.4,0.4,250,-0.125,0.125);
  // TH2D *h2_dpel_yptar = new TH2D("h2_dpel_yptar", "; yptar; p/p_{elastic}(#theta)-1",
  // 				 250,-0.15,0.15,250,-0.125,0.125);
  // TH2D *h2_dpel_ytar = new TH2D("h2_dpel_ytar", "; ytar; p/p_{elastic}(#theta)-1",
  // 				250,-0.15,0.15,250,-0.125,0.125);
  // TH2D *h2_W_xfp = new TH2D("h2_W_xfp", "; xfp (m); W",250,-0.75,0.75,250,0,2);
  // TH2D *h2_W_yfp = new TH2D("h2_W_yfp", "; yfp (m); W",250,-0.25,0.25,250,0,2);
  // TH2D *h2_W_xpfp = new TH2D("h2_W_xpfp", "; xpfp (m); W",250,-0.4,0.4,250,0,2);
  // TH2D *h2_W_ypfp = new TH2D("h2_W_ypfp", "; ypfp (m); W",250,-0.15,0.15,250,0,2);
  // TH2D *h2_W_xptar = new TH2D("h2_W_xptar", "; xptar; W",250,-0.4,0.4,250,0,2);
  // TH2D *h2_W_yptar = new TH2D("h2_W_yptar", "; yptar; W",250,-0.15,0.15,250,0,2);
  // TH2D *h2_W_ytar = new TH2D("h2_W_ytar", "; ytar; W",250,-0.15,0.15,250,0,2);

  
  TTree *Tout = new TTree("Tout","Tree containing variables for momentum calibration");
  double T_ebeam, T_etheta, T_ephi, T_precon, T_pelastic, T_thetabend, T_dpel, T_W2;
  double T_thetapq_n, T_thetapq_p;
  double T_pincident;
  double T_xfp, T_yfp, T_thfp, T_phfp;
  double T_thtgt, T_phtgt, T_ytgt, T_xtgt;
  double T_vx, T_vy, T_vz;
  double T_BBdist, T_BBtheta;
  double T_HCALdist, T_HCALtheta;
  double T_xHCAL, T_yHCAL, T_EHCAL, T_deltax, T_deltay;
  double T_pp_expect, T_ptheta_expect, T_pphi_expect;
  double T_EPS, T_ESH, T_Etot;
  double T_xSH, T_ySH;
  double T_Q2;
  int pcut, ncut;
  int BBcut;
  
  
  Tout->Branch( "pcut", &pcut, "pcut/I");
  Tout->Branch( "ncut", &ncut, "ncut/I");
  Tout->Branch( "BBcut", &BBcut, "BBcut/I");
  Tout->Branch( "Ebeam", &T_ebeam, "Ebeam/D" );
  Tout->Branch( "Q2", &T_Q2, "Q2/D"); 
  Tout->Branch( "etheta", &T_etheta, "etheta/D");
  Tout->Branch( "ephi", &T_ephi, "ephi/D");
  Tout->Branch( "ep_recon", &T_precon, "ep_recon/D");
  Tout->Branch( "ep_elastic", &T_pelastic, "ep_elastic/D");
  Tout->Branch( "ep_incident", &T_pincident, "ep_incident/D");
  Tout->Branch( "thetabend", &T_thetabend, "thetabend/D");
  Tout->Branch( "dpel", &T_dpel, "dpel/D");
  Tout->Branch( "W2", &T_W2, "W2/D");
  Tout->Branch( "xfp", &T_xfp, "xfp/D");
  Tout->Branch( "yfp", &T_yfp, "yfp/D");
  Tout->Branch( "thfp", &T_thfp, "thfp/D");
  Tout->Branch( "phfp", &T_phfp, "phfp/D");
  Tout->Branch( "thtgt", &T_thtgt, "thtgt/D");
  Tout->Branch( "phtgt", &T_phtgt, "phtgt/D");
  Tout->Branch( "ytgt", &T_ytgt, "ytgt/D");
  Tout->Branch( "xtgt", &T_xtgt, "xtgt/D");
  Tout->Branch( "vz", &T_vz, "vz/D");
  Tout->Branch( "BBdist", &T_BBdist, "BBdist/D");
  Tout->Branch( "BBtheta", &T_BBtheta, "BBtheta/D");
  Tout->Branch( "HCALdist", &T_HCALdist, "HCALdist/D");
  Tout->Branch( "HCALtheta", &T_HCALtheta, "HCALtheta/D");
  Tout->Branch( "xHCAL", &T_xHCAL, "xHCAL/D");
  Tout->Branch( "yHCAL", &T_yHCAL, "yHCAL/D");
  Tout->Branch( "EHCAL", &T_EHCAL, "EHCAL/D");
  Tout->Branch( "deltax", &T_deltax, "deltax/D");
  Tout->Branch( "deltay", &T_deltay, "deltay/D");
  Tout->Branch( "thetapq_n", &T_thetapq_n, "thetapq_n/D");
  Tout->Branch( "thetapq_p", &T_thetapq_p, "thetapq_p/D");
  Tout->Branch( "pp_expect", &T_pp_expect, "pp_expect/D");
  Tout->Branch( "ptheta_expect", &T_ptheta_expect, "ptheta_expect/D");
  Tout->Branch( "pphi_expect", &T_pphi_expect, "pphi_expect/D");
  Tout->Branch( "EPS", &T_EPS, "EPS/D");
  Tout->Branch( "ESH", &T_ESH, "ESH/D");
  Tout->Branch( "Etot", &T_Etot, "Etot/D");
  Tout->Branch( "xSH", &T_xSH, "xSH/D");
  Tout->Branch( "ySH", &T_ySH, "ySH/D");
  
  
  long nevent=0;
  long nevents=elist->GetN();

  // looping through all the events
  while( C->GetEntry(elist->GetEntry(nevent++)) ){
    if( nevent % 1000 == 0 ) cout << nevent << "/" << nevents << "\r";
    cout.flush();

    if( int(ntrack) == 1 ){

      // * ----
      T_xfp = xfp[0];
      T_yfp = yfp[0];
      T_thfp = thfp[0];
      T_phfp = phfp[0];
      T_thtgt = thtgt[0];
      T_phtgt = phtgt[0];
      T_ytgt = ytgt[0];
      //xtgt = -vy - vz * costheta * thtgt
      T_xtgt = -vz[0]*cos(bbtheta)*thtgt[0];
      T_vz = vz[0];
      T_BBdist = BBdist;
      T_BBtheta = bbtheta;
      T_HCALdist = hcaldist;
      T_HCALtheta = sbstheta;

      T_EPS = EPS;
      T_ESH = ESH;
      T_Etot = EPS + ESH;
      T_xSH = xSH;
      T_ySH = ySH;
      // * ----

      // *---- "true" trajectory bend angle for the electron from the reconstructed angles:
      TVector3 enhat_tgt( thtgt[0], phtgt[0], 1.0 );
      enhat_tgt = enhat_tgt.Unit();	
      TVector3 enhat_fp( thfp[0], phfp[0], 1.0 );
      enhat_fp = enhat_fp.Unit();
      TVector3 GEMzaxis(-sin(GEMpitch),0,cos(GEMpitch));
      TVector3 GEMyaxis(0,1,0);
      TVector3 GEMxaxis = (GEMyaxis.Cross(GEMzaxis)).Unit();	
      TVector3 enhat_fp_rot = enhat_fp.X() * GEMxaxis + enhat_fp.Y() * GEMyaxis + enhat_fp.Z() * GEMzaxis;
      double thetabend = acos( enhat_fp_rot.Dot( enhat_tgt ) );
      T_thetabend = thetabend;
      // *----

      //coincidence time analysis
      // double bbcal_time=0., hcal_time=0.;
      // for(int ihit=0; ihit<tdcElemN; ihit++){
      // 	if(tdcElem[ihit]==5) bbcal_time=tdcTrig[ihit];
      // 	if(tdcElem[ihit]==0) hcal_time=tdcTrig[ihit];
      // }
      // double coin_time = hcal_time-bbcal_time;
      // h_coin_time->Fill( coin_time );

      //if( fabs(coin_time-510.)>10. ) continue;
      
      //The first thing we want to do is to calculate the "true" electron momentum incident on BigBite:
      double Ebeam_corrected = ebeam - MeanEloss;
      T_ebeam = Ebeam_corrected;
      
      double etheta = acos(pz[0]/p[0]);
      double ephi = atan2( py[0], px[0] );
      T_etheta = etheta;
      T_ephi = ephi;
      
      // Calculate the expected momentum of an elastically scattered electron at the reconstructed scattering angle and then correct it for the mean energy loss of the
      // electron on its way out of the target:
      double pelastic = Ebeam_corrected/(1.+(Ebeam_corrected/Mp)*(1.0-cos(etheta))); 
      
      double precon = p[0] + MeanEloss_outgoing; //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)

      T_precon = precon;
      T_pelastic = pelastic; 
      T_pincident = pelastic - MeanEloss_outgoing;
      
      double nu_recon = Ebeam_corrected - precon;
      
      double Q2recon = 2.0*Ebeam_corrected*precon*(1.0-cos(etheta));
      double W2recon = pow(Mp,2) + 2.0*Mp*nu_recon - Q2recon;
      double Wrecon; //= sqrt( std::max(0.0,W2recon) );
      if(W2recon>0) Wrecon = sqrt(W2recon); 
      T_Q2 = Q2recon;
      T_W2 = W2recon;
      h_W->Fill( Wrecon );
      
      double dpel = precon/pelastic-1.0;
      T_dpel = dpel;
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
      T_pp_expect = pp_expect;
      T_pphi_expect = pphi_expect;
      
      double ptheta_expect;
      if(model==0)
	//will this give better resolution than methods based on electron angle only? Not entirely clear
	ptheta_expect = acos( (Ebeam_corrected-pz[0])/pp_expect ); 
      else if(model==1 || model==2)
	ptheta_expect = acos( (Ebeam_corrected-pelastic*cos(etheta))/pp_expect );
      
      T_ptheta_expect = ptheta_expect;
      
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
      
      h_dxHCAL->Fill( xHCAL - xexpect_HCAL );
      h_dyHCAL->Fill( yHCAL - yexpect_HCAL );
      h2_dxdyHCAL->Fill( yHCAL - yexpect_HCAL, xHCAL - xexpect_HCAL );
      h2_xyHCAL->Fill(yHCAL,xHCAL);

      // Calculation of thetapq for n using HCAL information for QE scattering
      TVector3 HCAL_pos = HCAL_origin + xHCAL*HCAL_xaxis + yHCAL*HCAL_yaxis;
      TVector3 n_dir = ( HCAL_pos - vertex ).Unit();
      T_thetapq_n = acos( n_dir.Dot(q_3v.Unit()) );

      // For p theta_pq calculation we need to estimate the deflection
      double BdL = sbsfield_frac*sbsmaxfield*dipGap;
      double proton_thetabend = 0.3*BdL/q_3v.Mag();
      double proton_deflection = tan(proton_thetabend)*( hcaldist - (sbsdist + dipGap/2.0) );
      TVector3 p_dir = ( HCAL_pos + proton_deflection*HCAL_xaxis - vertex ).Unit();
      T_thetapq_p = acos( p_dir.Dot(q_3v.Unit()) );

      h2_W_vs_dxHCAL->Fill( xHCAL - xexpect_HCAL, Wrecon );

      if( Wrecon >= Wmin && Wrecon <= Wmax ){
	//h_dxHCAL_cut->Fill( xHCAL - xexpect_HCAL );
	h_dyHCAL_cut->Fill( yHCAL - yexpect_HCAL );
	h2_dxdyHCAL_cut->Fill( yHCAL - yexpect_HCAL, xHCAL - xexpect_HCAL );
	h2_xyHCAL_W_cut->Fill(yHCAL,xHCAL);
      }

      // if( Wrecon >= 0.85 && Wrecon <= 0.9 ){
      // 	h_dxHCAL_cut_1->Fill( xHCAL - xexpect_HCAL );
      // }
      // if( Wrecon >= 0.9 && Wrecon <= 0.95 ){
      // 	h_dxHCAL_cut_2->Fill( xHCAL - xexpect_HCAL );
      // }
      // if( Wrecon >= 0.95 && Wrecon <= 1.0 ){
      // 	h_dxHCAL_cut_3->Fill( xHCAL - xexpect_HCAL );
      // }
      // if( Wrecon >= 1.0 && Wrecon <= 1.05 ){
      // 	h_dxHCAL_cut_4->Fill( xHCAL - xexpect_HCAL );
      // }
      // if( Wrecon >= 1.05 && Wrecon <= 1.1 ){
      // 	h_dxHCAL_cut_5->Fill( xHCAL - xexpect_HCAL );
      // }

      T_xHCAL = xHCAL;
      T_yHCAL = yHCAL;
      T_EHCAL = EHCAL;
      T_deltax = xHCAL - xexpect_HCAL;
      T_deltay = yHCAL - yexpect_HCAL;

      // h2_xyHCAL->Fill(yHCAL,xHCAL);
      
      // * ---- Defining various cuts
      bool passed_p_cut=true, passed_n_cut=true;
      //Cuts to choose only protons
      if( usepcut != 0 ){
	passed_p_cut = pow( (xHCAL-xexpect_HCAL - dx0_p)/dxsigma_p, 2 ) +
	  pow( (yHCAL-yexpect_HCAL - dy0_p)/dysigma_p, 2 ) <= pow(1.5,2); 
      }
      pcut = pow( (xHCAL-xexpect_HCAL - dx0_p)/dxsigma_p, 2 ) +
	pow( (yHCAL-yexpect_HCAL - dy0_p)/dysigma_p, 2 ) <= pow(1.5,2);      

      //Cuts to choose only neutrons
      if( usencut != 0 ){
	passed_n_cut = pow( (xHCAL-xexpect_HCAL - dx0_n)/dxsigma_n, 2 ) +
	  pow( (yHCAL-yexpect_HCAL - dy0_n)/dysigma_n, 2 ) <= pow(1.5,2); 
      }
      ncut = pow( (xHCAL-xexpect_HCAL - dx0_n)/dxsigma_n, 2 ) +
	pow( (yHCAL-yexpect_HCAL - dy0_n)/dysigma_n, 2 ) <= pow(1.5,2); 
      // * ----
	
      BBcut = Wrecon >= Wmin && Wrecon <= Wmax && dpel >= dpelmin_fit && dpel <= dpelmax_fit;

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
      
      if( Wrecon >= Wmin && Wrecon <= Wmax && HCAL_active){
	if(fidu_cut) 
	  {
	    h2_xyBB_p_exp->Fill(yexpect_HCAL,xexpect_HCAL-deltax_shift);  
	    h2_xyBB_n_exp->Fill(yexpect_HCAL,xexpect_HCAL);
	  }
	if( passed_p_cut && fidu_cut)
	  {
	    h2_xyHCAL_p_cut->Fill(yHCAL,xHCAL);
	  }
	if( passed_n_cut && fidu_cut)
	  {
	    h2_xyHCAL_n_cut->Fill(yHCAL,xHCAL);
	  }
      }

      if( Wrecon >= Wmin && Wrecon <= Wmax && fidu_cut && HCAL_active)
	{
	  h_dxHCAL_cut->Fill( xHCAL - xexpect_HCAL );
	}

      //filling W distributions for n and p
      if( passed_n_cut && fidu_cut && HCAL_active){
	h_dpel_n_cut->Fill( dpel );
	h_W_n_cut->Fill( Wrecon );
      }
            
      if( passed_p_cut && fidu_cut && HCAL_active){
	h_dpel_p_cut->Fill( dpel );
	h_W_p_cut->Fill( Wrecon );	
	// h2_dpel_xfp->Fill( xfp[0], dpel );
	// h2_dpel_yfp->Fill( yfp[0], dpel );
	// h2_dpel_xpfp->Fill( thfp[0], dpel );
	// h2_dpel_ypfp->Fill( phfp[0], dpel );
	// h2_dpel_xptar->Fill( thtgt[0], dpel );
	// h2_dpel_yptar->Fill( phtgt[0], dpel );
	// h2_dpel_ytar->Fill( ytgt[0], dpel );		
	// h2_W_xfp->Fill( xfp[0], Wrecon );
	// h2_W_yfp->Fill( yfp[0], Wrecon );
	// h2_W_xpfp->Fill( thfp[0], Wrecon );
	// h2_W_ypfp->Fill( phfp[0], Wrecon );
	// h2_W_xptar->Fill( thtgt[0], Wrecon );
	// h2_W_yptar->Fill( phtgt[0], Wrecon );
	// h2_W_ytar->Fill( ytgt[0], Wrecon );     	
      }

      Tout->Fill();     
    }
  }

  cout << endl;

  TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
  c1->Divide(2,2);

  //illustrating proton and neutron spot cut regions
  c1->cd(1)->SetGrid();
  h2_dxdyHCAL_cut->GetYaxis()->SetTitleOffset(1.3);
  h2_dxdyHCAL_cut->Draw("colz");

  TEllipse Ep;
  Ep.SetFillStyle(0);
  Ep.SetLineColor(2);
  Ep.SetLineWidth(4);
  Ep.DrawEllipse(dy0_p,dx0_p,1.5*dysigma_p,1.5*dxsigma_p,0,360,0);

  TEllipse En;
  En.SetFillStyle(0);
  En.SetLineColor(1);
  En.SetLineWidth(4);
  En.DrawEllipse(dy0_n,dx0_n,1.5*dysigma_n,1.5*dxsigma_n,0,360,0);
  // ******

  //fiducial cut visualization
  c1->cd(4);
  h2_xyBB_p_exp->SetStats(0);
  h2_xyBB_p_exp->GetYaxis()->SetTitleOffset(1.3);
  h2_xyBB_p_exp->Draw("colz");
  // h2_xyHCAL_p_cut->SetStats(0);
  // h2_xyHCAL_p_cut->GetYaxis()->SetTitleOffset(1.3);
  // h2_xyHCAL_p_cut->Draw("colz");

  //active area 
  TLine L1h_p;
  L1h_p.SetLineColor(2); L1h_p.SetLineWidth(4); L1h_p.SetLineStyle(9);
  L1h_p.DrawLine(-0.9,1.285,0.9,1.285);
  TLine L2h_p;
  L2h_p.SetLineColor(2); L2h_p.SetLineWidth(4); L2h_p.SetLineStyle(9);
  L2h_p.DrawLine(-0.9,-2.015,0.9,-2.015);
  
  TLine L1v_p;
  L1v_p.SetLineColor(2); L1v_p.SetLineWidth(4); L1v_p.SetLineStyle(9);
  L1v_p.DrawLine(-0.75,-2.165,-0.75,1.435);
  TLine L2v_p;
  L2v_p.SetLineColor(2); L2v_p.SetLineWidth(4); L2v_p.SetLineStyle(9);
  L2v_p.DrawLine(0.75,-2.165,0.75,1.435);

  //safety margin
  TLine L1h_p_s;
  L1h_p_s.SetLineColor(4); L1h_p_s.SetLineWidth(2); L1h_p_s.SetLineStyle(9);
  L1h_p_s.DrawLine(-0.9,1.098,0.9,1.098);
  TLine L2h_p_s;
  L2h_p_s.SetLineColor(4); L2h_p_s.SetLineWidth(2); L2h_p_s.SetLineStyle(9);
  L2h_p_s.DrawLine(-0.9,-1.84,0.9,-1.84);
  
  TLine L1v_p_s;
  L1v_p_s.SetLineColor(4); L1v_p_s.SetLineWidth(2); L1v_p_s.SetLineStyle(9);
  L1v_p_s.DrawLine(-0.49,-2.165,-0.49,1.435);
  TLine L2v_p_s;
  L2v_p_s.SetLineColor(4); L2v_p_s.SetLineWidth(2); L2v_p_s.SetLineStyle(9);
  L2v_p_s.DrawLine(0.49,-2.165,0.49,1.435);
  
  c1->cd(3);
  h2_xyBB_n_exp->SetStats(0);
  h2_xyBB_n_exp->GetYaxis()->SetTitleOffset(1.3);
  h2_xyBB_n_exp->Draw("colz");
  // h2_xyHCAL_n_cut->SetStats(0);
  // h2_xyHCAL_n_cut->GetYaxis()->SetTitleOffset(1.3);
  // h2_xyHCAL_n_cut->Draw("colz");

  TLine L1h_n;
  L1h_n.SetLineColor(2); L1h_n.SetLineWidth(4); L1h_n.SetLineStyle(9);
  L1h_n.DrawLine(-0.9,1.285,0.9,1.285);
  TLine L2h_n;
  L2h_n.SetLineColor(2); L2h_n.SetLineWidth(4); L2h_n.SetLineStyle(9);
  L2h_n.DrawLine(-0.9,-2.015,0.9,-2.015);
  
  TLine L1v_n;
  L1v_n.SetLineColor(2); L1v_n.SetLineWidth(4); L1v_n.SetLineStyle(9);
  L1v_n.DrawLine(-0.75,-2.165,-0.75,1.435);
  TLine L2v_n;
  L2v_n.SetLineColor(2); L2v_n.SetLineWidth(4); L2v_n.SetLineStyle(9);
  L2v_n.DrawLine(0.75,-2.165,0.75,1.435);

  //safety margin
  TLine L1h_n_s;
  L1h_n_s.SetLineColor(4); L1h_n_s.SetLineWidth(2); L1h_n_s.SetLineStyle(9);
  L1h_n_s.DrawLine(-0.9,1.098,0.9,1.098);
  TLine L2h_n_s;
  L2h_n_s.SetLineColor(4); L2h_n_s.SetLineWidth(2); L2h_n_s.SetLineStyle(9);
  L2h_n_s.DrawLine(-0.9,-1.84,0.9,-1.84);
  
  TLine L1v_n_s;
  L1v_n_s.SetLineColor(4); L1v_n_s.SetLineWidth(2); L1v_n_s.SetLineStyle(9);
  L1v_n_s.DrawLine(-0.49,-2.165,-0.49,1.435);
  TLine L2v_n_s;
  L2v_n_s.SetLineColor(4); L2v_n_s.SetLineWidth(2); L2v_n_s.SetLineStyle(9);
  L2v_n_s.DrawLine(0.49,-2.165,0.49,1.435);  
  // ******

  //fitting delta_x distribution
  c1->cd(2);
  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(1111);
  h_dxHCAL_cut->SetLineWidth(4); h_dxHCAL_cut->SetLineColor(1);

  double bgF_par0, bgF_par1, bgF_par2;
  double pF_low, pF_hi, pF_par0, pF_par1, pF_par2;
  double nF_low, nF_hi, nF_par0, nF_par1, nF_par2;

  //defining fit parameters and range
  if(SBSM==30){
    bgF_par0=480; bgF_par1=-0.3; bgF_par2=0.5;
    pF_low=-1.13; pF_hi=-0.2; pF_par0=1967; pF_par1=-0.58; pF_par2=0.16;
    nF_low=-0.2; nF_hi=0.41; nF_par0=662; nF_par1=-0.047; nF_par2=0.184;
  }else if(SBSM==50){
    bgF_par0=160; bgF_par1=-0.42; bgF_par2=0.7;
    pF_low=-1.4; pF_hi=-0.78; pF_par0=475; pF_par1=-1.09; pF_par2=0.184;
    nF_low=-0.22, nF_hi=0.26; nF_par0=182; nF_par1=0.009; nF_par2=0.169;
  }
    
  double par[9], par_f[9];
  TF1* bg_func = new TF1("bg_func", bg_fit, -2.5, 2.5, 3 );
  bg_func->SetParameters(bgF_par0,bgF_par1,bgF_par2);
  bg_func->SetLineColor(2);
  reject_bg = true;
  h_dxHCAL_cut->Fit(bg_func,"NR");
  bg_func->GetParameters(&par[0]);
  reject_bg = false;

  TF1* p_func = new TF1("p_func", "gaus", pF_low, pF_hi);
  p_func->SetParameters(pF_par0,pF_par1,pF_par2);
  p_func->SetLineColor(3);
  h_dxHCAL_cut->Fit(p_func,"NR+");
  p_func->GetParameters(&par[3]);

  TF1* n_func = new TF1("n_func", "gaus", nF_low, nF_hi);
  n_func->SetParameters(nF_par0,nF_par1,nF_par2);
  n_func->SetLineColor(4);
  h_dxHCAL_cut->Fit(n_func,"NR+");
  n_func->GetParameters(&par[6]);

  TF1* total = new TF1("total","gaus(0)+gaus(3)+gaus(6)",-2.5,2.5);
  total->SetLineColor(2); total->SetLineWidth(4);
  total->SetParameters(par);
  h_dxHCAL_cut->Fit(total,"R+");
  total->GetParameters(par_f);

  //visulaization and extraction of physics results
  double n_count, p_count;
  TF1 *backFcn = new TF1("backFcn","gaus",-2.5,2.5);
  backFcn->SetLineColor(6); backFcn->SetLineWidth(4);
  TF1 *pFcn = new TF1("pFcn","gaus",-2.5,2.5);
  pFcn->SetLineColor(kBlue); pFcn->SetLineWidth(4);
  TF1 *nFcn = new TF1("nFcn","gaus",-2.5,2.5);
  nFcn->SetLineColor(kGreen); nFcn->SetLineWidth(4);

  backFcn->SetParameters(&par_f[0]); backFcn->Draw("same");
  pFcn->SetParameters(&par_f[3]); pFcn->Draw("same");
  nFcn->SetParameters(&par_f[6]); nFcn->Draw("same");
  p_count = (int)pFcn->Integral(-2.5,2.5)/h_dxHCAL_cut->GetBinWidth(1);
  n_count = (int)nFcn->Integral(-2.5,2.5)/h_dxHCAL_cut->GetBinWidth(1);
  
  // draw the legend
  TLegend *legend=new TLegend(0.48,0.66,0.90,0.85);
  legend->SetTextFont(42);
  //legend->SetTextSize(0.02);
  legend->AddEntry(h_dxHCAL_cut,"Data","l");
  legend->AddEntry(total,"Global Fit","l");
  legend->AddEntry(pFcn,Form("p signal fit, #mu = %2.4f, #sigma = %2.4f",par_f[4],par_f[5]),"l");
  legend->AddEntry(nFcn,Form("n signal fit, #mu = %2.4f, #sigma = %2.4f",par_f[7],par_f[8]),"l");
  legend->AddEntry(backFcn,"Background fit","l");
  legend->Draw();  
 
  cout << endl << "------" << endl;
  cout << " n count = " << n_count << endl;
  cout << " p count = " << p_count << endl;
  cout << " R = " << n_count/p_count << endl;
  cout << " Sqrt(R) = " << sqrt(n_count/p_count) << endl;
  cout << "------" << endl;
  // *****

  //fitting W for neutron for SBS30%
  // double par_n[11];
  // TF1* total_n = new TF1("total_n",W_n_fit_total,0.07,1.75,11);
  // total_n->SetLineColor(2);
  // total_n->SetParameters(-2.031,54,-234,322,590,0.973,0.0967,-200,2220,-2413,670);
  // h_W_n_cut->Fit(total_n,"R+");
  // *****
 
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
