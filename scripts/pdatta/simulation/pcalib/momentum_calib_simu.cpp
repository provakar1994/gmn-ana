//#include "GEM_cosmic_tracks.C"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TEventList.h"
#include "TCut.h"
#include <iostream>
#include <fstream>
#include "TMinuit.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TEllipse.h"
#include <vector>
#include <set>
#include <map>
#include <sstream>

double PI = TMath::Pi();

double Mp = 0.938272;
double Mn = 0.939565;

void momentum_calib_simu( const char *configfilename, const char *outputfilename="NewMomentumFit.root" ){

  gStyle->SetOptFit();
  
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
  
  //int fitorder = 2;
  
  double ebeam=1.916;
  double bbtheta = 51.0*TMath::DegToRad();
  double sbstheta = 34.5*TMath::DegToRad();
  double hcaldist = 13.5; //meters
  double hcalvoff = 0.35; //meters
  double Ltgt = 15.0; //cm
  double rho_tgt = 0.0723; //g/cc
  double rho_Al = 2.7; //g/cc

  double sbsdist = 2.25;
  double sbsfield = 1.0; //fraction of max field:
  //double sbsmaxfield = 0.53/0.3; //About 1.77 Tesla.
  double sbsmaxfield = 1.26; //Tesla
  double Dgap = 48.0*2.54/100.0; //about 1.22 m
  
  double celldiameter = 1.6*2.54; //cm, right now this is a guess

  //Eventually we will grab HALLA:p from the EPICs tree for the beam energy:

  
  sbsmaxfield = 3.1 * atan( 0.85/(11.0 - 2.25 - Dgap/2.0 ))/0.3/Dgap/0.7;
  
  // 0.85 m deflection / 3.1 nucleon momentum @70% of max field
  // p * thetabend = 0.3 * BdL
  // deflection = tan(thetabend) * (11 - (2.25 + 0.6) )
    //BdL = p*thetabend/0.3 --> Bmax = p*thetabend/0.3/dL/0.7

    
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
  
  int usehcalcut = 0; //default to false:

  double dx0=0.0, dy0=0.0; //zero offsets for deltax and deltay cuts:

  double dxsigma = 0.06, dysigma = 0.1; //sigmas for deltax and deltay cuts:

  double GEMpitch = 10.0*TMath::DegToRad();

  double dpelmin_fit = -0.03;
  double dpelmax_fit = 0.03;
  double Wmin_fit = 0.86;
  double Wmax_fit = 1.02; 

  //To calculate a proper dx, dy for HCAL, we need
  double hcalheight = 0.0; //0.365; //m (we are guessing that this is the height of the center of HCAL above beam height:

  //The following are the positions of the "first" row and column from HCAL database (top right block as viewed from upstream)
  double xoff_hcal = 0.92835;
  double yoff_hcal = 0.47305;
  
  double blockspace_hcal = 0.15254; //don't need to make this configurable

  int nrows_hcal=24;
  int ncols_hcal=12;

  //By default we fix the zero-order coefficient for pth and the first-order pth coefficient for xfp:
  
  int fix_pth0_flag = 1;
  int fix_pthx_flag = 1;

  double pth0 = 0.275;
  double pthx = 0.102;
  
  int order = 2; 
  //int order_ptheta = 1; //default to first-order expansion for momentum, and 2nd-order for vertices, angles. 
  
  double A_pth0 = 0.275;
  double B_pth0 = 0.61;
  double C_pth0 = -0.074;

  double BBdist = 1.85;
  
  TString fname_oldcoeffs = "oldcoeffs.dat"; //mainly for merging angle/vertex reconstruction coefficients with momentum coefficients:

  double thtgtmin_fit = -0.15;
  double thtgtmax_fit = 0.15;

  //double Wmin_treefill = 0;
  //double Wmax_treefill = 

  int hcal_coordflag = 0; //0 = "new" HCAL coordinate system that matches transport.
                          //1 = "old" HCAL coordinate system that doesn't match transport
  
  
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endconfig") ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");

      int ntokens = tokens->GetEntries();

      if( ntokens >= 2 ){
	TString skey = ( (TObjString*) (*tokens)[0] )->GetString();

	if( skey == "oldcoeffs" ){ //this should be the name of a text file
	  fname_oldcoeffs = ( (TObjString*) (*tokens)[1] )->GetString();
	}
	  
	if ( skey == "order" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  order = stemp.Atoi();
	}

	// if( skey == "order_pth" ){
	//   TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	//   order_ptheta = stemp.Atoi();
	// }
	
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

	if( skey == "A_pth0" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  A_pth0 = stemp.Atof();
	}

	if( skey == "B_pth0" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  B_pth0 = stemp.Atof();
	}

	if( skey == "C_pth0" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  C_pth0 = stemp.Atof();
	}
	
	if( skey == "sbstheta" ){ //SBS/HCAL central angle (deg) (assumed to be on beam right)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sbstheta = stemp.Atof() * TMath::DegToRad();
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

	if( skey == "dx0" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dx0 = stemp.Atof();
	}

	if( skey == "dy0" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dy0 = stemp.Atof();
	}

	if( skey == "dxsigma" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dxsigma = stemp.Atof();
	}

	if( skey == "dysigma" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dysigma = stemp.Atof();
	}

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
	  Wmin_fit = stemp.Atof();
	}

	if( skey == "Wmax" ){ //Wmax to include in fit
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Wmax_fit = stemp.Atof();
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

	if( skey == "hcalvoff" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  hcalheight = stemp.Atof();
	}

	if( skey == "hcalxoff" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  xoff_hcal = stemp.Atof();
	}

	if( skey == "hcalyoff" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  yoff_hcal = stemp.Atof();
	}

	if( skey == "hcalcoordsys" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  hcal_coordflag = stemp.Atoi();
	}
	
	if( skey == "fit_thtgt_min" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  thtgtmin_fit = stemp.Atof();
	}

	if( skey == "fit_thtgt_max" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  thtgtmax_fit = stemp.Atof();
	}

	if( skey == "sbsdist" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sbsdist = stemp.Atof();
	}

	if( skey == "sbsfield" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sbsfield = stemp.Atof();
	}

	if( skey == "sbsmaxfield" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sbsmaxfield = stemp.Atof();
	}

	if( skey == "hcalvoff" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  hcalvoff = stemp.Atof();
	}
      }
      
      tokens->Delete();
    }
  }
  
  //Note that both of these calculations neglect the Aluminum end windows and cell walls:
  
  MeanEloss = Ltgt/2.0 * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //hopefully the user has provided these quantities with correct units
  
  MeanEloss_outgoing = celldiameter/2.0/sin(bbtheta) * rho_tgt * dEdx_tgt +
    cellthick_LH2/sin(bbtheta) * rho_Al * dEdx_Al;
  if( useAlshield != 0 ) MeanEloss_outgoing += Alshieldthick * rho_Al * dEdx_Al;

  cout << "mean eloss in aluminum shield = " << Alshieldthick * rho_Al * dEdx_Al << endl;
  
  cout << "use Aluminum shield = " << useAlshield << endl;
						     
						     
  
  cout << "Mean E loss (beam electron) = " << MeanEloss << endl;
  cout << "BigBite theta = " << bbtheta * TMath::RadToDeg() << endl;
  cout << "Mean E loss (outgoing electron) = " << MeanEloss_outgoing << endl;
  cout << "Ebeam corrected (GeV) = " << ebeam - MeanEloss << endl;
 
		   					   
  
  // In order to get 

  TEventList *elist = new TEventList("elist","Event list for BigBite momentum calibration");

  C->Draw(">>elist",globalcut);

  int MAXNTRACKS=1000;

  //Declare variables to hold branch addresses:
  
  double ntrack;
  double px[MAXNTRACKS], py[MAXNTRACKS], pz[MAXNTRACKS], p[MAXNTRACKS];
  double vx[MAXNTRACKS], vy[MAXNTRACKS], vz[MAXNTRACKS];

  //Use "rotated" versions of the focal plane variables:
  double xfp[MAXNTRACKS], yfp[MAXNTRACKS], thfp[MAXNTRACKS], phfp[MAXNTRACKS];
  //Also need target variables:
  double xtgt[MAXNTRACKS], ytgt[MAXNTRACKS], thtgt[MAXNTRACKS], phtgt[MAXNTRACKS];
  Double_t xTr[MAXNTRACKS], yTr[MAXNTRACKS];
  Double_t thTr[MAXNTRACKS], phTr[MAXNTRACKS];
  
  double xHCAL, yHCAL, EHCAL;
  double EPS, ESH, xSH, ySH, xPS;

  //Ignore hodo variables for meow:
  // int maxHODOclusters=100;

  // double nHODOclusters;
  
  // double xHODO[maxHODOclusters],yHODO[maxHODOclusters];
  // double tmeanHODO[maxHODOclusters], tdiffHODO[maxHODOclusters];
  
  C->SetBranchStatus("*",0);
  //Minimal set of HCAL variables:
  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);

  //BigBite track variables:
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.p",1);
  C->SetBranchStatus("bb.tr.vx",1);
  C->SetBranchStatus("bb.tr.vy",1);
  C->SetBranchStatus("bb.tr.vz",1);
  
  C->SetBranchStatus("bb.tr.r_x",1);
  C->SetBranchStatus("bb.tr.r_y",1);
  C->SetBranchStatus("bb.tr.r_th",1);
  C->SetBranchStatus("bb.tr.r_ph",1);
  
  C->SetBranchStatus("bb.tr.tg_x",1);
  C->SetBranchStatus("bb.tr.tg_y",1);
  C->SetBranchStatus("bb.tr.tg_th",1);
  C->SetBranchStatus("bb.tr.tg_ph",1);
  
  //Shower and preshower variables:
  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.ps.x",1);
  C->SetBranchStatus("bb.ps.y",1);
  C->SetBranchStatus("bb.sh.e",1);
  C->SetBranchStatus("bb.sh.x",1);
  C->SetBranchStatus("bb.sh.y",1);

  C->SetBranchAddress("bb.tr.n",&ntrack);
  C->SetBranchAddress("bb.tr.p",p);
  C->SetBranchAddress("bb.tr.px",px);
  C->SetBranchAddress("bb.tr.py",py);
  C->SetBranchAddress("bb.tr.pz",pz);
  C->SetBranchAddress("bb.tr.vx",vx);
  C->SetBranchAddress("bb.tr.vy",vy);
  C->SetBranchAddress("bb.tr.vz",vz);

  //Focal-plane track variables (use "rotated" versions):
  C->SetBranchAddress("bb.tr.r_x",xfp);
  C->SetBranchAddress("bb.tr.r_y",yfp);
  C->SetBranchAddress("bb.tr.r_th",thfp);
  C->SetBranchAddress("bb.tr.r_ph",phfp);

  //Target track variables (other than momentum):
  C->SetBranchAddress("bb.tr.tg_x",xtgt);
  C->SetBranchAddress("bb.tr.tg_y",ytgt);
  C->SetBranchAddress("bb.tr.tg_th",thtgt);
  C->SetBranchAddress("bb.tr.tg_ph",phtgt);

  C->SetBranchAddress("sbs.hcal.x",&xHCAL);
  C->SetBranchAddress("sbs.hcal.y",&yHCAL);
  C->SetBranchAddress("sbs.hcal.e",&EHCAL);

  C->SetBranchAddress("bb.ps.e",&EPS);
  C->SetBranchAddress("bb.sh.e",&ESH);
  C->SetBranchAddress("bb.ps.x",&xPS);
  C->SetBranchAddress("bb.sh.x",&xSH);
  C->SetBranchAddress("bb.sh.y",&ySH);

  // C->SetBranchStatus("bb.tr.x",1);
  // C->SetBranchAddress("bb.tr.x",xTr);
  // C->SetBranchStatus("bb.tr.y",1);
  // C->SetBranchAddress("bb.tr.y",yTr);
  // C->SetBranchStatus("bb.tr.th",1);
  // C->SetBranchAddress("bb.tr.th",thTr);
  // C->SetBranchStatus("bb.tr.ph",1);
  // C->SetBranchAddress("bb.tr.ph",phTr);
  
  
  int nparams = 0;

  vector<int> xtar_expon;
  vector<int> xfp_expon;
  vector<int> xpfp_expon;
  vector<int> yfp_expon;
  vector<int> ypfp_expon;

  //ifstream foldcoeffs(fname_oldcoeffs.Data());
  
  for(int i=0; i<=order; i++){
    for(int j=0; j<=order-i; j++){
      for(int k=0; k<=order-i-j; k++){
  	for(int l=0; l<=order-i-j-k; l++){
  	  for(int m=0; m<=order-i-j-k-l; m++){
  	    nparams++;
  	    xtar_expon.push_back( i );
  	    xfp_expon.push_back( m );
  	    yfp_expon.push_back( l );
  	    xpfp_expon.push_back( k );
  	    ypfp_expon.push_back( j );
  	  }
  	}
      }
    }
  }

  //For our purposes here we only need to set up the augmented matrix for the p*thetabend matrix elements:
  TMatrixD M(nparams,nparams);
  TVectorD b_pth(nparams);

  for( int i=0; i<nparams; i++ ){
    for( int j=0; j<nparams; j++ ){
      M(i,j) = 0.0;
    }
    b_pth(i) = 0.0;
  }


  double pcentral = ebeam/(1.+ebeam/Mp*(1.-cos(bbtheta)));

  TFile *fout = new TFile(outputfilename,"RECREATE");

  //After we get this working, we'll also want correlation plots of dpel and W vs. focal plane and target variables, new and old:
  TH1D *hdpel_old = new TH1D("hdpel_old","Old MEs;p/p_{elastic}(#theta)-1;",250,-0.25,0.25);
  TH1D *hW_old = new TH1D("hW_old","Old MEs;W (GeV);",250,0,2);
  TH1D *hp_old = new TH1D("hp_old","Old MEs;p (GeV);",250,0.5*pcentral,1.5*pcentral);
  TH1D *hpthetabend_old = new TH1D("hpthetabend_old","Old MEs; p#theta_{bend};",250,0.2,0.35);

  TH2D *hxpfp_vs_xfp_old = new TH2D("hxpfp_vs_xfp_old",";xfp (m); xpfp",250,-0.75,0.75,250,-0.4,0.4);
  
  TH2D *hdpel_xfp_old = new TH2D("hdpel_xfp_old", "Old MEs; xfp (m); p/p_{elastic}(#theta)-1",250,-0.75,0.75,250,-0.125,0.125);
  TH2D *hdpel_yfp_old = new TH2D("hdpel_yfp_old", "Old MEs; yfp (m); p/p_{elastic}(#theta)-1",250,-0.25,0.25,250,-0.125,0.125);
  TH2D *hdpel_xpfp_old = new TH2D("hdpel_xpfp_old", "Old MEs; xpfp; p/p_{elastic}(#theta)-1",250,-0.4,0.4,250,-0.125,0.125);
  TH2D *hdpel_ypfp_old = new TH2D("hdpel_ypfp_old", "Old MEs; ypfp; p/p_{elastic}(#theta)-1",250,-0.15,0.15,250,-0.125,0.125);

  TH2D *hpthetabend_xfp_old = new TH2D("hpthetabend_xfp_old", "Old MEs; xfp (m); p#theta_{bend}",250,-0.75,0.75,250,0.2,0.35);
  TH2D *hpthetabend_yfp_old = new TH2D("hpthetabend_yfp_old", "Old MEs; yfp (m); p#theta_{bend}",250,-0.25,0.25,250,0.2,0.35);
  TH2D *hpthetabend_xpfp_old = new TH2D("hpthetabend_xpfp_old", "Old MEs; xpfp; p#theta_{bend}",250,-0.4,0.4,250,0.2,0.35);
  TH2D *hpthetabend_ypfp_old = new TH2D("hpthetabend_ypfp_old", "Old MEs; ypfp; p#theta_{bend}",250,-0.15,0.15,250,0.2,0.35);

  TH2D *hpthetabend_xptar_old = new TH2D("hpthetabend_xptar_old", "Old MEs; xptar; p#theta_{bend}",250,-0.4,0.4,250,0.0,0.5);
  TH2D *hpthetabend_yptar_old = new TH2D("hpthetabend_yptar_old", "Old MEs; yptar; p#theta_{bend}",250,-0.15,0.15,250,0.2,0.35);
  TH2D *hpthetabend_ytar_old = new TH2D("hpthetabend_ytar_old", "Old MEs; ytar; p#theta_{bend}",250,-0.15,0.15,250,0.2,0.35);

  TH2D *hdpel_xptar_old = new TH2D("hdpel_xptar_old", "Old MEs; xptar; p/p_{elastic}(#theta)-1",250,-0.4,0.4,250,-0.125,0.125);
  TH2D *hdpel_yptar_old = new TH2D("hdpel_yptar_old", "Old MEs; yptar; p/p_{elastic}(#theta)-1",250,-0.15,0.15,250,-0.125,0.125);
  TH2D *hdpel_ytar_old = new TH2D("hdpel_ytar_old", "Old MEs; ytar; p/p_{elastic}(#theta)-1",250,-0.15,0.15,250,-0.125,0.125);

  TH1D *hdx_HCAL_old = new TH1D("hdx_HCAL_old","Old MEs; xHCAL - xBB (m);", 250, -2.5,2.5);
  TH1D *hdy_HCAL_old = new TH1D("hdy_HCAL_old","Old MEs; yHCAL - yBB (m);", 250, -1.25,1.25);
  
  TH2D *hdxdy_HCAL_old = new TH2D("hdxdy_HCAL_old","Old MEs; yHCAL - yBB (m); xHCAL - xBB (m)", 250, -1.25,1.25, 250,-2.5,2.5);
  TH2D *hdxdy_HCAL_new = new TH2D("hdxdy_HCAL_new","Old MEs; yHCAL - yBB (m); yHCAL - yBB (m)", 250, -1.25,1.25, 250,-2.5,2.5);
  
  TH1D *hW_new = new TH1D("hW_new","New MEs;W (GeV);",250,0,2);
  TH1D *hdpel_new = new TH1D("hdpel_new","New MEs;p/p_{elastic}(#theta)-1;",250,-0.25,0.25);

  TH1D *hW_new_nocut = new TH1D("hW_new_nocut","New MEs;W (GeV);",250,0,2);
  TH1D *hdpel_new_nocut = new TH1D("hdpel_new_nocut","New MEs;p/p_{elastic}(#theta)-1;",250,-0.25,0.25);
  
  TTree *Tout = new TTree("Tout","Tree containing variables for momentum calibration");

  double T_ebeam, T_etheta, T_ephi, T_precon, T_pelastic, T_thetabend, T_dpel, T_W2;
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
  double T_Q2_4vect, T_W2_4vect;
  double T_thetaq, T_phiq, T_px, T_py, T_pz, T_qx, T_qy, T_qz, T_qmag;
  double T_theta_recon_n, T_phi_recon_n, T_theta_recon_p, T_phi_recon_p; //Reconstructed nucleon angles under proton and neutron hypothesis
  double T_thetapq_n, T_phipq_n, T_thetapq_p, T_phipq_p;
  double T_dx_4vect, T_dy_4vect; //Here we want to use the 4-vector momentum transfer to calculate dx/dy
  int HCALcut;
  int BBcut;
  
  
  Tout->Branch( "HCALcut", &HCALcut, "HCALcut/I");
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
  Tout->Branch( "vx", &T_vx, "vx/D");
  Tout->Branch( "vy", &T_vy, "vy/D");
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
  Tout->Branch( "pp_expect", &T_pp_expect, "pp_expect/D");
  Tout->Branch( "ptheta_expect", &T_ptheta_expect, "ptheta_expect/D");
  Tout->Branch( "pphi_expect", &T_pphi_expect, "pphi_expect/D");
  Tout->Branch( "EPS", &T_EPS, "EPS/D");
  Tout->Branch( "ESH", &T_ESH, "ESH/D");
  Tout->Branch( "Etot", &T_Etot, "Etot/D");
  Tout->Branch( "xSH", &T_xSH, "xSH/D");
  Tout->Branch( "ySH", &T_ySH, "ySH/D");
  Tout->Branch( "Q2_4vect", &T_Q2_4vect, "Q2_4vect/D");
  Tout->Branch( "W2_4vect", &T_W2_4vect, "W2_4vect/D");
  Tout->Branch( "epx", &T_px, "epx/D");
  Tout->Branch( "epy", &T_py, "epy/D");
  Tout->Branch( "epz", &T_pz, "epz/D");
  Tout->Branch( "qx", &T_qx, "qx/D" );
  Tout->Branch( "qy", &T_qy, "qy/D" );
  Tout->Branch( "qz", &T_qz, "qz/D" );
  Tout->Branch( "qmag", &T_qmag, "qmag/D" );
  Tout->Branch( "thetaq", &T_thetaq, "thetaq/D");
  Tout->Branch( "phiq", &T_phiq, "phiq/D" );
  Tout->Branch( "thetarecon_n", &T_theta_recon_n, "thetarecon_n/D" );
  Tout->Branch( "phirecon_n", &T_phi_recon_n, "phirecon_n/D" );
  Tout->Branch( "thetarecon_p", &T_theta_recon_p, "thetarecon_p/D" );
  Tout->Branch( "phirecon_p", &T_phi_recon_p, "phirecon_p/D" );
  Tout->Branch( "thetapq_n", &T_thetapq_n, "thetapq_n/D" );
  Tout->Branch( "thetapq_p", &T_thetapq_p, "thetapp_n/D" );
  //Tout->Branch( "phipq_n", &T_phipq_n, "phipq_n/D" );
  //Tout->Branch( "phipq_p", &T_phipq_p, "phipp_n/D" );
  Tout->Branch( "deltax_4vect", &T_dx_4vect, "deltax_4vect/D" );
  Tout->Branch( "deltay_4vect", &T_dy_4vect, "deltay_4vect/D" );
  
  
  
  
  long nevent=0;

  //First pass: accumulate sums required for the fit: 

  while( C->GetEntry(elist->GetEntry(nevent++)) ){
    if( nevent % 100000 == 0 ) cout << nevent << endl;

    if( int(ntrack) == 1 ){
      //The first thing we want to do is to calculate the "true" electron momentum incident on BigBite:
      double Ebeam_corrected = ebeam - MeanEloss;

      T_ebeam = Ebeam_corrected;
      
      double etheta = acos(pz[0]/p[0]);

      T_etheta = etheta;
      
      // Calculate the expected momentum of an elastically scattered electron at the reconstructed scattering angle and then correct it for the mean energy loss of the
      // electron on its way out of the target:
      double pelastic = Ebeam_corrected/(1.+Ebeam_corrected/Mp*(1.0-cos(etheta))); 
      
      double precon = p[0] + MeanEloss_outgoing; //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)

      T_precon = precon;
      T_pelastic = pelastic; 
      T_pincident = pelastic - MeanEloss_outgoing;
      
      double nu_recon = Ebeam_corrected - precon;
      
      double Q2recon = 2.0*Ebeam_corrected*precon*(1.0-cos(etheta));
      double W2recon = pow(Mp,2) + 2.0*Mp*nu_recon - Q2recon;

      T_Q2 = Q2recon;
      T_W2 = W2recon;

      T_xfp = xfp[0];
      T_yfp = yfp[0];
      T_thfp = thfp[0];
      T_phfp = phfp[0];

      T_thtgt = thtgt[0];
      T_phtgt = phtgt[0];
      T_ytgt = ytgt[0];
      //T_xtgt = xtgt[0];

      //for now:
      T_vx = 0.0;
      T_vy = 0.0;

      //xtgt = -vy - vz * costheta * thtgt
      T_xtgt = -vz[0]*cos(bbtheta)*thtgt[0];

      T_vz = vz[0];

      T_BBdist = BBdist;
      T_BBtheta = bbtheta;
      T_HCALdist = hcaldist;
      T_HCALtheta = sbstheta;

      
      double Wrecon = sqrt( std::max(0.0,W2recon) );
      
      double dpel = precon/pelastic-1.0;

      T_dpel = dpel;
      
      bool passed_HCAL_cut = true;

      TVector3 vertex( 0, 0, vz[0] );
      TLorentzVector Pbeam(0,0,Ebeam_corrected,Ebeam_corrected);
      TLorentzVector Kprime(precon*px[0]/p[0],precon*py[0]/p[0],precon*pz[0]/p[0],precon);
      TLorentzVector Ptarg(0,0,0,Mp);
      
      TLorentzVector q = Pbeam - Kprime;

      T_Q2_4vect = -q.M2();
      T_W2_4vect = (Ptarg + q).M2();

      TVector3 qvect = q.Vect();
      T_thetaq = qvect.Theta();
      T_phiq = qvect.Phi();
      T_px = Kprime.Px();
      T_py = Kprime.Py();
      T_pz = Kprime.Pz();

      T_qx = q.Px();
      T_qy = q.Py();
      T_qz = q.Pz();

      T_qmag = qvect.Mag();
      
      double ephi = atan2( py[0], px[0] );
      //double etheta = acos( pz[0]/p[0] );

      T_ephi = ephi;
      
      //      double nu = Ebeam_corrected - p[0];

      // Usually when we are running this code, the angle reconstruction is already well calibrated, but the momentum reconstruction is
      // unreliable; use pel(theta) as electron momentum for kinematic correlation:
      double nu = Ebeam_corrected - pelastic; 
      
      double pp_expect = sqrt(pow(nu,2)+2.*Mp*nu); //sqrt(nu^2 + Q^2) = sqrt(nu^2 + 2Mnu)
      double pphi_expect = ephi + PI;
      //double ptheta_expect = acos( (Ebeam_corrected-pz[0])/pp_expect ); //will this give better resolution than methods based on electron angle only? Not entirely clear

      T_pp_expect = pp_expect;
      T_pphi_expect = pphi_expect;
      
      double ptheta_expect = acos( (Ebeam_corrected-pelastic*cos(etheta))/pp_expect ); //will this give better resolution than methods based on electron angle only? Not entirely clear

      T_ptheta_expect = ptheta_expect;

      T_EPS = EPS;
      T_ESH = ESH;
      T_Etot = EPS + ESH;
      T_xSH = xSH;
      T_ySH = ySH;
      
      TVector3 pNhat( sin(ptheta_expect)*cos(pphi_expect), sin(ptheta_expect)*sin(pphi_expect), cos(ptheta_expect) );

      TVector3 HCAL_zaxis(-sin(sbstheta),0,cos(sbstheta));
      TVector3 HCAL_xaxis(0,1,0);
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
      
      TVector3 HCAL_origin = hcaldist * HCAL_zaxis + hcalheight * HCAL_xaxis;
      
      TVector3 TopRightBlockPos_DB(xoff_hcal,yoff_hcal,0);
      
      TVector3 TopRightBlockPos_Hall( hcalheight + (nrows_hcal/2-0.5)*blockspace_hcal,
				      (ncols_hcal/2-0.5)*blockspace_hcal, 0 );
      
      
      //Assume that HCAL origin is at the vertical and horizontal midpoint of HCAL

      
      if( hcal_coordflag == 0 ){ //new HCAL coordinate system:
	HCAL_zaxis.SetXYZ( -sin(sbstheta),0,cos(sbstheta) );
	HCAL_xaxis.SetXYZ(0, -1, 0 );
	HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
	HCAL_origin = hcaldist * HCAL_zaxis - hcalvoff * HCAL_xaxis; //We no longer need to offset the origin. 
      }

      if( hcal_coordflag != 0 ){ //if we are still using the old HCAL coordinates, need to offset xHCAL, yHCAL:
	xHCAL += TopRightBlockPos_Hall.X() - TopRightBlockPos_DB.X();
	yHCAL += TopRightBlockPos_Hall.Y() - TopRightBlockPos_DB.Y();
      }
      
      double sintersect = (HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot(HCAL_zaxis));
      //Straight-line projection to the surface of HCAL:
      TVector3 HCAL_intersect = vertex + sintersect * pNhat;

      double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      double xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );

      TVector3 qunit = qvect.Unit();

      double s4vect = (HCAL_origin - vertex).Dot( HCAL_zaxis )/ (qunit.Dot( HCAL_zaxis ) );
      TVector3 HCAL_intersect4 = vertex + s4vect * qunit;

      T_dx_4vect = xHCAL - (HCAL_intersect4 - HCAL_origin).Dot( HCAL_xaxis );
      T_dy_4vect = yHCAL - (HCAL_intersect4 - HCAL_origin).Dot( HCAL_yaxis );

      TVector3 HCALpos = HCAL_origin + xHCAL * HCAL_xaxis + yHCAL * HCAL_yaxis;

      //Calculate expected proton deflection using crude model:
      double BdL = sbsfield * sbsmaxfield * Dgap;

      //thetabend = 0.3 * BdL: 
      double proton_deflection = tan( 0.3 * BdL / qvect.Mag() ) * (hcaldist - (sbsdist + Dgap/2.0) );

      //cout << "Expected proton deflection = " << proton_deflection << " meters" << endl;
      
      TVector3 NeutronDirection = (HCALpos - vertex).Unit();
      T_theta_recon_n = acos( NeutronDirection.Z() );
      T_phi_recon_n = TMath::ATan2( NeutronDirection.Y(), NeutronDirection.X() );

      TVector3 ProtonDirection = (HCALpos + proton_deflection * HCAL_xaxis - vertex).Unit();
      T_theta_recon_p = acos( ProtonDirection.Z() );
      T_phi_recon_p = TMath::ATan2( ProtonDirection.Y(), ProtonDirection.X() );

      T_thetapq_n = acos( NeutronDirection.Dot( qvect.Unit() ) );
      T_thetapq_p = acos( ProtonDirection.Dot( qvect.Unit() ) );
      
      //Now we need to calculate the "true" trajectory bend angle for the electron from the reconstructed angles:
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
      
      if( Wrecon >= Wmin_fit && Wrecon <= Wmax_fit ){
	hdx_HCAL_old->Fill( xHCAL - xexpect_HCAL );
	hdy_HCAL_old->Fill( yHCAL - yexpect_HCAL );
	hdxdy_HCAL_old->Fill( yHCAL - yexpect_HCAL,
			      xHCAL - xexpect_HCAL );
      }

      T_xHCAL = xHCAL;
      T_yHCAL = yHCAL;
      T_EHCAL = EHCAL;
      T_deltax = xHCAL - xexpect_HCAL;
      T_deltay = yHCAL - yexpect_HCAL;
      
      if( usehcalcut != 0 ){
	passed_HCAL_cut = pow( (xHCAL-xexpect_HCAL - dx0)/dxsigma, 2 ) +
	  pow( (yHCAL-yexpect_HCAL - dy0)/dysigma, 2 ) <= pow(2.5,2); 
      }

      HCALcut = pow( (xHCAL-xexpect_HCAL - dx0)/dxsigma, 2 ) +
	pow( (yHCAL-yexpect_HCAL - dy0)/dysigma, 2 ) <= pow(2.5,2); 

      BBcut = Wrecon >= Wmin_fit && Wrecon <= Wmax_fit && dpel >= dpelmin_fit && dpel <= dpelmax_fit;
      
      if( passed_HCAL_cut ){

	hdpel_old->Fill( dpel );
	hW_old->Fill( sqrt(W2recon) );

	hdpel_xfp_old->Fill( xfp[0], dpel );
	hdpel_yfp_old->Fill( yfp[0], dpel );
	hdpel_xpfp_old->Fill( thfp[0], dpel );
	hdpel_ypfp_old->Fill( phfp[0], dpel );

	
	hdpel_xptar_old->Fill( thtgt[0], dpel );
	hdpel_yptar_old->Fill( phtgt[0], dpel );
	hdpel_ytar_old->Fill( ytgt[0], dpel );
	
	
	double pincident = pelastic - MeanEloss_outgoing;
	
	
	
	// cout << "Reconstructed bend angle = " << thetabend * TMath::RadToDeg() << endl;
	// cout << "p*thetabend = " << pincident * thetabend << endl;

	//Increment sums only if dpel and W are within limits:

	hpthetabend_xptar_old->Fill( thtgt[0], pincident*thetabend );
	
	if( dpel >= dpelmin_fit && dpel <= dpelmax_fit &&
	    Wrecon >= Wmin_fit && Wrecon <= Wmax_fit ){
	  //cout << "dpel, W = " << dpel << ", " << Wrecon << endl;

	  hp_old->Fill( precon );
	  hpthetabend_old->Fill( pincident * thetabend );

	  hpthetabend_xfp_old->Fill( xfp[0], pincident*thetabend );
	  hpthetabend_yfp_old->Fill( yfp[0], pincident*thetabend );
	  hpthetabend_xpfp_old->Fill( thfp[0], pincident*thetabend );
	  hpthetabend_ypfp_old->Fill( phfp[0], pincident*thetabend );

	  
	  hpthetabend_yptar_old->Fill( phtgt[0], pincident*thetabend );
	  hpthetabend_ytar_old->Fill( ytgt[0], pincident*thetabend );
	  
	  hxpfp_vs_xfp_old->Fill( xfp[0], thfp[0] );
	  
	}
      }

      Tout->Fill();
      
    }
  }

  //Let's also automate the determination of the momentum coefficients:

  TCanvas *c1 = new TCanvas("c1","c1",1600,1200);

  c1->Divide(2,2,.001,.001);
  
  TProfile *hptheta_thtgt = new TProfile( "hptheta_thtgt", "; #theta_{tgt} ; p_{incident}#theta_{bend} (GeV*rad)", 250,-0.25,0.25);
  
  
  Tout->Project("hptheta_thtgt", "ep_incident*thetabend:thtgt", "HCALcut&&BBcut");

  double lmargin = 0.16, bmargin = 0.12, rmargin = 0.12, tmargin = 0.06;

  

  hptheta_thtgt->SetMarkerStyle(20);
  c1->cd(1)->SetGrid();
  gPad->SetLeftMargin(lmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetRightMargin(rmargin);
  gPad->SetTopMargin(tmargin);
  
  hpthetabend_xptar_old->SetTitle("");
  hpthetabend_xptar_old->SetXTitle("#theta_{tgt} (rad)");
  hpthetabend_xptar_old->SetYTitle("p_{incident}#theta_{bend} (GeV/c*rad)");
  hpthetabend_xptar_old->Draw("colz");
  hptheta_thtgt->Draw("SAME");

  
  //  hptheta_thgt->Draw("PSAME");
  
  TF1 *fitfunc = new TF1("fitfunc", "[0]*(1.0 + [1]*x)", -0.2,0.2 );

  hptheta_thtgt->Fit(fitfunc,"S","same",thtgtmin_fit,thtgtmax_fit);

  double A_pth_fit = fitfunc->GetParameter(0);
  double B_pth_fit = fitfunc->GetParameter(1);

  TString plotexpression;
  plotexpression.Form( "%g*(1.0+%g*thtgt)/thetabend/ep_incident-1.0", A_pth_fit, B_pth_fit );
  
  Tout->Project("hdpel_new", plotexpression.Data(), "HCALcut");
  Tout->Project("hdpel_new_nocut", plotexpression.Data(), "");
  
  //W^2 = M^2 + 2M nu - Q2
  //nu = E- Eprime
  //Q2 = 2E Eprime*(1-costheta)
  //-->W^2 = M^2 + 2M * (E-Eprime) - 2E Eprime * (1-pz/p)
  plotexpression.Form("sqrt(pow(0.938272,2)+2.*.938272*(%g-(%g*(1.0+%g*thtgt)/thetabend))-2.0*%g*(%g*(1.0+%g*thtgt)/thetabend)*(1.0-cos(etheta)))",ebeam-MeanEloss, A_pth_fit, B_pth_fit, ebeam-MeanEloss, A_pth_fit, B_pth_fit);

  Tout->Project("hW_new",plotexpression.Data(), "HCALcut");
  Tout->Project("hW_new_nocut",plotexpression.Data(), "");

  c1->cd(2)->SetGrid();
  hdxdy_HCAL_old->Draw("colz");

  TEllipse E;
  E.SetFillStyle(0);
  E.SetLineColor(6);
  E.DrawEllipse(dy0,dx0,2.5*dysigma,2.5*dxsigma,0,360,0);

  
  
  c1->cd(3)->SetGrid();
  gPad->SetLeftMargin(lmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetRightMargin(rmargin);
  gPad->SetTopMargin(tmargin);

  hdpel_new_nocut->SetLineColor( kGreen+2 );
  //hdpel_new_nocut->SetFillColorAlpha( kGreen+1, 0.2 );
  
  
  hdpel_new->SetLineColor(1);
  hdpel_new->SetLineWidth(1);
  hdpel_old->SetLineColor(4);
  hdpel_old->SetLineWidth(1);
  hdpel_new->SetFillColorAlpha(1,0.7);
  hdpel_new->Draw();
  hdpel_new->GetYaxis()->SetRangeUser(0,1.1*hdpel_new_nocut->GetMaximum());
  hdpel_new_nocut->Draw("same");
  
  
  hdpel_new->Fit("gaus","","",-0.02,0.03);
  hdpel_old->Draw("SAME");

  
  
  c1->cd(4)->SetGrid();
  gPad->SetLeftMargin(lmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetRightMargin(rmargin);
  gPad->SetTopMargin(tmargin);

  hW_new_nocut->SetLineColor( kGreen+1 );
  hW_new_nocut->SetFillColorAlpha( kGreen+1, 0.2 );
  hW_new_nocut->Draw();
  
  hW_new->SetLineColor(1);
  hW_new->SetLineWidth(1);
  hW_new->SetFillColorAlpha(1,0.7);
  hW_old->SetLineColor(4);
  hW_old->SetLineWidth(1);

  TLine L;
  L.SetLineColor(2);
  L.SetLineWidth(2);
  hW_new->Draw("same");
  L.DrawLine(0.938272,0,0.938272,hW_new->GetMaximum());
  //hW_new->GetXaxis()->SetRangeUser(0.

  // hW_new_nocut->SetLineColor(kGreen+1);
  // hW_new_nocut->Draw("SAME");
  
  hW_old->Draw("SAME");

  
  
  //File containing old optics coefficients: 
  ifstream foldcoeffs(fname_oldcoeffs.Data());

  TString newcoeffs_fname = outputfilename;
  newcoeffs_fname.ReplaceAll(".root",".dat");
  
  ofstream fnewcoeffs( newcoeffs_fname.Data() );

  fnewcoeffs << "bb.preconflag = 1" << endl;
  fnewcoeffs << "# NOTE on new momentum reconstruction formalism: 1st-order momentum is calculated from" << endl
	     << "#   p*thetabend = A_pth1*(1.0 + (B_pth1+C_pth1*bb.magdist)*thtgt)" << endl
	     << "# and momentum expansion coefficients are for delta = p(first order) * (1+delta)" << endl
	     << "# In most of the acceptance we don't need any corrections beyond the first-order model" << endl
	     << endl;
  
  currentline.Form("bb.A_pth1 = %15.9g", A_pth_fit );
  fnewcoeffs << currentline << endl;
  currentline.Form("bb.B_pth1 = %15.9g", B_pth_fit );
  fnewcoeffs << currentline << endl;
  fnewcoeffs << "bb.C_pth1 = 0.0" << endl << endl;

  fnewcoeffs << "#NOTE: angle and vertex reconstruction coefficients are from " << fname_oldcoeffs << endl << endl;
  fnewcoeffs << "bb.optics_parameters = " << endl;
  
  
  //TString currentline;
  while( currentline.ReadLine( foldcoeffs ) ){
    std::string thisline( currentline.Data() );
    std::istringstream sline( thisline );

    double coeffs[4];
    int expon[5];
    for( int i=0; i<4; i++ ){
      sline >> coeffs[i];
    }
    for( int i=0; i<5; i++ ){
      sline >> expon[i];
    }

    //Set all p*thetabend coefficients to zero: 
    coeffs[3] = 0.0; 
    TString newcoeffs;
    newcoeffs.Form(" %15.8g %15.8g %15.8g %15.8g    %d %d %d %d %d ",
		   coeffs[0], coeffs[1], coeffs[2], coeffs[3],
		   expon[0], expon[1], expon[2], expon[3], expon[4] );

    fnewcoeffs << newcoeffs << endl;
  }

  // TString plotsfilename = outputfilename;
  // plotsfilename.ReplaceAll(".root",".pdf");
  
  // c1->Print(plotsfilename.Data(),"pdf");
  // plotsfilename.ReplaceAll(".pdf",".png");
  // c1->Print(plotsfilename.Data(),"png");
 
  elist->Delete();

  c1->Write();
  fout->Write();
}
