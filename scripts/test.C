#include <iostream>
#include <vector>

#include "../src/ExpConstants.cpp"
#include "../src/KinematicVar.cpp"
#include "../src/SetROOTVar.cpp"

void test()
{
  //std::cout << kine::W(2,2,kine::Q2(2,2,0.2)) << std::endl;

  TChain *T = new TChain("T");
  T->Add("/lustre19/expphy/volatile/halla/sbs/pdbforce/gmn-replays/rootfiles/sbs4-sbs30p/e1209019_fullreplay_11616_stream0_seg5*");
  
  T->SetBranchStatus("*",0);
  //bbsh clus var
  double ESH, xSH, ySH;
  std::vector<std::string> shclvar = {"e","x","y"};
  std::vector<void*> shclvar_mem = {&ESH,&xSH,&ySH};
  setrootvar::setbranch(T,"bb.sh",shclvar,shclvar_mem);

  //bbps clus var
  double EPS, xPS, yPS;
  std::vector<std::string> psclvar = {"e","x","y"};
  std::vector<void*> psclvar_mem = {&EPS,&xPS,&yPS};
  setrootvar::setbranch(T,"bb.ps",psclvar,psclvar_mem);

  //hcal clus var
  double EHCAL, xHCAL, yHCAL;
  std::vector<std::string> hcalclvar = {"e","x","y"};
  std::vector<void*> hcalclvar_mem = {&EHCAL,&xHCAL,&yHCAL};
  setrootvar::setbranch(T,"sbs.hcal",hcalclvar,hcalclvar_mem);

  //track var
  double ntrack, p[1000],px[1000],py[1000],pz[1000];
  double vx[1000],vy[1000],vz[1000];
  double xfp[1000],yfp[1000],thfp[1000],phfp[1000];
  double xtgt[1000],ytgt[1000],thtgt[1000],phtgt[1000];
  std::vector<std::string> trvar = {"n","p","px","py","pz",
				    "vx","vy","vz",
				    "r_x","r_y","r_th","r_ph",
				    "tg_x","tg_y","tg_th","tg_ph"};
  std::vector<void*> trvar_mem = {&ntrack,&p,&px,&py,&pz,
				  &vx,&vy,&vz,&xfp,&yfp,&thfp,&phfp,
				  &xtgt,&ytgt,&thtgt,&phtgt};
  setrootvar::setbranch(T,"bb.tr",trvar,trvar_mem);

  //tdctrig variable
  int tdcElemN;
  double tdcTrig[1000], tdcElem[1000];
  std::vector<std::string> tdcvar = {"tdcelemID","tdcelemID","tdc"};
  std::vector<void*> tdcvar_mem = {&tdcElem,&tdcElemN,&tdcTrig};
  setrootvar::setbranch(T,"bb.tdctrig",tdcvar,tdcvar_mem,1);

  //T->GetEntry(5000);
  //std::cout << rblkSH << std::endl;
  //std::cout << px[0] << std::endl;
  //std::cout << p[0] << std::endl;
  
  TH1F *h_she = new TH1F("h_she","h_she",200,0,4);
  TH1D *h_coin_time = new TH1D("h_coin_time",";HCAL trigTime - BBCAL trigTime (ns);",150,450,600);

  Long64_t Nevents = T->GetEntries();  
  for(Long64_t nevent = 0; nevent<Nevents; nevent++){
    if( nevent%1000 == 0){
      cout << nevent << "/" << Nevents << endl;
    }
    T->GetEntry(nevent);

    //coincidence time analysis
    double bbcal_time=0., hcal_time=0.;
    for(int ihit=0; ihit<tdcElemN; ihit++){
      if(tdcElem[ihit]==5) bbcal_time=tdcTrig[ihit];
      if(tdcElem[ihit]==0) hcal_time=tdcTrig[ihit];
    }
    double coin_time = hcal_time-bbcal_time;
    h_coin_time->Fill( coin_time );

    if(ESH>0) h_she->Fill(p[0]);

  }

  TCanvas *c1 = new TCanvas("c1","",600,400);
  c1->cd();
  h_coin_time->Draw();

}
