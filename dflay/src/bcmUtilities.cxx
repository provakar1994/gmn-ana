#include "../include/bcmUtilities.h"
#include "Math.cxx"
// #include "CSVManager.cxx"
// #include "BCMManager.cxx"
namespace bcm_util {
  //______________________________________________________________________________
  TH1D *GetTH1D(std::vector<scalerData_t> data, std::string var, int NBin, double min, double max){
    // get a plot with data from a vector of scalerData 
    TH1D *h = new TH1D(var.c_str(), var.c_str(), NBin, min, max);
    const int N = data.size();
    for(int i=0;i<N;i++){
      h->Fill( data[i].getValue(var) ); 
    }
    return h;       
  }
  //______________________________________________________________________________
  TH2D *GetTH2D(std::vector<scalerData_t> data, std::string xVar, std::string yVar,
		int NXBin, double xMin, double xMax, int NYBin, double yMin, double yMax){
    // get a plot with data from a vector of scalerData 
    TH2D *h = new TH2D(Form("%s_vs_%s", yVar.c_str(), xVar.c_str()), 
		       Form("%s_vs_%s", yVar.c_str(), xVar.c_str()), 
		       NXBin, xMin, xMax, NYBin, yMin, yMax);
    const int N = data.size();
    for(int i=0;i<N;i++){
      h->Fill( data[i].getValue(xVar), data[i].getValue(yVar) ); 
    }
    return h;       
  }
  //______________________________________________________________________________
  TGraph *GetTGraph(std::string xAxis,std::string yAxis,std::vector<epicsData_t> data){
    // get a plot with data from a vector of scalerData 
    double ix=0, iy=0;
    std::vector<double> x,y; 
    const int N = data.size();
    for(int i=0;i<N;i++){
      ix = data[i].getValue(xAxis);
      iy = data[i].getValue(yAxis);
      x.push_back(ix); 
      y.push_back(iy); 
    }
    TGraph *g = graph_df::GetTGraph(x,y); 
    return g;       
  }
  //______________________________________________________________________________
  TGraph *GetTGraph(std::string xAxis,std::string yAxis,std::vector<scalerData_t> data){
    // get a plot with data from a vector of scalerData 
    double ix=0, iy=0;
    std::vector<double> x,y; 
    const int N = data.size();
    for(int i=0;i<N;i++){
      ix = data[i].getValue(xAxis);
      iy = data[i].getValue(yAxis);
      x.push_back(ix); 
      y.push_back(iy); 
    }
    TGraph *g = graph_df::GetTGraph(x,y); 
    return g;       
  }
  //______________________________________________________________________________
  TGraph *GetTGraph_timeStep(std::string xAxis, std::vector<scalerData_t> data){
    // get a plot with data from a vector of scalerData 
    double ix=0, iy=0, timeStep = 0.;
    std::vector<double> x,y; 
    const int N = data.size();
    for(int i=0;i<N;i++){
      if(i==0){
	timeStep = data[1].time103kHz - data[0].time103kHz;
      }else{
	timeStep = data[i].time103kHz - data[i-1].time103kHz;
      }
      ix = data[i].getValue(xAxis);
      iy = timeStep;
      //iy = data[i].time103kHz;
      x.push_back(ix); 
      y.push_back(iy); 
    }
    TGraph *g = graph_df::GetTGraph(x,y); 
    return g;       
  }
  //______________________________________________________________________________
  TGraph *GetTGraph_charge_pd(std::string var, std::string xAxis, std::vector<scalerData_t> data){
    // get a plot with data from a vector of scalerData 
    double ix=0, iy=0, timeStep = 0., current = 0., chargeSum = 0., gain = 0.;

   // extract charge from .cnt variable: charge = count / gain [Hz/uA]
    // most up-to-date gain coefficients for GMn can be found at:
    // https://sbs.jlab.org/DocDB/0001/000164/002/dflay_bcm-ana-update_02-21-22.pdf
    if (var.compare("dnew.cnt") == 0)   // dnew source is the best according to Mark J.
      gain = 3317.99; // +/- 31.69 [Hz/uA]
    else if (var.compare("u1.cnt") == 0)
      gain = 990.61; // +/- 12.00 [Hz/uA]
    else if (var.compare("unew.cnt") == 0)
      gain = 2956.77; // +/- 28.26 [Hz/uA]
    else if (var.compare("d1.cnt") == 0)
      gain = 913.73; // +/- 11.06 [Hz/uA]
    else if (var.compare("d3.cnt") == 0)
      gain = 3008.89; // +/- 36.45 [Hz/uA]
    else if (var.compare("d10.cnt") == 0)
      gain = 8606.94; // +/- 82.59 [Hz/uA]

    std::vector<double> x,y; 
    double MICROAMPS = 1E-6;
    const int N = data.size();
    for(int i=0;i<N;i++){

      chargeSum = (data[i].getValue("dnew.cnt") / 3318.)*MICROAMPS; // in C

      ix = data[i].getValue(xAxis);
      iy = chargeSum;
      x.push_back(ix); 
      y.push_back(iy); 
    }
    TGraph *g = graph_df::GetTGraph(x,y); 
    return g;       
  }
  // //______________________________________________________________________________
  // TGraph *GetTGraph_charge(std::string xAxis, std::string currSource, std::vector<scalerData_t> data){
  //   // get a plot with data from a vector of scalerData 
  //   double ix=0, iy=0, timeStep = 0., current = 0., chargeSum = 0.;
  //   std::vector<double> x,y; 
  //   double MICROAMPS = 1E-6;
  //   const int N = data.size();
  //   for(int i=0;i<N;i++){
  //     if(i==0){
  // 	timeStep = data[1].time103kHz - data[0].time103kHz;
  //     }else{
  // 	timeStep = data[i].time103kHz - data[i-1].time103kHz;
  //     }

  //     current = data[i].getValue(currSource)*MICROAMPS; // in A
  //     chargeSum += timeStep*current;

  //     ix = data[i].getValue(xAxis);
  //     iy = chargeSum;
  //     x.push_back(ix); 
  //     y.push_back(iy); 
  //   }
  //   TGraph *g = graph_df::GetTGraph(x,y); 
  //   return g;       
  // }
  //______________________________________________________________________________
  TGraph *GetTGraph_singleRun(int run,std::string xAxis,std::string yAxis,std::vector<scalerData_t> data){
    // get a plot with data from a single run  
    double ix=0,iy=0;
    std::vector<double> x,y; 
    const int N = data.size();
    for(int i=0;i<N;i++){
      if(run==data[i].runNumber){
	ix = data[i].getValue(xAxis);
	iy = data[i].getValue(yAxis);
	x.push_back(ix); 
	y.push_back(iy); 
      }
    }
    TGraph *g = graph_df::GetTGraph(x,y); 
    return g;       
  }
  //______________________________________________________________________________
  TGraphErrors *GetTGraphErrors_byRunByUnserCurrent(std::string var,std::vector<scalerData_t> data){
    // get a plot of BCM scaler rate vs Unser current (on a run-by-run basis) 
    // note: this is a sub-optimal way to do a BCM calibration; there is no toggling of on/off 
    // of the beam current here. we utilize the unser calibration against the 
    // precision current source
    // you must provide the calibration numbers to the BCMManager before calling this function  
    std::vector<double> run,mean,stdev;
    // first get the BCM scaler rate vs run number 
    GetStats_byRun(var,data,run,mean,stdev); 
    // now get the unser current 
    run.clear(); // clear this since we're grabbing them again 
    std::vector<double> mean_ui,stdev_ui; 
    GetStats_byRun("unser.current",data,run,mean_ui,stdev_ui); 
    // now make the TGraphErrors plot 
    TGraphErrors *g = graph_df::GetTGraphErrors(mean_ui,stdev_ui,mean,stdev); 
    return g; 
  }
  //______________________________________________________________________________
  double GetDAQLiveTime (std::vector<scalerData_t> runData) {
    // calculate the average DAQ live time during the run
    int N = runData.size();
    double liveTime = 0.;

    std::vector<double> frac;
    for (int i = 0; i < N; i++) {
      frac.push_back( runData[i].liveTime );
    }

    liveTime = math_df::GetMean<double>(frac); 
 
    // few sanity checks
    if (liveTime < 0.5)
      std::cout << " *!* For run " << runData[0].runNumber 
		<< " *!* Live Time seems very low! " << std::endl;
    else if (liveTime > 1. && liveTime <= 1.1) 
      liveTime = 1.;
    else if (liveTime > 1.1)
      std::cout << " **!** For run " << runData[0].runNumber 
		<< " **!** Live time unphysical! Needs attention!" << std::endl;
 
    return liveTime;
  }
  //______________________________________________________________________________
  int GetCharge_pd (std::string var, std::vector<scalerData_t> runData, charge_t &out) {
    // get the charge associated with the run
    int N = runData.size();
    double MICROAMPS = 1E-6; 
    double timeStep=0, chargeAtMaxEv=0, current=0, gain=0.;
 
    // extract charge from .cnt variable: charge = count / gain [Hz/uA]
    // most up-to-date gain coefficients for GMn can be found at:
    // https://sbs.jlab.org/DocDB/0001/000164/002/dflay_bcm-ana-update_02-21-22.pdf
    if (var.compare("dnew.cnt") == 0)   // dnew source is the best according to Mark J.
      gain = 3317.99; // +/- 31.69 [Hz/uA]
    else if (var.compare("u1.cnt") == 0)
      gain = 990.61; // +/- 12.00 [Hz/uA]
    else if (var.compare("unew.cnt") == 0)
      gain = 2956.77; // +/- 28.26 [Hz/uA]
    else if (var.compare("d1.cnt") == 0)
      gain = 913.73; // +/- 11.06 [Hz/uA]
    else if (var.compare("d3.cnt") == 0)
      gain = 3008.89; // +/- 36.45 [Hz/uA]
    else if (var.compare("d10.cnt") == 0)
      gain = 8606.94; // +/- 82.59 [Hz/uA]

    int maxChargeEv = 0;
    double a = 0., maxCharge = -1.;
    for(int i=0;i<N;i++){
      a = (runData[i].getValue("dnew.cnt") / gain) * MICROAMPS; // charge in C
      if (maxCharge < a) { 
	maxCharge = a;
	maxChargeEv = i;
      }
      if (i == N-1) chargeAtMaxEv = a;
    }
    
    // few sanity checks
    if (maxCharge != chargeAtMaxEv) {
      std::cout << " **!** For run " << runData[0].runNumber
		<< " **!** Scaler event corresponding to highest charge [" << maxChargeEv
		<< "] is not the last event [" << N << "], for which charge = "
		<< chargeAtMaxEv << std::endl;
    }
    if (maxCharge < 0.)
      std::cout << " **!** For run " << runData[0].runNumber
		<< " **!** Charge value unphysical! Needs Attention! " << std::endl;

    // calculate (statistical) uncertainty [Not necessary at the moment] 
    double chargeErr = 0.;   
    // save results
    out.runNumber = runData[0].runNumber; 
    out.totalTime = runData[N-1].time103kHz; // s
    out.value     = maxCharge;
    out.error     = chargeErr;
    return 0; 
  }
  //______________________________________________________________________________
  // int GetCharge(std::string var, std::vector<scalerData_t> runData, charge_t &out){
  //   // get the charge associated with the run

  //   int N = runData.size();
  //   double MICROAMPS = 1E-6; 
  //   double timeStep=0,chargeSum=0,current=0;
  //   double deltaTime = runData[N-1].time - runData[0].time;
  //   std::vector<double> I;
  //   // accumulate charge over the whole run for each variable  
  //   for(int i=0;i<N;i++){
  //     if(i==0){
  // 	timeStep = runData[1].time - runData[0].time;
  //     }else{
  // 	timeStep = runData[i].time - runData[i-1].time;
  //     }
  //     current = runData[i].getValue(var)*MICROAMPS; // in A
  //     // if(i<10){
  //     //    std::cout << Form("event %d, timeStep = %.3lf sec",i,timeStep) << std::endl;
  //     // }
  //     chargeSum += timeStep*current;  
  //     I.push_back(current);

  //     // std::cout << " timestep = " << timeStep << std::endl;
  //     //if (i == 2386) {
  // 	// std::cout << " evCont " << runData[i].event << " triggerEv " << runData[i].triggerEvent 
  // 	//  	  << " time " << runData[i].time << " timestep " << timeStep
  // 	// 	  << " time_num " << runData[i].time_num << " time_den " << runData[i].time_den
  // 	// 	  << " current " << current << " charge " << chargeSum << std::endl;
  //   }
  //   // calculate (statistical) uncertainty 
  //   double mean_i    = math_df::GetMean<double>(I); 
  //   double stdev_i   = math_df::GetStandardDeviation<double>(I); 
  //   double chargeErr = GetChargeError(chargeSum,mean_i,stdev_i,deltaTime,0);   
  //   // save results
  //   out.runNumber = runData[0].runNumber; 
  //   out.totalTime = deltaTime; 
  //   out.value     = chargeSum;
  //   out.error     = chargeErr;
  //   return 0; 
  // }
  //______________________________________________________________________________
  double GetChargeError(double Q,double I, double dI,double t,double dt){
    double T1=0,T2=0;
    if(I!=0) T1 = TMath::Power(dI/I,2.);
    if(t!=0) T2 = TMath::Power(dt/t,2.);
    double dQ = Q*TMath::Sqrt( T1 + T2 );
    return dQ;
  }
  //______________________________________________________________________________
  int GetStats_byRun(std::string var,std::vector<scalerData_t> data,
		     std::vector<double> &RUN,std::vector<double> &MEAN,std::vector<double> &STDEV){
    // get stats of var vs run number  
    int run_prev = data[0].runNumber;
    std::vector<double> v;
    double mean=0,stdev=0,theRun=0,theValue=0;
    const int NEV = data.size();
    for(int i=0;i<NEV;i++){
      theRun   = data[i].runNumber;
      theValue = data[i].getValue(var);  
      if(run_prev==theRun){
	v.push_back(theValue);
      }else{
	// new run! compute stats 
	mean  = math_df::GetMean<double>(v);
	stdev = math_df::GetStandardDeviation<double>(v);
	// save results
	RUN.push_back(run_prev);
	MEAN.push_back(mean);
	STDEV.push_back(stdev);
	// set up for next event 
	v.clear();
	// the current value counts for the "next" run that didn't match the previous 
	v.push_back(theValue);
      }
      run_prev = theRun;
    }
    // compute stats on the last run  
    mean  = math_df::GetMean<double>(v);
    stdev = math_df::GetStandardDeviation<double>(v);
    // save results
    RUN.push_back(run_prev);
    MEAN.push_back(mean);
    STDEV.push_back(stdev);
    // const int NR = RUN.size();
    // std::cout << "Found " << NR << " runs" << std::endl;
    // for(int i=0;i<NR;i++) std::cout << Form("%d",(int)RUN[i]) << std::endl;  
    return 0;
  }
  //______________________________________________________________________________
  int GetStats_byRun_epics(std::string var,std::vector<epicsData_t> data,
			   std::vector<double> &RUN,std::vector<double> &MEAN,std::vector<double> &STDEV){
    // get stats of var vs run number  
    int run_prev = data[0].runNumber;
    std::vector<double> v;
    double mean=0,stdev=0,theRun=0,theValue=0;
    const int NEV = data.size();
    for(int i=0;i<NEV;i++){
      theRun   = data[i].runNumber;
      theValue = data[i].getValue(var);  
      if(run_prev==theRun){
	v.push_back(theValue);
      }else{
	// new run! compute stats 
	mean  = math_df::GetMean<double>(v);
	stdev = math_df::GetStandardDeviation<double>(v);
	// save results
	RUN.push_back(run_prev);
	MEAN.push_back(mean);
	STDEV.push_back(stdev);
	// set up for next event 
	v.clear();
	// the current value counts for the "next" run that didn't match the previous 
	v.push_back(theValue);
      }
      run_prev = theRun;
    }
    // compute stats on the last run  
    mean  = math_df::GetMean<double>(v);
    stdev = math_df::GetStandardDeviation<double>(v);
    // save results
    RUN.push_back(run_prev);
    MEAN.push_back(mean);
    STDEV.push_back(stdev);
    // const int NR = RUN.size();
    // std::cout << "Found " << NR << " runs" << std::endl;
    // for(int i=0;i<NR;i++) std::cout << Form("%d",(int)RUN[i]) << std::endl;  
    return 0;
  }
  //______________________________________________________________________________
  int GetData(std::string var,std::vector<producedVariable_t> data,std::vector<double> &x,std::vector<double> &dx){
    // pull out the mean and stdev from the produced variable set
    const int N = data.size();
    for(int i=0;i<N;i++){
      if(var.compare(data[i].dev)==0){
	x.push_back(data[i].mean); 
	dx.push_back(data[i].stdev); 
      }
    }
    int NN = x.size();
    return NN; 
  }
  //______________________________________________________________________________
  int Print(producedVariable_t data){
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "dev:        " << data.dev                 << std::endl;
    std::cout << "beam_state: " << data.beam_state          << std::endl;
    std::cout << "group:      " << data.group               << std::endl;
    std::cout << "time:       " << Form("%.3lf",data.time)  << std::endl;
    std::cout << "mean:       " << Form("%.3lf",data.mean)  << std::endl;
    std::cout << "stdev:      " << Form("%.3lf",data.stdev) << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    return 0;
  }
  //______________________________________________________________________________
  int WriteToFile_cc(const char *outpath,std::vector<calibCoeff_t> data){
    // write the calibration coefficients to file
    const int NROW = data.size();
    const int NCOL = 9;
    std::vector<std::string> dev;
    std::vector<int> runMin,runMax;
    std::vector<double> pedestal,pedestalErr; 
    std::vector<double> offset,offsetErr; 
    std::vector<double> gain,gainErr; 
    for(int i=0;i<NROW;i++){
      dev.push_back( data[i].dev ); 
      runMin.push_back(data[i].runMin); 
      runMax.push_back(data[i].runMax); 
      pedestal.push_back( data[i].pedestal ); 
      pedestalErr.push_back( data[i].pedestalErr ); 
      offset.push_back( data[i].offset ); 
      offsetErr.push_back( data[i].offsetErr ); 
      gain.push_back( data[i].slope ); 
      gainErr.push_back( data[i].slopeErr ); 
    } 

    std::string header = "dev,runMin,runMax,pedestal,pedestalErr,offset,offsetErr,gain,gainErr";
 
    CSVManager *csv = new CSVManager(); 
    csv->InitTable(NROW,NCOL); 
    csv->SetColumn_str(0,dev); 
    csv->SetColumn<int>(1,runMin);  
    csv->SetColumn<int>(2,runMax);  
    csv->SetColumn<double>(3,pedestal);  
    csv->SetColumn<double>(4,pedestalErr);  
    csv->SetColumn<double>(5,offset);  
    csv->SetColumn<double>(6,offsetErr);  
    csv->SetColumn<double>(7,gain);  
    csv->SetColumn<double>(8,gainErr); 
    csv->SetHeader(header); 
    csv->WriteFile(outpath); 

    delete csv; 
    return 0;
  }
  //______________________________________________________________________________
  int WriteToFile(const char *outpath,std::vector<producedVariable_t> data){
    // write producedVariable data to file 
    std::vector<std::string> DEV,BEAM_STATE;
    std::vector<int> GRP;
    std::vector<double> MU,SIG,TIME;

    // write results to file
    const int N = data.size();
    for(int i=0;i<N;i++){
      MU.push_back(data[i].mean);
      SIG.push_back(data[i].stdev);
      TIME.push_back(data[i].time); 
      DEV.push_back(data[i].dev);
      BEAM_STATE.push_back(data[i].beam_state);
      GRP.push_back(data[i].group);
    }

    std::string header = "dev,beam_state,group,time,mean,stdev";

    const int NROW = DEV.size();
    const int NCOL = 6;
    CSVManager *csv = new CSVManager();
    csv->InitTable(NROW,NCOL);
    csv->SetColumn_str(0,DEV);
    csv->SetColumn_str(1,BEAM_STATE);
    csv->SetColumn<int>(2,GRP);
    csv->SetColumn<double>(3,TIME);
    csv->SetColumn<double>(4,MU);
    csv->SetColumn<double>(5,SIG);
    csv->SetHeader(header);
    csv->WriteFile(outpath);

    delete csv;
    return 0;
  }
  //______________________________________________________________________________
  int LoadCalibrationCoefficients(const char *inpath,std::vector<calibCoeff_t> &data){
    // load calibration coefficients from BCM or Unser calibration analysis  
    // note: this is a component of a calibCoeff data type,
    // so we fill a vector of calibCoeff structs
    CSVManager *csv = new CSVManager(); 
    int rc = csv->ReadFile(inpath,true); 
    if(rc!=0){
      delete csv;
      return 1;
    }
     
    // grab the data we want 
    std::vector<std::string> dev;
    std::vector<double> ped,pedErr,offset,offsetErr,gain,gainErr; 
    csv->GetColumn_byName_str("dev",dev); 
    csv->GetColumn_byName<double>("pedestal"   ,ped);  
    csv->GetColumn_byName<double>("pedestalErr",pedErr);  
    csv->GetColumn_byName<double>("offset"   ,offset);  
    csv->GetColumn_byName<double>("offsetErr",offsetErr);  
    csv->GetColumn_byName<double>("gain"     ,gain); 
    csv->GetColumn_byName<double>("gainErr"  ,gainErr);

    // fill a vector of type calibCoeff
    calibCoeff_t v;
    const int N = dev.size();
    for(int i=0;i<N;i++){
      v.dev         = dev[i];
      v.pedestal    = ped[i]; 
      v.pedestalErr = pedErr[i]; 
      v.offset      = offset[i];
      v.offsetErr   = offsetErr[i];
      v.slope       = gain[i];
      v.slopeErr    = gainErr[i];
      data.push_back(v);  
    }

    delete csv; 
    return 0; 
  }
  //______________________________________________________________________________
  int LoadFittedOffsetGainData(const char *inpath,std::vector<calibCoeff_t> &data){
    // load fitted offset and slope (gain) data from BCM or Unser calibration analysis  
    // note: this is a component of a calibCoeff data type,
    // so we fill a vector of calibCoeff structs
    CSVManager *csv = new CSVManager(); 
    int rc = csv->ReadFile(inpath,true); 
    if(rc!=0){
      delete csv;
      return 1;
    }
     
    // grab the data we want 
    std::vector<std::string> dev;
    std::vector<double> offset,offsetErr,gain,gainErr; 
    csv->GetColumn_byName_str("dev",dev); 
    csv->GetColumn_byName<double>("offset"   ,offset);  
    csv->GetColumn_byName<double>("offsetErr",offsetErr);  
    csv->GetColumn_byName<double>("gain"     ,gain); 
    csv->GetColumn_byName<double>("gainErr"  ,gainErr);

    // fill a vector of type calibCoeff
    calibCoeff_t v;
    const int N = dev.size();
    for(int i=0;i<N;i++){
      v.dev       = dev[i];
      v.offset    = offset[i];
      v.offsetErr = offsetErr[i];
      v.slope     = gain[i];
      v.slopeErr  = gainErr[i];
      data.push_back(v);  
    }

    delete csv; 
    return 0; 
  }
  //______________________________________________________________________________
  int LoadPedestalData(const char *inpath,std::vector<calibCoeff_t> &data){
    // load pedestal data
    // note: this is a component of a calibCoeff data type,
    // so we fill a vector of calibCoeff structs
    CSVManager *csv = new CSVManager(); 
    int rc = csv->ReadFile(inpath,true); 
    if(rc!=0){
      delete csv;
      return 1;
    }
     
    // grab the data we want 
    std::vector<std::string> dev;
    std::vector<int> runMin,runMax; 
    std::vector<double> ped,pedErr; 
    csv->GetColumn_byName_str("dev",dev); 
    csv->GetColumn_byName<int>("runMin",runMin);  
    csv->GetColumn_byName<int>("runMax",runMax);  
    csv->GetColumn_byName<double>("pedestal"   ,ped); 
    csv->GetColumn_byName<double>("pedestalErr",pedErr);

    // fill a vector of type calibCoeff
    calibCoeff_t v;
    const int N = ped.size();
    for(int i=0;i<N;i++){
      v.dev         = dev[i];
      v.runMin      = runMin[i];  
      v.runMax      = runMax[i];  
      v.pedestal    = ped[i];
      v.pedestalErr = pedErr[i];
      data.push_back(v);  
    }

    delete csv; 
    return 0; 
  }
  //______________________________________________________________________________
  int LoadProducedVariables(const char *inpath,std::vector<producedVariable_t> &data){
    // load variables that have been cut on => produced variable  
    CSVManager *csv = new CSVManager();
    int rc = csv->ReadFile(inpath,true);
    if(rc!=0){
      delete csv;
      return 1;
    }
   
    std::vector<std::string> dev,beam_state;
    csv->GetColumn_byIndex_str(0,dev); 
    csv->GetColumn_byIndex_str(1,beam_state); 

    std::vector<int> group;
    csv->GetColumn_byIndex<int>(2,group);
 
    std::vector<double> time,mean,stdev;  
    csv->GetColumn_byIndex<double>(3,time);
    csv->GetColumn_byIndex<double>(4,mean);
    csv->GetColumn_byIndex<double>(5,stdev);
     
    const int N = dev.size(); 
    producedVariable_t var; 
    for(int i=0;i<N;i++){
      var.dev        = dev[i];
      var.beam_state = beam_state[i];  
      var.group      = group[i];  
      var.time       = time[i];  
      var.mean       = mean[i];  
      var.stdev      = stdev[i]; 
      data.push_back(var);  
    }
 
    delete csv; 
    return 0;  
  }
  //______________________________________________________________________________
  int LoadConfigPaths(const char *inpath,std::vector<std::string> &label,std::vector<std::string> &path){
    // load configuration paths  
    CSVManager *csv = new CSVManager();
    int rc = csv->ReadFile(inpath);
    if(rc!=0){
      delete csv;
      return 1;
    }
     
    csv->GetColumn_byIndex_str(0,label);
    csv->GetColumn_byIndex_str(1,path );

    delete csv; 
    return 0;  
  }
  //______________________________________________________________________________
  int LoadRuns(const char *inpath,std::vector<codaRun_t> &runList){
    // load run list from a file 
    CSVManager *csv = new CSVManager();
    int rc = csv->ReadFile(inpath,true); // header is included 
    if(rc!=0){
      delete csv;
      return 1;
    }
 
    std::vector<std::string> data,STREAM,SEG_BEG,SEG_END,INFO;
    csv->GetColumn_byIndex_str(0,data);
    csv->GetColumn_byIndex_str(1,STREAM);
    csv->GetColumn_byIndex_str(2,SEG_BEG);
    csv->GetColumn_byIndex_str(3,SEG_END);
    csv->GetColumn_byIndex_str(4,INFO); 
  
    codaRun_t crun; 
    int aRun=0,aStr=0,aBeg=0,aEnd=0;
    std::string info;  
    const int N = data.size();
    for(int i=0;i<N;i++){
      aRun = std::atoi( data[i].c_str() );
      aStr = std::atoi( STREAM[i].c_str() );
      aBeg = std::atoi( SEG_BEG[i].c_str() );
      aEnd = std::atoi( SEG_END[i].c_str() );
      info = INFO[i]; 
      crun.runNumber    = aRun;
      crun.stream       = aStr; 
      // crun.segmentBegin = aBeg; 
      // crun.segmentEnd   = aEnd;
      crun.info         = info; 
      runList.push_back(crun);
    }

    delete csv; 
    return 0;  
  }
  //______________________________________________________________________________
  int LoadCuts(const char *inpath,std::vector<cut_t> &data){
    // load cuts from a file into a vector of type cut_t  
    CSVManager *csv = new CSVManager();
    int rc = csv->ReadFile(inpath,true);
    if(rc!=0){
      delete csv;
      return 1;
    }

    std::vector<std::string> arm,dev,beam_state,cut_var;
    csv->GetColumn_byName_str("arm"       ,arm);
    csv->GetColumn_byName_str("dev"       ,dev);
    csv->GetColumn_byName_str("cut_var"   ,cut_var);
    csv->GetColumn_byName_str("beam_state",beam_state);

    std::vector<double> low,high;
    csv->GetColumn_byName<double>("cut_low" ,low);
    csv->GetColumn_byName<double>("cut_high",high);

    std::vector<int> grp;
    csv->GetColumn_byName<int>("group",grp);

    cut_t aCut;
    const int N = dev.size();
    for(int i=0;i<N;i++){
      aCut.arm        = arm[i];
      aCut.dev        = dev[i];
      aCut.cut_var    = cut_var[i];
      aCut.beam_state = beam_state[i];
      aCut.low        = low[i];
      aCut.high       = high[i];
      aCut.group      = grp[i];
      data.push_back(aCut);
    }

    // delete CSV manager 
    delete csv;
  
    return 0;
  }
  //______________________________________________________________________________
  int ApplyCuts(double cutLo,double cutHi,std::vector<double> x,std::vector<double> y,std::vector<double> &out){
    // cuts are applied to the x variable. out has data that passed the cut 
    const int N = x.size();
    for(int i=0;i<N;i++){
      if(x[i]>cutLo&&x[i]<cutHi) out.push_back(y[i]);
    }
    return 0;
  }
  //______________________________________________________________________________
  int GetStatsWithCuts(std::vector<double> x,std::vector<double> y,
		       double cutLo,double cutHi,double &mean,double &stdev){
    // cuts are applied to the x variable. if true, compute stats on y 
    std::vector<double> Y;
    const int N = x.size();
    for(int i=0;i<N;i++){
      if(x[i]>cutLo&&x[i]<cutHi) Y.push_back(y[i]);
    }
    mean  = math_df::GetMean<double>(Y);
    stdev = math_df::GetStandardDeviation<double>(Y);
    return 0;
  }
  //______________________________________________________________________________
  int ConvertToCurrent(calibCoeff_t cc,std::vector<producedVariable_t> unser_ps,
		       std::vector<producedVariable_t> &unser_cur){
    double T1=0,T2=0,current=0,currentErr=0;
    producedVariable_t data;
    const int N = unser_ps.size();
    for(int i=0;i<N;i++){
      // compute current 
      current    = (unser_ps[i].mean - cc.offset)/cc.slope;
      // compute error
      T1         = cc.slope*cc.slope*(unser_ps[i].stdev*unser_ps[i].stdev + cc.offsetErr*cc.offsetErr);
      T2         = cc.slopeErr*cc.slopeErr*TMath::Power(unser_ps[i].mean - cc.offset,2.);
      currentErr = TMath::Power(1./cc.slope,2.)*TMath::Sqrt(T1 + T2);
      // save result
      data.mean  = current;
      data.stdev = currentErr;
      unser_cur.push_back(data);
    }
    return 0;
  }
  //______________________________________________________________________________
  int CalculateStatsForBeamState(std::string beamState,std::vector<producedVariable_t> data,
				 std::vector<producedVariable_t> &out,double &MEAN,double &STDEV,
				 std::string LOG_PATH){
    // compute the stats averaged over all cycles per group for a given beam state 
    // - input
    //   - vector of beam-on and beam-off data (alternates between on and off data; also has a given group/beam current) 
    //   - LOG_PATH: print screen messages named LOG_PATH (if not NONE) 
    // - output
    //   - vector of selected data, averaged over all cycles per group (index is group)  
    //   - mean and stdev of all the selected data averaged over all cycles and group 

    char msg[200];
    producedVariable_t aPt;
    std::vector<double> time,v,dv,w;
    double mean=0,err=0,stdev=0,argErr=0; 
   
    sprintf(msg,"Calculating stats for bcm: %s, beam state: %s",data[0].dev.c_str(),beamState.c_str()); 
    if(LOG_PATH.compare("NONE")!=0) util_df::LogMessage(LOG_PATH.c_str(),msg,'a'); 
 
    int M=0; 
    int grp_prev = data[0].group; 
    std::string dev_prev = data[0].dev;  
    const int N = data.size(); 
    for(int i=0;i<N;i++){
      // check the group 
      if(data[i].group==grp_prev){
	if(data[i].beam_state.compare(beamState)==0){
	  time.push_back(data[i].time);
	  v.push_back(data[i].mean); 
	  dv.push_back(data[i].stdev); 
	}
      }else{
	// new group
	sprintf(msg,"==== GROUP %d ====",grp_prev);
	if(LOG_PATH.compare("NONE")!=0) util_df::LogMessage(LOG_PATH.c_str(),msg,'a');
	M = v.size();
	for(int j=0;j<M;j++){
	  sprintf(msg,"time = %.3lf, value = %.3E, err = %.3E",time[j],v[j],dv[j]); 
	  if(LOG_PATH.compare("NONE")!=0) util_df::LogMessage(LOG_PATH.c_str(),msg,'a');
	}
	// compute stats
	for(int j=0;j<M;j++){
	  argErr = dv[j]*dv[j];
	  if(dv[j]!=0){
	    w.push_back(1./argErr);
	  }else{
	    w.push_back(1);
	  }
	}
	// compute the weighted mean
	math_df::GetWeightedMean<double>(v,w,mean,err);
	stdev = math_df::GetStandardDeviation<double>(v);
	// save result 
	aPt.dev   = dev_prev; 
	aPt.beam_state = beamState; 
	aPt.time  = time[0]; // keep first time stamp  
	aPt.mean  = mean; 
	aPt.stdev = stdev;
	out.push_back(aPt); 
	sprintf(msg,"--> time = %.3lf, mean = %.3E, stdev = %.3E",aPt.time,aPt.mean,aPt.stdev);
	if(LOG_PATH.compare("NONE")!=0) util_df::LogMessage(LOG_PATH.c_str(),msg,'a');
	// set up for next group 
	w.clear();
	time.clear();
	v.clear();
	dv.clear();  
	// save the current one since we'll need it 
	if(data[i].beam_state.compare(beamState)==0){
	  time.push_back(data[i].time);
	  v.push_back(data[i].mean); 
	  dv.push_back(data[i].stdev);
	} 
      }
      grp_prev = data[i].group;
      dev_prev = data[i].dev;  
    } 

    // get the last one
    sprintf(msg,"==== GROUP %d ====",grp_prev);
    if(LOG_PATH.compare("NONE")!=0) util_df::LogMessage(LOG_PATH.c_str(),msg,'a');
    M = v.size();
    for(int j=0;j<M;j++){
      sprintf(msg,"time = %.3lf, value = %.3E, err = %.3E",time[j],v[j],dv[j]); 
      if(LOG_PATH.compare("NONE")!=0) util_df::LogMessage(LOG_PATH.c_str(),msg,'a');
    }
    // compute stats
    for(int j=0;j<M;j++){
      argErr = dv[j]*dv[j];
      if(dv[j]!=0){
	w.push_back(1./argErr);
      }else{
	w.push_back(1);
      }
    }
    // compute the weighted mean
    math_df::GetWeightedMean<double>(v,w,mean,err);
    stdev = math_df::GetStandardDeviation<double>(v);
    // save result
    aPt.dev   = dev_prev; 
    aPt.beam_state = beamState;  
    aPt.time  = time[0]; // keep first time stamp  
    aPt.mean  = mean; 
    aPt.stdev = stdev;
    out.push_back(aPt); 
    sprintf(msg,"--> time = %.3lf, mean = %.3E, stdev = %.3E",aPt.time,aPt.mean,aPt.stdev);
    if(LOG_PATH.compare("NONE")!=0) util_df::LogMessage(LOG_PATH.c_str(),msg,'a');

    // now do the weighted mean on ALL groups 
    v.clear();
    w.clear();
    int NN = out.size();
    for(int i=0;i<NN;i++){
      v.push_back(out[i].mean);
      argErr = out[i].stdev*out[i].stdev;  
      if(argErr!=0){
	w.push_back(1./argErr); 
      }else{
	w.push_back(1);
      }
    }

    math_df::GetWeightedMean<double>(v,w,MEAN,err); 
    STDEV = math_df::GetStandardDeviation<double>(v); 

    return 0;
  }
  //______________________________________________________________________________
  int CalculatePedestalSubtraction(std::vector<producedVariable_t> data,std::vector<producedVariable_t> &out,
				   std::string LOG_PATH,std::string PLOT_PATH){
    // compute pedestal-subtracted rates 
    // utilizes ABA method to account for zero-point drift
    // - input:  
    //   - vector of beam-on and beam-off data (alternates between on and off data; also has a given group/beam current)
    //   - LOG_PATH:  print screen messages to file named LOG_PATH (if not NONE) 
    //   - PLOT_PATH: print plots to file named PLOT_PATH (if not NONE) 
    // - output: vector of (beam-on) - (beam-off) results (entry for each group/beam current)  
    ABA *myABA = new ABA();
    myABA->UseTimeWeight();
    // myABA->SetVerbosity(1);   

    char msg[200];
    double mean=0,err=0,stdev=0,argErr=0;
    std::vector<double> w,aba,abaErr;
    std::vector<double> timeOn,on,onErr;
    std::vector<double> timeOff,off,offErr;

    producedVariable_t aPt;

    // int NG = data.size()/2; // for each group, we have beam on and off => divide N by 2 to get number of groups 

    // find number of groups 
    const int N = data.size();
    std::vector<int> grp; 
    for(int i=0;i<N;i++) grp.push_back( data[i].group ); 
    std::sort(grp.begin(),grp.end()); 
    auto last = std::unique(grp.begin(),grp.end());
    grp.erase(last,grp.end());
    int NG = grp.size(); 

    int NROW=2;
    int NCOL=NG/2; 
    TCanvas *cp = new TCanvas("cp","Pedestal Subtraction",1200,800);
    cp->Divide(NCOL,NROW); 

    sprintf(msg,"[bcmUtilities::CalculatePedestalSubtraction]: Processing %d groups",NG);    
    util_df::LogMessage(LOG_PATH.c_str(),msg,'a'); 

    // create graph object for each group 
    TGraphErrors **gOff = new TGraphErrors*[NG]; 
    TGraphErrors **gOn  = new TGraphErrors*[NG];
    TMultiGraph **mg    = new TMultiGraph*[NG];  

    int M=0,k=0;
    int grp_prev = data[0].group; // effective beam current 
    std::string dev_prev = data[0].dev;  
    for(int i=0;i<N;i++){
      // check the group 
      if(data[i].group==grp_prev){
	// group match! store data based on beam state 
	if(data[i].beam_state.compare("on")==0){
	  timeOn.push_back( data[i].time );
	  on.push_back( data[i].mean );
	  onErr.push_back( data[i].stdev );
	}else if(data[i].beam_state.compare("off")==0){
	  timeOff.push_back( data[i].time );
	  off.push_back( data[i].mean );
	  offErr.push_back( data[i].stdev );
	}
      }else{
	// new group! compute ABA stats
	sprintf(msg,"==== GROUP %d ====",grp_prev);
	if(LOG_PATH.compare("NONE")!=0) util_df::LogMessage(LOG_PATH.c_str(),msg,'a');
	// check 
	M = timeOff.size();
	for(int j=0;j<M;j++){
	  sprintf(msg,"off: %.3lf, %.3lf ± %.3lf; on: %.3lf, %.3lf ± %.3lf",
		  timeOff[j],off[j],offErr[j],timeOn[j],on[j],onErr[j]);
	  if(LOG_PATH.compare("NONE")!=0) util_df::LogMessage(LOG_PATH.c_str(),msg,'a');
	}
	if(M==1){
	  sprintf(msg,"**** ONLY ONE CYCLE! GROUP %d",grp_prev);
	  if(LOG_PATH.compare("NONE")!=0) util_df::LogMessage(LOG_PATH.c_str(),msg,'a');
	  // account for one cycle
	  mean = on[0] - off[0];
	  err  = TMath::Sqrt( onErr[0]*onErr[0] + offErr[0]*offErr[0] );
	  stdev = err;
	}else{
	  // compute ABA stats
	  myABA->GetDifference(timeOff,off,offErr,timeOn,on,onErr,aba,abaErr);
	  M = aba.size();
	  for(int j=0;j<M;j++){
	    argErr = abaErr[j]*abaErr[j];
	    if(abaErr[j]!=0){
	      w.push_back(1./argErr);
	    }else{
	      w.push_back(1);
	    }
	  }
	  // compute the weighted mean
	  math_df::GetWeightedMean<double>(aba,w,mean,err);
	  stdev = math_df::GetStandardDeviation<double>(aba);
	  // we compute (A-B), but we actually want B-A
	  mean *= -1;
	}
	// store results
	aPt.mean  = mean;
	aPt.stdev = stdev;
	out.push_back(aPt);
	// make a plot
	gOff[k] = graph_df::GetTGraphErrors(timeOff,off,offErr);  
	gOn[k]  = graph_df::GetTGraphErrors(timeOn , on,onErr); 
	graph_df::SetParameters(gOff[k],20,kBlack);  
	graph_df::SetParameters(gOn[k] ,21,kRed);
	mg[k] = new TMultiGraph(); 
	mg[k]->Add(gOff[k],"lp");  
	mg[k]->Add(gOn[k] ,"lp"); 
	cp->cd(k+1); 
	mg[k]->Draw("a"); 
	graph_df::SetLabels(mg[k],Form("Group %d",grp_prev),"Time",Form("%s",dev_prev.c_str())); 
	mg[k]->Draw("a"); 
	cp->Update();  
	// increment plotter index 
	k++;
	// set up for next
	w.clear();
	aba.clear();
	abaErr.clear();
	timeOn.clear();
	on.clear();
	onErr.clear();
	timeOff.clear();
	off.clear();
	offErr.clear();
	// store this one since it's needed for the next set!  
	if(data[i].beam_state.compare("on")==0){
	  timeOn.push_back( data[i].time );
	  on.push_back( data[i].mean );
	  onErr.push_back( data[i].stdev );
	}else if(data[i].beam_state.compare("off")==0){
	  timeOff.push_back( data[i].time );
	  off.push_back( data[i].mean );
	  offErr.push_back( data[i].stdev );
	}
      }
      grp_prev = data[i].group;
      dev_prev = data[i].dev;  
    }
    // get the last index computed 
    sprintf(msg,"==== GROUP %d ====",grp_prev);
    if(LOG_PATH.compare("NONE")!=0) util_df::LogMessage(LOG_PATH.c_str(),msg,'a');

    M = timeOff.size();
    for(int j=0;j<M;j++){
      sprintf(msg,"off: %.3lf, %.3lf ± %.3lf; on: %.3lf, %.3lf ± %.3lf",
	      timeOff[j],off[j],offErr[j],timeOn[j],on[j],onErr[j]);
      if(LOG_PATH.compare("NONE")!=0) util_df::LogMessage(LOG_PATH.c_str(),msg,'a');
    }

    // make a plot 
    gOff[k] = graph_df::GetTGraphErrors(timeOff,off,offErr);  
    gOn[k]  = graph_df::GetTGraphErrors(timeOn , on,onErr); 
    graph_df::SetParameters(gOff[k],20,kBlack);  
    graph_df::SetParameters(gOn[k] ,21,kRed); 
    mg[k] = new TMultiGraph(); 
    mg[k]->Add(gOff[k],"lp");  
    mg[k]->Add(gOn[k] ,"lp"); 
    cp->cd(k+1); 
    mg[k]->Draw("a"); 
    graph_df::SetLabels(mg[k],Form("Group %d",grp_prev),"Time",Form("%s",dev_prev.c_str())); 
    mg[k]->Draw("a"); 
    cp->Update();  

    // save canvas if necessary 
    if(PLOT_PATH.compare("NONE")!=0){
      cp->Print(PLOT_PATH.c_str());
    }

    // compute ABA stats
    w.clear();
    aba.clear();
    abaErr.clear();
    if(M==1){
      sprintf(msg,"**** ONLY ONE CYCLE! GROUP %d",grp_prev);
      if(LOG_PATH.compare("NONE")!=0) util_df::LogMessage(LOG_PATH.c_str(),msg,'a');
      // account for one cycle
      mean = on[0] - off[0];
      err  = TMath::Sqrt( onErr[0]*onErr[0] + offErr[0]*offErr[0] );
      stdev = err;
    }else{
      // compute ABA stats
      myABA->GetDifference(timeOff,off,offErr,timeOn,on,onErr,aba,abaErr);
      M = aba.size();
      for(int j=0;j<M;j++){
	argErr = abaErr[j]*abaErr[j];
	if(abaErr[j]!=0){
	  w.push_back(1./argErr);
	}else{
	  w.push_back(1);
	}
      }
      // compute the weighted mean
      math_df::GetWeightedMean<double>(aba,w,mean,err);
      stdev = math_df::GetStandardDeviation<double>(aba);
      // we compute (A-B), but we actually want B-A
      mean *= -1;
    }
    // store results
    aPt.mean  = mean;
    aPt.stdev = stdev;
    out.push_back(aPt);

    delete myABA;
    delete cp; 

    return 0;
  }
  //______________________________________________________________________________
  int SubtractBaseline(std::vector<producedVariable_t> on,std::vector<producedVariable_t> off, 
		       std::vector<producedVariable_t> &diff,bool isDebug){
    // compute the quantity (beam-on) - (beam-off) for all variables
    // only compute this for producedVariables with matching group numbers and names  
    int cntr=0;
    double arg=0,argErr=0;
    producedVariable_t x;
    const int NOFF = off.size();
    const int NON  = on.size();
    for(int i=0;i<NON;i++){
      for(int j=0;j<NOFF;j++){
	if(on[i].group==off[j].group && on[i].dev.compare(off[j].dev)==0){
	  cntr++;
	  arg          = on[i].mean - off[j].mean;
	  argErr       = TMath::Sqrt( on[i].stdev*on[i].stdev + off[j].stdev*off[j].stdev );
	  x.mean       = arg;
	  x.stdev      = argErr;
	  x.dev        = on[i].dev;
	  x.beam_state = "DIFF";
	  x.group      = on[i].group;
	  diff.push_back(x); 
	  // print 
	  if(isDebug){
	    std::cout << Form("%s, group %d: on = %.3lf ± %.3lf, off = %.3lf ± %.3lf, on-off = %.3lf ± %.3lf",
			      on[i].dev.c_str(),on[i].group,on[i].mean,on[i].stdev,off[i].mean,off[i].stdev,arg,argErr) << std::endl;
	  }
	}
      }
    }
    std::cout << "Found " << cntr << " BCM matches" << std::endl;
    return cntr;
  }
  //______________________________________________________________________________
  TGraphErrors *GetTGraphErrors_byRun(std::string var,std::vector<scalerData_t> data){
    // make TGraphErrors of data as a function of run number 
    std::vector<double> run,mean,stdev;
    GetStats_byRun(var,data,run,mean,stdev); 
    TGraphErrors *g = graph_df::GetTGraphErrors(run,mean,stdev);
    return g;
  }
}
