#ifndef UTIL_GRAPH_H
#define UTIL_GRAPH_H

#include <cstdlib>
#include <vector>

#include "TString.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TGraphPolar.h"
#include "TLegend.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"

namespace graph_df {
   // ROOT only has the first four, so we add all of them here. 
   enum lineStyle{
      kSolidLine           = 1,
      kDashLine            = 2,
      kDotLine             = 3,
      kShortDashDotLine    = 4,
      kLongDashDotLine     = 5,
      kDashTripleDotLine   = 6,
      kLongDashLine        = 7,
      kDashDoubleDotLine   = 8,
      KLongLongDashLine    = 9,
      kLongLongDashDotLine = 10
   };

   TGraph *GetTGraph(std::vector<double>,std::vector<double>);
   TGraphErrors *GetTGraphErrors(std::vector<double> x,std::vector<double> y,std::vector<double> ey);
   TGraphErrors *GetTGraphErrors(std::vector<double> x,std::vector<double> ex,std::vector<double> y,std::vector<double> ey);
   TGraphAsymmErrors *GetTGraphAsymmErrors(std::vector<double> x,std::vector<double> y,std::vector<double> eyh);
   TGraphAsymmErrors *GetTGraphAsymmErrors(std::vector<double> x,std::vector<double> y,std::vector<double> eyl,std::vector<double> eyh);
   TGraphAsymmErrors *GetBand(std::vector<double> x,std::vector<double> min,std::vector<double> max,
	 int color=kBlue,int style=3001,double alpha=0.35);

   // TGraph *GetTGraph(CSVManager *data,std::string xAxis,std::string yAxis);
   // TGraph *GetTGraphErrors(CSVManager *data,std::string xAxis,std::string yAxis,std::string yAxisErr);
   // TGraph *GetTGraph(CSVManager *data,int xIndex,int yIndex);
   // TGraph *GetTGraphErrors(CSVManager *data,int xIndex,int yIndex,int eyIndex);

   TGraph *GetTGraphDifference(const int NPTS,TGraph *g1,TGraph *g2);

   void SetParameters(TGraph *g,int marker,int color,double mSize=1.0,int width=1);
   void SetParameters(TGraphErrors *g,int marker,int color,double mSize=1.0,int width=1);
   void SetParameters(TGraphAsymmErrors *g,int marker,int color,double mSize=1.0,int width=1);
   void SetLabels(TGraph *g,TString Title,TString xAxisTitle,TString yAxisTitle);
   void SetLabels(TGraphErrors *g,TString Title,TString xAxisTitle,TString yAxisTitle);
   void SetLabels(TGraphAsymmErrors *g,TString Title,TString xAxisTitle,TString yAxisTitle);
   void SetLabels(TMultiGraph *g,TString Title,TString xAxisTitle,TString yAxisTitle);
   void SetLabelSizes(TGraph *g,double xSize,double ySize,double offset=0.5);
   void SetLabelSizes(TGraphErrors *g,double xSize,double ySize,double offset=0.5);
   void SetLabelSizes(TGraphAsymmErrors *g,double xSize,double ySize,double offset=0.5);
   void SetLabelSizes(TMultiGraph *g,double xSize,double ySize,double offset=0.5);
   void UseTimeDisplay(TGraph *g);
   void UseTimeDisplay(TGraphErrors *g);
   void UseTimeDisplay(TGraphAsymmErrors *g);
   void UseTimeDisplay(TMultiGraph *g);

} // ::graph_df 

#endif 
