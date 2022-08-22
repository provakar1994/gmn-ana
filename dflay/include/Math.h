#ifndef UTIL_MATH_H
#define UTIL_MATH_H

#include <cstdlib>
#include <iostream>
#include <vector>
// #include <math.h>
#include <cmath>
#include <algorithm>

#include "Constants.h"

namespace math_df {
   // templated functions 
   //______________________________________________________________________________
   template < typename T >
   	T GetWeightedMean(std::vector<T> x,std::vector<T> weight,T &mean,T &meanErr){
   	   T num=0,den=0;
   
   	   int N = x.size();
   	   for(int i=0;i<N;i++){
   	   	den  += weight[i];
   	   	num  += weight[i]*x[i];
   	   }
   
   	   if( (num!=0)&&(den!=0) ){
   	   	mean    = num/den;
   	   	meanErr = sqrt(1.0/den);
   	   }else{
   	   	mean    = 0;
   	   	meanErr = 0;
   	   }
   
   	   return 0;
   	}
   //______________________________________________________________________________
   template < typename T >
   	T GetMean(std::vector<T> x){
   	   int N = x.size();
   	   T sum=0;
   	   for(int i=0;i<N;i++) sum += x[i];
   	   T mean = sum/( (T)N );
   	   return mean;
   	}
   
   template < typename T >
   	T GetMean(int N,T *x){
   	   T sum=0;
   	   for(int i=0;i<N;i++) sum += x[i];
   	   T mean = sum/( (T)N );
   	   return mean;
   	}
   //______________________________________________________________________________
   template < typename T >
   	T GetRMS(std::vector<T> x){
   	   T sum_sq = 0;
   	   const int N = x.size();
   	   for(int i=0;i<N;i++){
   	      sum_sq += pow(x[i],2.);
   	   }
   	   T arg = sum_sq/( (T)N );
   	   T rms = sqrt(arg);
   	   return rms;
   	}
   //______________________________________________________________________________
   template < typename T >
   	T GetRMS(int N,T *x){
   	   T sum_sq = 0;
   	   for(int i=0;i<N;i++){
   	   	sum_sq += pow(x[i],2.);
   	   }
   	   T arg = sum_sq/( (T)N );
   	   T rms = sqrt(arg);
   	   return rms;
   	}
      //______________________________________________________________________________
   template < typename T >
      T GetVariance(std::vector<T> x,int besselCor=Constants::kBesselsCorrectionDisabled){
         int N = x.size();
         // determine the denominator 
         int den = N;
         if(besselCor==Constants::kBesselsCorrectionEnabled){
                 den = N - 1;
         }
         T mean = GetMean<T>(x);
         T sum=0;
         for(int i=0;i<N;i++){
                 sum += pow(x[i]-mean,2);
         }
         T var   = sum/( (T)den );
         return var;
      }
   //______________________________________________________________________________
   template < typename T >
      T GetVariance(int N,T *x,int besselCor=Constants::kBesselsCorrectionDisabled){
         // determine the denominator 
         int den = N;
         if(besselCor==Constants::kBesselsCorrectionEnabled){
              den = N - 1;
         }
         T mean = GetMean<T>(N,x);
         T sum=0;
         for(int i=0;i<N;i++){
              sum += pow(x[i]-mean,2);
         }
         T var   = sum/( (T)den );
         return var;
      }
   //______________________________________________________________________________
   template < typename T >
      T GetStandardDeviation(std::vector<T> x,int besselCor=Constants::kBesselsCorrectionDisabled){
	 T var   = GetVariance<T>(x,besselCor);
	 T stdev = sqrt(var);
	 return stdev;
      }
   //______________________________________________________________________________
   template < typename T >
      T GetStandardErrorOfTheMean(std::vector<T> x, int besselCor=Constants::kBesselsCorrectionDisabled){
	 const int N = x.size();
	 // determine the denominator 
	 int den = N;
	 if(besselCor==Constants::kBesselsCorrectionEnabled){
	    den = N - 1;
	 }
	 T sd   = GetStandardDeviation<T>(x,besselCor);
	 T sdom = sd/sqrt( (T)den );
	 return sdom;
      }
   //______________________________________________________________________________
   template < typename T >
      T GetStandardErrorOfTheMean(int N,T *x,int besselCor=Constants::kBesselsCorrectionDisabled){
	 // determine the denominator 
	 int den = N;
	 if(besselCor==Constants::kBesselsCorrectionEnabled){
	    den = N - 1;
	 }
	 T sd   = GetStandardDeviation<T>(N,x,besselCor);
	 T sdom = sd/sqrt( (T)den );
	 return sdom;
      }
   //______________________________________________________________________________
   template < typename T >
      T GetCovariance(std::vector<T> x,std::vector<T> y,int besselCor=Constants::kBesselsCorrectionDisabled){
	 T mean_x = GetMean<T>(x);
	 T mean_y = GetMean<T>(y);
	 T sum=0,diff_x=0,diff_y=0;
	 const int N = x.size();
	 // determine the denominator 
	 int den = N;
	 if(besselCor==Constants::kBesselsCorrectionEnabled){
	    den = N - 1;
	 }
	 for (int i=0;i<N;i++) {
	    diff_x = x[i]-mean_x;
	    diff_y = y[i]-mean_y;
	    sum   += diff_x*diff_y;
	 }
	 T cov = sum/( (T)den );
	 return cov;
      }
   //______________________________________________________________________________
   template < typename T >
      T GetCovariance(int N,T *x,T *y,int besselCor=Constants::kBesselsCorrectionDisabled){
	 T mean_x = GetMean<T>(N,x);
	 T mean_y = GetMean<T>(N,y);
	 T sum=0,diff_x=0,diff_y=0;
	 // determine the denominator 
	 int den = N;
	 if(besselCor==Constants::kBesselsCorrectionEnabled){
	    den = N - 1;
	 }
	 for (int i=0;i<N;i++) {
	    diff_x = x[i]-mean_x;
	    diff_y = y[i]-mean_y;
	    sum   += diff_x*diff_y;
	 }
	 T cov = sum/( (T)den );
	 return cov;
      }
   //______________________________________________________________________________
   template < typename T >
      bool IsInfOrNaN(T v){
	 bool val = std::isinf(v) || std::isnan(v);
	 return val;
      }
   //______________________________________________________________________________
   template < typename T >
      int CheckVector(std::vector<T> &x){
	 // replace bad values  
	 std::replace_if(x.begin(),x.end(),IsInfOrNaN<T>,-1);
	 return 0;
      }
   // end of templated functions 
   int LeastSquaresFitting(std::vector<double> x,std::vector<double> y,double &intercept,double &slope,double &r);
   int LeastSquaresFitting(int N,double *x,double *y,double &intercept,double &slope,double &r);
   double LinearInterpolation(double x,double x0,double y0,double x1,double y1);
//   double BilinearInterpolation(double x0,double y0,
//	 std::vector<double> x,std::vector<double> y,std::vector<double> F,
//	 bool isDebug=false,double thr=1E-3);
   double AllanVariance(std::vector<double> data, const int M);
   double AllanDeviation(std::vector<double> data, const int M);
   double SimpsonIntegral(double (*f)(const double),double A,double B,double epsilon=1E-6,int Depth=10);
   double AdaptiveSimpsonAux(double (*f)(const double),double A,double B,double epsilon,
	 double S,double fa,double fb,double fc,int bottom);

}

#endif
