#include "../include/Math.h"
namespace math_df { 
   //______________________________________________________________________________
   double SimpsonIntegral(double (*f)(const double),
	 double A,double B,double epsilon,int Depth){
      // Adaptive Simpson's Rule
      double C   = (A + B)/2.0;
      double H   = B - A;
      double fa  = (*f)(A);
      double fb  = (*f)(B);
      double fc  = (*f)(C);
      double S   = (H/6.0)*(fa + 4.0*fc + fb);
      double ans = AdaptiveSimpsonAux(f,A,B,epsilon,S,fa,fb,fc,Depth);
      return ans;
   }
   //______________________________________________________________________________
   double AdaptiveSimpsonAux(double (*f)(const double),
	 double A,double B,double epsilon,
	 double S,double fa,double fb,double fc,int bottom){
      // Recursive auxiliary function for AdaptiveSimpson() function
      double C      = (A + B)/2.0;
      double H      = B - A;
      double D      = (A + C)/2.0;
      double E      = (C + B)/2.0;
      double fd     = (*f)(D);
      double fe     = (*f)(E);
      double Sleft  = (H/12.0)*(fa + 4.0*fd + fc);
      double Sright = (H/12.0)*(fc + 4.0*fe + fb);
      double S2     = Sleft + Sright;
      if( (bottom <= 0) || (fabs(S2 - S) <= 15.0*epsilon) ){
	 return S2 + (S2 - S)/15;
      }
      double arg = AdaptiveSimpsonAux(f,A,C,epsilon/2.0,Sleft, fa,fc,fd,bottom-1) +
	 AdaptiveSimpsonAux(f,C,B,epsilon/2.0,Sright,fc,fb,fe,bottom-1);
      return arg;
   }
   //______________________________________________________________________________
   int LeastSquaresFitting(std::vector<double> x,std::vector<double> y,double &intercept,double &slope,double &r){
      // we do this just in case we need to loop over N < x.size() 

      // linear regression to find slope b and y-intercept a of 
      // f(x) = a + bx 

      int rc=0;
      double num=0,rsq=0;

      const int N   = x.size();
      double mean_x = GetMean<double>(x);
      double mean_y = GetMean<double>(y);
      double var_x  = GetVariance<double>(x);
      double var_y  = GetVariance<double>(y);
      double cov_xy = GetCovariance<double>(x,y);

      double ss_xx = ( (double)N )*var_x;
      double ss_yy = ( (double)N )*var_y;
      double ss_xy = ( (double)N )*cov_xy;

      double den = ss_xx*ss_yy;
      if(den==0){
	 // singular matrix. can't solve the problem.
	 intercept = 0;
	 slope     = 0;
	 r         = 0;
	 rc        = 1;
      }else{
	 slope     = cov_xy/var_x;
	 intercept = mean_y - slope*mean_x;
	 num       = ss_xy*ss_xy;
	 rsq       = num/den;
	 r         = sqrt(rsq);
      }
      return rc;
   }

   //______________________________________________________________________________
   int LeastSquaresFitting(int N,double *x,double *y,double &intercept,double &slope,double &r){
      // we do this just in case we need to loop over N < x.size() 

      // linear regression to find slope b and y-intercept a of 
      // f(x) = a + bx 

      int rc=0;
      double num=0,rsq=0;

      double mean_x = GetMean<double>(N,x);
      double mean_y = GetMean<double>(N,y);
      double var_x  = GetVariance<double>(N,x);
      double var_y  = GetVariance<double>(N,y);
      double cov_xy = GetCovariance<double>(N,x,y);

      double ss_xx = ( (double)N )*var_x;
      double ss_yy = ( (double)N )*var_y;
      double ss_xy = ( (double)N )*cov_xy;

      double den = ss_xx*ss_yy;
      if(den==0){
	 // singular matrix. can't solve the problem.
	 intercept = 0;
	 slope     = 0;
	 r         = 0;
	 rc        = 1;
      }else{
	 slope     = cov_xy/var_x;
	 intercept = mean_y - slope*mean_x;
	 num       = ss_xy*ss_xy;
	 rsq       = num/den;
	 r         = sqrt(rsq);
      }
      return rc;
   }
   //______________________________________________________________________________
   double LinearInterpolation(double x,double x0,double y0,double x1,double y1){
      double b = (x-x0)/(x1-x0);
      double y = y0 + b*(y1-y0);
      return y;
   }
   // //______________________________________________________________________________
   // double BilinearInterpolation(double x0,double y0,
   //       std::vector<double> x,std::vector<double> y,std::vector<double> F,
   //       bool isDebug,double thr){
   //    // bilinear interpolation to estimate F(x0,y0)  
   //    // input 
   //    // - desired coordinate x0, y0 
   //    // - vectors x, y, F representing the grid F(x,y)  
   //    // optional input
   //    // - isDebug: print debug info (default = false)  
   //    // - thr: threshold for determining corners of F(x',y') given input x',y' (default = 1E-3) 
   //    // output 
   //    // - the value F(x0,y0)  

   //    // make copies of x and y vectors to find bounds
   //    // passing vector as argument to constructor => deep copy
   //    std::vector<double> xx(x),yy(y);

   //    // sort and strip out repeats
   //    Algorithm::SortedRemoveDuplicates<double>(xx);
   //    Algorithm::SortedRemoveDuplicates<double>(yy);

   //    // find bounding values for x0 and y0  
   //    int ixlo=0,ixhi=0,iylo=0,iyhi=0;
   //    Algorithm::BinarySearch<double>(xx,x0,ixlo,ixhi);
   //    Algorithm::BinarySearch<double>(yy,y0,iylo,iyhi);

   //    // get bounding (x,y) values
   //    double xl   = xx[ixlo];
   //    double xh   = xx[ixhi];
   //    double yl   = yy[iylo];
   //    double yh   = yy[iyhi];

   //    // weight factors for interpolation 
   //    double xwl  = (x0-xl)/(xh-xl);
   //    double xwh  = (xh-x0)/(xh-xl);
   //    double ywl  = (y0-yl)/(yh-yl);
   //    double ywh  = (yh-y0)/(yh-yl);
   //    // construct F(xlo,ylo), F(xlo,yhi), F(xhi,ylo), F(xhi,yhi)
   //    double fll=0,flh=0,fhl=0,fhh=0;
   //    const int NP = F.size();
   //    for(int i=0;i<NP;i++){
   //       if( abs(x[i]-xl)<thr && abs(y[i]-yl)<thr ) fll = F[i];
   //       if( abs(x[i]-xl)<thr && abs(y[i]-yh)<thr ) flh = F[i];
   //       if( abs(x[i]-xh)<thr && abs(y[i]-yl)<thr ) fhl = F[i];
   //       if( abs(x[i]-xh)<thr && abs(y[i]-yh)<thr ) fhh = F[i];
   //    }

   //    char msg[200];

   //    if(isDebug){
   //       std::cout << "---------------- BilinearInterpolation ----------------" << std::endl;
   //       sprintf(msg,"xl = %.3lf < x0 = %.3lf < xh = %.3lf",xl,x0,xh);
   //       std::cout << msg << std::endl;
   //       sprintf(msg,"yl = %.3lf < y0 = %.3lf < yh = %.3lf",yl,y0,yh);
   //       std::cout << msg << std::endl;
   //       sprintf(msg,"F(%.3lf,%.3lf) = %.3lf",xl,yl,fll);
   //       std::cout << msg << std::endl;
   //       sprintf(msg,"F(%.3lf,%.3lf) = %.3lf",xl,yh,flh);
   //       std::cout << msg << std::endl;
   //       sprintf(msg,"F(%.3lf,%.3lf) = %.3lf",xh,yl,fhl);
   //       std::cout << msg << std::endl;
   //       sprintf(msg,"F(%.3lf,%.3lf) = %.3lf",xh,yh,fhh);
   //       std::cout << msg << std::endl;
   //       std::cout << "-------------------------------------------------------" << std::endl;
   //    }

   //    // put everything together 
   //    double f00 = ywh*(xwh*fll + xwl*fhl) + ywl*(xwh*flh + xwl*fhh);
   //    return f00;
   // }
   //______________________________________________________________________________
   double AllanVariance(std::vector<double> x, const int M){
      // compute the allan variance for M points per group 
      // in the data vector x of size N points.

      const int N = x.size();

      double sum_sq=0,mean=0,mean_prev=0,diff=0;
      double v[M];

      int nAVG=0,j=0;
      for(int i=1;i<N+1;i++){
	 v[j] = x[i-1];
	 if(i%M==0){
	    // got the number of points we need, compute average 
	    nAVG++;
	    mean = GetMean<double>(M,v);
	    diff = mean-mean_prev;
	    // sum over squared differences for nAVG>1 
	    if(nAVG>1) sum_sq += pow(diff,2.);
	    // reset the j index and mean_prev now that we computed an average 
	    j = 0;
	    mean_prev = mean;
	 }else{
	    // keep going
	    j++;
	 }
      }

      // now get the variance 
      double allanVar = sum_sq/( 2.*(nAVG-1) );
      return allanVar;
   }
   //______________________________________________________________________________
   double AllanDeviation(std::vector<double> x, const int M){
      // compute the allan deviation for M points per group 
      // in the data vector x of size N points.
      double allanVar = AllanVariance(x,M);
      double allanDev = sqrt(allanVar);
      return allanDev;
   }
}
