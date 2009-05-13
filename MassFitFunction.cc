#include "MassFitFunction.h"
MassFitFunction::MassFitFunction(){
  
  
}

MassFitFunction::~MassFitFunction(){
 
}

Double_t MassFitFunction::fitG( Double_t* x, Double_t* par){

     Double_t A1 = (par[0]/x[0]) + (par[1]*x[0]) + par[2];
     Double_t A2 =  exp( A1 ) ;
     Double_t fitV = A2 ;

     return fitV;
}

// log-normal distribution
Double_t MassFitFunction::fitLG(Double_t *x, Double_t *par) {

     Double_t A0 =  log( x[0] ) - par[1] ;
     Double_t A1 =   (-1.*par[2]*A0*A0)  ;
     Double_t fitV = par[0]*exp( A1 ) / x[0]  ;
     return fitV;
}

Double_t MassFitFunction::fitTL(Double_t *x, Double_t *par) {

     //Double_t A1 = exp(par[1]*x[0]) + par[3]*log( x[0] + 100. );
     Double_t A0 = x[0] + 50 ;
     Double_t A2 = x[0]  ;
     Double_t A1 = exp(par[1]*x[0]) + par[3]/A0 
                 + par[4]*A2 + par[5]*A2*A2 + par[6]*A2*A2*A2;
     Double_t fitV = exp( par[0]/x[0] + par[2] ) * A1 ;
     return fitV;

}

Double_t MassFitFunction::fitLD(Double_t *x, Double_t *par) {

      Double_t ld_Value = TMath::Landau(x[0],par[1],par[2]) ;
      Double_t fitV = par[0]*ld_Value ;
      return fitV ;

}
Double_t MassFitFunction::fitParabola( Double_t *x, Double_t *par ){

     Double_t chi = (x[0] - par[0]) / (par[1]+1.);
     Double_t fV = chi*chi + par[2] ;
     return fV;
}

Double_t MassFitFunction::fitBW(Double_t *x, Double_t *par) {

     Double_t gm = par[0] + 0.001;
     Double_t chi = (x[0] - par[1]) / gm ;
     Double_t C1 = 1. + (chi*chi) ;
     Double_t fV = par[2]/C1 ;
     return fV;
}

Double_t MassFitFunction::fitGS(Double_t *x, Double_t *par) {

     Double_t gs_Value = TMath::Gaus(x[0],par[1],par[2]) ;
     Double_t fitV = par[0]*gs_Value ; 
     return fitV;
}

Double_t MassFitFunction::fitSG(Double_t *x, Double_t *par) {

     Double_t gs = TMath::Gaus(x[0],par[1],par[2]);
     Double_t cb_Val = gs + fitLG( x, &par[3] ) ;
     return cb_Val = par[0] * cb_Val ;
     //return fitGS(x,par) + fitLG(x,&par[3]);
     //return fitGS(x,par) + fitLG(x,&par[3]);
}

Double_t MassFitFunction::fitData(Double_t *x, Double_t *par) {

        //return fitBW(x,par) + fitGS(x, &par[3]) ;
        return fitSG(x,par) + fitLD(x, &par[6]) + fitLD(x, &par[9]) ;
        //return fitSG(x,par) + fitLD(x, &par[6]) ;
}

Double_t MassFitFunction::fitData1(Double_t *x, Double_t *par) {

     Double_t gs = TMath::Gaus(x[0],par[1],par[2]);

     Double_t A0 =  log( x[0] ) - par[4] ;
     Double_t A1 =   (-1.*par[5]*A0*A0)  ;
     Double_t LG_Val = par[3]*exp( A1 ) / x[0]  ;

     Double_t ld1_Val = TMath::Landau(x[0],par[6],par[7]) ;

     Double_t ld2_Val = TMath::Landau(x[0],par[8],par[9]) ;

     //Double_t sg_Val = gs + LG_Val + (17.14*ld1_Val) + (6.69*ld2_Val) ;
     Double_t sg_Val = gs + LG_Val + (par[10]*ld1_Val) + (par[11]*ld2_Val) ;
     Double_t fitV = par[0]*sg_Val ;

     return fitV ;
}


Double_t MassFitFunction::ConvBWGS(Double_t *x, Double_t *par) {

   Double_t np   = 150 ;
   Double_t sg   = 4.0 ;
   Double_t xlow = x[0] - sg*par[4];
   Double_t xup  = x[0] + sg*par[4];
   Double_t step = (xup - xlow) / np ;
   Double_t sum = 0. ;
   Double_t xx, fbw; 

   for (double i=1.; i<= np/2; i++  ) {
       xx = xlow + (i-.5) * step;
       //fbw = TMath::BreitWigner(xx,par[1],par[0]);
       fbw = fitBW(x,par);
       sum += fbw * TMath::Gaus(x[0],par[3],par[4]);

       xx = xup - (i-.5) * step;
       //fbw = TMath::BreitWigner(xx,par[1],par[0]);
       fbw = fitBW(x,par);
       sum += fbw * TMath::Gaus(x[0],par[3],par[4]);
   }
	
   //return par[2]*step*sum ;
   return step*sum ;

}

Double_t MassFitFunction::ConvSGGS(Double_t *x, Double_t *par) {

   Double_t np   = 200 ;
   Double_t sg   = 5.0 ;
   Double_t xlow = x[0] - sg*par[7];
   Double_t xup  = x[0] + sg*par[7];
   Double_t step = (xup - xlow) / np ;
   Double_t sum = 0. ;
   Double_t xx, fbw; 

   for (double i=1.; i<= np/2; i++  ) {
       xx = xlow + (i-.5) * step;
       fbw = fitSG(x,par);
       sum += fbw * TMath::Gaus(x[0],par[6],par[7]);

       xx = xup - (i-.5) * step;
       fbw = fitSG(x,par);
       sum += fbw * TMath::Gaus(x[0],par[6],par[7]);
   }
	
   return step*sum ;
}

