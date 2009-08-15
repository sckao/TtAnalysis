#include "TObject.h"
#include <vector>
#include <stdio.h>
#include <iostream>
#include <string>
#include <TMath.h>
#include <TF1.h>


class MassFitFunction : public TObject {

private:


public:

   MassFitFunction();     
   ~MassFitFunction();     
  
   static Double_t fitG( Double_t* x, Double_t* par);
   static Double_t fitLG( Double_t* x, Double_t* par);
   static Double_t fitGS( Double_t* x, Double_t* par);
   static Double_t fitTL( Double_t* x, Double_t* par);
   static Double_t fitLD( Double_t* x, Double_t* par);
   static Double_t fitSG( Double_t* x, Double_t* par);
   static Double_t fitBW( Double_t* x, Double_t* par);
   static Double_t fitData( Double_t* x, Double_t* par);
   static Double_t fitData1( Double_t* x, Double_t* par);
   static Double_t fitParabola( Double_t *x, Double_t *par);
   static Double_t ConvBWGS( Double_t* x, Double_t* par);
   static Double_t ConvSGGS( Double_t* x, Double_t* par);

   std::vector<bool> DataRejection( TF1* fitfunc, Double_t* x, Double_t* y, int N_data );
   bool DataRejection(double sigma, double deviation, int N_data );

   ClassDef(MassFitFunction, 1);

};

#if !defined(__CINT__)
    ClassImp(MassFitFunction);
#endif

