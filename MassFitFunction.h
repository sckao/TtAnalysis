#include "TObject.h"
#include <vector>
#include <stdio.h>
#include <iostream>
#include <string>

class MassFitFunction : public TObject {

private:


public:

   MassFitFunction();     
  
   Double_t fitG( Double_t* x, Double_t* par);
   Double_t fitLG( Double_t* x, Double_t* par);
   Double_t fitGS( Double_t* x, Double_t* par);
   Double_t fitTL( Double_t* x, Double_t* par);
   Double_t fitLD( Double_t* x, Double_t* par);
   Double_t fitSG( Double_t* x, Double_t* par);
   Double_t fitBW( Double_t* x, Double_t* par);
   Double_t fitData( Double_t* x, Double_t* par);
   Double_t fitData1( Double_t* x, Double_t* par);
   Double_t fitParabola( Double_t *x, Double_t *par);
   Double_t ConvBWGS( Double_t* x, Double_t* par);
   Double_t ConvSGGS( Double_t* x, Double_t* par);

   ClassDef(MassFitFunction, 1);

};

#if !defined(__CINT__)
    ClassImp(MassFitFunction);
#endif

