void TopMassFit( TString algo, TString channel , int step){

  //if (!TClass::GetDict("TemplateMassFit")) {
  //     gROOT->ProcessLine(".L TemplateMassFit.cc");
  //}
  
  //1. load the Class
  gROOT->ProcessLine(".L TemplateMassFit.cc");
  TemplateMassFit *theFitter = new TemplateMassFit(algo,channel);

  gROOT->ProcessLine(".L MassFitFunction.cc");
  //MassFitFunction *theFunctions = new MassFitFunction();

  if ( step == 1 ) theFitter->FitTester(algo, 10, 130, 330) ;
  
  if ( step == 2 ) {  
     int rb = 10;
     theFitter->FitTester("161", rb, 110, 330) ;
     theFitter->FitTester("163", rb, 110, 330) ;
     theFitter->FitTester("166", rb, 110, 330) ;
     theFitter->FitTester("168", rb, 110, 330) ;
     theFitter->FitTester("171", rb, 110, 330) ;
     theFitter->FitTester("173", rb, 110, 330) ;
     theFitter->FitTester("176", rb, 110, 330) ;
     theFitter->FitTester("178", rb, 110, 330) ;
     theFitter->FitTester("181", rb, 110, 330) ;
  }

  if ( step == 3 ) {
     //               signal   tt   wjets  QCD
     Bool_t comp[4] ={false, false, true, false}; 
     theFitter->TemplateDrawer(algo,10,110,330, comp) ;
  }
 
  if ( step == 4 ) theFitter->MoreCombinedFitting( 10, 110, 330, 9) ;
  if ( step == 5 ) theFitter->CombinedFitting(5,110,350) ;
  if ( step == 6 ) theFitter->MultiTemplatesFitting(5,100,330) ;
  if ( step == 7 ) theFitter->TemplateFitting(5,100,330) ;

  gROOT->Reset();

}
