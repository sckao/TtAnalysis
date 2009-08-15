void TopMassFit( TString algo, TString channel , int step){

  
  //1. load the Class
  gROOT->ProcessLine(".L MassFitFunction.cc+");
  gROOT->ProcessLine(".L TemplateMassFit.cc+");
  TemplateMassFit *theFitter = new TemplateMassFit(algo,channel);

  if ( step == 0 ) cout<<" Just compile "<<endl;  

  if ( step == 1 ) theFitter->FitSignal(algo, 10, 110, 330) ;
  
  if ( step == 2 ) {
     //               signal   tt   wjets  stt  sttw  QCD
     Bool_t comp[6] ={false, true,  true, true, true, false}; 
     theFitter->>FitBackground(algo,10,110,330, comp) ;
  }

  if ( step == 3 ) {  
     int rb = 10;
     //               signal   tt   wjets  stt  sttw  QCD
     Bool_t comp[6] ={true, true,  true, true, true, false}; 
     theFitter->GetAllCoeff("161", rb, 110, 330, comp) ;
     theFitter->GetAllCoeff("163", rb, 110, 330, comp) ;
     theFitter->GetAllCoeff("166", rb, 110, 330, comp) ;
     theFitter->GetAllCoeff("168", rb, 110, 330, comp) ;
     theFitter->GetAllCoeff("171", rb, 110, 330, comp) ;
     theFitter->GetAllCoeff("173", rb, 110, 330, comp) ;
     theFitter->GetAllCoeff("176", rb, 110, 330, comp) ;
     theFitter->GetAllCoeff("178", rb, 110, 330, comp) ;
     theFitter->GetAllCoeff("181", rb, 110, 330, comp) ;
     theFitter->GetAllCoeff("183", rb, 110, 330, comp) ;
  }

  if ( step == 4 ) {
     //               signal   tt   wjets  stt  sttw  QCD
     Bool_t comp[6] ={false, true,  true, true, true, false}; 
     theFitter->TemplateDrawer(algo,10,110,330, comp) ;
  }
  // Fit top mass with 12 parameters , separate tt-wrong combinatorics and other background
  if ( step == 5 ) theFitter->MoreCombinedFitting( 10, 110, 330 ) ;
  // Fit top mass with 9 parameters , combined tt-wrong combinatorics and other background
  if ( step == 6 ) theFitter->CombinedFitting( 10, 110, 330, 2 ) ;
  // TemplateFitting Methods...under developing
  if ( step == 7 ) theFitter->MultiTemplatesFitting(5,100,330) ;
  if ( step == 8 ) theFitter->TemplateFitting(5,100,330) ;


  gROOT->Reset();

}
