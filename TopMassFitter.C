void TopMassFitter( TString algo, TString channel ){

  //if (!TClass::GetDict("TemplateMassFit")) {
  //     gROOT->ProcessLine(".L TemplateMassFit.cc");
  //}
  
  //1. load the Class
  gROOT->ProcessLine(".L TemplateMassFit.cc");
  TemplateMassFit *theFitter = new TemplateMassFit(algo,channel);

  //theFitter->showAll(5,120,310) ;
  theFitter->MultiTemplatesFitting(5,110,300) ;
  //theFitter->TemplateFitting(5,120,310) ;
  //theFitter->CombinedFitting(5,120,310) ;
  gROOT->Reset();

}
