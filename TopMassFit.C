void TopMassFit( TString algo, TString channel ){

  //if (!TClass::GetDict("TemplateMassFit")) {
  //     gROOT->ProcessLine(".L TemplateMassFit.cc");
  //}
  
  //1. load the Class
  gROOT->ProcessLine(".L TemplateMassFit.cc");
  TemplateMassFit *theFitter = new TemplateMassFit(algo,channel);

  //theFitter->showAll(5,120,310) ;
  theFitter->MultiTemplatesFitting(5,100,330) ;
  //theFitter->TemplateFitting(5,100,330) ;
  //theFitter->CombinedFitting(5,110,350) ;
  //theFitter->MoreCombinedFitting(5,100,350) ;
  //Bool_t comp[4] ={true, true, true, true}; 
  //theFitter->TemplateDrawer(algo,5,100,330, comp) ;
  gROOT->Reset();

}
