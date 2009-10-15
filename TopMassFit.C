void TopMassFit( int type, int Nbtag, TString channel="had"  ){

  int rb =10 ;

  //1. load the Class
  gROOT->ProcessLine(".L MassFitFunction.cc+");
  gROOT->ProcessLine(".L MassAnaInput.cc+");
  gROOT->ProcessLine(".L MassAna.cc+");
  gROOT->ProcessLine(".L AlgoZero.cc+");
  gROOT->ProcessLine(".L AlgoKcon.cc+");
  gROOT->ProcessLine(".L MassAnaOutput.cc+");
  MassAna *theFitter = new MassAna( channel, Nbtag, 0, 480 );
  MassAnaOutput *FitResult = new MassAnaOutput( channel, Nbtag, 0, 480 );

  //               signal  tt   wjets   stt   sttw   QCD
  Bool_t comp[6] ={false, true, true, true, true, false};
  if ( type == 0 ) cout<<" Just Compile "<<endl;

  //theFitter->getJetPermutation();
  //theFitter->GetAllCoeff("171",rb, 110, 330) ;
  //theFitter->FitSignal("171",rb) ;
  //theFitter->FitTtbar("171",rb) ;
  if (type == 1) FitResult->CoeffCalib( rb, 110, 330 );
  
  if (type == 2) FitResult->MassCalib( rb, 110, 330, comp, Nbtag, 12 , false );
  if (type == 3) FitResult->MassCalib( rb, 110, 330, comp, Nbtag,  9 , true  );


  gROOT->Reset();

}
