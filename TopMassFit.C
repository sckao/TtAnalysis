void TopMassFit( int type,  TString channel="had"  ){

  int rb =10 ;
  vector<double> testx ;

  //1. load the Class
  gROOT->ProcessLine(".L MassFitFunction.cc+");
  gROOT->ProcessLine(".L MassAnaInput.cc+");
  gROOT->ProcessLine(".L MassAna.cc+");
  gROOT->ProcessLine(".L AlgoZero.cc+");
  gROOT->ProcessLine(".L AlgoKcon.cc+");
  gROOT->ProcessLine(".L HadWMassFitter.cc+");
  gROOT->ProcessLine(".L JES.cc+");
  gROOT->ProcessLine(".L MassAnaOutput.cc+");

  MassAna        *theFitter = new MassAna( channel, 0, 480 );
  MassAnaOutput  *FitResult = new MassAnaOutput( channel, 0, 480 );
  HadWMassFitter *WFitter   = new HadWMassFitter( 0, 480 );
  JES            *JESfromW  = new JES( 0, 480 );

  //               signal  tt   wjets   stt   sttw   QCD
  //Bool_t comp[6] ={false, true, true, true, true, false};

  if ( type == 0 ) cout<<" Just Compile "<<endl;

  if ( type == 99 )  {
     //FitResult->test();
     theFitter->FitBackground( "171", rb, 110, 330 );
     //theFitter->FitSignal("171",rb) ;
  }

  //theFitter->getJetPermutation();
  //theFitter->GetAllCoeff("171",rb, 110, 330) ;
  //theFitter->FitTtbar("171",rb) ;
  if (type == 1) FitResult->CoeffCalib( rb, 110, 330 );
  
  if (type == 2) FitResult->MassCalib( rb, 110, 330, -1, 12 , false );
  if (type == 3) FitResult->MassCalib( rb, 110, 330, -1,  9 , true  );

  if (type == 4) theFitter->FitSignal1("171",rb) ;

  if (type == 5) {
  
     int typeJES = 1;
     bool mcmatched =  true ;
     TString DrawOpt = "COLZ";
     WFitter->TWFitter( "tt171_336", typeJES, DrawOpt, mcmatched );
     //WFitter->TWFitter( "tt171_336", typeJES, DrawOpt );
     //WFitter->TWFitter( "QCD+",  typeJES, DrawOpt );
     //WFitter->TWFitter( "wj_336", typeJES, DrawOpt );
     
     /*
     WFitter->TWFitter( "zj_336", typeJES );
     WFitter->TWFitter( "ww_336", typeJES );
     WFitter->TWFitter( "wz_336", typeJES );
     
     WFitter->TWFitter( "tq_336", typeJES );
     WFitter->TWFitter( "tw_336", typeJES );
     */
  }

  if (type == 6) {
     JESfromW->GetJESFromW( "161" );
     JESfromW->GetJESFromW( "166" );
     JESfromW->GetJESFromW( "171" );
     JESfromW->GetJESFromW( "176" );
     JESfromW->GetJESFromW( "181" );
  }   

  if (type == 7) {
     JESfromW->CalibJES( "171" );
  }

  if (type == 8) {
     JESfromW->Smearing( "171" );
  }

  if (type == 9) {
     JESfromW->SmearAndMatch( "171" );
  }
 
  if (type == 10) {
     JESfromW->Matching( "171" );
  }

  gROOT->Reset();

}
