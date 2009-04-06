#include "TemplateMassFit.h"
TemplateMassFit::TemplateMassFit( TString aname, TString channel ){
  
  file0  = TFile::Open("ttFakeData1.root");
  file0a = TFile::Open("wjFakeData1.root");

  //file1  = TFile::Open("ttj_fall08_"+aname+".root");
  file1  = TFile::Open(aname+".root");
  //file2  = TFile::Open("wjets_fall08_"+aname+".root");
  file2  = TFile::Open("wjets_fall08_zero.root");
  //file3  = TFile::Open("qcd_fall08_"+aname+".root");
  file3  = TFile::Open("qcd_fall08_zero.root");

  hname =channel+"_selTmass_pt0";
  sname ="mcTtMass0";
  cname = channel ;
  hfolder = "tt_fall08";

  plot1 = hfolder+"/"+channel+"_TMass_"+ aname+".gif";
  plot2 = hfolder+"/"+channel+"_CombinedFit_"+aname+".gif";
  plot3 = hfolder+"/"+channel+"_Norm_"+aname+".gif";
  plot4 = hfolder+"/"+channel+"_TemplateFit_"+aname+".gif";
  plot5 = hfolder+"/"+channel+"_Norm2_"+aname+".gif";
  plot6 = hfolder+"/"+channel+"_TemplateFit2_"+aname+".gif";

  gSystem->mkdir(hfolder);
  
}

TemplateMassFit::~TemplateMassFit(){
 
  delete c1;
  delete c2;
  delete c3;
  delete c4;

  file0->Close();
  file0a->Close();
  file1->Close();
  file2->Close();
  file3->Close();
  
}

void TemplateMassFit::showAll( int rbin, int lowBound, int upBound ){

  int nbin = 200./rbin ;
 
  //TemplateFitting(rbin, lowBound, upBound );
  //MultiTemplatesFitting(rbin, lowBound, upBound );
  //CombinedFitting(rbin, lowBound, upBound );

  gStyle->SetOptFit(0111);
  c1 = new TCanvas("c1","", 900, 700);
  c1->SetFillColor(10);
  c1->SetFillColor(10);

  TH1D* sb = new TH1D("sb","", nbin, 0., 400.);
  TH1D* tt = new TH1D("tt","", nbin, 0., 400.);
  TH1D* wj = new TH1D("wj","", nbin, 0., 400.);
  TH1D* qm = new TH1D("qm","", nbin, 0., 400.);
  TH1D* sg = new TH1D("sg","", nbin, 0., 400.);

  getBackground( sb, 0, rbin );
  getBackground( tt, 1, rbin );
  getBackground( wj, 2, rbin );
  getBackground( qm, 3, rbin );
  getSignal( sg, rbin );

  //1. signal + bad combinatorics 
  gStyle->SetOptStat("nirm");
  gStyle->SetStatY(0.99);
  gStyle->SetStatX(0.55);
  gStyle->SetStatTextColor(3);
  sb->SetFillColor(3);
  sb->SetTitle(" mass spectrum ");
  //sb->SetAxisRange(1,380,"X");
  sb->SetAxisRange(1,600,"Y");
  sb->Draw();
  c1->Update();

  //2. bad combinatorics
  gStyle->SetStatY(0.99);
  gStyle->SetStatX(0.75);
  gStyle->SetStatTextColor(7);
  tt->SetFillColor(7);
  tt->Draw("sames");
  c1->Update();

  //3. w jets background
  gStyle->SetStatY(0.99);
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(2);
  wj->SetFillColor(2);
  wj->Draw("sames");
  c1->Update();

  //4. qcd background
  gStyle->SetStatY(0.99);
  gStyle->SetStatX(0.35);
  gStyle->SetStatTextColor(4);
  qm->SetFillColor(4);
  qm->Draw("sames");
  c1->Update();

  //5. signal only
  gStyle->SetStatY(0.70);
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(1);
  sg->SetLineColor(1);
  sg->SetLineWidth(2);
  sg->Draw("sames");
  c1->Update();
  
  TF1 *func2 = new TF1("func2", TemplateMassFit::fitBW, lowBound, upBound,3);
  func2->SetParLimits(0, 10.,70.);
  func2->SetParLimits(1, 150.,220.);
  func2->SetParLimits(2, 50.,120.);
  func2->SetLineColor(1);
  func2->SetLineStyle(2);
  func2->SetLineWidth(3);
  sg->Fit( func2 , "R","same",lowBound,upBound);
  c1->Update();
   
  c1->Print(plot1);

  delete c1;
  delete tt;
  delete wj;
  delete qm;
  delete sg;
  delete sb;
  delete func2;

}

void TemplateMassFit::CombinedFitting( int rbin, int lowBound, int upBound ){

  gStyle->SetOptStat("nirm");
  gStyle->SetOptFit(0111);
  gStyle->SetStatY(0.99);
  gStyle->SetStatX(0.99);
  gStyle->SetStatTextColor(1);
  c3 = new TCanvas("c3","", 900, 700);
  c3->SetFillColor(10);
  c3->SetFillColor(10);

  int nbin = 200./rbin ;

  THStack* ttstk = new THStack("ttstk", "Combined Fitting"); 
  TH1D*    fakedata = new TH1D("fakedata","", nbin, 0., 400.);

  getFakeData( fakedata, ttstk, rbin);
  ttstk->Draw();

  // pre-fit background
  TH1D* bgadd = new TH1D("bgadd","", nbin, 0., 400.);
  TH1D* tt = new TH1D("tt","", nbin, 0., 400.);
  TH1D* wj = new TH1D("wj","", nbin, 0., 400.);
  TH1D* qm = new TH1D("qm","", nbin, 0., 400.);

  getBackground( qm, 3, rbin );
  bgadd->Add(qm);
  getBackground( wj, 2, rbin );
  bgadd->Add(wj);
  getBackground( tt, 1, rbin );
  bgadd->Add(tt);

  // smooth the background channels and get parameters for background shap
  TF1 *func0 = new TF1("func0",TemplateMassFit::fitG, lowBound, upBound,3);
  bgadd->Fit( func0, "R0","",110,350);
  double bp0 = func0->GetParameter(0);
  double bp1 = func0->GetParameter(1);


  // Fit the data
  TF1 *func1 = new TF1("func1",TemplateMassFit::fitData, lowBound, upBound,6);
  func1->SetParLimits(0, 10.,70.);
  func1->SetParLimits(1, 150.,220.);
  func1->SetParLimits(2, 10.,120.);
  func1->FixParameter(3, bp0 );
  func1->FixParameter(4, bp1 );
  func1->SetLineColor(4);
  func1->SetLineStyle(2);
  func1->SetLineWidth(2);
  fakedata->Fit( func1, "R","sames",110,350);
  c3->Update();


  // Draw the expected signal  
  double width = func1->GetParameter(0);
  double peak = func1->GetParameter(1);
  double amp = func1->GetParameter(2);

  TF1 *func2 = new TF1("func2",TemplateMassFit::fitBW, 100, 360, 3 );
  func2->FixParameter(0, width);
  func2->FixParameter(1, peak );
  func2->FixParameter(2, amp  );
  func2->SetLineColor(1);
  func2->SetLineStyle(2);
  func2->SetLineWidth(3);
  func2->Draw("sames");
  c3->Update();

  // Draw the expected background
  double p0  = func1->GetParameter(3);
  double p1  = func1->GetParameter(4);
  double p2  = func1->GetParameter(5);

  TF1 *func3 = new TF1("func3",TemplateMassFit::fitG, 100, 360, 3 );
  func3->FixParameter(0, p0 );
  func3->FixParameter(1, p1 );
  func3->FixParameter(2, p2 );
  func3->SetLineColor(1);
  func3->SetLineStyle(2);
  func3->SetLineWidth(3);
  func3->Draw("sames");
  c3->Update();

  c3->Print(plot2);

  delete tt;
  delete wj;
  delete qm;
  delete func0;
  delete func1;
  delete func2;
  delete func3;
  delete ttstk;
  delete fakedata;
 
}

void TemplateMassFit::MultiTemplatesFitting( int rbin, int lowBound, int upBound ){

  gStyle->SetOptStat("nirm");
  gStyle->SetOptFit(0111);
  gStyle->SetStatY(0.95);
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(1);
  c5 = new TCanvas("c5","", 1000, 800);
  c5->SetFillColor(10);
  c5->SetFillColor(10);
  c5->Divide(2,2);

  int nbin = 200./rbin ;
  Int_t b1 = lowBound / (rbin*2);
  Int_t b2 =  upBound / (rbin*2);
  const Int_t sz = b2 - b1 + 1 ;
  cout<<" b1:"<<b1<<" b2:"<<b2<<" size:"<< sz <<endl;
 
  // tmp1: tt wrong combinatorics
  TH1D* tt = new TH1D("tt","", nbin, 0., 400.);
  getBackground( tt, 1, rbin );

  // tmp2: wjets and QCD 
  TH1D* wj = new TH1D("wj","", nbin, 0., 400.);
  TH1D* qm = new TH1D("qm","", nbin, 0., 400.);
  TH1D* bg = new TH1D("bg","", nbin, 0., 400.);
  getBackground( wj, 2, rbin );
  getBackground( qm, 3, rbin );
  bg->Add(wj,qm,1,2);

  // tmp3: tt correct combinatoric
  TH1D* sg = new TH1D("sg","", nbin, 0., 400.);
  getSignal( sg, rbin );

  // normalize signal and background shape
  ScaleTemplates( 5000. , bg, b1, b2 );
  ScaleTemplates( 5000. , tt, b1, b2 );
  ScaleTemplates( 5000. , sg, b1, b2 );

  // get data distribution
  TH1D* dt = new TH1D("dt","", nbin, 0., 400.);
  //getData( dt, rbin );
  getFakeData( dt, rbin );

  // smooth template distribution by parameterization
  TH1D* parabg = new TH1D("parabg","", nbin, 0., 400.);
  SmoothTemplate( 0 , bg, parabg, lowBound, upBound, 3);
  TH1D* paratt = new TH1D("paratt","", nbin, 0., 400.);
  SmoothTemplate( 0 , tt, paratt, lowBound, upBound, 3);
  TH1D* parasg = new TH1D("parasg","", nbin, 0., 400.);
  SmoothTemplate( 1 , sg, parasg, lowBound, upBound, 3);

  //1a display the signal distribution
  c5->cd(1);
  sg->SetLineColor(1);
  sg->SetLineWidth(2);
  sg->SetTitle(" Signal");
  sg->Draw();
  c5->Update();

  //1b display the parameterized signal distribution
  gStyle->SetStatX(0.32);
  gStyle->SetStatTextColor(2);
  parasg->SetLineColor(2);
  parasg->SetLineStyle(2);
  parasg->SetLineWidth(3);
  parasg->Draw("sames");
  c5->Update();
  gStyle->SetStatX(0.95);
  
  //2a display the background distribution
  c5->cd(2);
  gStyle->SetStatTextColor(1);
  bg->SetLineColor(1);
  bg->SetLineWidth(2);
  bg->SetTitle("background WJets and QCD");
  bg->Draw(); 
  c5->Update();

  //2b display the parameterized background distribution
  gStyle->SetStatX(0.32);
  gStyle->SetStatTextColor(2);
  parabg->SetLineColor(2);
  parabg->SetLineStyle(2);
  parabg->SetLineWidth(3);
  parabg->Draw("sames");
  c5->Update();
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(1);

  //3a display the background distribution
  c5->cd(3);
  gStyle->SetStatTextColor(1);
  tt->SetLineColor(1);
  tt->SetLineWidth(2);
  tt->SetTitle("tt wrong combinatorics");
  tt->Draw(); 
  c5->Update();

  //3b display the parameterized background distribution
  gStyle->SetStatX(0.32);
  gStyle->SetStatTextColor(2);
  paratt->SetLineColor(2);
  paratt->SetLineStyle(2);
  paratt->SetLineWidth(3);
  paratt->Draw("sames");
  c5->Update();
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(1);
  c5->Print(plot5);

  TObjArray *mctmp = new TObjArray(3);
  mctmp->Add(parasg);
  mctmp->Add(parabg);
  mctmp->Add(paratt);
 
  TFractionFitter* mfit = new TFractionFitter( dt, mctmp );
  mfit->Constrain( 0, 0.0, 1.0);
  mfit->Constrain( 1, 0.0, 1.0);
  mfit->Constrain( 2, 0.0, 1.0);
  mfit->SetRangeX(b1,b2);
  Int_t status = mfit->Fit();
  cout<<" fit status : "<< status << endl;
 
  if ( status == 0 ) {

     TH1D* result  = (TH1D*) mfit->GetPlot();
     double ratio[3];
     double err[3];
     mfit->GetResult(0,ratio[0],err[0]);
     mfit->GetResult(1,ratio[1],err[1]);
     mfit->GetResult(2,ratio[2],err[2]);

     double nData = dt->Integral(b1, b2);
     double nSig  = sg->Integral(b1, b2);
     double nBkgd = bg->Integral(b1, b2);
     double nBktt = tt->Integral(b1, b2);
     cout<<" nData= "<<nData<<"  nSignal= "<<nSig<<"  nBgd="<<nBkgd<<" nBadtt="<<nBktt<<endl;
     cout<<" r0:"<<ratio[0] <<"  r1:"<<ratio[1]<<" r2:"<<ratio[2]<<endl;

     double x[sz]     = {0.0};
     double sPred[sz] = {0.0};
     double bPred[sz] = {0.0};
     double tPred[sz] = {0.0};
     double aPred[sz] = {0.0};
     for(int k=0; k<sz ; k++) {
        int i = b1 + k ;
        x[k]  = 2.0*rbin*i - rbin ;

        double sN = sg->GetBinContent(i); 
        double bN = bg->GetBinContent(i); 
        double tN = tt->GetBinContent(i); 
        sPred[k] = (nData/nSig)*ratio[0]*sN ;
        bPred[k] = (nData/nBkgd)*ratio[1]*bN ;
        tPred[k] = (nData/nBktt)*ratio[2]*tN ;
        aPred[k] = sPred[k] + bPred[k] + tPred[k] ;
     }

    c6 = new TCanvas("c6","", 900, 700);
    c6->SetFillColor(10);
    c6->SetFillColor(10);
    
    result->SetLineColor(2);
    result->SetLineWidth(1);
    result->SetAxisRange(0,550,"Y");
    result->SetTitle(" Template Fitting Result ");
    result->Draw();
    
    hPred = new TGraph(sz,x,aPred);
    hPred->SetMarkerColor(6);
    hPred->SetMarkerStyle(20);
    hPred->SetLineWidth(1);
    hPred->SetLineColor(6);
    hPred->Draw("SAMEP");

    hPred0 = new TGraph(sz,x,sPred);
    hPred0->SetMarkerColor(1);
    hPred0->SetMarkerStyle(20);
    hPred0->SetMaximum(160.);
    hPred0->Draw("SAMECP");
    
    hPred1 = new TGraph(sz,x,bPred);
    hPred1->SetMarkerColor(4);
    hPred1->SetLineColor(4);
    hPred1->SetMarkerStyle(21);
    hPred1->Draw("SAMECP");
    
    hPred2 = new TGraph(sz,x,tPred);
    hPred2->SetMarkerColor(7);
    hPred2->SetLineColor(7);
    hPred2->SetMarkerStyle(22);
    hPred2->Draw("SAMECP");
    
    dt->SetLineWidth(2);
    dt->SetLineColor(8);
    dt->Draw("same");

    c6->Update();
  }

  c6->Print(plot6);

  delete bg;
  delete wj;
  delete qm;
  delete tt;
  delete sg;
  delete dt;
  delete mfit;
  delete hPred;
  delete hPred0;
  delete hPred1;
  delete hPred2;
  delete result;
  delete mctmp;
  delete parabg;
  delete parasg;
}

void TemplateMassFit::TemplateFitting( int rbin, int lowBound, int upBound ){

  gStyle->SetOptStat("nirm");
  gStyle->SetOptFit(0111);
  gStyle->SetStatY(0.95);
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(1);
  c2 = new TCanvas("c2","", 1000, 800);
  c2->SetFillColor(10);
  c2->SetFillColor(10);
  c2->Divide(1,2);

  int nbin = 200./rbin ;
  Int_t b1 = lowBound / (rbin*2);
  Int_t b2 =  upBound / (rbin*2);
  const Int_t sz = b2 - b1 + 1 ;
  cout<<" b1:"<<b1<<" b2:"<<b2<<" size:"<< sz <<endl;
 
  TH1D* bg = new TH1D("bg","", nbin, 0., 400.);
  TH1D* sg = new TH1D("sg","", nbin, 0., 400.);
  TH1D* dt = new TH1D("dt","", nbin, 0., 400.);

  // combine all background channel
  combineBG( bg, rbin );
  // get signal distribution
  getSignal( sg, rbin );
  // normalize signal and background shape
  ScaleTemplates( 5000., sg, b1, b2 );
  ScaleTemplates( 5000., bg, b1, b2 );
  // get data distribution
  getFakeData( dt, rbin );

  TH1D* parabg = new TH1D("parabg","", nbin, 0., 400.);
  SmoothTemplate( 0 , bg, parabg, lowBound, upBound, 3);
  TH1D* parasg = new TH1D("parasg","", nbin, 0., 400.);
  SmoothTemplate( 1 , sg, parasg, lowBound, upBound, 3);

  c2->cd(1);
  sg->SetLineColor(1);
  sg->SetLineWidth(2);
  sg->SetTitle(" Signal");
  sg->Draw();
  c2->Update();

  gStyle->SetStatX(0.32);
  gStyle->SetStatTextColor(2);
  parasg->SetLineColor(2);
  parasg->SetLineStyle(2);
  parasg->SetLineWidth(3);
  parasg->Draw("sames");
  c2->Update();
  gStyle->SetStatX(0.95);
  
  c2->cd(2);
  gStyle->SetStatTextColor(1);
  bg->SetLineColor(1);
  bg->SetLineWidth(2);
  bg->SetTitle("combined background");
  bg->Draw(); 
  c2->Update();

  gStyle->SetStatX(0.32);
  gStyle->SetStatTextColor(2);
  parabg->SetLineColor(2);
  parabg->SetLineStyle(2);
  parabg->SetLineWidth(3);
  parabg->Draw("sames");
  c2->Update();
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(1);
  c2->Print(plot3);

  TObjArray *mctmp = new TObjArray(2);
  mctmp->Add(parasg);
  mctmp->Add(parabg);
 
  TFractionFitter* mfit = new TFractionFitter( dt, mctmp );
  mfit->Constrain( 0, 0.0, 1.0);
  mfit->Constrain( 1, 0.0, 1.0);
  mfit->SetRangeX(b1,b2);
  Int_t status = mfit->Fit();
  cout<<" fit status : "<< status << endl;
 
  if ( status == 0 ) {

     TH1D* result  = (TH1D*) mfit->GetPlot();
     double ratio[2];
     double err[2];
     mfit->GetResult(0,ratio[0],err[0]);
     mfit->GetResult(1,ratio[1],err[1]);

     double nData = dt->Integral(b1, b2);
     double nSig  = sg->Integral(b1, b2);
     double nBkgd = bg->Integral(b1, b2);
     cout<<" nData= "<<nData<<"  nSignal= "<<nSig<<"  nBgd="<<nBkgd<<endl;
     cout<<" r0:"<<ratio[0] <<"  r1:"<<ratio[1]<<endl;

     double x[sz]     = {0.0};
     double sPred[sz] = {0.0};
     double bPred[sz] = {0.0};
     double aPred[sz] = {0.0};
     for(int k=0; k<sz ; k++) {
        int i = b1 + k ;
        x[k]  = 2.0*rbin*i - rbin ;

        double sN = sg->GetBinContent(i); 
        double bN = bg->GetBinContent(i); 
        sPred[k] = (nData/nSig)*ratio[0]*sN ;
        bPred[k] = (nData/nBkgd)*ratio[1]*bN ;
        aPred[k] = sPred[k] + bPred[k] ;
     }

    c4 = new TCanvas("c4","", 900, 700);
    c4->SetFillColor(10);
    c4->SetFillColor(10);
    
    result->SetLineColor(2);
    result->SetLineWidth(1);
    result->SetAxisRange(0,550,"Y");
    result->SetTitle(" Template Fitting Result ");
    result->Draw();
    
    hPred = new TGraph(sz,x,aPred);
    hPred->SetMarkerColor(6);
    hPred->SetMarkerStyle(20);
    hPred->SetLineWidth(1);
    hPred->SetLineColor(6);
    hPred->Draw("SAMEP");

    hPred0 = new TGraph(sz,x,sPred);
    hPred0->SetMarkerColor(1);
    hPred0->SetMarkerStyle(20);
    hPred0->SetMaximum(160.);
    hPred0->Draw("SAMECP");
    
    hPred1 = new TGraph(sz,x,bPred);
    hPred1->SetMarkerColor(4);
    hPred1->SetLineColor(4);
    hPred1->SetMarkerStyle(21);
    hPred1->Draw("SAMECP");
    
    dt->SetLineWidth(2);
    dt->SetLineColor(8);
    dt->Draw("same");

    c4->Update();
  }

  c4->Print(plot4);

  delete bg;
  delete sg;
  delete dt;
  delete mfit;
  delete hPred;
  delete hPred0;
  delete hPred1;
  delete result;
  delete mctmp;
  delete parabg;
  delete parasg;
}

void TemplateMassFit::ScaleTemplates( double factor, TH1D* tmp, int B1, int B2 ){

  double tmpN = tmp->Integral(B1,B2);
  double Scal = factor / tmpN ;
 
  tmp->Scale(Scal);
  
}

void TemplateMassFit::SmoothTemplate( int type, TH1D* ds, TH1D* ds1, int lowBound, int upBound, int npar ){

  TF1* func0;
  if (type == 0) func0 = new TF1("func0", TemplateMassFit::fitG  , lowBound, upBound, npar);
  if (type == 1) {
     func0 = new TF1("func0", TemplateMassFit::fitBW , lowBound, upBound, npar);
     func0->SetParLimits(0, 5.,80.);
     func0->SetParLimits(1, 150.,230.);
     func0->SetParLimits(2, 550.,1200.);
  }
  ds->Fit( func0, "R0","", lowBound, upBound );

  int nbin = ds->GetNbinsX() ;
  double bW = 400./nbin ;
  for (int i=0; i< nbin; i++) {
      double k = ( i - 1 )*bW + bW/2.  ;
      double theN = func0->Eval(k);
      ds1->SetBinContent(i,theN);
  }

  delete func0;

}

void TemplateMassFit::combineBG( TH1D* allbg, int rbin ) {

  int nbin = 200./rbin ;
  TH1D* tt = new TH1D("tt","", nbin, 0., 400.);
  TH1D* wj = new TH1D("wj","", nbin, 0., 400.);
  TH1D* qm = new TH1D("qm","", nbin, 0., 400.);

  getBackground( tt, 1, rbin );
  getBackground( wj, 2, rbin );
  getBackground( qm, 3, rbin );
  allbg->Add(tt,wj,1,1);
  allbg->Add(qm,2.);

  delete tt;
  delete wj;
  delete qm;
}

void TemplateMassFit::getFakeData( TH1D* ttadd, THStack* ttstk,  int rbin ){


  int nbin = 200./rbin ;
  hdata  = (TH2F *) file0->Get("Tops/"+hname);
  hdata->ProjectionX("hdata_pj",1,200,"");
  hdata_pj->Rebin(rbin);

  hMC   = (TH2D *) file0->Get("MObjs/"+sname);
  hMC->ProjectionX("hMC_px",1,200,"");
  hMC->ProjectionY("hMC_py",1,200,"");
  hMC_px->Rebin(rbin);
  hMC_py->Rebin(rbin);

  TH1D* sg1 = new TH1D("sg1","", nbin, 0., 400.);
  double theN = 0. ;
  for (int i=0; i< nbin; i++) {
      if ( cname == "Had" ) theN = hMC_py->GetBinContent(i);
      if ( cname == "Lep" ) theN = hMC_px->GetBinContent(i);
      sg1->SetBinContent(i,theN);
  }

  TH1D* qm1 = new TH1D("qm1","", nbin, 0., 400.);
  getBackground( qm1, 3, rbin );
  qm1->SetFillColor(4);
  ttstk->Add( qm1 );
  ttadd->Add( qm1 );

  TH2D* wj2d = (TH2D *) file0a->Get("Tops/"+hname);
  wj2d->ProjectionX("wj1",1,200,"");
  wj1->Rebin(rbin);

  //TH1D* wj1 = new TH1D("wj1","", nbin, 0., 400.);
  //getBackground( wj1, 2, rbin );
  wj1->SetFillColor(2);
  ttstk->Add( wj1 );
  ttadd->Add( wj1 );

  TH1D* tt1 = new TH1D("tt1","", nbin, 0., 400.);
  tt1->Add(hdata_pj,sg1,1,-1);
  tt1->SetFillColor(7);
  ttstk->Add( tt1 );
  ttadd->Add( tt1 );

  sg1->SetFillColor(3);
  ttstk->Add( sg1 );
  ttadd->Add( sg1 ); 


  delete hMC_px;
  delete hMC_py;
  delete hMC;
  delete hdata_pj;
  delete hdata;
  delete wj2d;
  
}

void TemplateMassFit::getFakeData( TH1D* ttadd,  int rbin ){


  int nbin = 200./rbin ;
  hdata  = (TH2F *) file0->Get("Tops/"+hname);
  hdata->ProjectionX("hdata_pj",1,200,"");
  hdata_pj->Rebin(rbin);

  hMC   = (TH2D *) file0->Get("MObjs/"+sname);
  hMC->ProjectionX("hMC_px",1,200,"");
  hMC->ProjectionY("hMC_py",1,200,"");
  hMC_px->Rebin(rbin);
  hMC_py->Rebin(rbin);

  TH1D* sg2 = new TH1D("sg2","", nbin, 0., 400.);
  double theN = 0. ;
  for (int i=0; i< nbin; i++) {
      if ( cname == "Had" ) theN = hMC_py->GetBinContent(i);
      if ( cname == "Lep" ) theN = hMC_px->GetBinContent(i);
      sg2->SetBinContent(i,theN);
  }

  qm2 = new TH1D("qm2","", nbin, 0., 400.);
  getBackground( qm2, 3, rbin );
  ttadd->Add( qm2 );

  wj2 = new TH1D("wj2","", nbin, 0., 400.);
  getBackground( wj2, 2, rbin );
  ttadd->Add( wj2 );

  tt2 = new TH1D("tt2","", nbin, 0., 400.);
  tt2->Add(hdata_pj,sg2,1,-1);
  ttadd->Add( tt2 );

  ttadd->Add( sg2 ); 
  
  delete sg2;
  delete tt2;
  delete wj2;
  delete qm2;
  delete hMC_px;
  delete hMC_py;
  delete hMC;
  delete hdata_pj;
  delete hdata;

}
void TemplateMassFit::getData( TH1D* h_data, int rbin ){

  int nbin = 200./rbin ;
  hdata  = (TH2F *) file0->Get("Tops/"+hname);

  hdata->ProjectionX("hdata_pj",1,200,"");
  hdata_pj->Rebin(rbin);

  double theN = 0. ;
  for (int i=0; i< nbin; i++) {
      theN  = hdata_pj->GetBinContent(i);
      h_data->SetBinContent(i,theN);
  }
  delete hdata_pj;
  delete hdata;
}

void TemplateMassFit::getSignal(TH1D* h_Sg, int rbin ) {

  int nbin = 200./rbin ;
  ttMC   = (TH2D *) file1->Get("MObjs/"+sname);
  
  ttMC->ProjectionX("ttMC_px",1,200,"");
  ttMC->ProjectionY("ttMC_py",1,200,"");
  ttMC_px->Rebin(rbin);
  ttMC_py->Rebin(rbin);

  double theN = 0. ;
  for (int i=0; i< nbin; i++) {
      if ( cname == "Had" ) theN = ttMC_py->GetBinContent(i);
      if ( cname == "Lep" ) theN = ttMC_px->GetBinContent(i);
      h_Sg->SetBinContent(i,theN);
  }
  delete ttMC_px;
  delete ttMC_py;
  delete ttMC;
}

void TemplateMassFit::getBackground(TH1D* h_Bg, int type, int rbin ){

  if (type == 0) ttmass = (TH2D *) file1->Get("Tops/"+hname); 
  if (type == 1) ttmass = (TH2D *) file1->Get("Tops/"+hname); 
  if (type == 2) ttmass = (TH2D *) file2->Get("Tops/"+hname); 
  if (type == 3) ttmass = (TH2D *) file3->Get("Tops/"+hname);

  int nbin = 200./rbin ;

  ttmass->ProjectionX("ttmass_pj",1,200,"");
  ttmass_pj->Rebin(rbin);

  if( type == 1 ) {
    TH1D* sgl = new TH1D("sgl","", nbin, 0., 400.);
    getSignal( sgl, rbin );
    ttmass_pj->Add(sgl, -1 );
    delete sgl;
  }

  for (int i=0; i< nbin; i++) {
      double theN = ttmass_pj->GetBinContent(i);
      h_Bg->SetBinContent(i,theN);
  }

  delete ttmass_pj;
  delete ttmass;
}

Double_t TemplateMassFit::fitG( Double_t* x, Double_t* par){

     Double_t A1 = (par[0]/x[0]) + (par[1]*x[0]) + par[2];
     Double_t A2 =  exp( A1 ) ;
     Double_t fitV = A2 ;

     return fitV;
}

Double_t TemplateMassFit::fitGS(Double_t *x, Double_t *par) {

     Double_t chi = (x[0] - par[1]) / par[2] ;
     Double_t A1 =  exp(-1.*chi*chi2  );
     Double_t fitV = par[0]*A1 ;
     return fitV;

}

Double_t TemplateMassFit::fitBG(Double_t *x, Double_t *par) {

        return fitG(x,par) + fitTail(x, &par[3] );

}

Double_t TemplateMassFit::fitBW(Double_t *x, Double_t *par) {

     Double_t gm = par[0] + 1.;
     Double_t chi = (x[0] - par[1]) / gm ;
     Double_t C1 = 1. + (chi*chi) ;
     Double_t fV = par[2]/C1 ;
     return fV;
}

Double_t TemplateMassFit::fitData(Double_t *x, Double_t *par) {

        return fitBW(x,par) + fitG(x, &par[3]) ;

}


