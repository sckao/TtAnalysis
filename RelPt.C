void RelPt( TString name ) {

 // name = "h" or "l"

 TFile *file  = TFile::Open("ttj_fall08_3j30.root");
 TFile *file1 = TFile::Open("wjets_fall08_3j30.root");
 TFile *file2 = TFile::Open("qcd_fall08_3j30.root");
 TString hfolder = "tt_fall08";

 TString name1 = name+"RelPtMC" ;
 TString name2 = name+"RelPtRC" ;
 TString name3 = name+"RelPtRC" ;
 TString name4 = name+"RelPtRC" ;

 TString plot1 = name+"RelPt_pfx.gif";
 TString plot2 = name+"RelPt_scl.gif";
 TString plot3 = name+"RelPt_py.gif";

 thMCPt     = (TH2F *) file->Get("Tops/"+name1); 
 thRelPt    = (TH2F *) file->Get("Tops/"+name2); 
 whRelPt    = (TH2F *) file1->Get("Tops/"+name3); 
 qhRelPt    = (TH2F *) file2->Get("Tops/"+name4); 


 // re-bin value for top mass plot
 int rbin = 8;

 thMCPt->RebinY(rbin);
 thRelPt->RebinY(rbin);
 whRelPt->RebinY(rbin);
 qhRelPt->RebinY(rbin);

 gSystem->mkdir(hfolder);
 gSystem->cd(hfolder);

 double bt[27]={0.0};
 double be[27]={0.0};
 double mMC[27]={0.0};
 double eMC[27]={0.0};
 double mTh[27]={0.0};
 double eTh[27]={0.0};
 double mWh[27]={0.0};
 double eWh[27]={0.0};
 double mQh[27]={0.0};
 double eQh[27]={0.0};

 for (int i=2; i< 29; i++){

     int j = i-2 ;
     bt[j]= (j*0.04) ;

     thMCPt->ProjectionY("hMC_py",i,i);
     double nMC = hMC_py->Integral();
     hMC_py->Fit("gaus","N0RQ","",10.,90.);
     mMC[j] = gaus->GetParameter(1);
     eMC[j] = gaus->GetParameter(2);
     if ( nMC < 8 || mMC[j] < 0. ) {
        mMC[j] =hMC_py->GetMean(1);
        eMC[j] =hMC_py->GetRMS(1);
     }

     thRelPt->ProjectionY("th_py",i,i,"");
     double nTh = th_py->Integral();
     th_py->Fit("gaus","N0RQ","",10.,80.);
     mTh[j] = gaus->GetParameter(1);
     eTh[j] = gaus->GetParameter(2);
     if ( nTh < 8 || mTh[j] < 0. ) {
        mTh[j] =th_py->GetMean(1);
        eTh[j] =th_py->GetRMS(1);
     }

     whRelPt->ProjectionY("wh_py",i,i,"");
     double nWh = wh_py->Integral();
     wh_py->Fit("gaus","N0RQ","",10.,60.);
     mWh[j] = gaus->GetParameter(1);
     eWh[j] = gaus->GetParameter(2);
     if ( nWh < 8 || mWh[j] < 0. ) {
        mWh[j] =wh_py->GetMean(1);
        eWh[j] =wh_py->GetRMS(1);
     }

     qhRelPt->ProjectionY("qh_py",i,i,"");
     double nQh = qh_py->Integral();
     qh_py->Fit("gaus","N0RQ","",10.,60.);
     mQh[j] = gaus->GetParameter(1);
     eQh[j] = gaus->GetParameter(2);
     if ( nQh < 8 || mQh[j] < 0.) {
        mQh[j] =qh_py->GetMean(1);
        eQh[j] =qh_py->GetRMS(1);
     }

 }

 TCanvas *c1 = new TCanvas("c1","", 800, 700);
 c1->SetFillColor(10);
 c1->SetFillColor(10);
 c1->Divide(2,2);

 c1->cd(1);
 c1_1->SetGridy();

 mcRelPt_pfx = new TGraphErrors(27,bt,mMC,be,eMC);
 mcRelPt_pfx->SetTitle(" MC RelPt of hadW  w.r.t. Top direction ");
 mcRelPt_pfx->SetMaximum(90);
 mcRelPt_pfx->SetMinimum(10);
 mcRelPt_pfx->SetMarkerColor(2);
 mcRelPt_pfx->SetMarkerStyle(21);
 mcRelPt_pfx->GetXaxis()->SetTitle(" #beta  ");
 mcRelPt_pfx->GetYaxis()->SetTitle(" RelPt  ");
 mcRelPt_pfx->Draw("AP");
 c1->Update();
 
 c1->cd(2);
 c1_2->SetGridy();

 tRelPt_pfx = new TGraphErrors(27,bt,mTh,be,eTh);
 tRelPt_pfx->SetTitle(" Tt RelPt of hadW  w.r.t. Top direction ");
 tRelPt_pfx->SetMaximum(90);
 tRelPt_pfx->SetMinimum(10);
 tRelPt_pfx->SetMarkerColor(4);
 tRelPt_pfx->SetMarkerStyle(20);
 tRelPt_pfx->GetXaxis()->SetTitle(" #beta  ");
 tRelPt_pfx->GetYaxis()->SetTitle(" RelPt  ");
 tRelPt_pfx->Draw("AP");
 c1->Update();
 
 c1->cd(3);
 c1_3->SetGridy();

 wRelPt_pfx = new TGraphErrors(27,bt,mWh,be,eWh);
 wRelPt_pfx->SetTitle(" Wj RelPt of hadW  w.r.t. Top direction ");
 wRelPt_pfx->SetMaximum(90);
 wRelPt_pfx->SetMinimum(10);
 wRelPt_pfx->SetMarkerColor(6);
 wRelPt_pfx->SetMarkerStyle(22);
 wRelPt_pfx->GetXaxis()->SetTitle(" #beta  ");
 wRelPt_pfx->GetYaxis()->SetTitle(" RelPt  ");
 wRelPt_pfx->Draw("AP");
 c1->Update();
 
 c1->cd(4);
 c1_4->SetGridy();

 qRelPt_pfx = new TGraphErrors(27,bt,mQh,be,eQh);
 qRelPt_pfx->SetTitle(" QCD RelPt of hadW  w.r.t. Top direction ");
 qRelPt_pfx->SetMaximum(90);
 qRelPt_pfx->SetMinimum(10);
 qRelPt_pfx->SetMarkerColor(8);
 qRelPt_pfx->SetMarkerStyle(23);
 qRelPt_pfx->GetXaxis()->SetTitle(" #beta  ");
 qRelPt_pfx->GetYaxis()->SetTitle(" RelPt  ");
 qRelPt_pfx->Draw("AP");
 c1->Update();

 
 c1->cd();

 c1->Update();
 c1->Print(plot1);


 TCanvas *c2 = new TCanvas("c2","", 800, 700);
 c2->SetFillColor(10);
 c2->SetFillColor(10);
 c2->Divide(2,2);

 c2->cd(1);
 c2_1->SetGridy();
 thMCPt->DrawCopy("BOX");
 c2->Update();
 c2->cd(2);
 c2_2->SetGridy();
 thRelPt->DrawCopy("BOX");
 c2->Update();
 c2->cd(3);
 c2_3->SetGridy();
 whRelPt->DrawCopy("BOX");
 c2->Update();
 c2->cd(4);
 c2_4->SetGridy();
 qhRelPt->DrawCopy("BOX");
 c2->Update();

 c2->Print(plot2);


 TCanvas *c3 = new TCanvas("c3","", 800, 700);
 c3->SetFillColor(10);
 c3->SetFillColor(10);
 c3->Divide(2,2);

 c3->cd(1);
 c3_1->SetGridx();
 thMCPt->ProjectionY("thMCPt_py",-1,-1,"");
 thMCPt_py->DrawCopy("");
 c3->Update();
 c3->cd(2);
 c3_2->SetGridx();
 thRelPt->ProjectionY("thRelPt_py",-1,-1,"");
 thRelPt_py->DrawCopy("");
 c3->Update();
 c3->cd(3);
 c3_3->SetGridx();
 whRelPt->ProjectionY("whRelPt_py",-1,-1,"");
 whRelPt_py->DrawCopy("");
 c3->Update();
 c3->cd(4);
 c3_4->SetGridx();
 qhRelPt->ProjectionY("qhRelPt_py",-1,-1,"");
 qhRelPt_py->DrawCopy("");
 c3->Update();

 c3->Print(plot3);

 gSystem->cd("../");

 file2->Close();
 file1->Close();
 file->Close();

}
