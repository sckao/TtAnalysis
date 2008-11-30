#include <vector>
#include <stdio.h>
#include <TFile.h>
#include <iostream>
#include <fstream>

// define the fitting function
Double_t fitf(Double_t *x, Double_t *par) {

         Double_t A1 = (x[0])/(par[1]+1.) ;
	 Double_t A2 = pow(x[0], par[2]) ;
         Double_t fitval =  par[0] / ( exp(A1) * A2 );
         return fitval;
}

Double_t fitG(Double_t *x, Double_t *par) {

	Double_t chi = (x[0]-par[1])/par[2] ;
	Double_t A3  = -0.5*chi*chi ;
	Double_t fitgau = par[0]*exp( A3 ) ;
	return fitgau ;
}

Double_t fitCb(Double_t *x, Double_t *par) {

        return fitf(x,par) + fitG(x, &par[3] );
}

void FCRatio() {

 TFile *file = TFile::Open("WJ_11pb_Et20.root");
 TString hfolder = "wjets_test";
 TFile *file1 = TFile::Open("tt_11pb_Et20.root");
 //TString hfolder = "tt_test";

 TString name1 = "eta_njets";
 TString name2 = "wEta_njets";
 TString name3 = "muNjets_dPhi0";
 TString name4 = "muNjets_dPhi1";
 TString name5 = "muNjets_dPhi2";
 TString name6 = "njets_dPhi0";
 TString name7 = "njets_dPhi1";
 TString name8 = "njets_dPhi2";
 //TString name9 = "thirdJetEt";
 TString name9 = "thirdCalEt";
 TString name10 = "m3_dR";

 TString plot1 = "njet_fcRatio.gif";
 TString plot2 = "eta_njets.gif";
 TString plot3 = "wEta_njets.gif";
 TString plot4 = "njet_WfcRatio.gif";
 TString plot5 = "dPhi_W_3j.gif";
 TString plot6 = "dPhi_W_1j.gif";
 TString plot7 = "dPhi_Mu_3j.gif";
 TString plot8 = "dPhi_Mu_1j.gif";
 TString plot9 = "3rdJetEt.gif";
 TString plot10 = "m3.gif";


 heta_njets  = (TH2F *) file->Get("Jets/"+name1); 
 wEta_njets  = (TH2F *) file->Get("Jets/"+name2); 
 
 muNjets_dPhi0 = (TH2F *) file->Get("Jets/"+name3); 
 muNjets_dPhi1 = (TH2F *) file->Get("Jets/"+name4); 
 muNjets_dPhi2 = (TH2F *) file->Get("Jets/"+name5); 
 njets_dPhi0  = (TH2F *) file->Get("Jets/"+name6); 
 njets_dPhi1  = (TH2F *) file->Get("Jets/"+name7); 
 njets_dPhi2  = (TH2F *) file->Get("Jets/"+name8); 
 thirdJetEt   = (TH1F *) file->Get("Jets/"+name9); 
 thirdJetEt1  = (TH1F *) file1->Get("Jets/"+name9); 
 m3_dR        = (TH2F *) file->Get("Jets/"+name10); 
 m3_dR1       = (TH2F *) file1->Get("Jets/"+name10); 


 gSystem->mkdir(hfolder);
 gSystem->cd(hfolder);

 // eta of iso mu and NJets 
 heta_njets->ProjectionY("hnjets_C",34,58,"");

 heta_njets->ProjectionY("heta_njets_pyf",1,33,"");
 heta_njets->ProjectionY("heta_njets_pyb",59,91,"");
 TH1F *hnjets_F = new TH1F("hnjets_F","",25 , -0.5, 24.5);

 hnjets_F->Add(heta_njets_pyf, heta_njets_pyb, 1., 1.);

 // eta of lep W and NJets
 wEta_njets->ProjectionY("wnjets_C",34,58,"");
 wEta_njets->ProjectionY("weta_njets_pyf",1,33,"");
 wEta_njets->ProjectionY("weta_njets_pyb",59,91,"");
 TH1F *wnjets_F = new TH1F("wnjets_F","",25 , -0.5, 24.5);
 wnjets_F->Add(weta_njets_pyf, weta_njets_pyb, 1., 1.);
 

 Float_t x1[8] = {0.0};
 Float_t ex1[8] = {0.0};
 Float_t y1[8] = {0.0};
 Float_t ey1[8] = {0.0};
 Float_t y2[8] = {0.0};
 Float_t ey2[8] = {0.0};

 for (int i=0; i<8; i++) {

     float nF = hnjets_F->GetBinContent(i+1);
     float nC = hnjets_C->GetBinContent(i+1);
     float nJ = nF + nC;

     x1[i]= i;
     if ( nJ != 0 ) {
        y1[i]= nC/nJ ;
     } else {
        y1[i]= 0. ;
     }

     if ( y1[i] < 0.999 && nJ != 0) { ey1[i] = sqrt( y1[i]*(1-y1[i])/nJ ); }
     else if ( nJ != 0 ) {  ey1[i] = 1.0/nJ ; }
     else { ey1[i] = 0. ; }

     float nFw = wnjets_F->GetBinContent(i+1);
     float nCw = wnjets_C->GetBinContent(i+1);
     float nJw = nFw + nCw;

     if ( nJw != 0 ) {
        y2[i]= nCw/nJw ;
     } else {
        y2[i]= 0. ;
     }

     if ( y2[i] < 0.999 && nJw != 0) { ey2[i] = sqrt( y2[i]*(1-y2[i])/nJw ); }
     else if ( nJw != 0 ) {  ey2[i] = 1.0/nJw ; }
     else { ey2[i] = 0. ; }
 }

 gStyle->SetOptStat(kFALSE);
 gStyle->SetOptFit(0111);  
 c1 = new TCanvas("c1","",200,8,800,600);
 c1->SetFillColor(10);
 c1->SetGrid(); 
 c1->GetFrame()->SetFillColor(21);
 c1->GetFrame()->SetBorderSize(12);
 c1->cd();

 hFCRatio = new TGraphErrors(8,x1,y1,ex1,ey1);
 hFCRatio->SetMaximum(1.02);
 hFCRatio->SetMinimum(0.0);
 hFCRatio->SetMarkerColor(4);
 hFCRatio->SetMarkerStyle(21);
 hFCRatio->SetTitle(" Central Event Ratio ");
 hFCRatio->GetXaxis()->SetTitle(" N of jets  ");
 hFCRatio->GetYaxis()->SetTitle(" Ratio  ");
 hFCRatio->Draw("ACP");

 c1->Update();
 c1->Print(plot1);

 //
 gStyle->SetOptStat(kFALSE);
 gStyle->SetOptFit(0111);  
 c4 = new TCanvas("c4","",200,8,800,600);
 c4->SetFillColor(10);
 c4->SetGrid(); 
 c4->GetFrame()->SetFillColor(21);
 c4->GetFrame()->SetBorderSize(12);
 c4->cd();

 wFCRatio = new TGraphErrors(8,x1,y2,ex1,ey2);
 wFCRatio->SetMaximum(1.02);
 wFCRatio->SetMinimum(0.0);
 wFCRatio->SetMarkerColor(4);
 wFCRatio->SetMarkerStyle(21);
 wFCRatio->SetTitle(" Central Event Ratio in term of W eta");
 wFCRatio->GetXaxis()->SetTitle(" N of jets  ");
 wFCRatio->GetYaxis()->SetTitle(" Ratio  ");
 wFCRatio->Draw("ACP");

 c4->Update();
 c4->Print(plot4);


 heta_njets->ProjectionX("heta_0jets_px",1,1,"");
 heta_njets->ProjectionX("heta_1jets_px",2,2,"");
 heta_njets->ProjectionX("heta_2jets_px",3,3,"");
 heta_njets->ProjectionX("heta_3jets_px",4,4,"");
 heta_njets->ProjectionX("heta_4jets_px",5,5,"");
 heta_njets->ProjectionX("heta_5jets_px",6,6,"");
 heta_njets->ProjectionX("heta_6jets_px",7,7,"");

 gStyle->SetOptStat("nimou");
 TCanvas *c2 = new TCanvas("c2","");
 c2->SetFillColor(10);
 c2->SetFillColor(10);
 c2->SetLogy();

 heta_0jets_px->SetTitle(" #eta of isoMuons ");
 heta_0jets_px->SetAxisRange(-4.55,4.55,"X");
 heta_0jets_px->SetAxisRange(1,3000,"Y");
 heta_0jets_px->DrawCopy();
 heta_0jets_px->GetXaxis()->SetTitle(" #eta  ");

 heta_1jets_px->SetLineColor(kRed);
 heta_1jets_px->DrawCopy("same");
 heta_2jets_px->SetLineColor(kBlue);
 heta_2jets_px->DrawCopy("same");
 heta_3jets_px->SetLineColor(kGreen);
 heta_3jets_px->DrawCopy("same");
 heta_4jets_px->SetLineColor(6);
 heta_4jets_px->DrawCopy("same");
 heta_5jets_px->SetLineColor(7);
 heta_5jets_px->DrawCopy("same");
 heta_6jets_px->SetLineColor(5);
 heta_6jets_px->DrawCopy("same");


 c2->Update();
 c2->Print(plot2);

 //***

 wEta_njets->ProjectionX("wEta_0jets_px",1,1,"");
 wEta_njets->ProjectionX("wEta_1jets_px",2,2,"");
 wEta_njets->ProjectionX("wEta_2jets_px",3,3,"");
 wEta_njets->ProjectionX("wEta_3jets_px",4,4,"");
 wEta_njets->ProjectionX("wEta_4jets_px",5,5,"");
 wEta_njets->ProjectionX("wEta_5jets_px",6,6,"");
 wEta_njets->ProjectionX("wEta_6jets_px",7,7,"");

 gStyle->SetOptStat("nimou");
 TCanvas *c3 = new TCanvas("c3","");
 c3->SetFillColor(10);
 c3->SetFillColor(10);
 c3->SetLogy();

 wEta_0jets_px->SetTitle(" Y of W ");
 wEta_0jets_px->SetAxisRange(-4.55,4.55,"X");
 wEta_0jets_px->SetAxisRange(1,3000,"Y");
 wEta_0jets_px->DrawCopy();
 wEta_0jets_px->GetXaxis()->SetTitle(" #eta  ");

 wEta_1jets_px->SetLineColor(kRed);
 wEta_1jets_px->DrawCopy("same");
 wEta_2jets_px->SetLineColor(kBlue);
 wEta_2jets_px->DrawCopy("same");
 wEta_3jets_px->SetLineColor(kGreen);
 wEta_3jets_px->DrawCopy("same");
 wEta_4jets_px->SetLineColor(6);
 wEta_4jets_px->DrawCopy("same");
 wEta_5jets_px->SetLineColor(7);
 wEta_5jets_px->DrawCopy("same");
 wEta_6jets_px->SetLineColor(5);
 wEta_6jets_px->DrawCopy("same");

 c3->Update();
 c3->Print(plot3);

 //***

 njets_dPhi2->ProjectionY("dphi2_0jets_py",1,1,"");
 njets_dPhi2->ProjectionY("dphi2_1jets_py",2,2,"");
 njets_dPhi2->ProjectionY("dphi2_2jets_py",3,3,"");
 njets_dPhi2->ProjectionY("dphi2_3jets_py",4,4,"");
 njets_dPhi2->ProjectionY("dphi2_4jets_py",5,5,"");
 njets_dPhi2->ProjectionY("dphi2_5jets_py",6,6,"");

 gStyle->SetOptStat("nimou");
 TCanvas *c5 = new TCanvas("c5","");
 c5->SetFillColor(10);
 c5->SetFillColor(10);
 c5->SetLogy();

 dphi2_0jets_py->SetTitle(" #Delta#phi of W & 3rdJet");
 dphi2_0jets_py->SetAxisRange(-0.5,3.15,"X");
 dphi2_0jets_py->SetAxisRange(1,1000,"Y");
 dphi2_0jets_py->DrawCopy();
 dphi2_0jets_py->GetXaxis()->SetTitle(" #Delta#phi  ");

 dphi2_1jets_py->SetLineColor(kRed);
 dphi2_1jets_py->DrawCopy("same");
 dphi2_2jets_py->SetLineColor(kBlue);
 dphi2_2jets_py->DrawCopy("same");
 dphi2_3jets_py->SetLineColor(kGreen);
 dphi2_3jets_py->DrawCopy("same");
 dphi2_4jets_py->SetLineColor(6);
 dphi2_4jets_py->DrawCopy("same");
 dphi2_5jets_py->SetLineColor(7);
 dphi2_5jets_py->DrawCopy("same");

 c5->Update();
 c5->Print(plot5);
 
 //***

 njets_dPhi0->ProjectionY("dphi0_0jets_py",1,1,"");
 njets_dPhi0->ProjectionY("dphi0_1jets_py",2,2,"");
 njets_dPhi0->ProjectionY("dphi0_2jets_py",3,3,"");
 njets_dPhi0->ProjectionY("dphi0_3jets_py",4,4,"");
 njets_dPhi0->ProjectionY("dphi0_4jets_py",5,5,"");
 njets_dPhi0->ProjectionY("dphi0_5jets_py",6,6,"");

 gStyle->SetOptStat("nimou");
 TCanvas *c6 = new TCanvas("c6","");
 c6->SetFillColor(10);
 c6->SetFillColor(10);
 c6->SetLogy();

 dphi0_0jets_py->SetTitle(" #Delta#phi of W & leading Jet");
 dphi0_0jets_py->SetAxisRange(-0.5,3.15,"X");
 dphi0_0jets_py->SetAxisRange(1,1000,"Y");
 dphi0_0jets_py->DrawCopy();
 dphi0_0jets_py->GetXaxis()->SetTitle(" #Delta#phi  ");

 dphi0_1jets_py->SetLineColor(kRed);
 dphi0_1jets_py->DrawCopy("same");
 dphi0_2jets_py->SetLineColor(kBlue);
 dphi0_2jets_py->DrawCopy("same");
 dphi0_3jets_py->SetLineColor(kGreen);
 dphi0_3jets_py->DrawCopy("same");
 dphi0_4jets_py->SetLineColor(6);
 dphi0_4jets_py->DrawCopy("same");
 dphi0_5jets_py->SetLineColor(7);
 dphi0_5jets_py->DrawCopy("same");

 c6->Update();
 c6->Print(plot6);
 
 //***

 muNjets_dPhi2->ProjectionY("mdphi2_0jets_py",1,1,"");
 muNjets_dPhi2->ProjectionY("mdphi2_1jets_py",2,2,"");
 muNjets_dPhi2->ProjectionY("mdphi2_2jets_py",3,3,"");
 muNjets_dPhi2->ProjectionY("mdphi2_3jets_py",4,4,"");
 muNjets_dPhi2->ProjectionY("mdphi2_4jets_py",5,5,"");
 muNjets_dPhi2->ProjectionY("mdphi2_5jets_py",6,6,"");

 gStyle->SetOptStat("nimou");
 TCanvas *c7 = new TCanvas("c7","");
 c7->SetFillColor(10);
 c7->SetFillColor(10);
 c7->SetLogy();

 mdphi2_0jets_py->SetTitle(" #Delta#phi of isoMu & 3rdJet");
 mdphi2_0jets_py->SetAxisRange(-0.5,3.15,"X");
 mdphi2_0jets_py->SetAxisRange(1,1000,"Y");
 mdphi2_0jets_py->DrawCopy();
 mdphi2_0jets_py->GetXaxis()->SetTitle(" #Delta#phi  ");

 mdphi2_1jets_py->SetLineColor(kRed);
 mdphi2_1jets_py->DrawCopy("same");
 mdphi2_2jets_py->SetLineColor(kBlue);
 mdphi2_2jets_py->DrawCopy("same");
 mdphi2_3jets_py->SetLineColor(kGreen);
 mdphi2_3jets_py->DrawCopy("same");
 mdphi2_4jets_py->SetLineColor(6);
 mdphi2_4jets_py->DrawCopy("same");
 mdphi2_5jets_py->SetLineColor(7);
 mdphi2_5jets_py->DrawCopy("same");

 c7->Update();
 c7->Print(plot7);
 
 //***

 muNjets_dPhi0->ProjectionY("mdphi0_0jets_py",1,1,"");
 muNjets_dPhi0->ProjectionY("mdphi0_1jets_py",2,2,"");
 muNjets_dPhi0->ProjectionY("mdphi0_2jets_py",3,3,"");
 muNjets_dPhi0->ProjectionY("mdphi0_3jets_py",4,4,"");
 muNjets_dPhi0->ProjectionY("mdphi0_4jets_py",5,5,"");
 muNjets_dPhi0->ProjectionY("mdphi0_5jets_py",6,6,"");

 gStyle->SetOptStat("nimou");
 TCanvas *c8 = new TCanvas("c8","");
 c8->SetFillColor(10);
 c8->SetFillColor(10);
 c8->SetLogy();

 mdphi0_0jets_py->SetTitle(" #Delta#phi of isoMu & leading Jet");
 mdphi0_0jets_py->SetAxisRange(-0.5,3.15,"X");
 mdphi0_0jets_py->SetAxisRange(1,1000,"Y");
 mdphi0_0jets_py->DrawCopy();
 mdphi0_0jets_py->GetXaxis()->SetTitle(" #Delta#phi  ");

 mdphi0_1jets_py->SetLineColor(kRed);
 mdphi0_1jets_py->DrawCopy("same");
 mdphi0_2jets_py->SetLineColor(kBlue);
 mdphi0_2jets_py->DrawCopy("same");
 mdphi0_3jets_py->SetLineColor(kGreen);
 mdphi0_3jets_py->DrawCopy("same");
 mdphi0_4jets_py->SetLineColor(6);
 mdphi0_4jets_py->DrawCopy("same");
 mdphi0_5jets_py->SetLineColor(7);
 mdphi0_5jets_py->DrawCopy("same");

 c8->Update();
 c8->Print(plot8);
 
 gStyle->SetOptStat("nimou");
 gStyle->SetOptFit(0111);  
 TCanvas *c9 = new TCanvas("c9","");
 c9->SetFillColor(10);
 c9->SetFillColor(10);
 c9->SetLogy();

 thirdJetEt->SetTitle(" 3rd Jet Et ");
 thirdJetEt->SetAxisRange(0.0, 100, "X");
 thirdJetEt->SetAxisRange(0.1, 5000, "Y");
 
 TF1 *func0 = new TF1("fitf",fitf,8,48,3);
 func0->SetParLimits(1,0.1,1000000.);
 func0->SetParLimits(2,0.1,1000000.);
 thirdJetEt->Fit("fitf","R","",8,48);
 thirdJetEt->Draw(); 

 float nSoup[100] ={0.0} ;
 float et[100] = {0.0};
 for (i=0; i<100; i++) {
    float nWj = thirdJetEt->GetBinContent(i+1);
    float nTt = thirdJetEt1->GetBinContent(i+1);
    nSoup[i] = nWj+nTt ;
    et[i]= i;
 }	 
 
 SoupEt = new TGraph(100, et, nSoup);
 SoupEt->SetMarkerSize(0.5);
 SoupEt->SetMarkerColor(4);
 SoupEt->SetMarkerStyle(21);

 func0->SetLineColor(2);
 SoupEt->Fit("fitf","R","",8,15);
 func0->Draw("SAME");
 
 TF1 *func1 = new TF1("fitc",fitCb,8,60,6);
 func1->SetParLimits(1,0.1,1000000.);
 func1->SetParLimits(2,1.1,10.);
 func1->SetParLimits(3,21.,38.);
 func1->SetParLimits(4,41.,58.);
 func1->SetParLimits(5,10.,999.);
 SoupEt->Draw("P");

 //func1->SetLineColor(3);
 //SoupEt->Fit("fitc","R","",8,60);
 //func1->Draw("SAME");

 c9->Update();
 c9->Print(plot9);


 gStyle->SetOptStat("nimou");
 gStyle->SetOptFit(0111);  
 TCanvas *c10 = new TCanvas("c10","");
 c10->SetFillColor(10);
 c10->SetFillColor(10);
 c10->SetLogy();

 TF1 *func2 = new TF1("fitm",fitf,50,450,3);
 func2->SetParLimits(1,0.1,1000000.);
 func2->SetParameter(2,0.);

 m3_dR->ProjectionX("m3_wj0");
 m3_wj0->Rebin(5,"m3_wj");
 m3_dR1->ProjectionX("m3_tt0");
 m3_tt0->Rebin(5,"m3_tt");
 m3_wj->SetAxisRange(0.5, 500, "Y");
 m3_wj->Fit("fitm","R","",50,450);
 m3_wj->Draw(); 

 double mp0 = func2->GetParameter(0);
 double mp1 = func2->GetParameter(1);
 double mp2 = func2->GetParameter(2);

 float mSoup[100] ={0.0} ;
 float m3[100] = {0.0};
 for (i=0; i<100; i++) {
    float nWj = m3_wj->GetBinContent(i+1);
    float nTt = m3_tt->GetBinContent(i+1);
    mSoup[i] = nWj+nTt ;
    m3[i]= (i*5.) + 2.5;
 }	 
 
 SoupM3 = new TGraph(100, m3, mSoup);
 SoupM3->SetMarkerSize(0.5);
 SoupM3->SetMarkerColor(4);
 SoupM3->SetMarkerStyle(21);

 func2->SetLineColor(3);
 func2->SetParameter(0,mp0);
 func2->SetParameter(1,mp1);
 func2->SetParameter(2,mp2);
 SoupM3->Fit("fitm","R","",50,120);
 func2->Draw("SAME");
 SoupM3->Draw("P");
 
 c10->Update();
 c10->Print(plot10);

 gSystem->cd("../");
 file1->Close();
 file->Close();

}
