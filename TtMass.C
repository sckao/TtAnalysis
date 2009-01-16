void TtMass( TString name1 ) {

 TFile *file  = TFile::Open("ttj_1Jskim_NoBTag.root");
 TFile *file1 = TFile::Open("wjets_1Jskim_NoBTag.root");
 TFile *file2 = TFile::Open("qcd_1Jskim_NoBTag.root");
 //TFile *file  = TFile::Open("ttj_1Jskim.root");
 //TFile *file1 = TFile::Open("wjets_1Jskim.root");
 //TFile *file2 = TFile::Open("qcd_1Jskim.root");
 TString hfolder = "tt_test";

 //TString name1 = "hRecott0";

 TString plot1 = name1+".gif";
 TString plot2 = "TtMassFit.gif";

 ttmass   = (TH2F *) file->Get("Tops/"+name1); 
 ttmass1  = (TH2F *) file1->Get("Tops/"+name1); 
 ttmass2  = (TH2F *) file2->Get("Tops/"+name1);

 gSystem->mkdir(hfolder);
 gSystem->cd(hfolder);

 TCanvas *c1 = new TCanvas("c1","", 1000, 800);
 c1->SetFillColor(10);
 c1->SetFillColor(10);
 c1->Divide(2,2);

 THStack *ttstkx = new THStack("ttstk", " Leptonic Top (projectX) ");
 THStack *ttstky = new THStack("ttstk", " Hadronic Top (projectY) ");
 TH1F * ttaddx = new TH1F("ttaddx","", 100, 0., 400.);
 TH1F * ttaddy = new TH1F("ttaddy","", 100, 0., 400.);

 ttmass1->ProjectionY("ttmass1_py",1,200,"");
 ttmass1->ProjectionX("ttmass1_px",1,200,"");
 ttmass1_px->Rebin(2);
 ttmass1_py->Rebin(2);
 ttmass1_px->SetFillColor(2);
 ttmass1_py->SetFillColor(2);
 //ttmass1->Scale(1);
 ttstkx->Add(ttmass1_px);
 ttstky->Add(ttmass1_py);
 ttaddx->Add(ttmass1_px);
 ttaddy->Add(ttmass1_py);

 ttmass2->ProjectionY("ttmass2_py",1,200,"");
 ttmass2->ProjectionX("ttmass2_px",1,200,"");
 ttmass2_px->Rebin(2);
 ttmass2_py->Rebin(2);
 ttmass2_px->SetFillColor(4);
 ttmass2_py->SetFillColor(4);
 //ttmass2->Scale(1);
 ttstkx->Add(ttmass2_px);
 ttstky->Add(ttmass2_py);
 ttaddx->Add(ttmass2_px);
 ttaddy->Add(ttmass2_py);

 ttmass->ProjectionY("ttmass_py",1,200,"");
 ttmass->ProjectionX("ttmass_px",1,200,"");
 ttmass_px->Rebin(2);
 ttmass_py->Rebin(2);
 ttmass_px->SetFillColor(3);
 ttmass_py->SetFillColor(3);
 //ttmass_py->Rebin(2);
 ttstkx->Add(ttmass_px);
 ttstky->Add(ttmass_py);
 ttaddx->Add(ttmass_px);
 ttaddy->Add(ttmass_py);

 c1->cd(1);
 gStyle->SetStatY(0.95); 
 gStyle->SetStatTextColor(3);
 ttmass->RebinX(4);
 ttmass->RebinY(4);
 ttmass->SetFillColor(3);
 ttmass->Draw("BOX");
 c1->Update();
 gStyle->SetStatY(0.7); 
 gStyle->SetStatTextColor(2);
 ttmass1->RebinX(4);
 ttmass1->RebinY(4);
 ttmass1->SetLineColor(2);
 ttmass1->Draw("BOXSAMES");
 c1->Update();
 gStyle->SetStatY(0.45); 
 gStyle->SetStatTextColor(4);
 ttmass2->RebinX(4);
 ttmass2->RebinY(4);
 ttmass2->SetLineColor(4);
 ttmass2->Draw("BOXSAMES");
 c1->Update();

 gStyle->SetStatY(0.95); 
 gStyle->SetStatTextColor(1);

 c1->cd(2);

 gStyle->SetOptStat(kTRUE);
 gStyle->SetOptFit(111);

 ttstky->Draw();

 ttaddy->SetMarkerStyle(21);
 ttaddy->SetMarkerColor(1);
 ttaddy->SetMarkerSize(0.7);
 ttaddy->SetLineColor(1);
 ttaddy->SetTitle(" Fitting of Top Mass ");
 ttaddy->Fit("gaus","N0R","",120,220);
 gaus->SetLineColor(1);
 ttaddy->Fit("gaus","R","sames",100,250);

 c1->cd(3);

 gStyle->SetOptStat(kTRUE);
 gStyle->SetOptFit(111);

 ttstkx->Draw();

 ttaddx->SetMarkerStyle(21);
 ttaddx->SetMarkerColor(1);
 ttaddx->SetMarkerSize(0.7);
 ttaddx->SetLineColor(1);
 ttaddx->SetTitle(" Fitting of Top Mass ");
 ttaddx->Fit("gaus","N0R","",120,220);
 gaus->SetLineColor(1);
 ttaddx->Fit("gaus","R","sames",100,250);
 //ttaddx->DrawCopy("samesP");


 //c1->cd(4);

 /*

 c1->Update();
 gStyle->SetStatY(0.99); 
 gStyle->SetStatTextColor(1);

 gStyle->SetOptStat(kTRUE);
 gStyle->SetOptFit(111);
 ttaddx->SetMarkerStyle(21);
 ttaddx->SetMarkerColor(1);
 ttaddx->SetMarkerSize(0.7);
 ttaddx->SetLineColor(1);
 ttaddx->SetTitle(" Fitting of Top Mass ");
 ttaddx->Fit("gaus","N0R","",120,220);
 gaus->SetLineColor(1);
 ttaddx->Fit("gaus","R","",120,220);
 ttaddx->DrawCopy("samesP");

 c1->Update();
 gStyle->SetStatY(0.62); 
 gStyle->SetStatTextColor(14);

 gStyle->SetOptStat(kTRUE);
 ttaddy->SetMarkerStyle(22);
 ttaddy->SetMarkerColor(14);
 ttaddy->SetMarkerSize(0.7);
 ttaddy->SetLineColor(14);
 ttaddy->SetTitle(" Reco Hadronic Top Mass ");
 ttaddy->Fit("gaus","N0R","",120,220);
 gaus->SetLineColor(14);
 ttaddy->Fit("gaus","R","sames",120,220);

 ttaddy->DrawCopy("sameP");

 */

 c1->Update();
 c1->Print(plot1);

 gStyle->SetStatY(0.95); 
 gStyle->SetStatTextColor(1);

 gSystem->cd("../");

 file2->Close();
 file1->Close();
 file->Close();

}
