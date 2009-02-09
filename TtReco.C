void TtReco() {

 TFile *file  = TFile::Open("ttj_fall08_3j30.root");
 TFile *file1 = TFile::Open("wjets_fall08_3j30.root");
 TFile *file2 = TFile::Open("qcd_fall08_3j30.root");
 TString hfolder = "tt_fall08";

 //TString name1 = "hRecott0";

 TString plot1 = "recoLepT.gif";
 TString plot2 = "recoHadT.gif";

 tlepTa  = (TH1F *) file->Get("Tops/alepRCTMass"); 
 wlepTa  = (TH1F *) file1->Get("Tops/alepRCTMass"); 
 qlepTa  = (TH1F *) file2->Get("Tops/alepRCTMass"); 

 tlepTp  = (TH1F *) file->Get("Tops/lepRCTMass"); 
 wlepTp  = (TH1F *) file1->Get("Tops/lepRCTMass"); 
 qlepTp  = (TH1F *) file2->Get("Tops/lepRCTMass"); 

 thadTa  = (TH1F *) file->Get("Tops/ahadRCTMass"); 
 whadTa  = (TH1F *) file1->Get("Tops/ahadRCTMass"); 
 qhadTa  = (TH1F *) file2->Get("Tops/ahadRCTMass"); 

 thadTp  = (TH1F *) file->Get("Tops/hadRCTMass"); 
 whadTp  = (TH1F *) file1->Get("Tops/hadRCTMass"); 
 qhadTp  = (TH1F *) file2->Get("Tops/hadRCTMass"); 

 tSelT  = (TH2F *) file->Get("Tops/hRecott0"); 
 wSelT  = (TH2F *) file1->Get("Tops/hRecott0"); 
 qSelT  = (TH2F *) file2->Get("Tops/hRecott0"); 

 tMCT  = (TH2F *) file->Get("Tops/hMCtt"); 
 // re-bin value for top mass plot
 int rbin = 4;

 gSystem->mkdir(hfolder);
 gSystem->cd(hfolder);

 gStyle->SetOptStat(kTRUE);
 gStyle->SetOptStat(1000101);

 TCanvas *c1 = new TCanvas("c1","", 1000, 800);
 c1->SetFillColor(10);
 c1->SetFillColor(10);
 c1->Divide(2,2);

 c1->cd(1);
 c1_1->SetLogy();
 c1_1->SetGridx();
 gStyle->SetStatY(0.95); 
 gStyle->SetStatX(0.89); 
 gStyle->SetStatTextColor(1);

 tlepTa->Rebin(rbin);
 tlepTa->SetLineColor(1);
 tlepTa->DrawCopy();
 c1->Update();

 gStyle->SetStatX(0.69); 
 gStyle->SetStatTextColor(2);
 tSelT->ProjectionX("tSelT_px",0,-1,"");
 tSelT_px->Rebin(rbin);
 tSelT_px->SetLineWidth(2);
 tSelT_px->SetLineColor(2);
 tSelT_px->DrawCopy("sames");
 c1->Update();

 gStyle->SetStatX(0.49); 
 gStyle->SetStatTextColor(4);
 tlepTp->Rebin(rbin);
 tlepTp->SetLineColor(4);
 tlepTp->DrawCopy("sames");
 c1->Update();

 c1->cd(2);
 c1_2->SetLogy();
 c1_2->SetGridx();
 gStyle->SetStatX(0.89); 
 gStyle->SetStatTextColor(1);

 wlepTa->Rebin(rbin);
 wlepTa->SetLineColor(1);
 wlepTa->DrawCopy();
 c1->Update();

 gStyle->SetStatX(0.69); 
 gStyle->SetStatTextColor(2);
 wSelT->ProjectionX("wSelT_px",0,-1,"");
 wSelT_px->Rebin(rbin);
 wSelT_px->SetLineWidth(2);
 wSelT_px->SetLineColor(2);
 wSelT_px->DrawCopy("sames");
 c1->Update();

 gStyle->SetStatX(0.49); 
 gStyle->SetStatTextColor(4);
 wlepTp->Rebin(rbin);
 wlepTp->SetLineColor(4);
 wlepTp->DrawCopy("sames");
 c1->Update();

 c1->cd(3);
 c1_3->SetLogy();
 c1_3->SetGridx();
 gStyle->SetOptStat(1000101);
 gStyle->SetStatX(0.89); 
 gStyle->SetStatTextColor(1);

 qlepTa->Rebin(rbin);
 qlepTa->SetLineColor(1);
 qlepTa->DrawCopy();
 c1->Update();

 gStyle->SetStatX(0.69); 
 gStyle->SetStatTextColor(2);
 qSelT->ProjectionX("qSelT_px",0,-1,"");
 qSelT_px->Rebin(rbin);
 qSelT_px->SetLineWidth(2);
 qSelT_px->SetLineColor(2);
 qSelT_px->DrawCopy("sames");
 c1->Update();

 gStyle->SetStatX(0.49); 
 gStyle->SetStatTextColor(4);
 qlepTp->Rebin(rbin);
 qlepTp->SetLineColor(4);
 qlepTp->DrawCopy("sames");
 c1->Update();

 c1->cd(4);
 c1_4->SetLogy();
 c1_4->SetGridx();
 gStyle->SetStatY(0.95); 
 gStyle->SetStatX(0.89); 
 gStyle->SetStatTextColor(1);

 tMCT->ProjectionX("tMCT_px",0,-1,"");
 tMCT_px->Rebin(rbin);
 tMCT_px->SetLineColor(1);
 tMCT_px->SetAxisRange(1,5000,"Y");
 tMCT_px->DrawCopy();
 c1->Update();

 c1->Print(plot1);

 TCanvas *c2 = new TCanvas("c2","", 1000, 800);
 c2->SetFillColor(10);
 c2->SetFillColor(10);
 c2->Divide(2,2);

 c2->cd(1);
 c2_1->SetLogy();
 c2_1->SetGridx();
 gStyle->SetStatX(0.89); 
 gStyle->SetStatTextColor(1);

 thadTa->Rebin(rbin);
 thadTa->SetLineColor(1);
 thadTa->DrawCopy();
 c2->Update();

 gStyle->SetStatX(0.69); 
 gStyle->SetStatTextColor(2);
 tSelT->ProjectionY("tSelT_py",0,-1,"");
 tSelT_py->Rebin(rbin);
 tSelT_py->SetLineWidth(2);
 tSelT_py->SetLineColor(2);
 tSelT_py->DrawCopy("sames");
 c2->Update();

 gStyle->SetStatX(0.49); 
 gStyle->SetStatTextColor(4);
 thadTp->Rebin(rbin);
 thadTp->SetLineColor(4);
 thadTp->DrawCopy("sames");
 c2->Update();

 c2->cd(2);
 c2_2->SetLogy();
 c2_2->SetGridx();
 gStyle->SetStatX(0.89); 
 gStyle->SetStatTextColor(1);

 whadTa->Rebin(rbin);
 whadTa->SetLineColor(1);
 whadTa->DrawCopy();
 c2->Update();

 gStyle->SetStatX(0.69); 
 gStyle->SetStatTextColor(2);
 wSelT->ProjectionY("wSelT_py",0,-1,"");
 wSelT_py->Rebin(rbin);
 wSelT_py->SetLineWidth(2);
 wSelT_py->SetLineColor(2);
 wSelT_py->DrawCopy("sames");
 c2->Update();

 gStyle->SetStatX(0.49); 
 gStyle->SetStatTextColor(4);
 whadTp->Rebin(rbin);
 whadTp->SetLineColor(4);
 whadTp->DrawCopy("sames");
 c2->Update();

 c2->cd(3);
 c2_3->SetLogy();
 c2_3->SetGridx();
 gStyle->SetOptStat(1000101);
 gStyle->SetStatX(0.89); 
 gStyle->SetStatTextColor(1);

 qhadTa->Rebin(rbin);
 qhadTa->SetLineColor(1);
 qhadTa->DrawCopy();
 c2->Update();

 gStyle->SetStatX(0.69); 
 gStyle->SetStatTextColor(2);
 qSelT->ProjectionY("qSelT_py",0,-1,"");
 qSelT_py->Rebin(rbin);
 qSelT_py->SetLineWidth(2);
 qSelT_py->SetLineColor(2);
 qSelT_py->DrawCopy("sames");
 c2->Update();

 gStyle->SetStatX(0.49); 
 gStyle->SetStatTextColor(4);
 qhadTp->Rebin(rbin);
 qhadTp->SetLineColor(4);
 qhadTp->DrawCopy("sames");
 c2->Update();

 c2->cd(4);
 c2_4->SetLogy();
 c2_4->SetGridx();
 gStyle->SetStatY(0.95); 
 gStyle->SetStatX(0.89); 
 gStyle->SetStatTextColor(1);

 tMCT->ProjectionY("tMCT_py",0,-1,"");
 tMCT_py->Rebin(rbin);
 tMCT_py->SetLineColor(1);
 tMCT_py->SetAxisRange(1,20000,"Y");
 tMCT_py->DrawCopy();
 c2->Update();

 c2->Print(plot2);

 // Reset gStyle
 gStyle->SetStatX(0.95); 
 gStyle->SetStatY(0.99); 
 gStyle->SetStatTextColor(1);
 gStyle->SetOptStat(1111);

 gSystem->cd("../");

 file2->Close();
 file1->Close();
 file->Close();

}
