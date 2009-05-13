#include <stdio.h>

void TtMassFit4( TString idx, const Int_t nM ) {

    TString plot0 = "p1Scan.gif";
    TString plot1 = "p2Scan.gif";
    TString plot2 = "p3Scan.gif";
    TString plot3 = "p4Scan.gif";
    TString plot4 = "p5Scan.gif";
    TString plot5 = "p6Scan.gif";
    TString plot6 = "p7Scan.gif";
    TString plot7 = "p8Scan.gif";
    TString plot8 = "p9Scan.gif";

    TString plot = "p"+idx+"Scan.gif" ;
    TString title_plot = "para"+idx+"_InputMass";
    TString hfolder = "tt_fall08j10";

    FILE *pfile = fopen("paraf.log","r");
    FILE *ffile = fopen("extrf.log","a");

    Double_t pa0[nM] ;
    Double_t pa1[nM] ;
    Double_t pa2[nM] ;
    Double_t pa3[nM] ;
    Double_t pa4[nM] ;
    Double_t pa5[nM] ;
    Double_t pa6[nM] ;
    Double_t pa7[nM] ;
    Double_t pa8[nM] ;
    Double_t pa9[nM] ;
    Double_t pa10[nM] ;
    Double_t pa11[nM] ;
    Double_t pr[nM] ;
    Double_t inputM[nM] ;

    double inputM5[5] = { 161.2, 166.2, 171.2, 176.2, 181.2 };
    double inputM9[9] = { 161.2, 163.7, 166.2, 168.7, 171.2, 173.7, 176.2, 178.7, 181.2 };

    float f1;
    for ( int i=0 ; i < nM;  i++) {
        for ( int j=0 ; j< 12; j++) {

            fscanf(pfile, "%f" , &f1 );
            
            if ( j == 0 ) pa0[i] = f1;
            if ( j == 1 ) pa1[i] = f1;
            if ( j == 2 ) pa2[i] = f1;
            if ( j == 3 ) pa3[i] = f1;
            if ( j == 4 ) pa4[i] = f1;
            if ( j == 5 ) pa5[i] = f1;
            if ( j == 6 ) pa6[i] = f1;
            if ( j == 7 ) pa7[i] = f1;
            if ( j == 8 ) pa8[i] = f1;
            if ( j == 9 ) pa9[i] = f1;
            if ( j == 10 ) pa10[i] = f1;
            if ( j == 11 ) pa11[i] = f1;

            //cout<<" para => "<< pa1[i] <<endl;
        }
    }
    if ( nM == 5 ) FillPlotArray( inputM5, inputM, nM);
    if ( nM == 9 ) FillPlotArray( inputM9, inputM, nM);

    if ( idx == "0" ) FillPlotArray( pa0, pr, nM);
    if ( idx == "1" ) FillPlotArray( pa1, pr, nM);
    if ( idx == "2" ) FillPlotArray( pa2, pr, nM);
    if ( idx == "3" ) FillPlotArray( pa3, pr, nM);
    if ( idx == "4" ) FillPlotArray( pa4, pr, nM);
    if ( idx == "5" ) FillPlotArray( pa5, pr, nM);
    if ( idx == "6" ) FillPlotArray( pa6, pr, nM);
    if ( idx == "7" ) FillPlotArray( pa7, pr, nM);
    if ( idx == "8" ) FillPlotArray( pa8, pr, nM);
    if ( idx == "9" ) FillPlotArray( pa9, pr, nM);
    if ( idx == "10" ) FillPlotArray( pa10, pr, nM);
    if ( idx == "11" ) FillPlotArray( pa11, pr, nM);
  


    gSystem->mkdir(hfolder);
    gSystem->cd(hfolder);

  
    TCanvas *c1 = new TCanvas("c1","", 700, 600);
    c1->SetFillColor(10);
    c1->SetFillColor(10);
    gPad->SetGridx();
    gPad->SetGridy();
    //c1->Divide(1,2);
    c1->cd();

    gStyle->SetOptStat(kTRUE);
    gStyle->SetOptFit(111);

    hMass = new TGraph( nM, inputM, pr );   
    hMass->SetMarkerColor(4);
    hMass->SetMarkerStyle(20);
    hMass->SetLineWidth(2);
    hMass->SetLineColor(4);
    hMass->SetTitle(title_plot);
    hMass->Draw("ACP");

    double a0 = 0;
    double a1 = 1;
    int fstatus = 999;
    hMass->LeastSquareLinearFit(nM, a0, a1, fstatus, 150., 185);
    cout<<" a0= "<<a0 <<"   a1= "<<a1<<"  fstatus:"<<fstatus<<endl;
    fprintf(ffile,"  %.3f  %.3f \n", a0, a1 );


    TF1 *fith = new TF1("fith",fitf,150,185,2);
    fith->FixParameter(0, a0);
    fith->FixParameter(1, a1);
    fith->SetLineStyle(2);
    fith->SetLineColor(7);
    fith->Draw("same");

    c1->Update();

    c1->Print(plot);

    gSystem->cd("../");

    fclose(pfile);
    fclose(ffile);
}

Double_t fitf(Double_t *x, Double_t *par) {

         Double_t fitval =  par[0] + par[1]*x[0] ;
         return fitval;
}

void FillPlotArray(Double_t* a1, Double_t* a2, const Int_t nM) {

    for ( int k=0; k<nM; k++ ) {
        a2[k] = a1[k] ;
    }
}
