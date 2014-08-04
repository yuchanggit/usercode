#include "TH1.h"
#include "TFile.h"
#include "fit_chi2_err.C"
#include "TStyle.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TCanvas.h"

void Fit_new_samples(){

  TFile *fGetData = new TFile("skimmed_data.root");
//  TH1F *h_Data_TIGHT_Quark_Barrel = new TH1F("Data_TIGHT_Quark_Barrel","",50,0,0.05);
TH1F  *h_Data_TIGHT_Quark_Barrel = (TH1F*)fGetData->Get("Inclusive_Barrel_Sigma_Ieta_Ieta_TIGHT_Quark");

  TFile *fSignal = new TFile("Tprime_M950_sample.root");
//  TH1F *h_Signal_TIGHT_Quark_Barrel = new TH1F("Signal_TIGHT_Quark_Barrel","",50,0,0.05);
TH1F  *h_Signal_TIGHT_Quark_Barrel = (TH1F*)fSignal->Get("Inclusive_Barrel_Sigma_Ieta_Ieta_TIGHT_Quark");

  TFile *fBackground = new TFile("QCD_PtHat80to170.root");
//  TH1F *h_Background_TIGHT_Quark_Barrel = new TH1F("Background_TIGHT_Quark_Barrel","",50,0,0.05);
TH1F  *h_Background_TIGHT_Quark_Barrel = (TH1F*)fBackground->Get("Inclusive_Barrel_Sigma_Ieta_Ieta_TIGHT_Quark");

 gStyle->SetOptStat();
TCanvas *c1 = new TCanvas("c1","c1",640,640);
  h_Data_TIGHT_Quark_Barrel->GetXaxis()->SetRangeUser(0,0.05);
  h_Data_TIGHT_Quark_Barrel->Draw();

 gStyle->SetOptStat();
TCanvas *c2 = new TCanvas("c2","c2",640,640);
  h_Signal_TIGHT_Quark_Barrel->GetXaxis()->SetRangeUser(0,0.05);
  h_Signal_TIGHT_Quark_Barrel->Draw();

 gStyle->SetOptStat();
TCanvas *c3 = new TCanvas("c3","c3",640,640);
  h_Background_TIGHT_Quark_Barrel->GetXaxis()->SetRangeUser(0,0.05);
  h_Background_TIGHT_Quark_Barrel->Draw();


TCanvas *c4 = new TCanvas("c4","c4",640,640);

      double sig, sigerr,chi2;
     fit_chi2_err(h_Data_TIGHT_Quark_Barrel,h_Signal_TIGHT_Quark_Barrel,h_Background_TIGHT_Quark_Barrel,sig,sigerr,chi2);
      cout << "Fitted result = " << sig << " +- " << sigerr << endl;



}
