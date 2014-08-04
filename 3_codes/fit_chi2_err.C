#include <TH1.h>
#include <TMinuit.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <iostream>
#include <TH2.h>

TH1D* signal_pos;
TH1D* background_pos;
TH1D* data;
TH1D* fit_result;

using namespace std;


// -- returns prediction and error on the prediction based on stat unc in templates
Double_t ftotal_pos(Double_t *x, Double_t *par, Double_t& err) {  
//Double_t ftotal_pos(Double_t *x, Double_t *par) {  
  Double_t xx = x[0];
  Int_t bin_background = background_pos->GetXaxis()->FindBin(xx);
  Int_t bin_signal = signal_pos->GetXaxis()->FindBin(xx);
  Double_t br = (1-par[0])*background_pos->GetBinContent(bin_background);
  Double_t sr = par[0]*signal_pos->GetBinContent(bin_signal);

  double bkgerr = background_pos->GetBinError(bin_background);
  double sigerr = signal_pos->GetBinError(bin_signal);
  
 
  err = sqrt((1-par[0])*(1-par[0])*bkgerr*bkgerr + 
	     par[0]*par[0]*sigerr*sigerr);
 

  return sr + br;


}


Double_t ferr_sigtemplate(Double_t *x, Double_t *par){

  Double_t xx = x[0];

  Int_t bin_signal = signal_pos->GetXaxis()->FindBin(xx);
  
  Double_t sigerr = signal_pos->GetBinError(bin_signal);

  return sigerr;
}

Double_t ferr_bkgtemplate(Double_t *x, Double_t *par){

  Double_t xx = x[0];

  Int_t bin_background = background_pos->GetXaxis()->FindBin(xx);
  
  Double_t bkgerr = background_pos->GetBinError(bin_background);
  
  return bkgerr;
}


// -- Defines Chi2 to minimize taking into account error on the prediction from templates stat uncertainties
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  Double_t delta2;
  Double_t sum=0;

  for (int i=1; i<= data->GetNbinsX(); i++) {
    Double_t data_value = data->GetBinContent(i);
    if(data_value==0)continue;
    Double_t data_err = data->GetBinError(i);
    Double_t temp_x[2];
    Double_t center = data->GetBinCenter(i);
    temp_x[0]=center;
    Double_t sr, br, serr, berr;
    Double_t theory_err;
    Double_t theory_value = ftotal_pos(temp_x,par,theory_err);
    
    // delta2 = (theory_value-data_value)*(theory_value-data_value)/
    //   ( pow(data_err,2) + pow(theory_err,2));
    delta2 = (theory_value-data_value)*(theory_value-data_value)/
      ( pow(data_err,2));

    sum +=  delta2;
    
    // for this round, set the fit content
    fit_result->SetBinContent(i,theory_value);
    fit_result->SetBinError(i,theory_err);
   

  }
  f=sum  ;

}



// -- main function
//void fit_chi2_err(TH1F* dataInput, TH1F* sigTemplate, TH1F* bkgTemplate, std::string prefix, Double_t& sigFrac, Double_t& sigFrac_err)
//void fit_chi2_err(TH1F* dataInput, TH1F* sigTemplate, TH1F* bkgTemplate, Double_t& sigFrac, Double_t& sigFrac_err)
void fit_chi2_err(TH1F* dataInput, TH1F* sigTemplate, TH1F* bkgTemplate, Double_t& sigFrac, Double_t& sigFrac_err, Double_t& FitChi2)
{
  gStyle->SetOptStat(kFALSE);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);

  Double_t scale=1.;




  data = (TH1D*)dataInput->Clone();
  data->SetName("data");
  data->SetLineColor(1);
  data->SetMarkerColor(1);
  data->SetXTitle("SigmaIetaIeta");
  data->Sumw2();
//  scale = 1.0/(Double_t)data->Integral(0,1000); 
  scale = 1.0/(Double_t)data->Integral(0,2000);
//  scale = 1.0/(Double_t)data->Integral();//changed by Yu-hsiang
  cout << "scale for data = " << scale << endl;
  data->Scale(scale);


  fit_result = (TH1D*)dataInput->Clone();
  fit_result->SetName("fit_result");
  fit_result->SetLineColor(8); 
  fit_result->SetMarkerColor(8);
  fit_result->SetLineStyle(2);
  fit_result->Sumw2();
  fit_result->Reset();


  signal_pos = (TH1D*)sigTemplate->Clone();
  signal_pos->SetName("signal_pos");
  signal_pos->SetLineColor(2);
  signal_pos->SetMarkerColor(2);
  signal_pos->SetFillColor(2);
  signal_pos->SetXTitle("SigmaIetaIeta");
  signal_pos->Sumw2();
//  scale = 1.0/(Double_t)signal_pos->Integral(0,1000); 
  scale = 1.0/(Double_t)signal_pos->Integral(0,2000);
//  scale = 1.0/(Double_t)signal_pos->Integral();//changed by Yu-hsiang
  cout << "scale for signal template = " << scale << endl;
  signal_pos->Scale(scale);


  background_pos = (TH1D*)bkgTemplate->Clone();
  background_pos->SetName("background_pos");
  background_pos->SetLineColor(4);
  background_pos->SetMarkerColor(4);
  background_pos->SetFillColor(4);
  background_pos->SetXTitle("SigmaIetaIeta");
  background_pos->Sumw2();
//  scale = 1.0/(Double_t)background_pos->Integral(0,1000); 
  scale = 1.0/(Double_t)background_pos->Integral(0,2000);
//  scale = 1.0/(Double_t)background_pos->Integral();//changed by Yu-hsiang
  cout << "scale for background template = " << scale << endl;
  background_pos->Scale(scale);



  TMinuit *gMinuit = new TMinuit(1);  //initialize TMinuit with a maximum of 5 (1param??) params
  gMinuit->SetFCN(fcn); // sets function to minimize: fcn is Chi2 with errors on templates

  Double_t arglist[10];
  Int_t ierflg = 0; // status flag, it is 0 when ereything goes fine

  // -- sets error
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);


  Double_t vstart = 0.5;
  Double_t step = 0.001;
  gMinuit->mnparm(0, "fsig", vstart, step, 0,1,ierflg);



  // Now ready for minimization step
  arglist[0] = 1000;
  arglist[1] = 0.01;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

  Double_t fsig=0;
  Double_t fsigerr=0;
  Double_t chi2 = 0;


  if ( ierflg == 0 ) 
    {

      // Print results
      Double_t amin,edm,errdef;
      Int_t nvpar,nparx,icstat;
      gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
      gMinuit->mnprin(3,amin);
      chi2 = gMinuit->fAmin/9;
      gMinuit->GetParameter(0, fsig, fsigerr );  
      cout << "Fsig = " << fsig << " +- " << fsigerr << endl;
      cout << "Chi2/degree of freedom = " << chi2 <<endl;

//      TCanvas* c1 = new TCanvas("c1","",500,500);

  TH2F *htmp2 = new TH2F("htmp2","",100, 0., 0.025, 100, 0., data->GetMaximum()*1.25);

  htmp2->SetNdivisions(505,"XY");
  htmp2->SetXTitle("SigmaIetaIeta");
  htmp2->SetYTitle("A.U.");
  htmp2->SetLineColor(1); 

//     htmp2->Draw();

  data->Draw("e");
  TH1D* signal_display = (TH1D*)signal_pos->Clone();
  signal_display->SetName("signal_display");
//  signal_display->Scale(fsig/signal_display->Integral(0,1000));
  signal_display->Scale(fsig/signal_display->Integral());
  signal_display->SetFillStyle(3001);
  signal_display->Draw("histsame");
      
      
  TH1D* background_display = (TH1D*)background_pos->Clone();
  background_display->SetName("background_display");
//  background_display->Scale((1-fsig)/background_display->Integral(0,1000));
  background_display->Scale((1-fsig)/background_display->Integral());
  background_display->SetFillStyle(3001);
  background_display->Draw("histsame");
//      background_display->Draw("same");
  data->Draw("same");

  fit_result->Draw("histesame");

  char result[300];
  sprintf(result,"fsig = %.3lf #pm %.3lf",fsig,fsigerr );
  sigFrac = fsig;
  sigFrac_err = fsigerr ;
  FitChi2 = chi2;

  TLegend* leg = new TLegend(0.6,0.6,0.4,0.9);
  leg->SetHeader(result);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->SetBorderSize(0);
  leg->AddEntry(data,"data");
  leg->AddEntry(fit_result,"fit");
  leg->AddEntry(signal_display,"signal template");
  leg->AddEntry(background_display,"background template");
  leg->Draw("same");

 /*     std::string outputFile = prefix + ".eps";
      c1->Print(outputFile.data());

      outputFile = prefix + ".gif";
      c1->Print(outputFile.data());

      outputFile = prefix + ".C";
      c1->Print(outputFile.data());*/

    }

    
  else{
    cout << "Fit failed!\n";
    sigFrac = 0;
    sigFrac_err = 0;
  }

  
  return;

}



