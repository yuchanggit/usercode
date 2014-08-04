#include "TH1.h"
#include "TFile.h"
#include "fit_chi2_err.C"
#include "TStyle.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TCanvas.h"

void Fillrandom_times(){


  const int NTOYS=10;//1000
  const int NBINS=12;//12
  const int NRUNBINS=12;//12
  char nsigma[500];
  TH1F *Pull[NBINS];
  char npull[500];
  TH1F *Pull2[NBINS];
  char npull2[500];
  TH1F *Mean[NBINS];
  char nmean[500];

  int Nsig[NBINS] ={0};
  int Nbkg[NBINS] ={0};
  int Numsig[NBINS]={0};
  int Numbkg[NBINS]={0};

  int ptBoundary[NBINS+1]={25,30,35,40,45,50,55,60,65,70,75,80,200};

  double InputP[NBINS]={0.7,0.3,0.4,0.4,0.4,0.4,0.5,0.5,0.6,0.6,0.7,0.8};
  double Sigma[NBINS] = {(0.006/0.279)*0.3,(0.005/0.333)*0.3,(0.004/0.385)*0.4,(0.004/0.437)*0.4,(0.006/0.467)*0.4,(0.007/0.519)*0.4,(0.009/0.553)*0.5,(0.010/0.590)*0.5,(0.012/0.615)*0.6,(0.013/0.642)*0.6,(0.016/0.660)*0.7,(0.008/0.742)*0.8};
  TFile *fGetData = new TFile("12FakeData.root");
  TH1F *GethData[NBINS];
  char numentries[500];
  int GetNum[NBINS]={0};

  TFile *fSignal = new TFile("MC_SignalCorrect_barrel.root");
  TH1F *hsignal[NBINS];
  char osignal[5000];
  TH1F *hfillSig[NBINS];
  char signal[5000];

  TFile *fBackground = new TFile("MC_BkgCorrect_barrel.root");
  TH1F *hbkg[NBINS];
  char obkg[5000];
  TH1F *hfillBag[NBINS];
  char bkgrond[5000];

  char name[5000];
  char bname[500];
  char histogram[500];

  TH1F *hSig[NBINS];
  char fillsig[500];
  TH1F *hBkg[NBINS];
  char fillbkg[500];

cout<<"test1"<<endl;


TCanvas *checkpdf = new TCanvas("checkpdf","",640,640);
//    checkpdf->cd();



  for (int f = 0; f < NRUNBINS; f++){
    cout<<"f:"<<f<<endl;
    sprintf(name,"h_%d",f);
    sprintf(npull,"pull_%d",f);
    sprintf(npull2,"pull2_%d",f);
    Pull[f] = new TH1F(npull,"",100,-5,5);
    Pull2[f] = new TH1F(npull2,"",40,-1,1);

    sprintf(nsigma,"hsigmai_%d",f);
    //Sigma[f] = new TH1F(nsigma,"",10,0,2);
    sprintf(nmean,"hmean_%d",f);
    Mean[f] = new TH1F(nmean,"",100,-5,5);

    sprintf(numentries,"hdata_%d",f);
    GethData[f] = new TH1F(numentries,"",30,0,0.03);
    GethData[f] = (TH1F*)fGetData->Get(name);
    GetNum[f] = GethData[f]->GetEntries();

    sprintf(osignal,"hosig_%d",f);
    hsignal[f] = new TH1F(osignal,"",30,0,0.03);
    hsignal[f] = (TH1F*)fSignal->Get(name);
    Numsig[f] = hsignal[f]->GetEntries();
    sprintf(signal,"hsig_%d",f);
    hfillSig[f] = new TH1F(signal,"",30,0,0.03);

    sprintf(obkg,"hobkg_%d",f);
    hbkg[f] = new TH1F(obkg,"",30,0,0.03);
    hbkg[f] = (TH1F*)fBackground->Get(name);
    Numbkg[f] = hbkg[f]->GetEntries();
    sprintf(bkgrond,"hbkg_%d",f);
    hfillBag[f] = new TH1F(bkgrond,"",30,0,0.03);

    sprintf(fillsig,"hsig2_%d",f);
    hSig[f] = new TH1F(fillsig,"",30,0,0.03);
    sprintf(fillbkg,"hbkg2_%d",f);
    hBkg[f] = new TH1F(fillbkg,"",30,0,0.03);

cout<<"test2"<<endl;

    for (int i = 0; i< NTOYS; i++){
      cout<<"f:"<<f<<" i:"<<i<<endl;
      double sig, sigerr,chi2;
      double Purity = InputP[f];
      Nsig[f] = GetNum[f]*InputP[f];
      Nbkg[f] = GetNum[f]*(1.0-InputP[f]);
      // int temp_nsig = gRandom->Poisson(Nsig[f]);
      // int temp_nbkg = gRandom->Poisson(Nbkg[f]);
      int temp_nsig = gRandom->Poisson(GetNum[f]*InputP[f]);
      int temp_nbkg = gRandom->Poisson(GetNum[f]*(1.0-InputP[f]));
     

      cout << "For this toy: nsig = " << temp_nsig << "\t, nbkg = " << temp_nbkg << endl;

      hfillSig[f]->Reset();
      hfillBag[f]->Reset();
      hSig[f]->Reset();
      hBkg[f]->Reset();

      // pseudo-data
      hfillSig[f]->FillRandom(hsignal[f], temp_nsig);
      hfillBag[f]->FillRandom(hbkg[f], temp_nbkg);
      hfillSig[f]->Add(hfillBag[f]);

      // // template
      // fluctuate bin by bin

      // for(int ib=1; ib <= hSig[f]->GetNbinsX(); ib++){

      // 	double mean = hsignal[f]->GetBinContent(ib);
      // 	double err  = hsignal[f]->GetBinError(ib);

      // 	double mean_temp = gRandom->Gaus(mean,err);
      // 	hSig[f]->SetBinContent(ib,mean_temp);
      // 	hSig[f]->SetBinError(ib,err);

      // }

      // for(int ib=1; ib <= hbkg[f]->GetNbinsX(); ib++){

      // 	double mean = hbkg[f]->GetBinContent(ib);
      // 	double err  = hbkg[f]->GetBinError(ib);

      // 	double mean_temp = gRandom->Poisson(mean);
      // 	hBkg[f]->SetBinContent(ib,mean_temp);
      // 	hBkg[f]->SetBinError(ib,err);

      // }


      // // //sprintf(histogram,"testpull_%d",f);
      // fit_chi2_err(hfillSig[f],hSig[f],hBkg[f],sig,sigerr,chi2);
      
      fit_chi2_err(hfillSig[f],hsignal[f],hbkg[f],sig,sigerr,chi2);

      cout << "Fitted result = " << sig << " +- " << sigerr << endl;
      cout << "Purity = " << Purity << endl;
      float pull_value = (sig-Purity)/(sigerr);
      if(sigerr>1e-6)Pull[f]->Fill(pull_value);
      if(Purity>1e-6)Pull2[f]->Fill((sig-Purity)/Purity);


    }//i


//TCanvas *checkpdf = new TCanvas("checkpdf","",640,640);
//    checkpdf->cd();

// gStyle->SetOptStat(1);
// gStyle->SetOptFit(1);
//    Pull[f]->Fit("gaus");


//TCanvas *checkpdf1 = new TCanvas("checkpdf1","",640,640);
//    checkpdf1->cd();

// gStyle->SetOptStat(1);
// gStyle->SetOptFit(1);
//    Pull2[f]->Fit("gaus");

// gStyle->SetOptStat(1);
// gStyle->SetOptFit(1);


//    Pull[f]->Fit("gaus");

/*            if(f==0){
                checkpdf->Print("checkpdf_comv2.pdf(");
            }else if(f==NRUNBINS-1){
                checkpdf->Print("checkpdf_comv2.pdf)");
            }else{
                checkpdf->Print("checkpdf_comv2.pdf");
            }   
*/
  }//f

  for (int f = 0; f < NRUNBINS; f++){

 gStyle->SetOptStat();
 gStyle->SetOptFit();

    Pull[f]->Fit("gaus");
    Pull2[f]->Fit("gaus");

}//f


cout<<"test3"<<endl;

  TFile *myFile = new TFile("histograms_fillrandom_poisson.root","recreate");

  for(int eiko=0; eiko < NRUNBINS; eiko++){
    std::string title = Form("Pt = %d -- %d GeV",
			      ptBoundary[eiko], ptBoundary[eiko+1]);
    Pull[eiko]->SetTitle(title.data());
    Pull[eiko]->SetXTitle("Pull: (fit-input)/fit_err");

    Pull2[eiko]->SetTitle(title.data());
    Pull2[eiko]->SetXTitle("Bias: (fit-input)/input");

    Pull[eiko]->Write();
    Pull2[eiko]->Write();
  }

}//end

