#ifndef ins1345023_table1
#define ins1345023_table1

#include "TGraphAsymmErrors.h"
#include <string>
#include "TFile.h"
#include "TH1.h"

namespace ins1345023table1 {
   // metadata
   string location = "Figure 4";
   string dscomment = "Differential production cross section in rapidity for prompt J/PSI     mesons (assuming no polarisation) and from J/PSI from b-hadron decays.";
   string reackey = "P P --> J/PSI X";
   string obskey = "DSIG/DYRAP";
   string qual = "RE : P P --> J/PSI < MU+ MU- > X";
   double sqrts = 7000.0; // GEV
   string yheader = "d^{2}#sigma/dp_{T}dy  [nb]"; // "B*D2(SIG)/DPT/DYRAP";
   string xheader = "p_{T} [GeV]"; // "";


   // Get the histos from the files
   TGraphAsymmErrors* graph_table1_y1(int irap, bool stat = false, bool correl=true) {
      TFile *f = new TFile("expdata/ins1345023/HEPData-ins1345023-v1-Table1.root");
      TH1F *hy = (TH1F*) f->Get(Form("Table 1/Hist1D_y%d",irap));
      TH1F *he1 = (TH1F*) f->Get(Form("Table 1/Hist1D_y%d_e1",irap));
      TH1F *he2 = (TH1F*) f->Get(Form("Table 1/Hist1D_y%d_e2",irap));
      TH1F *he3 = (TH1F*) f->Get(Form("Table 1/Hist1D_y%d_e3",irap));

      int nbins = hy->GetNbinsX();
      double *x = new double[nbins];
      double *ex = new double[nbins];
      double *y = new double[nbins];
      double *ey = new double[nbins];
      for (int i=0; i<nbins; i++) {
         x[i] = hy->GetBinCenter(i+1);
         ex[i] = hy->GetBinWidth(i+1)/2.;
         if (stat) ey[i] = he1->GetBinContent(i+1);
         else {
            if (!correl) ey[i] = he2->GetBinContent(i+1);
            else ey[i] = sqrt(pow(he2->GetBinContent(i+1),2)+pow(he3->GetBinContent(i+1),2));
         }
         y[i] = hy->GetBinContent(i+1);

      }

      TGraphAsymmErrors *ans = new TGraphAsymmErrors(nbins,x,y,ex,ex,ey,ey);
      ans->SetTitle(hy->GetTitle());
      ans->SetName(hy->GetName());
      delete[] x;
      delete[] y;
      delete[] ex;
      delete[] ey;
      delete f;

      return ans;
   };
}

#endif // #ifndef ins1345023_table1
