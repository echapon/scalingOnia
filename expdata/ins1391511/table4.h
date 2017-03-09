#ifndef ins1391511_table4
#define ins1391511_table4

#include "TGraphAsymmErrors.h"
#include <string>
#include "TFile.h"
#include "TH1.h"

namespace ins1391511table4 {
   // metadata
   string location = "Figure 2";
   string dscomment = "Differential cross-sections as a function of $p_\\perp$ integrated over $y$ for prompt $J/\\psi$ mesons. The first uncertainties are statistical and the second (third) are uncorrelated (correlated) systematic uncertainties amongst bins.";
   string reackey = "P P --> J/PSI X";
   string obskey = "DSIG/DPT";
   string qual = "RE : P P --> J/PSI < MU+ MU- > X";
   double sqrts = 13000.0; // GEV
   string yheader = "d#sigma/dp_{T}  [nb/GeV]"; // "DSIG/DPT [NB/(GEV/c)]";
   string xheader = "p_{T} [GeV]"; // "PT [GEV] ";


   // Get the histos from the files
   TGraphAsymmErrors* graph_table4_y1(bool stat = false, bool correl=true) {
      TFile *f = new TFile("expdata/ins1391511/HEPData-ins1391511-v1-Table4.root");
      TH1F *hy = (TH1F*) f->Get("Table 4/Hist1D_y1");
      TH1F *he1 = (TH1F*) f->Get("Table 4/Hist1D_y1_e1");
      TH1F *he2 = (TH1F*) f->Get("Table 4/Hist1D_y1_e2");
      TH1F *he3 = (TH1F*) f->Get("Table 4/Hist1D_y1_e3");

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

#endif // #ifndef ins1391511_table4
