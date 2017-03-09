#ifndef ins1391511_table6
#define ins1391511_table6

#include "TGraphAsymmErrors.h"
#include <string>
#include "TFile.h"
#include "TH1.h"

namespace ins1391511table6 {
   // metadata
   string location = "Figure 6a";
   string dscomment = "Differential cross-sections as a function of $y$ integrated over $p_\\perp$ for prompt $J/\\psi$ mesons. The first uncertainties are statistical and the second (third) are the correlated (uncorrelated) systematic uncertainties.";
   string reackey = "P P --> J/PSI X";
   string obskey = "DSIG/DYRAP";
   string qual = "RE : P P --> J/PSI < MU+ MU- > X";
   double sqrts = 13000.0; // GEV
   string yheader = "d#sigma/y  [#mu b]"; // "DSIG/DYRAP [MUB]";
   string xheader = "p_{T} [GeV]"; // "PT [GEV] ";


   // Get the histos from the files
   TGraphAsymmErrors* graph_table6_y1(bool stat = false, bool correl=true) {
      TFile *f = new TFile("expdata/ins1391511/HEPData-ins1391511-v1-Table6.root");
      TH1F *hy = (TH1F*) f->Get("Table 6/Hist1D_y1");
      TH1F *he1 = (TH1F*) f->Get("Table 6/Hist1D_y1_e1");
      TH1F *he2 = (TH1F*) f->Get("Table 6/Hist1D_y1_e2");
      TH1F *he3 = (TH1F*) f->Get("Table 6/Hist1D_y1_e3");

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

#endif // #ifndef ins1391511_table6
