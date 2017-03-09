#ifndef ins1230344_table7
#define ins1230344_table7

#include "TGraphAsymmErrors.h"
#include <string>
#include "TFile.h"
#include "TH1.h"

namespace ins1230344table7 {
   // metadata
   string location = "Figure 3,5";
   string dscomment = "The double-differential cross sections for prompt J/PSI production (assuming no polarisation) and         production of J/PSI from b-hadron decays as a function of transverse momentum for the rapidity range 4.0-4.5.         Also shown in the final column is the fraction (in %) of J/PSIs from the latter.";
   string reackey = "P P --> J/PSI X";
   string obskey = "DSIG/DYRAP";
   string qual = "RE : P P --> J/PSI < MU+ MU- > X";
   double sqrts = 8000.0; // GEV
   string yheader = "d^{2}#sigma/dp_{T}/dy  [nb/GeV]"; // "D2(SIG)/DPT/DYRAP [NB/GEV]";
   string xheader = "p_{T} [GeV]"; // "PT [GEV] ";


   // Get the histos from the files
   TGraphAsymmErrors* graph_table7_y1(bool stat = false, bool correl=true) {
      TFile *f = new TFile("expdata/ins1230344/HEPData-ins1230344-v1-Table7.root");
      TH1F *hy = (TH1F*) f->Get("Table 7/Hist1D_y1");
      TH1F *he1 = (TH1F*) f->Get("Table 7/Hist1D_y1_e1");
      TH1F *he2 = (TH1F*) f->Get("Table 7/Hist1D_y1_e2");
      TH1F *he3 = (TH1F*) f->Get("Table 7/Hist1D_y1_e3");

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

#endif // #ifndef ins1230344_table7
