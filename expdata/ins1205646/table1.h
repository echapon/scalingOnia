#ifndef ins1205646_table1
#define ins1205646_table1

#include "TGraphAsymmErrors.h"
#include <string>
#include "TFile.h"
#include "TH1.h"

namespace ins1205646table1 {
   // metadata
   string location = "Figure 1";
   string dscomment = "he production of J/psi mesons is studied with the LHCb detector using data from pp collisions at sqrt(s)=2.76 TeV corresponding to an integrated luminosity of 71 nb^-1. The differential cross-section for inclusive J/psi production is measured as a function of its transverse momentum pT. The cross-section in the fiducial region 0<pT<12 GeV/c and rapidity 2.0<y<4.5 is measured to be 5.6 +/- 0.1(stat) +/- 0.4 (syst) mub, with the assumption of unpolarised J/psi production. The fraction of J/psi production from b-hadron decays is measured to be (7.1 +/- 0.6 (stat) +/- 0.7 (syst)) %";
   string reackey = "P P --> J/PSI X";
   string obskey = "DSIG/DPT";
   string qual = "RE : P P --> J/PSI < MU+ MU- > X";
   double sqrts = 2760.0; // GEV
   string yheader = "d#sigma/dp_{T}  [nb/GeV]"; // "DSIG/DPT [NB/(GEV)]";
   string xheader = "p_{T} [GeV]"; // "YRAP ";


   // Get the histos from the files
   TGraphAsymmErrors* graph_table1_y1(bool stat = false) {
      TFile *f = new TFile("expdata/ins1205646/HEPData-ins1205646-v1-Table1.root");
      TH1F *hy = (TH1F*) f->Get("Table 1/Hist1D_y1");
      TH1F *he1 = (TH1F*) f->Get("Table 1/Hist1D_y1_e1");
      TH1F *he2 = (TH1F*) f->Get("Table 1/Hist1D_y1_e2");

      int nbins = hy->GetNbinsX();
      double *x = new double[nbins];
      double *ex = new double[nbins];
      double *y = new double[nbins];
      double *ey = new double[nbins];
      for (int i=0; i<nbins; i++) {
         x[i] = hy->GetBinCenter(i+1);
         ex[i] = hy->GetBinWidth(i+1)/2.;
         if (stat) ey[i] = he1->GetBinContent(i+1);
         else ey[i] = he2->GetBinContent(i+1);
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
};

#endif // #ifndef ins1205646_table1
