#ifndef ins1391511_table1
#define ins1391511_table1

#include "TGraphAsymmErrors.h"
#include <string>
#include "TFile.h"
#include "TH1.h"

namespace ins1391511table1 {
   // metadata
   string location = "Figure 2";
   string dscomment = "Double differential cross-section for prompt $J/\\psi$ mesons as a function of $p_\\perp$ in bins of $y$. The first uncertainties are statistical, the second are the correlated systematic\nuncertainties shared between bins and the last are the uncorrelated systematic uncertainties.";
   string reackey = "P P --> J/PSI X";
   string obskey = "D2SIG/DYRAP/DPT";
   string qual = "RE : P P --> J/PSI < MU+ MU- > X";
   double sqrts = 13000.0; // GEV
   // rapidities: 2-2.5-3-3.5-4-4.5
   string yheader = "d^{2}#sigma/dp_{T}/dy  [nb/GeV]"; // "D2SIG/DYRAP/DPT [NB/(GEV/c)]";
   string xheader = "p_{T} [GeV]"; // "PT [GEV] ";


   // Get the histos from the files
   TGraphAsymmErrors* graph_table1_y1(bool stat = false, bool correl=true) {
      TFile *f = new TFile("expdata/ins1391511/HEPData-ins1391511-v1-Table1.root");
      TH1F *hy = (TH1F*) f->Get("Table 1/Hist1D_y1");
      TH1F *he1 = (TH1F*) f->Get("Table 1/Hist1D_y1_e1");
      TH1F *he2 = (TH1F*) f->Get("Table 1/Hist1D_y1_e2");
      TH1F *he3 = (TH1F*) f->Get("Table 1/Hist1D_y1_e3");

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

   TGraphAsymmErrors* graph_table1_y2(bool stat = false, bool correl=true) {
      TFile *f = new TFile("expdata/ins1391511/HEPData-ins1391511-v1-Table1.root");
      TH1F *hy = (TH1F*) f->Get("Table 1/Hist1D_y2");
      TH1F *he1 = (TH1F*) f->Get("Table 1/Hist1D_y2_e1");
      TH1F *he2 = (TH1F*) f->Get("Table 1/Hist1D_y2_e2");
      TH1F *he3 = (TH1F*) f->Get("Table 1/Hist1D_y2_e3");

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

   TGraphAsymmErrors* graph_table1_y3(bool stat = false, bool correl=true) {
      TFile *f = new TFile("expdata/ins1391511/HEPData-ins1391511-v1-Table1.root");
      TH1F *hy = (TH1F*) f->Get("Table 1/Hist1D_y3");
      TH1F *he1 = (TH1F*) f->Get("Table 1/Hist1D_y3_e1");
      TH1F *he2 = (TH1F*) f->Get("Table 1/Hist1D_y3_e2");
      TH1F *he3 = (TH1F*) f->Get("Table 1/Hist1D_y3_e3");

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

   TGraphAsymmErrors* graph_table1_y4(bool stat = false, bool correl=true) {
      TFile *f = new TFile("expdata/ins1391511/HEPData-ins1391511-v1-Table1.root");
      TH1F *hy = (TH1F*) f->Get("Table 1/Hist1D_y4");
      TH1F *he1 = (TH1F*) f->Get("Table 1/Hist1D_y4_e1");
      TH1F *he2 = (TH1F*) f->Get("Table 1/Hist1D_y4_e2");
      TH1F *he3 = (TH1F*) f->Get("Table 1/Hist1D_y4_e3");

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

   TGraphAsymmErrors* graph_table1_y5(bool stat = false, bool correl=true) {
      TFile *f = new TFile("expdata/ins1391511/HEPData-ins1391511-v1-Table1.root");
      TH1F *hy = (TH1F*) f->Get("Table 1/Hist1D_y5");
      TH1F *he1 = (TH1F*) f->Get("Table 1/Hist1D_y5_e1");
      TH1F *he2 = (TH1F*) f->Get("Table 1/Hist1D_y5_e2");
      TH1F *he3 = (TH1F*) f->Get("Table 1/Hist1D_y5_e3");

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

#endif // #ifndef ins1391511_table1
