#ifndef cmsbph15005_tableA1
#define cmsbph15005_tableA1

#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"

#include <fstream>
#include <vector>
#include <string>

using namespace std;

namespace cmsbph15005tableA1 {
   // metadata
   string location = "Table A.1.";
   string dscomment = "J/psi cross section";
   string reackey = "P P --> J/PSI X";
   string obskey = "DSIG/DPT";
   string qual = "RE : P P --> J/PSI < MU+ MU- > X";
   double sqrts = 13000.0; // GEV
   string yheader = "d#sigma/dp_{T}  [nb/GeV]"; // "DSIG/DPT [NB/(GEV)]";
   string xheader = "p_{T} [GeV]"; // "YRAP ";


   // Get the histos from the files
   TGraphAsymmErrors* graph_tableA1(int irap, bool stat = false) {
      ifstream file("expdata/cmsbph15005/jpsixsec_all.dat");
      vector<double> x, y, dx, dy;

      int icol = 3+3*irap;
      int icol_dy = 4+3*irap;
      if (!stat) icol_dy++;

      string theline;
      TString thelineT;
      while (file.good()) {
         getline(file,theline);
         thelineT = TString(theline.c_str());
         TString tok;
         Ssiz_t from = 0;
         int ccol = -2;
         double xl=-1, xh=-1;
         while (thelineT.Tokenize(tok, from, " ")) {
            ccol++;
            double val = atof(tok.Data());
            if (ccol==0 && irap<4 && val>100) break; // ignore the last line for irap<4
            if (ccol==0) xl = val;
            else if (ccol==1) xh = val;
            else if (ccol==icol) y.push_back(val*1e-3); // 1e-3 for pb -> nb
            else if (ccol==icol_dy) dy.push_back(val*0.01*y.back()); 
         }
         if (xl>0) {
            x.push_back((xl+xh)/2.);
            dx.push_back((xh-xl)/2.);
         }
      }


      const int n = x.size();

      if (y.size() != x.size() || dx.size() != x.size() || dy.size() != x.size()) {
         cout << "Error in " << __FILE__ << "): inconsistent vector sizes." << endl;
         return NULL;
      }

      // for (int i=0; i<n; i++)
      //    cout << x[i] << " " << y[i] << " " << dy[i] << endl;

      TGraphAsymmErrors *ans = new TGraphAsymmErrors(n,x.data(),y.data(),dx.data(),dx.data(),dy.data(),dy.data());
      ans->SetTitle(Form("CMS-BPH-15-005_%d_%d",irap,(int) stat));
      ans->SetName(Form("CMS-BPH-15-005_%d_%d",irap,(int) stat));

      return ans;
   };
};

#endif // #ifndef cmsbph15005_tableA1
