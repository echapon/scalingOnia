#ifndef hwu_reader_h
#define hwu_reader_h

#include "TGraphAsymmErrors.h"
#include "TString.h"

#include <fstream>
#include <vector>
#include <string>

using namespace std;

const int icol_dy = 1; // hard-coding where the error column is

TGraphAsymmErrors* read_hwu(const char* fname, int ihist, int icol) {
   ifstream file(fname);

   vector<double> x, y, dx, dy;

   bool isInHist = false;
   string theline;
   TString thelineT;
   int chist=0;
   TString title;

   while (file.good()) {
      getline(file,theline);
      thelineT = TString(theline.c_str());
      if (!isInHist && thelineT.Contains("<histogram>")) { // entering the histo
         isInHist = true;

         if (chist==ihist) { // get the title
            TString tok;
            Ssiz_t from = 0;
            int cnt=0;
            while (thelineT.Tokenize(tok, from, "\"")) {
               if (cnt==1) title = tok;
               cnt++;
            }
         }

         continue;
      }
      if (isInHist && thelineT.Contains("<\\histogram>")) { // leaving the histo
         isInHist = false;
         chist++;
         continue;
      }
      if (!isInHist || ihist!=chist) continue;

      // we are in an histogram, and it is the histo we want. Read it!
      TString tok;
      Ssiz_t from = 0;
      int ccol = -4;
      double xl, xh;
      while (thelineT.Tokenize(tok, from, " ")) {
         ccol++;
         double val = atof(tok.Data());
         if (ccol==-2) xl = val;
         else if (ccol==-1) xh = val;
         else if (ccol==icol) y.push_back(val);
         else if (ccol==icol_dy) dy.push_back(val);
      }
      x.push_back((xl+xh)/2.);
      dx.push_back((xh-xl)/2.);
   }

   const int n = x.size();

   // for (int i=0; i<n; i++)
      // cout << x[i] << " " << y[i] << " " << dy[i] << endl;

   TGraphAsymmErrors *ans = new TGraphAsymmErrors(n,x.data(),y.data(),dx.data(),dx.data(),dy.data(),dy.data());
   ans->SetName(Form("gHwU_%i%i",ihist,icol));
   ans->SetTitle(title);
   return ans;

}

#endif // ifndef hwu_reader_h
