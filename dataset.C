#ifndef dataset_C
#define dataset_C

#include "include/dataset.h"

#include <fstream>
#include "TFile.h"
#include "TH1.h"

void dataset::set_graph(TGraphAsymmErrors *gstat, TGraphAsymmErrors *gsyst) {
   if (gstat) fstat = gstat;
   if (gsyst) fsyst = gsyst;
   if (gstat && gsyst) ftot = combgraph(gstat,gsyst);
   else if (gsyst) ftot = gsyst;
   else ftot = gstat;
}

void dataset::set_graph(TH1F *hist) {
   int nbins = hist->GetNbinsX();
   double *x = new double[nbins];
   double *ex = new double[nbins];
   double *y = new double[nbins];
   double *ey = new double[nbins];
   for (int i=0; i<nbins; i++) {
      x[i] = hist->GetBinCenter(i+1);
      ex[i] = hist->GetBinWidth(i+1)/2.;
      ey[i] = hist->GetBinError(i+1);
      y[i] = hist->GetBinContent(i+1);

   }
   fstat = new TGraphAsymmErrors(nbins,x,y,ex,ex,ey,ey);
   ftot = fstat;
   delete[] x;
   delete[] y;
   delete[] ex;
   delete[] ey;
}

void dataset::set_graph(const char* file_theory) {
   TFile *f = TFile::Open(file_theory);
   if (!f || !f->IsOpen()) return;
   TH1F *h = (TH1F*) f->Get("id2");
   if (!h) return;
   set_graph(h);
}

void dataset::interpolate(int n, Einterpolation type, bool doLW) {
   if (fstat) fstat = interpolgraph(doLW ? xlw(fstat) : fstat, n, type);
   if (fsyst) fsyst = interpolgraph(doLW ? xlw(fsyst) : fsyst, n, type);
   if (ftot) ftot = interpolgraph(doLW ? xlw(ftot) : ftot, n, type);
}

#endif // ifndef dataset_C
