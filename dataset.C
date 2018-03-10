#ifndef dataset_C
#define dataset_C

#include "include/dataset.h"

#include <fstream>
#include "TFile.h"
#include "TH1.h"

using namespace std;

// dataset::dataset (const dataset &d, const char* suffix) {
//    fsqrts = d.fsqrts;
//    fraprange = d.fraprange;
//    if (d.fstat) fstat = (TGraphAsymmErrors*) d.fstat->Clone(Form("%s%s",d.fstat->GetName(),suffix));
//    if (d.fsyst) fsyst = (TGraphAsymmErrors*) d.fsyst->Clone(Form("%s%s",d.fsyst->GetName(),suffix));
//    for (int i=0; i<d.fsyst_all.size(); i++)
//       fsyst_all.push_back((TGraphAsymmErrors*) d.fsyst_all[i]->Clone(Form("%s%s",d.fsyst_all[i]->GetName(),suffix)));
//    if (d.ftot) ftot = (TGraphAsymmErrors*) d.ftot->Clone(Form("%s%s",d.ftot->GetName(),suffix));
//    fname = d.fname + "_copy";
//    fexpname = d.fexpname;
//    flegend = d.flegend;
//    flocation = d.flocation;
//    fcomment = d.fcomment;
//    fxheader = d.fxheader;
//    fyheader = d.fyheader;
//    freaction = d.freaction;
// }

dataset dataset::operator=(const dataset &d) {
   if (fstat) delete fstat;
   if (fsyst) delete fsyst;
   for (int i=0; i<fsyst_all.size(); i++) delete fsyst_all[i];
   fsyst_all.clear();
   if (ftot) delete ftot;

   const char* suffix="_copy";
   fsqrts = d.fsqrts;
   fraprange = d.fraprange;
   if (d.fstat) fstat = (TGraphAsymmErrors*) d.fstat->Clone(Form("%s%s",d.fstat->GetName(),suffix));
   if (d.fsyst) fsyst = (TGraphAsymmErrors*) d.fsyst->Clone(Form("%s%s",d.fsyst->GetName(),suffix));
   for (int i=0; i<d.fsyst_all.size(); i++)
      fsyst_all.push_back((TGraphAsymmErrors*) d.fsyst_all[i]->Clone(Form("%s%s",d.fsyst_all[i]->GetName(),suffix)));
   if (d.ftot) ftot = (TGraphAsymmErrors*) d.ftot->Clone(Form("%s%s",d.ftot->GetName(),suffix));
   fname = d.fname + "_copy";
   fexpname = d.fexpname;
   flegend = d.flegend;
   flocation = d.flocation;
   fcomment = d.fcomment;
   fxheader = d.fxheader;
   fyheader = d.fyheader;
   freaction = d.freaction;

   return *this;
}

void dataset::set_graph(TGraphAsymmErrors *gstat, TGraphAsymmErrors *gsyst) {
   if (gstat) fstat = gstat;
   if (gsyst) fsyst = gsyst;
   if (gstat && gsyst) ftot = combgraph(gstat,gsyst);
   else if (gsyst) ftot = gsyst;
   else ftot = gstat;
}

void dataset::set_graph(TGraphAsymmErrors *gstat, vector<TGraphAsymmErrors*> gsyst) {
   if (gstat) fstat = gstat;
   if (gsyst.size()>0) fsyst_all = gsyst;
   get_graphsyst(); // to compute the total syst;
   if (gstat && fsyst) ftot = combgraph(gstat,fsyst);
   else if (fsyst) ftot = fsyst;
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

void dataset::set_graphHwU(const char* file_theory, int nsysts, int ihist) {
   fstat = read_hwu(file_theory,ihist,0);
   for (int i=0; i<nsysts; i++) {
      fsyst_all.push_back(read_hwu(file_theory,ihist,i+4));
      TGraphAsymmErrors *glast = fsyst_all.back();
      for (int j=0; j<glast->GetN(); j++) {
         glast->SetPoint(j,glast->GetX()[j],glast->GetY()[j]);
         glast->SetPointEYhigh(j,fstat->GetEYhigh()[j]);
         glast->SetPointEYlow(j,fstat->GetEYlow()[j]);
      }
   }
   get_graphsyst(); // to compute the total syst
}

void dataset::interpolate(int n, Einterpolation type, bool doLW) {
   if (fstat) fstat = interpolgraph(doLW ? xlw(fstat) : fstat, n, type);
   if (fsyst) fsyst = interpolgraph(doLW ? xlw(fsyst) : fsyst, n, type);
   if (fsyst_all.size()>0) 
      for (int i=0; i<fsyst_all.size(); i++)
         fsyst_all[i] = interpolgraph(doLW ? xlw(fsyst_all[i]) : fsyst_all[i], n, type);
   if (ftot) ftot = interpolgraph(doLW ? xlw(ftot) : ftot, n, type);
}

TGraphAsymmErrors* dataset::get_graphsyst() {
   if (fsyst) return fsyst;
   if (fsyst_all.size()==0) return NULL;

   // we don't have a fsyst yet: compute it from the vector of systs
   fsyst = combo(fstat,fsyst_all,takeMaxSyst);

   // compute the graph for stat+syst
   if (fstat && fsyst) ftot = combgraph(fstat,fsyst);

   return fsyst;
}
dataset dataset::operator+(dataset d) {
   add(fstat,d.fstat);
   add(fsyst,d.fsyst);
   for (int i=0; i<fsyst_all.size(); i++) add(fsyst_all[i],d.fsyst_all[i]);
   add(ftot,d.ftot);

   return *this;
}

dataset dataset::operator*(double s) {
   scale(fstat,s);
   scale(fsyst,s);
   for (int i=0; i<fsyst_all.size(); i++) scale(fsyst_all[i],s);
   scale(ftot,s);

   return *this;
}

#endif // ifndef dataset_C
