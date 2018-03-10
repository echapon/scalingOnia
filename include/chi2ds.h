#ifndef chi2ds_h
#define chi2ds_h

#include "dataset.h"
#include "parameters.h"
#include "TGraphAsymmErrors.h"
#include <vector>

using namespace std;

class chi2ds {
   private: 
      vector<dataset> fdata;
      vector< vector<dataset> > ftheory;
      vector<TGraphAsymmErrors*> fdatan;
      vector< vector<TGraphAsymmErrors*> > ftheoryn;
      bool ffitxsec;
      bool fchecked;
      bool fok;
      Einterpolation ftype;

   public:
      chi2ds(): 
         fdata(), ftheory(), 
         fdatan(), ftheoryn(), 
         ffitxsec(false), fchecked(false), fok(false),
         ftype(cspline) {};
      ~chi2ds() {};

      // setters
      void setdata(vector<dataset> data) {fdata = data;};
      void addtheory(vector<dataset> theory) {ftheory.push_back(theory);};
      void setmode(bool fitxsec) {ffitxsec = fitxsec;};
      void setinterpolation(Einterpolation type) {ftype = type;};

      // check all OK + create neff graphs if needed
      void sanitycheck() {
         bool fchecked=true;

         if (fdata.size()==0) {fok=false; return;}
         if (ftheory.size()==0) {fok=false; return;}
         for (int i=0; i<ftheory.size(); i++) {
            if (fdata.size()!=ftheory[i].size()) {fok=false; return;}
         }

         // create the neff graphs if needed
         if (!ffitxsec) {
            if (fdata.size()<=1) {fok=false; return;}
            for (int i=0; i<ftheory.size(); i++) if (ftheory[i].size()<=1) {fok=false; return;}
            for (int i=0; i<ftheoryn.size(); i++) {
               if (fdatan.size()!=ftheoryn[i].size()) {fok=false; return;}
            }

            if (fdatan.size()==0) {
               // data
               TGraphAsymmErrors *g0tot = fdata[0].get_graphtot();
               for (int i=1; i<fdata.size(); i++) {
                  fdatan.push_back(ngraph(g0tot,fdata[i].get_graphtot(),fdata[i].get_sqrts(),fdata[0].get_sqrts(),ginterpolation));
               } // data loop

               // theory
               for (int j=0; j<ftheory.size(); j++) {
                  ftheoryn.push_back(vector<TGraphAsymmErrors*>());
                  TGraphAsymmErrors *g0tot = ftheory[j][0].get_graphtot();
                  for (int i=1; i<ftheory[j].size(); i++) {
                     ftheoryn[j].push_back(ngraph(g0tot,ftheory[j][i].get_graphtot(),ftheory[j][i].get_sqrts(),ftheory[j][0].get_sqrts(),ginterpolation));
                  } // theory loop
               } // theory loop 2
            } // if fdatan.size() == 0
         } // create ngraphs if needed

         fok = true;
      };

      // the main chi2 function
      double chi2(const double *xx) {
         if (!fchecked) sanitycheck();
         if (!fok) return -1;

         double chi2=0.;

         if (ffitxsec) {
            for (int i=0; i<fdata.size(); i++) {
               TGraphAsymmErrors *gd = fdata[i].get_graphtot();
               vector<TGraphAsymmErrors*> gt;
               vector<TGraph*> gth;
               vector<TGraph*> gtl;
               for (int k=0; k<ftheory.size(); k++) {
                  gt.push_back(ftheory[k][i].get_graphtot());
                  gth.push_back(graphhigh(gt.back()));
                  gtl.push_back(graphlow(gt.back()));
                  if (ftype==loglin || ftype==logcspline) {
                     gt.back() = lngraph(gt.back());
                     gth.back() = lngraph(gth.back());
                     gtl.back() = lngraph(gtl.back());
                  }
               } // k in theory.size()

               for (int j=0; j<gd->GetN(); j++) {
                  double x = gd->GetX()[j];
                  double num = gd->GetY()[j];
                  double denl = pow(gd->GetEYlow()[j],2);
                  double denh = pow(gd->GetEYhigh()[j],2);
                  for (int k=0; k<gt.size(); k++) {
                     double y0,yl,yh;
                     if (ftype==lin) {
                       y0 = gt[k]->Eval(x);
                       yl = gtl[k]->Eval(x);
                       yh = gth[k]->Eval(x);
                     } else if (ftype==cspline) {
                       y0 = gt[k]->Eval(x,0,"S");
                       yl = gtl[k]->Eval(x,0,"S");
                       yh = gth[k]->Eval(x,0,"S");
                     } else if (ftype==loglin) {
                       y0 = gt[k]->Eval(log(x));
                       yl = gtl[k]->Eval(log(x));
                       yh = gth[k]->Eval(log(x));
                     } else if (ftype==logcspline) {
                       y0 = gt[k]->Eval(log(x),0,"S");
                       yl = gtl[k]->Eval(log(x),0,"S");
                       yh = gth[k]->Eval(log(x),0,"S");
                     }
                     y0 *= xx[k];
                     yl *= xx[k];
                     yh *= xx[k];
                     num += -y0;
                     denl += pow(y0-yl,2);
                     denh += pow(yh-y0,2);
                  }
                  chi2 += pow(num,2)/((denl+denh)/2.);
               } // bin loop
            } // data size loop
         } else { // if (ffitxsec)
            for (int i=0; i<fdatan.size(); i++) {
               TGraphAsymmErrors *gd = fdatan[i];
               vector<TGraph*> gt;
               vector<TGraph*> gth;
               vector<TGraph*> gtl;
               for (int k=0; k<ftheory.size(); k++) {
                  gt.push_back(ftheoryn[k][i]);
                  gth.push_back(graphhigh(ftheoryn[k][i]));
                  gtl.push_back(graphlow(ftheoryn[k][i]));
                  if (ftype==loglin || ftype==logcspline) {
                     gt.back() = lngraph(gt.back());
                     gth.back() = lngraph(gth.back());
                     gtl.back() = lngraph(gtl.back());
                  }
               } // k in theory.size()

               for (int j=0; j<gd->GetN(); j++) {
                  double x = gd->GetX()[j];
                  double num = gd->GetY()[j];
                  double denl = pow(gd->GetEYlow()[j],2);
                  double denh = pow(gd->GetEYhigh()[j],2);
                  for (int k=0; k<gt.size(); k++) {
                     double a=1.;
                     // the xx are fractions, constrain the sum to 1
                     if (k<gt.size()-1) a=xx[k];
                     else for (int k2=0; k2<gt.size()-1; k2++) a+=-xx[k2];

                     double y0,yl,yh;
                     if (ftype==lin) {
                       y0 = a*gt[k]->Eval(x);
                       yl = a*gtl[k]->Eval(x);
                       yh = a*gth[k]->Eval(x);
                     } else if (ftype==cspline) {
                       y0 = a*gt[k]->Eval(x,0,"S");
                       yl = a*gtl[k]->Eval(x,0,"S");
                       yh = a*gth[k]->Eval(x,0,"S");
                     } else if (ftype==loglin) {
                       y0 = a*gt[k]->Eval(log(x));
                       yl = a*gtl[k]->Eval(log(x));
                       yh = a*gth[k]->Eval(log(x));
                     } else if (ftype==logcspline) {
                       y0 = a*gt[k]->Eval(log(x),0,"S");
                       yl = a*gtl[k]->Eval(log(x),0,"S");
                       yh = a*gth[k]->Eval(log(x),0,"S");
                     }
                     num += -y0;
                     denl += pow(y0-yl,2);
                     denh += pow(yh-y0,2);
                  }
                  chi2 += pow(num,2)/((denl+denh)/2.);
               } // bin loop
            } // data size loop
         }

         return chi2;
      };

      // number of degrees of freedom
      int ndf() {
         if (!fchecked) sanitycheck();
         if (!fok) return -1;

         int ans=0;

         if (ffitxsec) {
            for (int i=0; i<fdata.size(); i++) ans+=fdata[i].get_graphtot()->GetN();
            ans += ftheory.size();
         } else {
            for (int i=0; i<fdatan.size(); i++) ans+=fdatan[i]->GetN();
            ans += ftheory.size()-1;
         }

         return ans;
      }
};

#endif // #ifndef chi2ds_h
