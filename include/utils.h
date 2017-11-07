#ifndef range_h
#define range_h

#include "TGraphAsymmErrors.h"
#include <math.h>
#include <fstream>

using namespace std;

class range {
   public:
      double min;
      double max;

      range(): min(0), max(0) {};
      range(double mmin, double mmax): min(mmin), max(mmax) {};
};

typedef enum Escaling {
   pt2,
   mtpt
} Escaling;

typedef enum Einterpolation {
   lin,
   cspline,
   loglin,
   logcspline
} Einterpolation;

typedef enum Elwmode {
   expo,
   powlaw
} Elwmode;

const double mjpsi = 3.096916;
const double mjpsi2 = pow(mjpsi,2);

TGraphAsymmErrors *combgraph(TGraphAsymmErrors *gstat, TGraphAsymmErrors *gsyst) {
   // Combine two graphs with same values but different uncertainties.
   // These uncertainties will be added in quadrature.

   if (gstat->GetN() != gsyst->GetN()) return NULL;
   TGraphAsymmErrors *ans = new TGraphAsymmErrors(gstat->GetN());
   ans->SetName(TString(gstat->GetName()) + "_statsyst");
   ans->SetTitle(TString(gstat->GetTitle()) + " (stat+syst)");

   for (int i=0; i<ans->GetN(); i++) {
      ans->SetPoint(i, gstat->GetX()[i], gstat->GetY()[i]);
      ans->SetPointError(i, gstat->GetEXlow()[i], gstat->GetEXhigh()[i],
            sqrt(pow(gstat->GetEYlow()[i],2)+pow(gsyst->GetEYlow()[i],2)),
            sqrt(pow(gstat->GetEYhigh()[i],2)+pow(gsyst->GetEYhigh()[i],2)));
   }
   return ans;
};

TGraph *lngraph(TGraph* g) {
   // take the ln of graph coordinates

   TGraph *ans = new TGraph(g->GetN());
   ans->SetName(TString(g->GetName()) + "_ln");
   ans->SetTitle(TString(g->GetTitle()) + " (log scale)");
   for (int i=0; i<g->GetN(); i++) {
      double lnx = g->GetX()[i]>0 ? log(g->GetX()[i]) : 0;
      double lny = g->GetY()[i]>0 ? log(g->GetY()[i]) : 0;
      ans->SetPoint(i,lnx,lny);
   }
   return ans;
};

TGraphAsymmErrors *lngraph(TGraphAsymmErrors* g) {
   // take the ln of graph coordinates

   TGraphAsymmErrors *ans = new TGraphAsymmErrors(g->GetN());
   ans->SetName(TString(g->GetName()) + "_ln");
   ans->SetTitle(TString(g->GetTitle()) + " (log scale)");
   for (int i=0; i<g->GetN(); i++) {
      double lnx = g->GetX()[i]>0 ? log(g->GetX()[i]) : 0;
      double lny = g->GetY()[i]>0 ? log(g->GetY()[i]) : 0;
      double lnexl = g->GetEXlow()[i]!=0 ? log(fabs(g->GetEXlow()[i])) : 0;
      double lnexh = g->GetEXhigh()[i]!=0 ? log(fabs(g->GetEXhigh()[i])) : 0;
      double lneyl = g->GetEYlow()[i]!=0 ? log(fabs(g->GetEYlow()[i])) : 0;
      double lneyh = g->GetEYhigh()[i]!=0 ? log(fabs(g->GetEYhigh()[i])) : 0;
      ans->SetPoint(i,lnx,lny);
      ans->SetPointError(i,lnexl,lnexh,lneyl,lneyh);
   }
   return ans;
};

TGraphAsymmErrors *xtgraph(TGraphAsymmErrors* g, double sqrts, Escaling type=pt2) {
   // take the ln of graph coordinates

   TGraphAsymmErrors *ans = new TGraphAsymmErrors(g->GetN());
   ans->SetName(TString(g->GetName()) + "_scal");
   ans->SetTitle(TString(g->GetTitle()) + " (xt scaled)");
   for (int i=0; i<g->GetN(); i++) {
      double x = g->GetX()[i];
      double y = g->GetY()[i];
      double exl = g->GetEXlow()[i];
      double exh = g->GetEXhigh()[i];
      double eyl = g->GetEYlow()[i];
      double eyh = g->GetEYhigh()[i];
      if (type==pt2) {
         double newx = 2*x/sqrts;
         double tmpx = x-exl;
         exl = -2*tmpx/sqrts+newx;
         tmpx = x+exh;
         exh = 2*tmpx/sqrts-newx;
         x = newx;
      } else if (type==mtpt) {
         double newx = (x+sqrt(mjpsi2+x*x))/sqrts;
         double tmpx = x-exl;
         exl = -(tmpx+sqrt(mjpsi2+tmpx*tmpx))/sqrts + newx;
         tmpx = x+exh;
         exh = (tmpx+sqrt(mjpsi2+tmpx*tmpx))/sqrts - newx;
         x = newx;
      }
      ans->SetPoint(i,x,y);
      ans->SetPointError(i,exl,exh,eyl,eyh);
   }
   return ans;
};

TGraph *graphlow(TGraphAsymmErrors *g) {
   // return the graph of y-dy

   TGraph *ans = new TGraph(g->GetN());
   ans->SetName(TString(g->GetName()) + "_low");
   for (int i=0; i<g->GetN(); i++) {
      ans->SetPoint(i,g->GetX()[i],g->GetY()[i]-g->GetEYlow()[i]);
   }
   return ans;
};

TGraph *graphhigh(TGraphAsymmErrors *g) {
   // return the graph of y+dy

   TGraph *ans = new TGraph(g->GetN());
   ans->SetName(TString(g->GetName()) + "_high");
   for (int i=0; i<g->GetN(); i++) {
      ans->SetPoint(i,g->GetX()[i],g->GetY()[i]+g->GetEYhigh()[i]);
   }
   return ans;
};

TGraphAsymmErrors *ngraph(TGraphAsymmErrors* g1, TGraphAsymmErrors *g2, double sqrts1, double sqrts2, Einterpolation type=loglin, bool correl=false) {
   // derive the n exponent for each pt

   TGraphAsymmErrors *ans = new TGraphAsymmErrors(g1->GetN());
   ans->SetName(TString(g1->GetName()) + "_" + TString(g2->GetName()) + "_n");
   ans->SetTitle(TString("n exponent of (") + TString(g1->GetTitle()) + ") vs (" + TString(g2->GetTitle()) + ")");

   TGraph *g1low = graphlow(g1);
   TGraph *g1high = graphhigh(g1);
   TGraph *g2low = correl ? graphlow(g2) : graphhigh(g2);
   TGraph *g2high = correl ? graphhigh(g2) : graphlow(g2);

   TGraphAsymmErrors *lng1=NULL, *lng2=NULL;
   TGraph *lng1low=NULL, *lng1high=NULL, *lng2low=NULL, *lng2high=NULL;
   if (type==loglin || type==logcspline) {
      lng1 = lngraph(g1);
      lng2 = lngraph(g2);
      lng1low = lngraph(g1low);
      lng2low = lngraph(g2low);
      lng1high = lngraph(g1high);
      lng2high = lngraph(g2high);
   }

   for (int i=0; i<g1->GetN(); i++) {
      double x = g1->GetX()[i];
      double exl = g1->GetEXlow()[i];
      double exh = g1->GetEXhigh()[i];
      double y=0,eyl=0,eyh=0,ya=0,yb=0;
      if (type==lin) {
         y = -(log(g2->Eval(x))-log(g1->Eval(x))) / (log(sqrts2) - log(sqrts1));
         ya = -(log(g2low->Eval(x))-log(g1low->Eval(x))) / (log(sqrts2) - log(sqrts1));
         yb = -(log(g2high->Eval(x))-log(g1high->Eval(x))) / (log(sqrts2) - log(sqrts1));
      } else if (type==cspline) {
         y = -(log(g2->Eval(x,0,"S"))-log(g1->Eval(x,0,"S"))) / (log(sqrts2) - log(sqrts1));
         ya = -(log(g2low->Eval(x,0,"S"))-log(g1low->Eval(x,0,"S"))) / (log(sqrts2) - log(sqrts1));
         yb = -(log(g2high->Eval(x,0,"S"))-log(g1high->Eval(x,0,"S"))) / (log(sqrts2) - log(sqrts1));
      } else if (type==loglin) {
         y = -(lng2->Eval(log(x))-lng1->Eval(log(x))) / (log(sqrts2) - log(sqrts1));
         ya = -(lng2low->Eval(log(x))-lng1low->Eval(log(x))) / (log(sqrts2) - log(sqrts1));
         yb = -(lng2high->Eval(log(x))-lng1high->Eval(log(x))) / (log(sqrts2) - log(sqrts1));
      } else if (type==logcspline) {
         y = -(lng2->Eval(log(x),0,"S")-lng1->Eval(log(x),0,"S")) / (log(sqrts2) - log(sqrts1));
         ya = -(lng2low->Eval(log(x),0,"S")-lng1low->Eval(log(x),0,"S")) / (log(sqrts2) - log(sqrts1));
         yb = -(lng2high->Eval(log(x),0,"S")-lng1high->Eval(log(x),0,"S")) / (log(sqrts2) - log(sqrts1));
      }
      // add one unit because we use dsigma/dpt instead of dsigma/ptdpt :)
      y+=1.;
      ya+=1.;
      yb+=1.;
      eyl = y - min(ya,yb);
      eyh = max(ya,yb) - y;
      ans->SetPoint(i,x,y);
      // cout << i << " " << x << " " << y << " " << exl << " " << exh << " " << eyl << " " << eyh << endl;
      ans->SetPointError(i,exl,exh,eyl,eyh);
   }

   if (lng1) {delete lng1; lng1=0;}
   if (lng2) {delete lng2; lng2=0;}
   if (lng1low) {delete lng1low; lng1low=0;}
   if (lng2low) {delete lng2low; lng2low=0;}
   if (lng1high) {delete lng1high; lng1high=0;}
   if (lng2high) {delete lng2high; lng2high=0;}

   return ans;
};

TGraphAsymmErrors *ratiograph(TGraphAsymmErrors* gnum, TGraphAsymmErrors *gden, Einterpolation type=loglin, bool correl=false) {
   // take the ratio of two graphs

   TGraph *gnumlow = graphlow(gnum);
   TGraph *gnumhigh = graphhigh(gnum);
   TGraph *gdenlow = correl ? graphlow(gden) : graphhigh(gden);
   TGraph *gdenhigh = correl ? graphhigh(gden) : graphlow(gden);

   TGraphAsymmErrors *lngnum=NULL, *lngden=NULL;
   TGraph *lngnumlow=NULL, *lngnumhigh=NULL, *lngdenlow=NULL, *lngdenhigh=NULL;
   if (type==loglin || type==logcspline) {
      lngnum = lngraph(gnum);
      lngden = lngraph(gden);
      lngnumlow = lngraph(gnumlow);
      lngdenlow = lngraph(gdenlow);
      lngnumhigh = lngraph(gnumhigh);
      lngdenhigh = lngraph(gdenhigh);
   }

   // find overlap in x between the numerator and denominator. We will only compute the ratio where they overlap.
   int inumfirst=-1, inumlast=-1;
   int idenfirst=0, idenlast=-1;
   int iden=0;
   for (int i=0; i<gnum->GetN(); i++) {
      double xnum = gnum->GetX()[i];
      for (int iden=idenfirst; iden<gden->GetN()-1; iden++) {
         double xlden = gden->GetX()[iden];
         double xhden = gden->GetX()[iden+1];
         if (xlden<xnum && xhden>xnum) {
            if (inumfirst<0) {
               inumfirst = i;
               idenfirst = iden;
            } else {
               inumlast = i;
               idenlast = iden;
            }
            break;
         }
      }
   }
   // cout << inumfirst << " " << idenfirst << ", " << inumlast << " " << idenlast << endl;

   TGraphAsymmErrors *ans = new TGraphAsymmErrors(inumlast-inumfirst+1);
   ans->SetName(TString(gnum->GetName()) + "_" + TString(gden->GetName()) + "_ratio");
   ans->SetTitle(TString("ratio of (") + TString(gnum->GetTitle()) + ") / (" + TString(gden->GetTitle()) + ")");


   for (int i=inumfirst; i<=inumlast; i++) {
      double x = gnum->GetX()[i];
      double exl = gnum->GetEXlow()[i];
      double exh = gnum->GetEXhigh()[i];
      double y=0,eyl=0,eyh=0,ya=0,yb=0;
      if (type==lin || (type==loglin && x<=0)) {
         y = gnum->Eval(x)/gden->Eval(x);
         // if (y<0) cout << x << ", " << y << " = " <<  gnum->Eval(x) << "/" << gden->Eval(x) << endl;
         ya = gnumhigh->Eval(x)/gdenlow->Eval(x);
         yb = gnumlow->Eval(x)/gdenhigh->Eval(x);
      } else if (type==cspline || (type==logcspline && x<=0)) {
         y = gnum->Eval(x,0,"S")/gden->Eval(x,0,"S");
         // if (y<0) cout << x << ", " << y << " = " << gnum->Eval(x,0,"S") << "/" << gden->Eval(x,0,"S") << endl;
         ya = gnumhigh->Eval(x,0,"S")/gdenlow->Eval(x,0,"S");
         yb = gnumlow->Eval(x,0,"S")/gdenhigh->Eval(x,0,"S");
      } else if (type==loglin) {
         y = exp(lngnum->Eval(log(x)))/exp(lngden->Eval(log(x),0,"S"));
         // if (y<0) cout << x << ", " << y << " = " << exp(lngnum->Eval(log(x))) << "/" << exp(lngden->Eval(log(x),0,"S")) << endl;
         ya = exp(lngnumhigh->Eval(log(x)))/exp(lngdenlow->Eval(log(x),0,"S"));
         yb = exp(lngnumlow->Eval(log(x)))/exp(lngdenhigh->Eval(log(x),0,"S"));
      } else if (type==logcspline) {
         y = exp(lngnum->Eval(log(x),0,"S"))/exp(lngden->Eval(log(x),0,"S"));
         // if (y<0) cout << x << ", " << y  << " = " << exp(lngnum->Eval(log(x),0,"S")) << "/" << exp(lngden->Eval(log(x),0,"S")) << endl;
         ya = exp(lngnumhigh->Eval(log(x),0,"S"))/exp(lngdenlow->Eval(log(x),0,"S"));
         yb = exp(lngnumlow->Eval(log(x),0,"S"))/exp(lngdenhigh->Eval(log(x),0,"S"));
      }
      // cout << y << " " << y-ya << " " << y-yb << endl;
      eyl = y - min(ya,yb);
      eyh = max(ya,yb) - y;
      int iout = i-inumfirst;
      ans->SetPoint(iout,x,y);
      ans->SetPointError(iout,exl,exh,eyl,eyh);
   }

   if (lngnum) {delete lngnum; lngnum=0;}
   if (lngden) {delete lngden; lngden=0;}
   if (lngnumlow) {delete lngnumlow; lngnumlow=0;}
   if (lngdenlow) {delete lngdenlow; lngdenlow=0;}
   if (lngnumhigh) {delete lngnumhigh; lngnumhigh=0;}
   if (lngdenhigh) {delete lngdenhigh; lngdenhigh=0;}

   return ans;
};

TGraphAsymmErrors *xlw(TGraphAsymmErrors *g, Elwmode mode=powlaw) {
   // // change the x of each bin, according to the Lafferty and Wyatt prescription (http://www.sciencedirect.com/science/article/pii/0168900294011125)

   TGraphAsymmErrors *ans = new TGraphAsymmErrors(g->GetN());
   ans->SetName(TString(g->GetName()) + "_xlw");
   ans->SetTitle(TString(g->GetTitle()) + " with x_{LW}");
   TGraphAsymmErrors *tmp = (TGraphAsymmErrors*) g->Clone("tmp");

   // we iterate 5 times, so that we are sure that we get the correct spectrum
   for (int j=0; j<5; j++) {
      for (int i=0; i<tmp->GetN(); i++) {
         double x = tmp->GetX()[i];
         double exl = tmp->GetEXlow()[i];
         double exh = tmp->GetEXhigh()[i];
         double y = tmp->GetY()[i];
         double eyl = tmp->GetEYlow()[i];
         double eyh = tmp->GetEYhigh()[i];

         // Lafferty and Wyatt: f(xLW) = \int f(x)dx
         // here we assume an exponentially falling function in each bin
         // first let's get the slope. 
         double b;
         if (i>0) {
            if (mode==expo) b = -log(tmp->GetY()[i]/tmp->GetY()[i-1])/(tmp->GetX()[i]-tmp->GetX()[i-1]);
            else b = log(tmp->GetY()[i]/tmp->GetY()[i-1])/log(tmp->GetX()[i]/tmp->GetX()[i-1]);
         } else {
            if (mode==expo) b = -log(tmp->GetY()[i+1]/tmp->GetY()[i])/(tmp->GetX()[i+1]-tmp->GetX()[i]);
            else b = log(tmp->GetY()[i+1]/tmp->GetY()[i])/log(tmp->GetX()[i+1]/tmp->GetX()[i]);
         }
         double xlw;
         if (mode==expo) xlw = (b>0 && exl+exh>0) ? (x-exl) + (1./b) * (log(b*(exl+exh)) - log(1-exp(-b*(exl+exh)))) : x;
         else {
            double intgx = (pow(x+exh,b+1)-pow(x-exl,b+1))/((b+1)*(exl+exh));
            xlw = pow(intgx,1./b);
            // cout << x-exl << " " << x+exh << " -> " << xlw << " (" << b << ", " << intgx << ")" << endl;
         }
         double exl_lw = exl - x + xlw;
         double exh_lw = exh + x - xlw;

         ans->SetPoint(i,xlw,y);
         ans->SetPointError(i,exl_lw,exh_lw,eyl,eyh);
      }
      delete tmp;
      tmp = (TGraphAsymmErrors*) ans->Clone("tmp");
   }

   if (tmp) {delete tmp; tmp=0;}

   return ans;
}

TGraphAsymmErrors *interpolgraph(TGraphAsymmErrors* gr, int n, Einterpolation type=loglin) {
   // interpolate a graph by making it with more points

   TGraphAsymmErrors *ans = new TGraphAsymmErrors(n);
   ans->SetName(TString(gr->GetName()) + "_interpol");
   ans->SetTitle(TString(gr->GetTitle()) + " with interpolation");

   TGraph *grlow = graphlow(gr);
   TGraph *grhigh = graphhigh(gr);

   TGraphAsymmErrors *lngr=NULL;
   TGraph *lngrlow=NULL, *lngrhigh=NULL;
   if (type==loglin || type==logcspline) {
      lngr = lngraph(gr);
      lngrlow = lngraph(grlow);
      lngrhigh = lngraph(grhigh);
   }

   double xmin = gr->GetX()[0]-gr->GetEXlow()[0];
   if (xmin<=0) xmin=1e-2;
   double xmax = gr->GetX()[gr->GetN()-1]+gr->GetEXhigh()[gr->GetN()-1];
   double dr = (xmax-xmin)/n;

   for (int i=0; i<n; i++) {
      double x = xmin+dr*i;
      double exl = 0;
      double exh = 0;
      double y=0,eyl=0,eyh=0,ya=0,yb=0;
      if (type==lin) {
         y = gr->Eval(x);
         ya = grhigh->Eval(x);
         yb = grlow->Eval(x);
      } else if (type==cspline) {
         y = gr->Eval(x,0,"S");
         ya = grhigh->Eval(x,0,"S");
         yb = grlow->Eval(x,0,"S");
      } else if (type==loglin) {
         y = exp(lngr->Eval(log(x)));
         ya = exp(lngrhigh->Eval(log(x)));
         yb = exp(lngrlow->Eval(log(x)));
      } else if (type==logcspline) {
         y = exp(lngr->Eval(log(x),0,"S"));
         ya = exp(lngrhigh->Eval(log(x),0,"S"));
         yb = exp(lngrlow->Eval(log(x),0,"S"));
      }
      eyl = y - min(ya,yb);
      eyh = max(ya,yb) - y;
      ans->SetPoint(i,x,y);
      // ans->SetPointError(i,exl,exh,eyl,eyh);
      ans->SetPointError(i,exl,exh,0,0);
   }

   if (lngr) {delete lngr; lngr=0;}
   if (lngrlow) {delete lngrlow; lngrlow=0;}
   if (lngrhigh) {delete lngrhigh; lngrhigh=0;}

   return ans;
};

TGraphAsymmErrors *combo(TGraphAsymmErrors* g0, vector<TGraphAsymmErrors*> gall, bool dominmax=true) {
   if (gall.size()==0) return NULL;
   TGraphAsymmErrors *ans = (TGraphAsymmErrors*) gall[0]->Clone(TString(gall[0]->GetName())+"_all");
   int n = ans->GetN();
   if (dominmax) { // min max
      double *thevalsmin = new double[n];
      double *thevalsmax = new double[n];
      for (int i=0; i<n; i++) {thevalsmin[i]=g0->GetY()[i]; thevalsmax[i]=g0->GetY()[i];}

      for (unsigned int i=0; i<gall.size(); i++) 
         for (int j=0; j<n; j++) {
            thevalsmin[j] = min(gall[i]->GetY()[j],thevalsmin[j]);
            thevalsmax[j] = max(gall[i]->GetY()[j],thevalsmax[j]);
         }
      for (int i=0; i<n; i++) ans->SetPoint(i,ans->GetX()[i],g0->GetY()[i]);
      for (int i=0; i<n; i++) ans->SetPointError(i,
            ans->GetEXlow()[i],ans->GetEXhigh()[i],
            g0->GetY()[i]-thevalsmin[i],thevalsmax[i]-g0->GetY()[i]);
   } else { // quadratic sum
      double *thevals = new double[n];
      for (int i=0; i<n; i++) thevals[i]=0;

      for (unsigned int i=0; i<gall.size(); i++) 
         for (int j=0; j<n; j++) {
            thevals[j] += pow(gall[i]->GetY()[j]-g0->GetY()[j],2);
         }
      for (int i=0; i<n; i++) ans->SetPoint(i,ans->GetX()[i],g0->GetY()[i]);
      for (int i=0; i<n; i++) ans->SetPointError(i,ans->GetEXlow()[i],ans->GetEXhigh()[i],sqrt(thevals[i]),sqrt(thevals[i]));
   }

   return ans;
};

TGraphAsymmErrors* ratiograph(TGraphAsymmErrors* g0num, TGraphAsymmErrors* g0den,
      vector<TGraphAsymmErrors*> gnum, vector<TGraphAsymmErrors*> gden, 
      bool dominmax=true, Einterpolation type=loglin, bool correl=false) {
   if (gnum.size() != gden.size()) {
      cout << "ratiograph: ERROR different vector sizes" << endl;
      return NULL;
   }
   vector<TGraphAsymmErrors*> gall;
   for (int i=0; i<gnum.size(); i++) {
      gall.push_back(ratiograph(gnum[i],gden[i],type,correl));
   }

   TGraphAsymmErrors *g0 = ratiograph(g0num,g0den,type,correl);

   return combo(g0, gall, dominmax);
};

void printgraph(ofstream &f, TGraphAsymmErrors* g) {
   if (!g) {
      cout << "Error, cannot print null graph" << endl;
      return;
   }
   f << "# " << g->GetName() << " ; " << g->GetTitle() << endl;
   f << "# x_low, x_high, x, y, dy_low, dy_up" << endl;

   double *x = g->GetX();
   double *y = g->GetY();
   double *exl = g->GetEXlow();
   double *exh = g->GetEXhigh();
   double *eyl = g->GetEYlow();
   double *eyh = g->GetEYhigh();
   for (int i=0; i<g->GetN(); i++) {
      f << x[i]-exl[i] << ", " << x[i]+exh[i] << ", " << x[i] << ", " << y[i] << ", " << eyl[i] << ", " << eyh[i] << endl;
   }
   f << endl;
}

void setUncert(TGraphAsymmErrors *g, double xerr=-1, double yerr=-1) {
   if (!g) return;

   for (int i=0; i<g->GetN(); i++) {
      double exl = (xerr>=0) ? xerr : g->GetEXlow()[i];
      double exh = (xerr>=0) ? xerr : g->GetEXhigh()[i];
      double eyl = (yerr>=0) ? yerr : g->GetEYlow()[i];
      double eyh = (yerr>=0) ? yerr : g->GetEYhigh()[i];
      g->SetPointError(i,exl,exh,eyl,eyh);
   }
}

#endif // ifndef range_h
