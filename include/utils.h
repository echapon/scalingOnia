#ifndef range_h
#define range_h

#include "TGraphAsymmErrors.h"
#include <math.h>

typedef struct range {
   double min;
   double max;
} range;

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
   ans->SetTitle(TString("n exponent of ") + TString(g1->GetTitle()) + " vs " + TString(g2->GetTitle()));

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
      eyl = y - min(ya,yb);
      eyh = max(ya,yb) - y;
      ans->SetPoint(i,x,y);
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

TGraphAsymmErrors *xlw(TGraphAsymmErrors *g) {
   // change the x of each bin, according to the Lafferty and Wyatt prescription (http://www.sciencedirect.com/science/article/pii/0168900294011125)

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
         if (i>0) b = -log(tmp->GetY()[i]/tmp->GetY()[i-1])/(tmp->GetX()[i]-tmp->GetX()[i-1]);
         else b = -log(tmp->GetY()[i+1]/tmp->GetY()[i])/(tmp->GetX()[i+1]-tmp->GetX()[i]);
         double xlw = (x-exl) + (1./b) * (log(b*(exl+exh)) - log(1-exp(-b*(exl+exh))));
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

#endif // ifndef range_h
