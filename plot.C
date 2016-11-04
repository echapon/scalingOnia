#ifndef plot_C
#define plot_C

#include "include/dataset.h"
#include "include/utils.h"

#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TH1F.h"
#include "Rtypes.h"

#include <vector>
#include <iostream>

using namespace std;

// global settings
const Escaling       gscaling = mtpt;
const Einterpolation ginterpolation = loglin;
float                gTextSize = 0.04;
// should we use the Lafferty & Wyatt prescription to change the x position of the points? (assuming a locally exponentially falling spectrum)
bool                 doxLW = true;
bool                 plotxt = true;

// declarations
int   mycolor(int i);
TH1F* haxes(TGraphAsymmErrors *graph, dataset data, float textSize, int ndiv=503, bool log=false);
void  plot(vector<dataset> data, vector<dataset> theory);

// main function
void plot(vector<dataset> data, vector<dataset> theory) {
   if (data.size()==0) return;
   if (theory.size()!=0 && theory.size()!=data.size()) {
      cout << "Error, inconsistent data and theory vector sizes." << endl;
      return;
   }

   float lTextSize = gTextSize;
   float lTextSize2 = gTextSize;
   float lTextSize3 = gTextSize;
   int npads=1;
   if (data.size()>1) npads=2;
   if (data.size()>1 && plotxt) npads=3;

   TCanvas *c1 = new TCanvas("c1","c1",600,600);
   c1->cd();
   float ysep = (npads>=2)*0.3*0.96 + (npads==3)*0.15;
	TPad *pad1 = new TPad("pad1","pad1",0,ysep,1,1);
	c1->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,(npads==3)*0.3*0.96,1, ysep);
	TPad *pad3 = new TPad("pad3","pad3",0,0,1,0.3*0.96);
   if (npads>=2)
   {
      pad1->SetBottomMargin(0.); // 0.04
      pad2->SetTopMargin(0);
      pad3->SetTopMargin(0);
      pad2->SetFillColor(0);
      pad2->SetFillStyle(0);
      pad3->SetFillColor(0);
      pad3->SetFillStyle(0);
      if (npads==2) {
         pad2->SetBottomMargin(gStyle->GetPadBottomMargin()/0.3);
      } else {
         pad2->SetBottomMargin(0);
         pad3->SetBottomMargin(gStyle->GetPadBottomMargin()/0.3);
         pad3->Draw();
      }
      pad1->SetTopMargin(gStyle->GetPadTopMargin()/0.7);
      // pad2->SetGridy();
      pad2->Draw();
      pad1->Draw();
      pad1->cd();

      // adapt text size
      lTextSize *= 1./(1.-ysep);
      lTextSize2 *= 1./(ysep-(npads==3)*0.3*0.96);
      lTextSize3 *= 1./(0.3*0.96);
   } else {
      pad1 = (TPad*) c1;
   }

   pad1->cd();
   pad1->SetLogy();

   TLegend *tleg = new TLegend(0.4,0.5,0.9,0.9);
   tleg->SetBorderSize(0);

   // loop on datasets...
   bool isfirst=true;
   for (unsigned int i=0; i<data.size(); i++) {
      pad1->cd();
      TGraphAsymmErrors *gstat = data[i].get_graphstat();
      if (gstat) {
         if (doxLW) gstat = xlw(gstat);
         if (plotxt) gstat = xtgraph(gstat, data[i].get_sqrts(), gscaling);

         // gstat->SetMarkerStyle(20+i);
         gstat->SetMarkerStyle(kBlack);
         gstat->SetMarkerColor(mycolor(i));
         gstat->SetLineColor(mycolor(i));
         if (i==0) {
            // draw axes
            TH1F *axes = haxes(gstat,data[i],lTextSize,503,true);
            axes->Draw();
         }
         gstat->Draw("P");

         tleg->AddEntry(gstat,data[i].get_legend().c_str(),"LP");
      }

      TGraphAsymmErrors *gtot  = data[i].get_graphtot();
      if (gtot) {
         if (doxLW) gtot = xlw(gtot);
         if (plotxt) gtot = xtgraph(gtot, data[i].get_sqrts(), gscaling);

         gtot->SetMarkerStyle(20+i);
         gtot->SetMarkerColor(mycolor(i));
         gtot->SetLineColor(mycolor(i));
         if (!gstat && i==0) {
            // draw axes
            TH1F *axes = haxes(gtot,data[i],lTextSize,503,true);
            axes->Draw();
         }
         gtot->Draw("P");

         if (!gstat) tleg->AddEntry(gstat,data[i].get_legend().c_str(),"LP");
      }

      if (i>0) {
         // make ratios
         pad2->cd();
         TGraphAsymmErrors *g0stat = data[0].get_graphstat();
         TGraphAsymmErrors *g0tot = data[0].get_graphtot();
         TGraphAsymmErrors *gratiostat=NULL, *gratiotot=NULL;
         if (gstat && g0stat) {
            if (doxLW) g0stat = xlw(g0stat);
            if (plotxt) g0stat = xtgraph(g0stat, data[0].get_sqrts(), gscaling);
            gratiostat = ratiograph(gstat,g0stat,ginterpolation);
            gratiostat->SetMarkerStyle(20+i);
            gratiostat->SetMarkerColor(mycolor(i));
            gratiostat->SetLineColor(mycolor(i));
            if (i==1) {
               // draw axes
               TH1F *axes = haxes(gratiostat,data[i],lTextSize2,505);
               axes->GetYaxis()->SetTitle("ratio");
               axes->Draw();
            }
            gratiostat->Draw("P");
         }
         if (gtot && g0tot) {
            if (doxLW) g0tot = xlw(g0tot);
            if (plotxt) g0tot = xtgraph(g0tot, data[0].get_sqrts(), gscaling);
            gratiotot = ratiograph(gtot,g0tot,ginterpolation);
            gratiotot->SetMarkerStyle(20+i);
            gratiotot->SetMarkerColor(mycolor(i));
            gratiotot->SetLineColor(mycolor(i));
            if (i==1 && !gratiostat) {
               // draw axes
               TH1F *axes = haxes(gratiotot,data[i],lTextSize2,505);
               axes->GetYaxis()->SetTitle("ratio");
               axes->Draw();
            }
            gratiotot->Draw("P");
         }

         if (plotxt) {
            // graph of n
            pad3->cd();
            TGraphAsymmErrors *gnstat=NULL, *gntot=NULL;
            if (gstat && g0stat) {
               gnstat = ngraph(gstat,g0stat,data[i].get_sqrts(),data[0].get_sqrts(),ginterpolation);
               gnstat->SetMarkerStyle(20+i);
               gnstat->SetMarkerColor(mycolor(i));
               gnstat->SetLineColor(mycolor(i));
               if (i==1) {
                  // draw axes
                  TH1F *axes = haxes(gnstat,data[i],lTextSize3,505);
                  axes->GetYaxis()->SetTitle("n_{eff}");
                  axes->Draw();
               }
               gnstat->Draw("P");
            }
            if (gtot && g0tot) {
               gntot = ngraph(gtot,g0tot,data[i].get_sqrts(),data[0].get_sqrts(),ginterpolation);
               gntot->SetMarkerStyle(20+i);
               gntot->SetMarkerColor(mycolor(i));
               gntot->SetLineColor(mycolor(i));
               if (i==1 && !gnstat) {
                  // draw axes
                  TH1F *axes = haxes(gntot,data[i],lTextSize3,505);
                  axes->GetYaxis()->SetTitle("n_{eff}");
                  axes->Draw();
               }
            }
         } // if (plotxt)
      } // if (i>1) 

      // take care of theory here
      // [not implemented yet] //FIXME

   } // for (int i=0; i<data.size(); i++) 

   pad1->cd();
   tleg->Draw();
}

// auxiliary functions
int mycolor(int i) {
   if (i==0) return kBlue+2;
   else if (i==1) return kRed+2;
   else if (i==2) return kGreen+2;
   else if (i==3) return kMagenta+2;
   else if (i==4) return kYellow+2;
   else if (i==5) return kCyan+2;
   else if (i==6) return kViolet+2;
   else if (i==7) return kOrange+2;
   else if (i==8) return kTeal+2;
   else if (i==9) return kPink+2;
   else if (i==10) return kSpring+2;
   else if (i==11) return kAzure+2;
   else return kBlack;
}

TH1F* haxes(TGraphAsymmErrors *graph, dataset data, float textSize, int ndiv, bool log) {
   if (!graph || graph->GetN()==0) return new TH1F("haxes","haxes",1,0,1);

   float xmin=0, xmax=-1e99, ymin= log ? 1e99 : 0, ymax=-1e99;
   for (int i=0; i<graph->GetN(); i++) {
      float x = graph->GetX()[i];
      float y = graph->GetY()[i];
      float ylow = y - graph->GetEYlow()[i];
      float xhigh = x + graph->GetEXhigh()[i];
      float yhigh = y + graph->GetEYhigh()[i];
      if (xhigh>xmax) xmax=xhigh;
      if (yhigh>ymax) ymax=yhigh;
      if (log && ylow<ymin) {
         if (ylow>0) ymin=ylow;
         else if (y>0) ymin=y;
         else ymin=1e-10;
      }
   }

   if (log) {
      ymin=0.5*ymin;
      ymax=1.5*ymax;
   } else {
      ymin=0.8*ymin;
      ymax=1.2*ymax;
   }

   TH1F *ans = new TH1F(Form("haxes_%f_%f%f_%f%f",textSize,xmin,xmax,ymin,ymax),"haxes",1,xmin,xmax);
   ans->SetMinimum(ymin);
   ans->SetMaximum(ymin);
   ans->GetXaxis()->SetTitle(plotxt ? "x_{T}" : data.get_xheader().c_str());
   ans->GetYaxis()->SetTitle(data.get_yheader().c_str());
   ans->GetXaxis()->SetLabelSize(textSize);
   ans->GetXaxis()->SetTitleSize(textSize);
   ans->GetYaxis()->SetLabelSize(textSize);
   ans->GetYaxis()->SetTitleSize(textSize);
   ans->GetYaxis()->SetNdivisions(ndiv);
   ans->GetYaxis()->SetTitleOffset(1.7*(gTextSize/textSize));
   ans->GetYaxis()->SetRangeUser(ymin,ymax);

   return ans;
}

#endif // ifndef plot_C
