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
#include "TString.h"

#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

// global settings
const Escaling       gscaling = pt2; // pt2 or mtpt
const Einterpolation ginterpolation = logcspline; // lin, cspline, loglin, logcspline
float                gTextSize = 0.04;
// should we use the Lafferty & Wyatt prescription to change the x position of the points? (assuming a locally exponentially falling spectrum)
bool                 doxLW = true;//true;
Elwmode              lwmode = powlaw; // expo or powlaw
bool                 plotxt = true;//true;

// declarations
int   mycolor(int i);
TH1F* haxes(TGraphAsymmErrors *graph, dataset data, float textSize, int ndiv=505, bool log=false, double xmin_=-1, double xmax_=-1);
void  plot(vector<dataset> data, vector<dataset> theory, TString cname="graphs");

// main function
void plot(vector<dataset> data, vector<dataset> theory, TString cname) {
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

   // grid style
   gStyle->SetGridStyle(3);
   gStyle->SetGridWidth(2);

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
         pad3->SetGridy();
         pad3->Draw();
      }
      pad1->SetTopMargin(gStyle->GetPadTopMargin()/0.7);
      pad2->SetGridy();
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

   ofstream flog(TString(cname+".dat").Data());

   double xmin=0, xmax=0;

   // loop on datasets...
   for (unsigned int i=0; i<data.size()+theory.size(); i++) {
      pad1->cd();
      int idata = i - theory.size();
      bool isth = (idata<0);
      dataset *di = (!isth) ? &(data[idata]) : &(theory[i]);
      TGraphAsymmErrors *gstat = di->get_graphstat();
      if (isth) setUncert(gstat,0);
      if (gstat) {
         if (doxLW && !isth) gstat = xlw(gstat,lwmode);
         if (plotxt) gstat = xtgraph(gstat, di->get_sqrts(), gscaling);

         // gstat->SetMarkerStyle(20+i);
         gstat->SetMarkerStyle(kBlack);
         gstat->SetMarkerColor(mycolor(i));
         if (isth) gstat->SetMarkerSize(0);
         gstat->SetLineColor(mycolor(i));
         if (i==0) {
            // draw axes
            TH1F *axes = haxes(gstat,*di,lTextSize,503,true);
            xmin = axes->GetXaxis()->GetXmin();
            xmax = axes->GetXaxis()->GetXmax();
            if (plotxt) axes->GetYaxis()->SetTitle("d#sigma/dp_{T} (x_{T}) [nb/GeV]");
            axes->Draw();
         }
         if (!isth) gstat->Draw("PL");
         else gstat->Draw("L");

         printgraph(flog,gstat);

         tleg->AddEntry(gstat,di->get_legend().c_str(),isth ? "L" : "LP");
      }

      TGraphAsymmErrors *gtot  = di->get_graphtot();
      if (isth) setUncert(gtot,0);
      if (gtot) {
         if (doxLW && !isth) gtot = xlw(gtot,lwmode);
         if (plotxt) gtot = xtgraph(gtot, di->get_sqrts(), gscaling);

         gtot->SetMarkerStyle(20+i);
         gtot->SetMarkerColor(mycolor(i));
         if (isth) gtot->SetMarkerSize(0);
         gtot->SetLineColor(mycolor(i));
         if (!gstat && i==0) {
            // draw axes
            TH1F *axes = haxes(gtot,*di,lTextSize,503,true);
            xmin = axes->GetXaxis()->GetXmin();
            xmax = axes->GetXaxis()->GetXmax();
            if (plotxt) axes->GetYaxis()->SetTitle("d#sigma/dp_{T} (x_{T}) [nb/GeV]");
            axes->Draw();
         }
         gtot->SetFillColor(mycolor(i));
         gtot->SetFillStyle(0);
         gtot->Draw(isth ? "3" : "P");
         if (isth) gstat->Draw("P");

         printgraph(flog,gtot);

         if (!gstat) tleg->AddEntry(gstat,di->get_legend().c_str(),isth ? "L" : "LP");
      }

      if (i>0 && idata != 0) {
         // make ratios
         pad2->cd();
         TGraphAsymmErrors *g0stat = isth ? theory[0].get_graphstat() : data[0].get_graphstat();
         TGraphAsymmErrors *g0tot = isth ? theory[0].get_graphtot() : data[0].get_graphtot();
         TGraphAsymmErrors *gratiostat=NULL, *gratiotot=NULL;
         if (gstat && g0stat) {
            if (doxLW && !isth) g0stat = xlw(g0stat,lwmode);
            if (plotxt) g0stat = xtgraph(g0stat, data[0].get_sqrts(), gscaling);
            gratiostat = ratiograph(gstat,g0stat,ginterpolation);
            gratiostat->SetMarkerStyle(20+i);
            gratiostat->SetMarkerColor(mycolor(i));
            if (isth) gratiostat->SetMarkerSize(0);
            gratiostat->SetLineColor(mycolor(i));
            if (i==1) {
               // draw axes
               TH1F *axes = haxes(gratiostat,*di,lTextSize2,505,false,xmin,xmax);
               axes->GetYaxis()->SetTitle("ratio");
               axes->Draw();
            }
            gratiostat->Draw(isth ? "L" : "P");

            printgraph(flog,gratiostat);
         }
         if (gtot && g0tot) {
            if (doxLW && !isth) g0tot = xlw(g0tot,lwmode);
            if (plotxt) g0tot = xtgraph(g0tot, data[0].get_sqrts(), gscaling);
            gratiotot = ratiograph(gtot,g0tot,ginterpolation);
            gratiotot->SetMarkerStyle(20+i);
            gratiotot->SetMarkerColor(mycolor(i));
            if (isth) gratiotot->SetMarkerSize(0);
            gratiotot->SetLineColor(mycolor(i));
            if (i==1 && !gratiostat) {
               // draw axes
               TH1F *axes = haxes(gratiotot,*di,lTextSize2,505,false,xmin,xmax);
               axes->GetYaxis()->SetTitle("ratio");
               axes->Draw();
            }
            gratiotot->SetFillColor(mycolor(i));
            gratiotot->SetFillStyle(0);
            gratiotot->Draw(isth ? "L" : "P");
            if (isth) gratiostat->Draw("P");
            printgraph(flog,gratiotot);
         }

         if (plotxt && i!=0) {
            // graph of n
            pad3->cd();
            TGraphAsymmErrors *gnstat=NULL, *gntot=NULL;
            if (gstat && g0stat) {
               if (!isth) gnstat = ngraph(gstat,g0stat,di->get_sqrts(),data[0].get_sqrts(),ginterpolation);
               else gnstat = ngraph(gstat,g0stat,di->get_sqrts(),theory[0].get_sqrts(),ginterpolation);
               gnstat->SetMarkerStyle(20+i);
               gnstat->SetMarkerColor(mycolor(i));
               if (isth) gnstat->SetMarkerSize(0);
               gnstat->SetLineColor(mycolor(i));
               if (i==1) {
                  // draw axes
                  TH1F *axes = haxes(gnstat,*di,lTextSize3,505,false,xmin,xmax);
                  axes->GetYaxis()->SetTitle("n_{eff}");
                  axes->Draw();
               }
               gnstat->Draw(isth ? "L" : "P");
               printgraph(flog,gnstat);
            }
            if (gtot && g0tot) {
               if (!isth) gntot = ngraph(gtot,g0tot,di->get_sqrts(),data[0].get_sqrts(),ginterpolation);
               else gntot = ngraph(gtot,g0tot,di->get_sqrts(),theory[0].get_sqrts(),ginterpolation);
               gntot->SetMarkerStyle(20+i);
               if (isth) gntot->SetMarkerSize(0);
               gntot->SetMarkerColor(mycolor(i));
               gntot->SetLineColor(mycolor(i));
               if (i==1 && !gnstat) {
                  // draw axes
                  TH1F *axes = haxes(gntot,*di,lTextSize3,505,false,xmin,xmax);
                  axes->GetYaxis()->SetTitle("n_{eff}");
                  axes->Draw();
               }
               gntot->Draw(isth ? "L" : "P");
               printgraph(flog,gntot);
            }
         } // if (plotxt)
      } // if (i>1) 
   } // for (unsigned int i=0; i<data.size()+theory.size(); i++) 

   flog.close();

   pad1->cd();
   tleg->Draw();

   c1->SaveAs(cname+".pdf");
   c1->SaveAs(cname+".png");
   c1->SaveAs(cname+".C");
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

TH1F* haxes(TGraphAsymmErrors *graph, dataset data, float textSize, int ndiv, bool log, double xmin_, double xmax_) {
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
   if (xmin_>0) xmin = xmin_;
   if (xmax_>0) xmax = xmax_;

   if (log) {
      // ymin=0.5*ymin;
      ymin=0.01*ymin;
      // ymax=1.5*ymax;
      ymax=100.*ymax;
   } else {
      ymin=0.8*ymin;
      ymax=1.2*ymax;
   }

   TH1F *ans = new TH1F(Form("haxes_%f_%f%f_%f%f",textSize,xmin,xmax,ymin,ymax),"haxes",1000,xmin,xmax);
   ans->SetMinimum(ymin);
   ans->SetMaximum(ymin);
   ans->GetXaxis()->SetTitle(plotxt ? "x_{T}" : data.get_xheader().c_str());
   ans->GetYaxis()->SetTitle(data.get_yheader().c_str());
   ans->GetXaxis()->SetLabelSize(textSize);
   ans->GetXaxis()->SetTitleSize(textSize);
   ans->GetYaxis()->SetLabelSize(textSize);
   ans->GetYaxis()->SetTitleSize(textSize);
   ans->GetYaxis()->SetNdivisions(ndiv);
   ans->GetXaxis()->SetNdivisions(505);
   ans->GetYaxis()->SetTitleOffset(1.7*(gTextSize/textSize));
   ans->GetYaxis()->SetRangeUser(ymin,ymax);

   return ans;
}

#endif // ifndef plot_C
