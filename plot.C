#ifndef plot_C
#define plot_C

#include "include/dataset.h"
#include "include/utils.h"

#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLegend.h"
#include "Rtypes.h"

#include <vector>
#include <iostream>

using namespace std;

// global settings
const Escaling       gscaling = mtpt;
const Einterpolation ginterpolation = loglin;
float                gTextSize = 0.05;
// should we use the Lafferty & Wyatt prescription to change the x position of the points? (assuming a locally exponentially falling spectrum)
bool                 doxLW = true;

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

void plot(vector<dataset> data, vector<dataset> theory, bool plotxt = false) {
   if (data.size()==0) return;
   if (theory.size()!=0 && theory.size()!=data.size()) {
      cout << "Error, inconsistent data and theory vector sizes." << endl;
      return;
   }

   float lTextSize = gTextSize;
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
      pad1->SetBottomMargin(0.04);
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
   } else {
      pad1 = (TPad*) c1;
   }

   pad1->cd();
   pad1->SetLogy();

   TLegend *tleg = new TLegend(0.6,0.6,0.9,0.9);
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
         if (i==0) gstat->Draw("AP");
         else gstat->Draw("P");

         tleg->AddEntry(gstat,data[i].get_legend().c_str(),"LP");
      }

      TGraphAsymmErrors *gtot  = data[i].get_graphtot();
      if (gtot) {
         if (doxLW) gtot = xlw(gtot);
         if (plotxt) gtot = xtgraph(gtot, data[i].get_sqrts(), gscaling);

         gtot->SetMarkerStyle(20+i);
         gtot->SetMarkerColor(mycolor(i));
         gtot->SetLineColor(mycolor(i));
         if (!gstat && i==0) gtot->Draw("AP");
         else gtot->Draw("P");

         if (!gstat) tleg->AddEntry(gstat,data[i].get_legend().c_str(),"LP");
      }

      if (i>0) {
         // make ratios
         pad2->cd();
         TGraphAsymmErrors *g0stat = data[0].get_graphstat();
         TGraphAsymmErrors *g0tot = data[0].get_graphtot();
         TGraphAsymmErrors *gratiostat, *gratiotot;
         if (gstat && g0stat) {
            if (doxLW) g0stat = xlw(g0stat);
            if (plotxt) g0stat = xtgraph(g0stat, data[0].get_sqrts(), gscaling);
            gratiostat = ratiograph(gstat,g0stat,ginterpolation);
            gratiostat->SetMarkerStyle(20+i);
            gratiostat->SetMarkerColor(mycolor(i));
            gratiostat->SetLineColor(mycolor(i));
            if (i==1) gratiostat->Draw("AP");
            else gratiostat->Draw("P");
         }
         if (gtot && g0tot) {
            if (doxLW) g0tot = xlw(g0tot);
            if (plotxt) g0tot = xtgraph(g0tot, data[0].get_sqrts(), gscaling);
            gratiotot = ratiograph(gtot,g0tot,ginterpolation);
            gratiotot->SetMarkerStyle(20+i);
            gratiotot->SetMarkerColor(mycolor(i));
            gratiotot->SetLineColor(mycolor(i));
            if (i==1) gratiotot->Draw("AP");
            else gratiotot->Draw("P");
         }

         if (plotxt) {
            // graph of n
            pad3->cd();
            TGraphAsymmErrors *gnstat, *gntot;
            if (gstat && g0stat) {
               gnstat = ngraph(gstat,g0stat,data[i].get_sqrts(),data[0].get_sqrts(),ginterpolation);
               gnstat->SetMarkerStyle(20+i);
               gnstat->SetMarkerColor(mycolor(i));
               gnstat->SetLineColor(mycolor(i));
               if (i==1) gnstat->Draw("AP");
               else gnstat->Draw("P");
            }
            if (gtot && g0tot) {
               gntot = ngraph(gtot,g0tot,data[i].get_sqrts(),data[0].get_sqrts(),ginterpolation);
               gntot->SetMarkerStyle(20+i);
               gntot->SetMarkerColor(mycolor(i));
               gntot->SetLineColor(mycolor(i));
               if (i==1) gntot->Draw("AP");
               else gntot->Draw("P");
            }
         } // if (plotxt)
      } // if (i>1) 

      // take care of theory here
      // [not implemented yet] //FIXME

   } // for (int i=0; i<data.size(); i++) 

   pad1->cd();
   tleg->Draw();
};

#endif // ifndef plot_C
