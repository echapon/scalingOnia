#ifndef dummy_input_h
#define dummy_input_h

#include "../include/utils.h"
#include "../include/dataset.h"

#include <sstream>

using namespace std;

dataset dummy_dataset(int n, double ptmin, double ptmax, double sqrts, double neff=-4., int binned=false) {
   dataset ans;
   ans.set_sqrts(sqrts);
   ans.set_expname("dummy");
   ans.set_comment("dummy");
   ans.set_xheader("p_T [GeV]");
   ans.set_yheader("dummy");

   TGraphAsymmErrors *g = new TGraphAsymmErrors(n);
   double dpt = (ptmax-ptmin)/n;
   for (int i=0; i<n; i++) {
      double pt = ptmin + i*dpt;
      double y = pow(pt,neff);
      if (binned) { // emulate data: integrate the xsec in a bin
         pt += dpt/2.;
         y = (pow(pt+dpt/2.,neff+1)-pow(pt-dpt/2.,neff+1))/((neff+1)*dpt);
      }
      g->SetPoint(i,pt,y);
      g->SetPointError(i,dpt/2.,dpt/2.,0,0);
   }

   ans.set_graph(g,g);

   ostringstream oss;
   oss << "dummy, neff=" << neff <<", " << sqrts/1000. << "TeV";
   ans.set_legend(oss.str());

   return ans;
}

#endif // #ifndef dummy_input_h
