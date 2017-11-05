#include "all.C"

void plotMidRap() {
   vector<dataset> data;
   dataset data7;
   data7 = dataset_atlas_7tev(0.1);
   data.push_back(data7); 
   dataset data8;
   data8 = dataset_atlas_8tev(0.1);
   data.push_back(data8); 
   // FIXME put CMS 13TeV data here
   dataset data13;
   data13 = dataset_atlas_8tev(0.1);
   data13.set_sqrts(13);
   data.push_back(data13); 

   vector<dataset> theory;
   dataset th7;
   th7.set_sqrts(7000);
   th7.set_expname("theory");
   th7.set_legend("J/#psi 1S08 7TeV |y|<0.75");
   th7.set_graphHwU("th_inputs/ForEmilien/LHC7/direct_psi1S_1S08.HwU",0,0);
   theory.push_back(th7);
   dataset th8;
   th8.set_sqrts(8000);
   th8.set_expname("theory");
   th8.set_legend("J/#psi 1S08 8TeV |y|<0.75");
   th8.set_graphHwU("th_inputs/ForEmilien/LHC8/direct_psi1S_1S08.HwU",0,0);
   theory.push_back(th8);
   dataset th13;
   th13.set_sqrts(13000);
   th13.set_expname("theory");
   th13.set_legend("J/#psi 1S08 13TeV |y|<0.75");
   th13.set_graphHwU("th_inputs/ForEmilien/LHC13/direct_psi1S_1S08.HwU",0,0);
   theory.push_back(th13);

   plot(data,theory);
   // cout << th7.get_graphtot()->GetN() << endl;
   // th7.get_graphtot()->Draw("AP");
   // th8.get_graphtot()->Draw("P");
   // th13.get_graphtot()->Draw("P");
}
