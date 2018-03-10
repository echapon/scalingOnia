#include "all.C"

void plotMidRap_wLDMEs() {
   vector<dataset> data;
   // dataset data7;
   // data7 = dataset_atlas_7tev(0.1);
   // data.push_back(data7); 
   dataset data72;
   data72 = dataset_cms_7tev(0.1);
   data.push_back(data72); 
   // dataset data8;
   // data8 = dataset_atlas_8tev(0.1);
   // data.push_back(data8); 
   dataset data13;
   data13 = dataset_cms_13tev(0.1);
   data.push_back(data13); 

   vector<dataset> theory;
   // available: 7, 8, 13 TeV; 1S 08, 3P J9, 3S 11, 3S 18
   dataset th7;
   th7.set_sqrts(7000);
   th7.set_expname("theory");
   th7.set_legend("J/#psi PKU 7TeV |y|<0.75");
   th7.set_graphHwU("th_inputs/ForEmilien_wLDMEs/LHC7/prompt_jpsi_PKU_cen.HwU",0,0);
   theory.push_back(th7);
   // dataset th8;
   // th8.set_sqrts(8000);
   // th8.set_expname("theory");
   // th8.set_legend("J/#psi 3PJ8 8TeV |y|<0.75");
   // th8.set_graphHwU("th_inputs/ForEmilien/LHC8/direct_psi1S_3PJ8.HwU",0,0);
   // theory.push_back(th8);
   dataset th13;
   th13.set_sqrts(13000);
   th13.set_expname("theory");
   th13.set_legend("J/#psi PKU 13TeV |y|<0.75");
   th13.set_graphHwU("th_inputs/ForEmilien_wLDMEs/LHC13/prompt_jpsi_PKU_cen.HwU",0,0);
   theory.push_back(th13);

   plot(data,theory,"plot_wLDMEs_PKU");
   // cout << th7.get_graphtot()->GetN() << endl;
   // th7.get_graphtot()->Draw("AP");
   // th8.get_graphtot()->Draw("P");
   // th13.get_graphtot()->Draw("P");
}
