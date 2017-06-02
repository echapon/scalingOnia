#include "all.C"

void plotLHCb() {
   vector<dataset> data;
   dataset data1;
   data1 = dataset_lhcb_8tev(2.7);
   data.push_back(data1); 
   dataset data2;
   data2 = dataset_lhcb_13tev(2.7);
   data.push_back(data2); 

   vector<dataset> theory;
   dataset th1;
   th1.set_sqrts(8000);
   th1.set_expname("theory");
   // th1.set_legend("COM 8TeV");
   // th1.set_graph("th_inputs/pp2psi8X_8TeV.root");
   th1.set_legend("CSM LO 8TeV");
   th1.set_graphHwU("th_inputs/ggpsi1g_8000/results.HwU",8);
   theory.push_back(th1);
   dataset th2;
   th2.set_sqrts(13000);
   th2.set_expname("theory");
   // th2.set_legend("COM 13TeV");
   // th2.set_graph("th_inputs/pp2psi8X_13TeV.root");
   th2.set_legend("CSM LO 13TeV");
   th2.set_graphHwU("th_inputs/ggpsi1g_8000/results.HwU",8);
   theory.push_back(th2);

   plot(data,theory);
}
