#ifndef ins1409298_table1
#define ins1409298_table1

#include "TGraphAsymmErrors.h"
#include <string>

namespace ins1409298table1 {
   // metadata
   string location = "Figure 7";
   string dscomment = "Summary of results for cross-section of prompt $J/psi$ decaying to a muon pair for 7 TeV data in nb/GeV. Uncertainties are statistical and systematic, respectively.";
   string reackey = "P P --> J/PSI X";
   string obskey = "D2SIG/DPT/DYRAP";
   string qual = "RE : P P --> J/PSI < MU+ MU- > X";
   double sqrts = 7000.0; // GEV
   // qual: YRAP : 0.0 TO 0.25 : 0.25 TO 0.5 : 0.5 TO 0.75 : 0.75 TO 1.0 : 1.0 TO 1.25 : 1.25 TO 1.5 : 1.5 TO 1.75 : 1.75 TO 2.0
   string yheader = "d^{2}#sigma/dp_{T}/dy #times BR [nb/GeV]"; // "D2SIG/DPT/DYRAP*BR [NB/GEV]";
   string xheader = "p_{T} [GeV]"; // "PT IN GEV ";


   // Plot: p9066_d1x1y8
   TGraphAsymmErrors* graph_9066_d1x1y1(bool stat = false) {
      double p9066_d1x1y1_xval[] = { 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.5, 
         13.5, 14.5, 15.5, 16.5, 17.5, 19.0, 21.0, 23.0, 25.0, 28.0, 
         35.0, 50.0, 80.0 };
      double p9066_d1x1y1_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 
         0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 
         5.0, 10.0, 20.0 };
      double p9066_d1x1y1_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 
         0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 
         5.0, 10.0, 20.0 };
      double p9066_d1x1y1_yval[] = { 3.12, 2.526, 1.903, 1.465, 1.106, 0.8471, 0.6552, 0.5117, 0.3546, 
         0.2315, 0.1513, 0.1039, 0.0723, 0.05171, 0.03326, 0.01871, 0.0111, 0.006976, 0.00369, 
         0.001198, 1.5E-4, 9.043E-6 };
      double p9066_d1x1y1_yerrminus[] = { 0.24318717071424636, 0.17075128110793197, 0.11941942890501528, 0.08935882720806042, 0.06527633568147036, 0.04898754943860736, 0.037315010384562405, 0.02824340631014609, 0.017791290003819286, 
         0.012180722474467595, 0.007679192665899196, 0.0055731499172371095, 0.0038989485762189785, 0.0028580587817607947, 0.0019006577808748212, 0.0010988175462741756, 6.868041933477111E-4, 4.819647289999549E-4, 2.7432280255202995E-4, 
         1.0793516572461452E-4, 2.2585836269662454E-5, 3.34952190618303E-6 };
      double p9066_d1x1y1_yerrplus[] = { 0.24318717071424636, 0.17075128110793197, 0.11941942890501528, 0.08935882720806042, 0.06527633568147036, 0.04898754943860736, 0.037315010384562405, 0.02824340631014609, 0.017791290003819286, 
         0.012180722474467595, 0.007679192665899196, 0.0055731499172371095, 0.0038989485762189785, 0.0028580587817607947, 0.0019006577808748212, 0.0010988175462741756, 6.868041933477111E-4, 4.819647289999549E-4, 2.7432280255202995E-4, 
         1.0793516572461452E-4, 2.2585836269662454E-5, 3.34952190618303E-6 };
      double p9066_d1x1y1_ystatminus[] = { 0.024, 0.016, 0.01, 0.008, 0.006, 0.0053, 0.004, 0.0037, 0.0018, 
         0.0014, 0.0011, 9.0E-4, 7.3E-4, 6.2E-4, 3.4E-4, 2.5E-4, 1.9E-4, 1.47E-4, 7.8E-5, 
         2.5E-5, 6.6E-6, 1.216E-6 };
      double p9066_d1x1y1_ystatplus[] = { 0.024, 0.016, 0.01, 0.008, 0.006, 0.0053, 0.004, 0.0037, 0.0018, 
         0.0014, 0.0011, 9.0E-4, 7.3E-4, 6.2E-4, 3.4E-4, 2.5E-4, 1.9E-4, 1.47E-4, 7.8E-5, 
         2.5E-5, 6.6E-6, 1.216E-6 };
      int p9066_d1x1y1_numpoints = 22;
      TGraphAsymmErrors* p9066_d1x1y1 = new TGraphAsymmErrors(p9066_d1x1y1_numpoints, p9066_d1x1y1_xval, p9066_d1x1y1_yval, p9066_d1x1y1_xerrminus, p9066_d1x1y1_xerrplus, p9066_d1x1y1_yerrminus, p9066_d1x1y1_yerrplus);
      if (stat) p9066_d1x1y1 = new TGraphAsymmErrors(p9066_d1x1y1_numpoints, p9066_d1x1y1_xval, p9066_d1x1y1_yval, p9066_d1x1y1_xerrminus, p9066_d1x1y1_xerrplus, p9066_d1x1y1_ystatminus, p9066_d1x1y1_ystatplus);
      p9066_d1x1y1->SetName("/HepData/9066/d1x1y1");
      p9066_d1x1y1->SetTitle("/HepData/9066/d1x1y1");
      // p9066_d1x1y1.Draw("AP");
      return p9066_d1x1y1;
   };

   TGraphAsymmErrors* graph_9066_d1x1y2(bool stat = false) {
      double p9066_d1x1y2_xval[] = { 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.5, 
         13.5, 14.5, 15.5, 16.5, 17.5, 19.0, 21.0, 23.0, 25.0, 28.0, 
         35.0, 50.0, 80.0 };
      double p9066_d1x1y2_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 
         0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 
         5.0, 10.0, 20.0 };
      double p9066_d1x1y2_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 
         0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 
         5.0, 10.0, 20.0 };
      double p9066_d1x1y2_yval[] = { 3.178, 2.544, 2.015, 1.566, 1.213, 0.9467, 0.7513, 0.5943, 0.4293, 
         0.28, 0.1871, 0.1292, 0.09165, 0.06731, 0.04278, 0.02392, 0.01467, 0.008978, 0.00471, 
         0.001463, 2.009E-4, 8.783E-6 };
      double p9066_d1x1y2_yerrminus[] = { 0.20850419660045214, 0.14975313018431366, 0.11543396380615195, 0.08438009243891595, 0.06537583651472462, 0.049822183814040105, 0.03710552519504339, 0.02948525055006316, 0.02077811348510735, 
         0.01367186892857008, 0.008967719888578144, 0.006462971452822611, 0.004547361872558638, 0.003483819742753634, 0.0023052114870440845, 0.001341640786499874, 8.612200647918045E-4, 5.697174738412014E-4, 3.2136116753584277E-4, 
         1.2041594578792295E-4, 2.9231831964486935E-5, 3.252642771655074E-6 };
      double p9066_d1x1y2_yerrplus[] = { 0.20850419660045214, 0.14975313018431366, 0.11543396380615195, 0.08438009243891595, 0.06537583651472462, 0.049822183814040105, 0.03710552519504339, 0.02948525055006316, 0.02077811348510735, 
         0.01367186892857008, 0.008967719888578144, 0.006462971452822611, 0.004547361872558638, 0.003483819742753634, 0.0023052114870440845, 0.001341640786499874, 8.612200647918045E-4, 5.697174738412014E-4, 3.2136116753584277E-4, 
         1.2041594578792295E-4, 2.9231831964486935E-5, 3.252642771655074E-6 };
      double p9066_d1x1y2_ystatminus[] = { 0.025, 0.015, 0.01, 0.008, 0.007, 0.0047, 0.0039, 0.0033, 0.0018, 
         0.0014, 0.0011, 9.0E-4, 7.2E-4, 6.1E-4, 3.4E-4, 2.4E-4, 1.9E-4, 1.37E-4, 7.7E-5, 
         2.4E-5, 6.5E-6, 1.066E-6 };
      double p9066_d1x1y2_ystatplus[] = { 0.025, 0.015, 0.01, 0.008, 0.007, 0.0047, 0.0039, 0.0033, 0.0018, 
         0.0014, 0.0011, 9.0E-4, 7.2E-4, 6.1E-4, 3.4E-4, 2.4E-4, 1.9E-4, 1.37E-4, 7.7E-5, 
         2.4E-5, 6.5E-6, 1.066E-6 };
      int p9066_d1x1y2_numpoints = 22;
      TGraphAsymmErrors* p9066_d1x1y2 = new TGraphAsymmErrors(p9066_d1x1y2_numpoints, p9066_d1x1y2_xval, p9066_d1x1y2_yval, p9066_d1x1y2_xerrminus, p9066_d1x1y2_xerrplus, p9066_d1x1y2_yerrminus, p9066_d1x1y2_yerrplus);
      if (stat) p9066_d1x1y2 = new TGraphAsymmErrors(p9066_d1x1y2_numpoints, p9066_d1x1y2_xval, p9066_d1x1y2_yval, p9066_d1x1y2_xerrminus, p9066_d1x1y2_xerrplus, p9066_d1x1y2_ystatminus, p9066_d1x1y2_ystatplus);
      p9066_d1x1y2->SetName("/HepData/9066/d1x1y2");
      p9066_d1x1y2->SetTitle("/HepData/9066/d1x1y2");
      // p9066_d1x1y2.Draw("AP");
      return p9066_d1x1y2;
   };

   TGraphAsymmErrors* graph_9066_d1x1y3(bool stat = false) {
      double p9066_d1x1y3_xval[] = { 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.5, 
         13.5, 14.5, 15.5, 16.5, 17.5, 19.0, 21.0, 23.0, 25.0, 28.0, 
         35.0, 50.0, 80.0 };
      double p9066_d1x1y3_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 
         0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 
         5.0, 10.0, 20.0 };
      double p9066_d1x1y3_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 
         0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 
         5.0, 10.0, 20.0 };
      double p9066_d1x1y3_yval[] = { 3.491, 2.838, 2.179, 1.681, 1.284, 0.9892, 0.7687, 0.6059, 0.432, 
         0.2756, 0.1835, 0.1262, 0.0889, 0.06556, 0.04143, 0.02371, 0.0143, 0.008743, 0.004748, 
         0.001419, 1.792E-4, 8.063E-6 };
      double p9066_d1x1y3_yerrminus[] = { 0.21557365330670628, 0.16001562423713506, 0.11950732195141851, 0.08537564055396599, 0.06637017402418048, 0.04796133859683235, 0.03847596652457219, 0.02830494656416083, 0.02019925741209315, 
         0.01288759093081403, 0.008570880934886449, 0.006067124524847005, 0.004219016473065731, 0.003202577087284551, 0.0020780038498520643, 0.0012649505919204908, 7.907591289387686E-4, 5.338773267333985E-4, 3.322845166419886E-4, 
         1.1618089343777658E-4, 2.5733441277839233E-5, 2.8785352177800434E-6 };
      double p9066_d1x1y3_yerrplus[] = { 0.21557365330670628, 0.16001562423713506, 0.11950732195141851, 0.08537564055396599, 0.06637017402418048, 0.04796133859683235, 0.03847596652457219, 0.02830494656416083, 0.02019925741209315, 
         0.01288759093081403, 0.008570880934886449, 0.006067124524847005, 0.004219016473065731, 0.003202577087284551, 0.0020780038498520643, 0.0012649505919204908, 7.907591289387686E-4, 5.338773267333985E-4, 3.322845166419886E-4, 
         1.1618089343777658E-4, 2.5733441277839233E-5, 2.8785352177800434E-6 };
      double p9066_d1x1y3_ystatminus[] = { 0.026, 0.018, 0.011, 0.008, 0.007, 0.005, 0.0046, 0.0034, 0.002, 
         0.0015, 0.0011, 9.0E-4, 7.6E-4, 6.3E-4, 3.4E-4, 2.5E-4, 1.8E-4, 1.37E-4, 7.8E-5, 
         2.7E-5, 6.1E-6, 1.022E-6 };
      double p9066_d1x1y3_ystatplus[] = { 0.026, 0.018, 0.011, 0.008, 0.007, 0.005, 0.0046, 0.0034, 0.002, 
         0.0015, 0.0011, 9.0E-4, 7.6E-4, 6.3E-4, 3.4E-4, 2.5E-4, 1.8E-4, 1.37E-4, 7.8E-5, 
         2.7E-5, 6.1E-6, 1.022E-6 };
      int p9066_d1x1y3_numpoints = 22;
      TGraphAsymmErrors* p9066_d1x1y3 = new TGraphAsymmErrors(p9066_d1x1y3_numpoints, p9066_d1x1y3_xval, p9066_d1x1y3_yval, p9066_d1x1y3_xerrminus, p9066_d1x1y3_xerrplus, p9066_d1x1y3_yerrminus, p9066_d1x1y3_yerrplus);
      if (stat) p9066_d1x1y3 = new TGraphAsymmErrors(p9066_d1x1y3_numpoints, p9066_d1x1y3_xval, p9066_d1x1y3_yval, p9066_d1x1y3_xerrminus, p9066_d1x1y3_xerrplus, p9066_d1x1y3_ystatminus, p9066_d1x1y3_ystatplus);
      p9066_d1x1y3->SetName("/HepData/9066/d1x1y3");
      p9066_d1x1y3->SetTitle("/HepData/9066/d1x1y3");
      // p9066_d1x1y3.Draw("AP");
      return p9066_d1x1y3;
   };

   TGraphAsymmErrors* graph_9066_d1x1y4(bool stat = false) {
      double p9066_d1x1y4_xval[] = { 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.5, 
         13.5, 14.5, 15.5, 16.5, 17.5, 19.0, 21.0, 23.0, 25.0, 28.0, 
         35.0, 50.0, 80.0 };
      double p9066_d1x1y4_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 
         0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 
         5.0, 10.0, 20.0 };
      double p9066_d1x1y4_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 
         0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 
         5.0, 10.0, 20.0 };
      double p9066_d1x1y4_yval[] = { 2.917, 2.305, 1.801, 1.408, 1.09, 0.842, 0.6615, 0.5097, 0.3628, 
         0.2366, 0.1558, 0.1065, 0.07393, 0.0544, 0.03516, 0.02043, 0.01219, 0.007452, 0.004237, 
         0.00129, 1.753E-4, 1.023E-5 };
      double p9066_d1x1y4_yerrminus[] = { 0.21191035840656774, 0.14057738082636198, 0.09687620966986683, 0.07467261881037787, 0.0565685424949238, 0.0431335136523794, 0.033391765451979326, 0.02476307735318856, 0.01691862878604528, 
         0.011088733020503288, 0.0074148499647666505, 0.004883646178829912, 0.00356202189774291, 0.002549078264785136, 0.001735655495770978, 0.0010623088063270493, 6.992138442565337E-4, 4.6675903847702833E-4, 2.864349838968697E-4, 
         1.038123306741545E-4, 2.5551907952244977E-5, 3.874596753211875E-6 };
      double p9066_d1x1y4_yerrplus[] = { 0.21191035840656774, 0.14057738082636198, 0.09687620966986683, 0.07467261881037787, 0.0565685424949238, 0.0431335136523794, 0.033391765451979326, 0.02476307735318856, 0.01691862878604528, 
         0.011088733020503288, 0.0074148499647666505, 0.004883646178829912, 0.00356202189774291, 0.002549078264785136, 0.001735655495770978, 0.0010623088063270493, 6.992138442565337E-4, 4.6675903847702833E-4, 2.864349838968697E-4, 
         1.038123306741545E-4, 2.5551907952244977E-5, 3.874596753211875E-6 };
      double p9066_d1x1y4_ystatminus[] = { 0.035, 0.021, 0.013, 0.01, 0.008, 0.0061, 0.0051, 0.0036, 0.002, 
         0.0014, 0.0013, 9.0E-4, 7.6E-4, 6.3E-4, 3.5E-4, 2.6E-4, 2.0E-4, 1.5E-4, 6.9E-5, 
         2.4E-5, 7.9E-6, 1.15E-6 };
      double p9066_d1x1y4_ystatplus[] = { 0.035, 0.021, 0.013, 0.01, 0.008, 0.0061, 0.0051, 0.0036, 0.002, 
         0.0014, 0.0013, 9.0E-4, 7.6E-4, 6.3E-4, 3.5E-4, 2.6E-4, 2.0E-4, 1.5E-4, 6.9E-5, 
         2.4E-5, 7.9E-6, 1.15E-6 };
      int p9066_d1x1y4_numpoints = 22;
      TGraphAsymmErrors* p9066_d1x1y4 = new TGraphAsymmErrors(p9066_d1x1y4_numpoints, p9066_d1x1y4_xval, p9066_d1x1y4_yval, p9066_d1x1y4_xerrminus, p9066_d1x1y4_xerrplus, p9066_d1x1y4_yerrminus, p9066_d1x1y4_yerrplus);
      if (stat) p9066_d1x1y4 = new TGraphAsymmErrors(p9066_d1x1y4_numpoints, p9066_d1x1y4_xval, p9066_d1x1y4_yval, p9066_d1x1y4_xerrminus, p9066_d1x1y4_xerrplus, p9066_d1x1y4_ystatminus, p9066_d1x1y4_ystatplus);
      p9066_d1x1y4->SetName("/HepData/9066/d1x1y4");
      p9066_d1x1y4->SetTitle("/HepData/9066/d1x1y4");
      // p9066_d1x1y4.Draw("AP");
      return p9066_d1x1y4;
   };

   TGraphAsymmErrors* graph_9066_d1x1y5(bool stat = false) {
      double p9066_d1x1y5_xval[] = { 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.5, 
         13.5, 14.5, 15.5, 16.5, 17.5, 19.0, 21.0, 23.0, 25.0, 28.0, 
         35.0, 50.0, 80.0 };
      double p9066_d1x1y5_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 
         0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 
         5.0, 10.0, 20.0 };
      double p9066_d1x1y5_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 
         0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 
         5.0, 10.0, 20.0 };
      double p9066_d1x1y5_yval[] = { 3.051, 2.518, 1.963, 1.539, 1.19, 0.9163, 0.7295, 0.5689, 0.4017, 
         0.2563, 0.1685, 0.1108, 0.07616, 0.05382, 0.03426, 0.01938, 0.01164, 0.007111, 0.003883, 
         0.001122, 1.583E-4, 9.517E-6 };
      double p9066_d1x1y5_yerrminus[] = { 0.35875479090877654, 0.25068306683938585, 0.1856906028855526, 0.13861096637712328, 0.09940824915468535, 0.07334446127690897, 0.05939377071713834, 0.044658705758228154, 0.028229417280560363, 
         0.017990275150758534, 0.012180722474467595, 0.00795110055275369, 0.005576710499927354, 0.003952480233979671, 0.002596497641054195, 0.0015339491516996252, 9.610411021387171E-4, 6.211682541791716E-4, 3.5993888370110835E-4, 
         1.1618089343777658E-4, 2.4961570463414358E-5, 2.425217722184959E-6 };
      double p9066_d1x1y5_yerrplus[] = { 0.35875479090877654, 0.25068306683938585, 0.1856906028855526, 0.13861096637712328, 0.09940824915468535, 0.07334446127690897, 0.05939377071713834, 0.044658705758228154, 0.028229417280560363, 
         0.017990275150758534, 0.012180722474467595, 0.00795110055275369, 0.005576710499927354, 0.003952480233979671, 0.002596497641054195, 0.0015339491516996252, 9.610411021387171E-4, 6.211682541791716E-4, 3.5993888370110835E-4, 
         1.1618089343777658E-4, 2.4961570463414358E-5, 2.425217722184959E-6 };
      double p9066_d1x1y5_ystatminus[] = { 0.064, 0.029, 0.016, 0.013, 0.009, 0.0071, 0.0059, 0.0048, 0.0027, 
         0.0018, 0.0014, 9.0E-4, 8.6E-4, 7.0E-4, 3.7E-4, 2.7E-4, 2.0E-4, 1.57E-4, 8.4E-5, 
         2.7E-5, 1.52E-5, 1.185E-6 };
      double p9066_d1x1y5_ystatplus[] = { 0.064, 0.029, 0.016, 0.013, 0.009, 0.0071, 0.0059, 0.0048, 0.0027, 
         0.0018, 0.0014, 9.0E-4, 8.6E-4, 7.0E-4, 3.7E-4, 2.7E-4, 2.0E-4, 1.57E-4, 8.4E-5, 
         2.7E-5, 1.52E-5, 1.185E-6 };
      int p9066_d1x1y5_numpoints = 22;
      TGraphAsymmErrors* p9066_d1x1y5 = new TGraphAsymmErrors(p9066_d1x1y5_numpoints, p9066_d1x1y5_xval, p9066_d1x1y5_yval, p9066_d1x1y5_xerrminus, p9066_d1x1y5_xerrplus, p9066_d1x1y5_yerrminus, p9066_d1x1y5_yerrplus);
      if (stat) p9066_d1x1y5 = new TGraphAsymmErrors(p9066_d1x1y5_numpoints, p9066_d1x1y5_xval, p9066_d1x1y5_yval, p9066_d1x1y5_xerrminus, p9066_d1x1y5_xerrplus, p9066_d1x1y5_ystatminus, p9066_d1x1y5_ystatplus);
      p9066_d1x1y5->SetName("/HepData/9066/d1x1y5");
      p9066_d1x1y5->SetTitle("/HepData/9066/d1x1y5");
      // p9066_d1x1y5.Draw("AP");
      return p9066_d1x1y5;
   };

   TGraphAsymmErrors* graph_9066_d1x1y6(bool stat = false) {
      double p9066_d1x1y6_xval[] = { 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.5, 
         13.5, 14.5, 15.5, 16.5, 17.5, 19.0, 21.0, 23.0, 25.0, 28.0, 
         35.0, 50.0, 80.0 };
      double p9066_d1x1y6_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 
         0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 
         5.0, 10.0, 20.0 };
      double p9066_d1x1y6_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 
         0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 
         5.0, 10.0, 20.0 };
      double p9066_d1x1y6_yval[] = { 3.035, 2.38, 1.829, 1.39, 1.065, 0.8144, 0.6345, 0.5058, 0.3607, 
         0.2302, 0.1538, 0.1022, 0.06964, 0.05065, 0.03185, 0.01759, 0.01037, 0.006521, 0.003606, 
         0.001036, 1.413E-4, 6.244E-6 };
      double p9066_d1x1y6_yerrminus[] = { 0.2609310253687744, 0.17129214809792073, 0.12317873193047572, 0.08868483523128404, 0.06562011886609168, 0.04733465960583217, 0.03630385654444993, 0.028016602220826133, 0.018531055015837603, 
         0.011511733144926529, 0.007808969202141855, 0.005080354318352215, 0.0035027132340515687, 0.0026998148084637213, 0.0016735889578985636, 9.604686356149274E-4, 6.08769250208977E-4, 4.0925786492137205E-4, 2.490401574043833E-4, 
         8.612200647918045E-5, 2.0593445559206454E-5, 2.414164244619657E-6 };
      double p9066_d1x1y6_yerrplus[] = { 0.2609310253687744, 0.17129214809792073, 0.12317873193047572, 0.08868483523128404, 0.06562011886609168, 0.04733465960583217, 0.03630385654444993, 0.028016602220826133, 0.018531055015837603, 
         0.011511733144926529, 0.007808969202141855, 0.005080354318352215, 0.0035027132340515687, 0.0026998148084637213, 0.0016735889578985636, 9.604686356149274E-4, 6.08769250208977E-4, 4.0925786492137205E-4, 2.490401574043833E-4, 
         8.612200647918045E-5, 2.0593445559206454E-5, 2.414164244619657E-6 };
      double p9066_d1x1y6_ystatminus[] = { 0.039, 0.021, 0.017, 0.011, 0.009, 0.0064, 0.0054, 0.0042, 0.0022, 
         0.0016, 0.0013, 9.0E-4, 7.1E-4, 6.1E-4, 2.8E-4, 2.4E-4, 1.5E-4, 1.36E-4, 7.0E-5, 
         1.9E-5, 6.0E-6, 8.3E-7 };
      double p9066_d1x1y6_ystatplus[] = { 0.039, 0.021, 0.017, 0.011, 0.009, 0.0064, 0.0054, 0.0042, 0.0022, 
         0.0016, 0.0013, 9.0E-4, 7.1E-4, 6.1E-4, 2.8E-4, 2.4E-4, 1.5E-4, 1.36E-4, 7.0E-5, 
         1.9E-5, 6.0E-6, 8.3E-7 };
      int p9066_d1x1y6_numpoints = 22;
      TGraphAsymmErrors* p9066_d1x1y6 = new TGraphAsymmErrors(p9066_d1x1y6_numpoints, p9066_d1x1y6_xval, p9066_d1x1y6_yval, p9066_d1x1y6_xerrminus, p9066_d1x1y6_xerrplus, p9066_d1x1y6_yerrminus, p9066_d1x1y6_yerrplus);
      if (stat) p9066_d1x1y6 = new TGraphAsymmErrors(p9066_d1x1y6_numpoints, p9066_d1x1y6_xval, p9066_d1x1y6_yval, p9066_d1x1y6_xerrminus, p9066_d1x1y6_xerrplus, p9066_d1x1y6_ystatminus, p9066_d1x1y6_ystatplus);
      p9066_d1x1y6->SetName("/HepData/9066/d1x1y6");
      p9066_d1x1y6->SetTitle("/HepData/9066/d1x1y6");
      // p9066_d1x1y6.Draw("AP");
      return p9066_d1x1y6;
   };

   TGraphAsymmErrors* graph_9066_d1x1y7(bool stat = false) {
      double p9066_d1x1y7_xval[] = { 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.5, 
         13.5, 14.5, 15.5, 16.5, 17.5, 19.0, 21.0, 23.0, 25.0, 28.0, 
         35.0, 50.0, 80.0 };
      double p9066_d1x1y7_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 
         0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 
         5.0, 10.0, 20.0 };
      double p9066_d1x1y7_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 
         0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 
         5.0, 10.0, 20.0 };
      double p9066_d1x1y7_yval[] = { 3.403, 2.567, 1.954, 1.491, 1.132, 0.8648, 0.6742, 0.5289, 0.366, 
         0.2316, 0.1492, 0.1014, 0.06867, 0.04895, 0.03047, 0.01707, 0.01, 0.0062, 0.003277, 
         9.401E-4, 1.161E-4, 4.894E-6 };
      double p9066_d1x1y7_yerrminus[] = { 0.26780029873022926, 0.18388311504866345, 0.13055267136294071, 0.09533624704172071, 0.06835934464285041, 0.04997369307945932, 0.03809107507015259, 0.029347061181658377, 0.019204426573058618, 
         0.011684177335182823, 0.007454528824815154, 0.005080354318352215, 0.0034595086356302104, 0.002495956730394179, 0.001560160248179654, 9.374433316206372E-4, 5.920304046246273E-4, 4.007056276120913E-4, 2.332402195162747E-4, 
         8.001262400396577E-5, 1.7477127910500624E-5, 1.7108389754737295E-6 };
      double p9066_d1x1y7_yerrplus[] = { 0.26780029873022926, 0.18388311504866345, 0.13055267136294071, 0.09533624704172071, 0.06835934464285041, 0.04997369307945932, 0.03809107507015259, 0.029347061181658377, 0.019204426573058618, 
         0.011684177335182823, 0.007454528824815154, 0.005080354318352215, 0.0034595086356302104, 0.002495956730394179, 0.001560160248179654, 9.374433316206372E-4, 5.920304046246273E-4, 4.007056276120913E-4, 2.332402195162747E-4, 
         8.001262400396577E-5, 1.7477127910500624E-5, 1.7108389754737295E-6 };
      double p9066_d1x1y7_ystatminus[] = { 0.031, 0.018, 0.012, 0.008, 0.007, 0.0061, 0.0047, 0.0038, 0.002, 
         0.0014, 9.0E-4, 9.0E-4, 6.9E-4, 5.7E-4, 2.5E-4, 1.8E-4, 1.6E-4, 1.21E-4, 6.5E-5, 
         2.21E-5, 6.8E-6, 7.39E-7 };
      double p9066_d1x1y7_ystatplus[] = { 0.031, 0.018, 0.012, 0.008, 0.007, 0.0061, 0.0047, 0.0038, 0.002, 
         0.0014, 9.0E-4, 9.0E-4, 6.9E-4, 5.7E-4, 2.5E-4, 1.8E-4, 1.6E-4, 1.21E-4, 6.5E-5, 
         2.21E-5, 6.8E-6, 7.39E-7 };
      int p9066_d1x1y7_numpoints = 22;
      TGraphAsymmErrors* p9066_d1x1y7 = new TGraphAsymmErrors(p9066_d1x1y7_numpoints, p9066_d1x1y7_xval, p9066_d1x1y7_yval, p9066_d1x1y7_xerrminus, p9066_d1x1y7_xerrplus, p9066_d1x1y7_yerrminus, p9066_d1x1y7_yerrplus);
      if (stat) p9066_d1x1y7 = new TGraphAsymmErrors(p9066_d1x1y7_numpoints, p9066_d1x1y7_xval, p9066_d1x1y7_yval, p9066_d1x1y7_xerrminus, p9066_d1x1y7_xerrplus, p9066_d1x1y7_ystatminus, p9066_d1x1y7_ystatplus);
      p9066_d1x1y7->SetName("/HepData/9066/d1x1y7");
      p9066_d1x1y7->SetTitle("/HepData/9066/d1x1y7");
      // p9066_d1x1y7.Draw("AP");
      return p9066_d1x1y7;
   };

   TGraphAsymmErrors* graph_9066_d1x1y8(bool stat = false) {
      double p9066_d1x1y8_xval[] = { 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.5, 
         13.5, 14.5, 15.5, 16.5, 17.5, 19.0, 21.0, 23.0, 25.0, 28.0, 
         35.0, 50.0, 80.0 };
      double p9066_d1x1y8_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 
         0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 
         5.0, 10.0, 20.0 };
      double p9066_d1x1y8_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 
         0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 
         5.0, 10.0, 20.0 };
      double p9066_d1x1y8_yval[] = { 3.565, 2.702, 2.016, 1.535, 1.175, 0.8908, 0.6991, 0.5493, 0.3794, 
         0.2387, 0.1543, 0.1038, 0.07227, 0.05123, 0.0311, 0.01764, 0.01017, 0.006351, 0.003232, 
         9.264E-4, 1.17E-4, 6.537E-6 };
      double p9066_d1x1y8_yerrminus[] = { 0.2982633064927699, 0.2008083663595718, 0.14051334456200237, 0.10259142264341595, 0.07355270219373317, 0.05453897688809353, 0.04124196891517183, 0.031468237955119126, 0.020479257799051215, 
         0.012403628501370072, 0.0077646635471216655, 0.005277309920783505, 0.0036808966298987533, 0.0026181863951980195, 0.0016650825805346714, 0.001008959860450355, 6.063827174318213E-4, 4.3232510914819645E-4, 3.088462400612965E-4, 
         1.0745287339108248E-4, 1.9313207915827966E-5, 2.946512005745098E-6 };
      double p9066_d1x1y8_yerrplus[] = { 0.2982633064927699, 0.2008083663595718, 0.14051334456200237, 0.10259142264341595, 0.07355270219373317, 0.05453897688809353, 0.04124196891517183, 0.031468237955119126, 0.020479257799051215, 
         0.012403628501370072, 0.0077646635471216655, 0.005277309920783505, 0.0036808966298987533, 0.0026181863951980195, 0.0016650825805346714, 0.001008959860450355, 6.063827174318213E-4, 4.3232510914819645E-4, 3.088462400612965E-4, 
         1.0745287339108248E-4, 1.9313207915827966E-5, 2.946512005745098E-6 };
      double p9066_d1x1y8_ystatminus[] = { 0.044, 0.018, 0.012, 0.011, 0.009, 0.0051, 0.0053, 0.0041, 0.0018, 
         0.0016, 0.001, 9.0E-4, 6.1E-4, 5.0E-4, 3.4E-4, 2.4E-4, 1.4E-4, 1.31E-4, 2.05E-4, 
         7.24E-5, 1.02E-5, 2.147E-6 };
      double p9066_d1x1y8_ystatplus[] = { 0.044, 0.018, 0.012, 0.011, 0.009, 0.0051, 0.0053, 0.0041, 0.0018, 
         0.0016, 0.001, 9.0E-4, 6.1E-4, 5.0E-4, 3.4E-4, 2.4E-4, 1.4E-4, 1.31E-4, 2.05E-4, 
         7.24E-5, 1.02E-5, 2.147E-6 };
      int p9066_d1x1y8_numpoints = 22;
      TGraphAsymmErrors* p9066_d1x1y8 = new TGraphAsymmErrors(p9066_d1x1y8_numpoints, p9066_d1x1y8_xval, p9066_d1x1y8_yval, p9066_d1x1y8_xerrminus, p9066_d1x1y8_xerrplus, p9066_d1x1y8_yerrminus, p9066_d1x1y8_yerrplus);
      if (stat) p9066_d1x1y8 = new TGraphAsymmErrors(p9066_d1x1y8_numpoints, p9066_d1x1y8_xval, p9066_d1x1y8_yval, p9066_d1x1y8_xerrminus, p9066_d1x1y8_xerrplus, p9066_d1x1y8_ystatminus, p9066_d1x1y8_ystatplus);
      p9066_d1x1y8->SetName("/HepData/9066/d1x1y8");
      p9066_d1x1y8->SetTitle("/HepData/9066/d1x1y8");
      // p9066_d1x1y8.Draw("AP");
      return p9066_d1x1y8;
   };
}

#endif // #ifndef ins1409298_table1
