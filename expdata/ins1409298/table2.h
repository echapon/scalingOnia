#ifndef ins1409298_table2
#define ins1409298_table2

#include "TGraphAsymmErrors.h"
#include <string>

namespace ins1409298table2 {
   // metadata
   string location = "Figure 7";
   string dscomment = "Summary of results for cross-section of prompt $J/psi$ decaying to a muon pair for 8 TeV data in nb/GeV. Uncertainties are statistical and systematic, respectively.";
   string reackey = "P P --> J/PSI X";
   string obskey = "D2SIG/DPT/DYRAP";
   string qual = "RE : P P --> J/PSI < MU+ MU- > X";
   double sqrts = 8000.0; // GEV
   // string qual = YRAP : 0.0 TO 0.25 : 0.25 TO 0.5 : 0.5 TO 0.75 : 0.75 TO 1.0 : 1.0 TO 1.25 : 1.25 TO 1.5 : 1.5 TO 1.75 : 1.75 TO 2.0
   string yheader = "D2SIG/DPT/DYRAP*BR [NB/GEV]";
   string xheader = "PT IN GEV ";


   // Plot: p9066_d2x1y8
   TGraphAsymmErrors* graph_9066_d2x1y1(bool stat = false) {
      double p9066_d2x1y1_xval[] = { 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.25, 
         12.75, 13.5, 14.5, 15.5, 16.5, 17.5, 19.0, 21.0, 23.0, 25.0, 
         28.0, 32.5, 37.5, 50.0, 85.0 };
      double p9066_d2x1y1_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
         0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 
         2.0, 2.5, 2.5, 10.0, 25.0 };
      double p9066_d2x1y1_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
         0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 
         2.0, 2.5, 2.5, 10.0, 25.0 };
      double p9066_d2x1y1_yval[] = { 3.249, 2.594, 1.997, 1.491, 1.135, 0.8753, 0.6843, 0.5386, 0.4291, 
         0.3446, 0.2536, 0.1728, 0.1213, 0.08658, 0.06268, 0.04095, 0.0238, 0.0144, 0.009095, 
         0.004966, 0.002201, 9.747E-4, 2.221E-4, 1.147E-5 };
      double p9066_d2x1y1_yerrminus[] = { 0.6170518616777686, 0.43402880088768303, 0.3000149996250187, 0.21102132593650338, 0.14501379244747722, 0.09901459488378468, 0.07041391907854583, 0.051513978685401494, 0.038512984823303426, 
         0.029613679271579884, 0.02040612653102004, 0.012806248474865698, 0.008309632964216891, 0.005707889277132135, 0.004017785459678005, 0.0024539967400141346, 0.0013146102083887831, 7.343704787094861E-4, 4.4803794482164115E-4, 
         2.3858122306669483E-4, 1.0707940978544849E-4, 5.334285331700958E-5, 1.6846661390317074E-5, 2.2875751353780706E-6 };
      double p9066_d2x1y1_yerrplus[] = { 0.6170518616777686, 0.43402880088768303, 0.3000149996250187, 0.21102132593650338, 0.14501379244747722, 0.09901459488378468, 0.07041391907854583, 0.051513978685401494, 0.038512984823303426, 
         0.029613679271579884, 0.02040612653102004, 0.012806248474865698, 0.008309632964216891, 0.005707889277132135, 0.004017785459678005, 0.0024539967400141346, 0.0013146102083887831, 7.343704787094861E-4, 4.4803794482164115E-4, 
         2.3858122306669483E-4, 1.0707940978544849E-4, 5.334285331700958E-5, 1.6846661390317074E-5, 2.2875751353780706E-6 };
      double p9066_d2x1y1_ystatminus[] = { 0.008, 0.005, 0.003, 0.003, 0.002, 0.0017, 0.0014, 0.0012, 0.001, 
         9.0E-4, 5.0E-4, 4.0E-4, 4.0E-4, 3.0E-4, 2.5E-4, 1.4E-4, 1.1E-4, 8.0E-5, 6.7E-5, 
         3.5E-5, 2.1E-5, 1.39E-5, 3.4E-6, 5.1E-7 };
      double p9066_d2x1y1_ystatplus[] = { 0.008, 0.005, 0.003, 0.003, 0.002, 0.0017, 0.0014, 0.0012, 0.001, 
         9.0E-4, 5.0E-4, 4.0E-4, 4.0E-4, 3.0E-4, 2.5E-4, 1.4E-4, 1.1E-4, 8.0E-5, 6.7E-5, 
         3.5E-5, 2.1E-5, 1.39E-5, 3.4E-6, 5.1E-7 };
      int p9066_d2x1y1_numpoints = 24;
      TGraphAsymmErrors* p9066_d2x1y1 = new TGraphAsymmErrors(p9066_d2x1y1_numpoints, p9066_d2x1y1_xval, p9066_d2x1y1_yval, p9066_d2x1y1_xerrminus, p9066_d2x1y1_xerrplus, p9066_d2x1y1_yerrminus, p9066_d2x1y1_yerrplus);
      if (stat) p9066_d2x1y1 = new TGraphAsymmErrors(p9066_d2x1y1_numpoints, p9066_d2x1y1_xval, p9066_d2x1y1_yval, p9066_d2x1y1_xerrminus, p9066_d2x1y1_xerrplus, p9066_d2x1y1_ystatminus, p9066_d2x1y1_ystatplus);
      p9066_d2x1y1->SetName("/HepData/9066/d2x1y1");
      p9066_d2x1y1->SetTitle("/HepData/9066/d2x1y1");
      // p9066_d2x1y1.Draw("AP");
      return p9066_d2x1y1;
   };

   TGraphAsymmErrors* graph_9066_d2x1y2(bool stat = false) {
      double p9066_d2x1y2_xval[] = { 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.25, 
         12.75, 13.5, 14.5, 15.5, 16.5, 17.5, 19.0, 21.0, 23.0, 25.0, 
         28.0, 32.5, 37.5, 50.0, 85.0 };
      double p9066_d2x1y2_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
         0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 
         2.0, 2.5, 2.5, 10.0, 25.0 };
      double p9066_d2x1y2_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
         0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 
         2.0, 2.5, 2.5, 10.0, 25.0 };
      double p9066_d2x1y2_yval[] = { 3.194, 2.524, 1.958, 1.497, 1.162, 0.9132, 0.7224, 0.5784, 0.4647, 
         0.3787, 0.2809, 0.1936, 0.1353, 0.09739, 0.07058, 0.04546, 0.02642, 0.01593, 0.01013, 
         0.005533, 0.002445, 0.001075, 2.603E-4, 1.29E-5 };
      double p9066_d2x1y2_yerrminus[] = { 0.612029411057998, 0.4290291365396993, 0.299015049788468, 0.21400934559032694, 0.1530130713370593, 0.11121623082985684, 0.08301018009858792, 0.06300960244280232, 0.04740854353383998, 
         0.03640879014743555, 0.02490501957437496, 0.01520526224699857, 0.009304837451562492, 0.006344927107540323, 0.004564832965180653, 0.002892490276561012, 0.0016324827717314506, 9.525754563287888E-4, 5.821511831131154E-4, 
         3.103578579639961E-4, 1.3705838172107534E-4, 6.413267497929585E-5, 2.111421322237701E-5, 2.6305322655310656E-6 };
      double p9066_d2x1y2_yerrplus[] = { 0.612029411057998, 0.4290291365396993, 0.299015049788468, 0.21400934559032694, 0.1530130713370593, 0.11121623082985684, 0.08301018009858792, 0.06300960244280232, 0.04740854353383998, 
         0.03640879014743555, 0.02490501957437496, 0.01520526224699857, 0.009304837451562492, 0.006344927107540323, 0.004564832965180653, 0.002892490276561012, 0.0016324827717314506, 9.525754563287888E-4, 5.821511831131154E-4, 
         3.103578579639961E-4, 1.3705838172107534E-4, 6.413267497929585E-5, 2.111421322237701E-5, 2.6305322655310656E-6 };
      double p9066_d2x1y2_ystatminus[] = { 0.006, 0.005, 0.003, 0.002, 0.002, 0.0019, 0.0013, 0.0011, 9.0E-4, 
         8.0E-4, 5.0E-4, 4.0E-4, 3.0E-4, 2.5E-4, 2.1E-4, 1.2E-4, 9.0E-5, 7.0E-5, 5.0E-5, 
         2.9E-5, 1.7E-5, 1.2E-5, 3.0E-6, 4.6E-7 };
      double p9066_d2x1y2_ystatplus[] = { 0.006, 0.005, 0.003, 0.002, 0.002, 0.0019, 0.0013, 0.0011, 9.0E-4, 
         8.0E-4, 5.0E-4, 4.0E-4, 3.0E-4, 2.5E-4, 2.1E-4, 1.2E-4, 9.0E-5, 7.0E-5, 5.0E-5, 
         2.9E-5, 1.7E-5, 1.2E-5, 3.0E-6, 4.6E-7 };
      int p9066_d2x1y2_numpoints = 24;
      TGraphAsymmErrors* p9066_d2x1y2 = new TGraphAsymmErrors(p9066_d2x1y2_numpoints, p9066_d2x1y2_xval, p9066_d2x1y2_yval, p9066_d2x1y2_xerrminus, p9066_d2x1y2_xerrplus, p9066_d2x1y2_yerrminus, p9066_d2x1y2_yerrplus);
      if (stat) p9066_d2x1y2 = new TGraphAsymmErrors(p9066_d2x1y2_numpoints, p9066_d2x1y2_xval, p9066_d2x1y2_yval, p9066_d2x1y2_xerrminus, p9066_d2x1y2_xerrplus, p9066_d2x1y2_ystatminus, p9066_d2x1y2_ystatplus);
      p9066_d2x1y2->SetName("/HepData/9066/d2x1y2");
      p9066_d2x1y2->SetTitle("/HepData/9066/d2x1y2");
      // p9066_d2x1y2.Draw("AP");
      return p9066_d2x1y2;
   };

   TGraphAsymmErrors* graph_9066_d2x1y3(bool stat = false) {
      double p9066_d2x1y3_xval[] = { 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.25, 
         12.75, 13.5, 14.5, 15.5, 16.5, 17.5, 19.0, 21.0, 23.0, 25.0, 
         28.0, 32.5, 37.5, 50.0, 85.0 };
      double p9066_d2x1y3_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
         0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 
         2.0, 2.5, 2.5, 10.0, 25.0 };
      double p9066_d2x1y3_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
         0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 
         2.0, 2.5, 2.5, 10.0, 25.0 };
      double p9066_d2x1y3_yval[] = { 3.379, 2.704, 2.029, 1.553, 1.189, 0.9159, 0.7126, 0.5609, 0.449, 
         0.3636, 0.2701, 0.1842, 0.1296, 0.09281, 0.06711, 0.04369, 0.02501, 0.01524, 0.009775, 
         0.005314, 0.002403, 0.001067, 2.541E-4, 1.225E-5 };
      double p9066_d2x1y3_yerrminus[] = { 0.46805234749972147, 0.32505537989702615, 0.23601906702637396, 0.16901183390520322, 0.12601587201618691, 0.09691160921169352, 0.07631107390149873, 0.06131174438882, 0.048208401757370054, 
         0.037808464660707926, 0.026704681237565822, 0.016604818577750254, 0.010604244433244642, 0.0072043389703705645, 0.005094330181682377, 0.0031822633454822685, 0.0017223530416264838, 9.82496819333274E-4, 5.973558403497869E-4, 
         3.0528675044947496E-4, 1.3110682667199293E-4, 6.001666435249463E-5, 1.9022355269524328E-5, 2.280898945591409E-6 };
      double p9066_d2x1y3_yerrplus[] = { 0.46805234749972147, 0.32505537989702615, 0.23601906702637396, 0.16901183390520322, 0.12601587201618691, 0.09691160921169352, 0.07631107390149873, 0.06131174438882, 0.048208401757370054, 
         0.037808464660707926, 0.026704681237565822, 0.016604818577750254, 0.010604244433244642, 0.0072043389703705645, 0.005094330181682377, 0.0031822633454822685, 0.0017223530416264838, 9.82496819333274E-4, 5.973558403497869E-4, 
         3.0528675044947496E-4, 1.3110682667199293E-4, 6.001666435249463E-5, 1.9022355269524328E-5, 2.280898945591409E-6 };
      double p9066_d2x1y3_ystatminus[] = { 0.007, 0.006, 0.003, 0.002, 0.002, 0.0015, 0.0013, 0.0012, 9.0E-4, 
         8.0E-4, 5.0E-4, 4.0E-4, 3.0E-4, 2.5E-4, 2.1E-4, 1.2E-4, 9.0E-5, 7.0E-5, 5.3E-5, 
         2.8E-5, 1.7E-5, 1.1E-5, 2.9E-6, 4.3E-7 };
      double p9066_d2x1y3_ystatplus[] = { 0.007, 0.006, 0.003, 0.002, 0.002, 0.0015, 0.0013, 0.0012, 9.0E-4, 
         8.0E-4, 5.0E-4, 4.0E-4, 3.0E-4, 2.5E-4, 2.1E-4, 1.2E-4, 9.0E-5, 7.0E-5, 5.3E-5, 
         2.8E-5, 1.7E-5, 1.1E-5, 2.9E-6, 4.3E-7 };
      int p9066_d2x1y3_numpoints = 24;
      TGraphAsymmErrors* p9066_d2x1y3 = new TGraphAsymmErrors(p9066_d2x1y3_numpoints, p9066_d2x1y3_xval, p9066_d2x1y3_yval, p9066_d2x1y3_xerrminus, p9066_d2x1y3_xerrplus, p9066_d2x1y3_yerrminus, p9066_d2x1y3_yerrplus);
      if (stat) p9066_d2x1y3 = new TGraphAsymmErrors(p9066_d2x1y3_numpoints, p9066_d2x1y3_xval, p9066_d2x1y3_yval, p9066_d2x1y3_xerrminus, p9066_d2x1y3_xerrplus, p9066_d2x1y3_ystatminus, p9066_d2x1y3_ystatplus);
      p9066_d2x1y3->SetName("/HepData/9066/d2x1y3");
      p9066_d2x1y3->SetTitle("/HepData/9066/d2x1y3");
      // p9066_d2x1y3.Draw("AP");
      return p9066_d2x1y3;
   };

   TGraphAsymmErrors* graph_9066_d2x1y4(bool stat = false) {
      double p9066_d2x1y4_xval[] = { 10.25, 10.75, 11.25, 11.75, 12.25, 
         12.75, 13.5, 14.5, 15.5, 16.5, 17.5, 19.0, 21.0, 23.0, 25.0, 
         28.0, 32.5, 37.5, 50.0, 85.0 };
      double p9066_d2x1y4_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 
         0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 
         2.0, 2.5, 2.5, 10.0, 25.0 };
      double p9066_d2x1y4_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 
         0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 
         2.0, 2.5, 2.5, 10.0, 25.0 };
      double p9066_d2x1y4_yval[] = { 1.072, 0.8291, 0.6527, 0.5139, 0.4098, 
         0.3329, 0.2448, 0.166, 0.1154, 0.08154, 0.05863, 0.0373, 0.02168, 0.0131, 0.008336, 
         0.004604, 0.002041, 9.231E-4, 2.313E-4, 1.33E-5 };
      double p9066_d2x1y4_yerrminus[] = { 0.06603029607687672, 0.0535302718095098, 0.04482510457321879, 0.037019454344979205, 0.028917295862511072, 
         0.022418073066166947, 0.015408114745159449, 0.010259142264341596, 0.006906518659932803, 0.005056688639811631, 0.003935606687665829, 0.002692675249635574, 0.0016524527224704494, 0.0010323759005323595, 6.640391554720249E-4, 
         3.699864862397004E-4, 1.6876018487783188E-4, 8.124807690031808E-5, 2.5465663156493685E-5, 3.0464733709651887E-6 };
      double p9066_d2x1y4_yerrplus[] = { 0.06603029607687672, 0.0535302718095098, 0.04482510457321879, 0.037019454344979205, 0.028917295862511072, 
         0.022418073066166947, 0.015408114745159449, 0.010259142264341596, 0.006906518659932803, 0.005056688639811631, 0.003935606687665829, 0.002692675249635574, 0.0016524527224704494, 0.0010323759005323595, 6.640391554720249E-4, 
         3.699864862397004E-4, 1.6876018487783188E-4, 8.124807690031808E-5, 2.5465663156493685E-5, 3.0464733709651887E-6 };
      double p9066_d2x1y4_ystatminus[] = { 0.002, 0.0018, 0.0015, 0.0012, 0.001, 
         9.0E-4, 5.0E-4, 0.0018, 3.0E-4, 2.6E-4, 2.1E-4, 1.2E-4, 9.0E-5, 7.0E-5, 5.2E-5, 
         2.7E-5, 1.6E-5, 1.1E-5, 2.9E-6, 4.7E-7 };
      double p9066_d2x1y4_ystatplus[] = { 0.002, 0.0018, 0.0015, 0.0012, 0.001, 
         9.0E-4, 5.0E-4, 0.0018, 3.0E-4, 2.6E-4, 2.1E-4, 1.2E-4, 9.0E-5, 7.0E-5, 5.2E-5, 
         2.7E-5, 1.6E-5, 1.1E-5, 2.9E-6, 4.7E-7 };
      int p9066_d2x1y4_numpoints = 24;
      TGraphAsymmErrors* p9066_d2x1y4 = new TGraphAsymmErrors(p9066_d2x1y4_numpoints, p9066_d2x1y4_xval, p9066_d2x1y4_yval, p9066_d2x1y4_xerrminus, p9066_d2x1y4_xerrplus, p9066_d2x1y4_yerrminus, p9066_d2x1y4_yerrplus);
      if (stat) p9066_d2x1y4 = new TGraphAsymmErrors(p9066_d2x1y4_numpoints, p9066_d2x1y4_xval, p9066_d2x1y4_yval, p9066_d2x1y4_xerrminus, p9066_d2x1y4_xerrplus, p9066_d2x1y4_ystatminus, p9066_d2x1y4_ystatplus);
      p9066_d2x1y4->SetName("/HepData/9066/d2x1y4");
      p9066_d2x1y4->SetTitle("/HepData/9066/d2x1y4");
      // p9066_d2x1y4.Draw("AP");
      return p9066_d2x1y4;
   };

   TGraphAsymmErrors* graph_9066_d2x1y5(bool stat = false) {
      double p9066_d2x1y5_xval[] = { 10.25, 10.75, 11.25, 11.75, 12.25, 
         12.75, 13.5, 14.5, 15.5, 16.5, 17.5, 19.0, 21.0, 23.0, 25.0, 
         28.0, 32.5, 37.5, 50.0, 85.0 };
      double p9066_d2x1y5_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 
         0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 
         2.0, 2.5, 2.5, 10.0, 25.0 };
      double p9066_d2x1y5_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 
         0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 
         2.0, 2.5, 2.5, 10.0, 25.0 };
      double p9066_d2x1y5_yval[] = { 1.095, 0.8632, 0.6884, 0.5475, 0.44, 
         0.3533, 0.2603, 0.1743, 0.1193, 0.08404, 0.05997, 0.03809, 0.02192, 0.0131, 0.008188, 
         0.004414, 0.001941, 8.71E-4, 2.158E-4, 1.177E-5 };
      double p9066_d2x1y5_yerrminus[] = { 0.08002499609497023, 0.05913383126434478, 0.0446323873437216, 0.03472434304634142, 0.027231048455760937, 
         0.021718655575334308, 0.015811388300841896, 0.010807404868885036, 0.007705841939723394, 0.005797257972524597, 0.004447032718566393, 0.0029924070578716392, 0.0018521878954360976, 0.0011621101496846157, 7.507343071952953E-4, 
         4.118215633013891E-4, 1.8959430371189953E-4, 8.702367493963926E-5, 2.2848413511664216E-5, 2.0729930052945187E-6 };
      double p9066_d2x1y5_yerrplus[] = { 0.08002499609497023, 0.05913383126434478, 0.0446323873437216, 0.03472434304634142, 0.027231048455760937, 
         0.021718655575334308, 0.015811388300841896, 0.010807404868885036, 0.007705841939723394, 0.005797257972524597, 0.004447032718566393, 0.0029924070578716392, 0.0018521878954360976, 0.0011621101496846157, 7.507343071952953E-4, 
         4.118215633013891E-4, 1.8959430371189953E-4, 8.702367493963926E-5, 2.2848413511664216E-5, 2.0729930052945187E-6 };
      double p9066_d2x1y5_ystatminus[] = { 0.002, 0.002, 0.0017, 0.0013, 0.0013, 
         9.0E-4, 6.0E-4, 4.0E-4, 3.0E-4, 2.9E-4, 2.5E-4, 1.2E-4, 9.0E-5, 7.0E-5, 5.1E-5, 
         2.6E-5, 1.5E-5, 1.04E-5, 2.6E-6, 4.2E-7 };
      double p9066_d2x1y5_ystatplus[] = { 0.002, 0.002, 0.0017, 0.0013, 0.0013, 
         9.0E-4, 6.0E-4, 4.0E-4, 3.0E-4, 2.9E-4, 2.5E-4, 1.2E-4, 9.0E-5, 7.0E-5, 5.1E-5, 
         2.6E-5, 1.5E-5, 1.04E-5, 2.6E-6, 4.2E-7 };
      int p9066_d2x1y5_numpoints = 24;
      TGraphAsymmErrors* p9066_d2x1y5 = new TGraphAsymmErrors(p9066_d2x1y5_numpoints, p9066_d2x1y5_xval, p9066_d2x1y5_yval, p9066_d2x1y5_xerrminus, p9066_d2x1y5_xerrplus, p9066_d2x1y5_yerrminus, p9066_d2x1y5_yerrplus);
      if (stat) p9066_d2x1y5 = new TGraphAsymmErrors(p9066_d2x1y5_numpoints, p9066_d2x1y5_xval, p9066_d2x1y5_yval, p9066_d2x1y5_xerrminus, p9066_d2x1y5_xerrplus, p9066_d2x1y5_ystatminus, p9066_d2x1y5_ystatplus);
      p9066_d2x1y5->SetName("/HepData/9066/d2x1y5");
      p9066_d2x1y5->SetTitle("/HepData/9066/d2x1y5");
      // p9066_d2x1y5.Draw("AP");
      return p9066_d2x1y5;
   };

   TGraphAsymmErrors* graph_9066_d2x1y6(bool stat = false) {
      double p9066_d2x1y6_xval[] = { 10.25, 10.75, 11.25, 11.75, 12.25, 
         12.75, 13.5, 14.5, 15.5, 16.5, 17.5, 19.0, 21.0, 23.0, 25.0, 
         28.0, 32.5, 37.5, 50.0, 85.0 };
      double p9066_d2x1y6_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 
         0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 
         2.0, 2.5, 2.5, 10.0, 25.0 };
      double p9066_d2x1y6_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 
         0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 
         2.0, 2.5, 2.5, 10.0, 25.0 };
      double p9066_d2x1y6_yval[] = { 1.036, 0.8231, 0.6446, 0.5153, 0.411, 
         0.3315, 0.2474, 0.1667, 0.1151, 0.08055, 0.05782, 0.03689, 0.02104, 0.01251, 0.0078, 
         0.004182, 0.001817, 7.826E-4, 1.843E-4, 9.98E-6 };
      double p9066_d2x1y6_yerrminus[] = { 0.06503076195155644, 0.047537984812147854, 0.035132036661713766, 0.028429737951659, 0.023920911353876133, 
         0.02051974658712919, 0.016407620180879372, 0.012506398362438326, 0.009604686356149273, 0.00714370352128362, 0.005223830012548264, 0.0033817894671312704, 0.001951640335717624, 0.0011515641536623134, 7.213431915530915E-4, 
         3.906200199682551E-4, 1.7548219282878817E-4, 8.14343907695023E-5, 2.3394443784796424E-5, 2.381537318624254E-6 };
      double p9066_d2x1y6_yerrplus[] = { 0.06503076195155644, 0.047537984812147854, 0.035132036661713766, 0.028429737951659, 0.023920911353876133, 
         0.02051974658712919, 0.016407620180879372, 0.012506398362438326, 0.009604686356149273, 0.00714370352128362, 0.005223830012548264, 0.0033817894671312704, 0.001951640335717624, 0.0011515641536623134, 7.213431915530915E-4, 
         3.906200199682551E-4, 1.7548219282878817E-4, 8.14343907695023E-5, 2.3394443784796424E-5, 2.381537318624254E-6 };
      double p9066_d2x1y6_ystatminus[] = { 0.002, 0.0019, 0.0015, 0.0013, 0.001, 
         9.0E-4, 5.0E-4, 4.0E-4, 3.0E-4, 2.3E-4, 2.0E-4, 1.1E-4, 8.0E-5, 6.0E-5, 4.4E-5, 
         2.2E-5, 1.3E-5, 8.4E-6, 2.1E-6, 3.34E-7 };
      double p9066_d2x1y6_ystatplus[] = { 0.002, 0.0019, 0.0015, 0.0013, 0.001, 
         9.0E-4, 5.0E-4, 4.0E-4, 3.0E-4, 2.3E-4, 2.0E-4, 1.1E-4, 8.0E-5, 6.0E-5, 4.4E-5, 
         2.2E-5, 1.3E-5, 8.4E-6, 2.1E-6, 3.34E-7 };
      int p9066_d2x1y6_numpoints = 24;
      TGraphAsymmErrors* p9066_d2x1y6 = new TGraphAsymmErrors(p9066_d2x1y6_numpoints, p9066_d2x1y6_xval, p9066_d2x1y6_yval, p9066_d2x1y6_xerrminus, p9066_d2x1y6_xerrplus, p9066_d2x1y6_yerrminus, p9066_d2x1y6_yerrplus);
      if (stat) p9066_d2x1y6 = new TGraphAsymmErrors(p9066_d2x1y6_numpoints, p9066_d2x1y6_xval, p9066_d2x1y6_yval, p9066_d2x1y6_xerrminus, p9066_d2x1y6_xerrplus, p9066_d2x1y6_ystatminus, p9066_d2x1y6_ystatplus);
      p9066_d2x1y6->SetName("/HepData/9066/d2x1y6");
      p9066_d2x1y6->SetTitle("/HepData/9066/d2x1y6");
      // p9066_d2x1y6.Draw("AP");
      return p9066_d2x1y6;
   };

   TGraphAsymmErrors* graph_9066_d2x1y7(bool stat = false) {
      double p9066_d2x1y7_xval[] = { 10.25, 10.75, 11.25, 11.75, 12.25, 
         12.75, 13.5, 14.5, 15.5, 16.5, 17.5, 19.0, 21.0, 23.0, 25.0, 
         28.0, 32.5, 37.5, 50.0, 85.0 };
      double p9066_d2x1y7_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 
         0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 
         2.0, 2.5, 2.5, 10.0, 25.0 };
      double p9066_d2x1y7_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 
         0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 
         2.0, 2.5, 2.5, 10.0, 25.0 };
      double p9066_d2x1y7_yval[] = { 1.036, 0.8132, 0.6416, 0.5072, 0.4073, 
         0.325, 0.2401, 0.1609, 0.1103, 0.07714, 0.0554, 0.03535, 0.01992, 0.01185, 0.007341, 
         0.00405, 0.00172, 7.485E-4, 1.741E-4, 8.565E-6 };
      double p9066_d2x1y7_yerrminus[] = { 0.06702984409947557, 0.054423524325423835, 0.04521869082580786, 0.03751612986436634, 0.031212978070027216, 
         0.025412595302329904, 0.019104188022525323, 0.013003461077728498, 0.008905054744357275, 0.006213549710109352, 0.00448361461323339, 0.0028817355881482257, 0.0016115210206509875, 9.718538984847465E-4, 6.044609168507092E-4, 
         3.3272811723688157E-4, 1.4558159224297557E-4, 6.802065862662607E-5, 2.0044201156444224E-5, 1.951279067688679E-6 };
      double p9066_d2x1y7_yerrplus[] = { 0.06702984409947557, 0.054423524325423835, 0.04521869082580786, 0.03751612986436634, 0.031212978070027216, 
         0.025412595302329904, 0.019104188022525323, 0.013003461077728498, 0.008905054744357275, 0.006213549710109352, 0.00448361461323339, 0.0028817355881482257, 0.0016115210206509875, 9.718538984847465E-4, 6.044609168507092E-4, 
         3.3272811723688157E-4, 1.4558159224297557E-4, 6.802065862662607E-5, 2.0044201156444224E-5, 1.951279067688679E-6 };
      double p9066_d2x1y7_ystatminus[] = { 0.002, 0.0016, 0.0013, 0.0011, 9.0E-4, 
         8.0E-4, 4.0E-4, 3.0E-4, 3.0E-4, 2.1E-4, 1.8E-4, 1.0E-4, 7.0E-5, 6.0E-5, 4.2E-5, 
         2.2E-5, 1.3E-5, 8.4E-6, 2.4E-6, 3.31E-7 };
      double p9066_d2x1y7_ystatplus[] = { 0.002, 0.0016, 0.0013, 0.0011, 9.0E-4, 
         8.0E-4, 4.0E-4, 3.0E-4, 3.0E-4, 2.1E-4, 1.8E-4, 1.0E-4, 7.0E-5, 6.0E-5, 4.2E-5, 
         2.2E-5, 1.3E-5, 8.4E-6, 2.4E-6, 3.31E-7 };
      int p9066_d2x1y7_numpoints = 24;
      TGraphAsymmErrors* p9066_d2x1y7 = new TGraphAsymmErrors(p9066_d2x1y7_numpoints, p9066_d2x1y7_xval, p9066_d2x1y7_yval, p9066_d2x1y7_xerrminus, p9066_d2x1y7_xerrplus, p9066_d2x1y7_yerrminus, p9066_d2x1y7_yerrplus);
      if (stat) p9066_d2x1y7 = new TGraphAsymmErrors(p9066_d2x1y7_numpoints, p9066_d2x1y7_xval, p9066_d2x1y7_yval, p9066_d2x1y7_xerrminus, p9066_d2x1y7_xerrplus, p9066_d2x1y7_ystatminus, p9066_d2x1y7_ystatplus);
      p9066_d2x1y7->SetName("/HepData/9066/d2x1y7");
      p9066_d2x1y7->SetTitle("/HepData/9066/d2x1y7");
      // p9066_d2x1y7.Draw("AP");
      return p9066_d2x1y7;
   };

   TGraphAsymmErrors* graph_9066_d2x1y8(bool stat = false) {
      double p9066_d2x1y8_xval[] = { 10.25, 10.75, 11.25, 11.75, 12.25, 
         12.75, 13.5, 14.5, 15.5, 16.5, 17.5, 19.0, 21.0, 23.0, 25.0, 
         28.0, 32.5, 37.5, 50.0, 85.0 };
      double p9066_d2x1y8_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 
         0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 
         2.0, 2.5, 2.5, 10.0, 25.0 };
      double p9066_d2x1y8_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 
         0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 
         2.0, 2.5, 2.5, 10.0, 25.0 };
      double p9066_d2x1y8_yval[] = { 1.048, 0.819, 0.6363, 0.5018, 0.4006, 
         0.3251, 0.236, 0.1567, 0.1057, 0.07464, 0.05331, 0.03374, 0.01886, 0.0112, 0.006931, 
         0.003691, 0.00157, 6.451E-4, 1.51E-4, 6.975E-6 };
      double p9066_d2x1y8_yerrminus[] = { 0.09502105029939419, 0.06621933252457321, 0.049417102302745354, 0.03791595970036892, 0.02951372562046344, 
         0.023413671220037235, 0.016604818577750254, 0.010604244433244642, 0.0071063352017759485, 0.00505478980769725, 0.0036544356609468443, 0.0023721087664776252, 0.001381774221788784, 8.521150157109074E-4, 5.397156658834353E-4, 
         2.888390555309306E-4, 1.2467958934805648E-4, 5.11775341336411E-5, 1.284562182223967E-5, 1.1861757879842264E-6 };
      double p9066_d2x1y8_yerrplus[] = { 0.09502105029939419, 0.06621933252457321, 0.049417102302745354, 0.03791595970036892, 0.02951372562046344, 
         0.023413671220037235, 0.016604818577750254, 0.010604244433244642, 0.0071063352017759485, 0.00505478980769725, 0.0036544356609468443, 0.0023721087664776252, 0.001381774221788784, 8.521150157109074E-4, 5.397156658834353E-4, 
         2.888390555309306E-4, 1.2467958934805648E-4, 5.11775341336411E-5, 1.284562182223967E-5, 1.1861757879842264E-6 };
      double p9066_d2x1y8_ystatminus[] = { 0.002, 0.0016, 0.0013, 0.0011, 9.0E-4, 
         8.0E-4, 4.0E-4, 3.0E-4, 3.0E-4, 2.2E-4, 1.8E-4, 1.0E-4, 7.0E-5, 6.0E-5, 4.3E-5, 
         2.2E-5, 1.3E-5, 8.3E-6, 2.5E-6, 3.38E-7 };
      double p9066_d2x1y8_ystatplus[] = { 0.002, 0.0016, 0.0013, 0.0011, 9.0E-4, 
         8.0E-4, 4.0E-4, 3.0E-4, 3.0E-4, 2.2E-4, 1.8E-4, 1.0E-4, 7.0E-5, 6.0E-5, 4.3E-5, 
         2.2E-5, 1.3E-5, 8.3E-6, 2.5E-6, 3.38E-7 };
      int p9066_d2x1y8_numpoints = 24;
      TGraphAsymmErrors* p9066_d2x1y8 = new TGraphAsymmErrors(p9066_d2x1y8_numpoints, p9066_d2x1y8_xval, p9066_d2x1y8_yval, p9066_d2x1y8_xerrminus, p9066_d2x1y8_xerrplus, p9066_d2x1y8_yerrminus, p9066_d2x1y8_yerrplus);
      if (stat) p9066_d2x1y8 = new TGraphAsymmErrors(p9066_d2x1y8_numpoints, p9066_d2x1y8_xval, p9066_d2x1y8_yval, p9066_d2x1y8_xerrminus, p9066_d2x1y8_xerrplus, p9066_d2x1y8_ystatminus, p9066_d2x1y8_ystatplus);
      p9066_d2x1y8->SetName("/HepData/9066/d2x1y8");
      p9066_d2x1y8->SetTitle("/HepData/9066/d2x1y8");
      // p9066_d2x1y8.Draw("AP");
      return p9066_d2x1y8;
   };
}

#endif // #ifndef ins1409298_table2
