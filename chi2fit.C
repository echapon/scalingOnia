#include "all.C"
#include "include/chi2ds.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>

using namespace std;

// fit cross sections or neff?
// true = fit cross sections
// false = fit neff
bool fitxsec = true;

chi2ds thechi2;

double myfunc(const double *xx) {
   return thechi2.chi2(xx);
}

void chi2fit(const char * minName = "Minuit2",
                          const char *algoName = "" ,
                          int randomSeed = -1) {

   vector<dataset> data;
   dataset data7;
   data7 = dataset_cms_7tev(0.1);
   data.push_back(data7); 
   dataset data13;
   data13 = dataset_cms_13tev(0.1);
   data.push_back(data13); 

   thechi2.setdata(data);

   const int ntheories=4;
   const char* tags[ntheories] = {"1S08", "3PJ8", "3S11", "3S18"};

   for (int i=0; i<ntheories; i++) {
      vector<dataset> theory;
      dataset th7;
      th7.set_sqrts(7000);
      th7.set_expname("theory");
      th7.set_legend(Form("J/#psi %s 7TeV |y|<0.75",tags[i]));
      th7.set_graphHwU(Form("th_inputs/ForEmilien/LHC7/direct_psi1S_%s.HwU",tags[i]),0,0);
      theory.push_back(th7);
      dataset th13;
      th13.set_sqrts(13000);
      th13.set_expname("theory");
      th13.set_legend(Form("J/#psi %s 13TeV |y|<0.75",tags[i]));
      th13.set_graphHwU(Form("th_inputs/ForEmilien/LHC13/direct_psi1S_%s.HwU",tags[i]),0,0);
      theory.push_back(th13);

      thechi2.addtheory(theory);
   }

   // now do the minimisation
   ROOT::Math::Minimizer* minimum =
      ROOT::Math::Factory::CreateMinimizer(minName, algoName);

   // set tolerance , etc...
   minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
   minimum->SetMaxIterations(10000);  // for GSL
   minimum->SetTolerance(0.001);
   minimum->SetPrintLevel(1);

   // create function wrapper for minimizer
   // a IMultiGenFunction type
   int npars;
   ROOT::Math::Functor f(&myfunc,fitxsec ? ntheories : ntheories-1);
   double step[ntheories] = {0.01,0.01,0.01,0.01};
   // starting point

   double variable[ntheories] = {1./ntheories,1./ntheories,1./ntheories,1./ntheories};

   minimum->SetFunction(f);

   // Set the free variables to be minimized !
   for (int i=0; i<ntheories; i++) {
      if (!fitxsec && i==ntheories-1) break;
      minimum->SetVariable(i,tags[i],variable[i], step[i]);
   }

   // do the minimization
   minimum->Minimize();

   const double *xs = minimum->X();
   const double *err = minimum->Errors();
   std::cout << "Minimum: chi2/ndf: " << minimum->MinValue() << "/" << thechi2.ndf()  << std::endl;
   cout << "Best fit parameters: " << endl;
   for (int i=0; i<ntheories; i++) {
      if (!fitxsec && i==ntheories-1) break;
      cout << tags[i] << ": " << xs[i] << "+/-" << err[i] << endl;
   }

   // expected minimum is 0
   if ( minimum->MinValue()  < 1.E-4  && f(xs) < 1.E-4)
      std::cout << "Minimizer " << minName << " - " << algoName
                << "   converged to the right minimum" << std::endl;
   else {
      std::cout << "Minimizer " << minName << " - " << algoName
                << "   failed to converge !!!" << std::endl;
      Error("NumericalMinimization","fail to converge");
   }

   return 0;
}
