#ifndef dataset_h
#define dataset_h

#include <iostream>
#include <string>
#include "TGraphAsymmErrors.h"
#include "utils.h"

class dataset {
   private:
      double         fsqrts;
      range          fraprange;
      TGraphAsymmErrors *fstat;
      TGraphAsymmErrors *fsyst;
      TGraphAsymmErrors *ftot;
      string         fname;
      string         fexpname;
      string         flegend;
      string         flocation;
      string         freaction;

   public:
      // constructor
      dataset() {};

      // destructor
      ~dataset() {
         if (fstat) {
            delete fstat;
            fstat=0;
         }
         if (fsyst) {
            delete fsyst;
            fsyst=0;
         }
         if (ftot) {
            delete ftot;
            ftot=0;
         }
      };
      
      // setters
      void set_sqrts(double sqrts) {fsqrts = sqrts;};
      void set_raprange(range raprange) {fraprange.min = raprangemin; fraprange.max = raprange.max;};
      void set_graph(TGraphAsymmErrors *gstat, TGraphAsymmErrors *gsyst);
      void set_name(string name) {fname = name;};
      void set_expname(string expname) {fexpname = expname;};
      void set_legend(string legend) {flegend = legend;};
      void set_location(string location) {flocation = location;};
      void set_reaction(string reaction) {freaction = reaction;};
};

#endif // ifndef dataset_h
