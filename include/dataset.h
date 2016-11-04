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
      string         fcomment;
      string         fxheader;
      string         fyheader;
      string         freaction;

   public:
      // constructor
      dataset() {};

      // destructor
      ~dataset() {
         // if (fstat) {
         //    delete fstat;
         //    fstat=0;
         // }
         // if (fsyst) {
         //    delete fsyst;
         //    fsyst=0;
         // }
         // if (ftot) {
         //    delete ftot;
         //    ftot=0;
         // }
      };

      // getters
      double             get_sqrts() const {return fsqrts;};
      range              get_raprange() const {return fraprange;};
      TGraphAsymmErrors* get_graphstat() const {return fstat;};
      TGraphAsymmErrors* get_graphsyst() const {return fsyst;};
      TGraphAsymmErrors* get_graphtot() const {return ftot;};
      string             get_name() const {return fname;};
      string             get_expname() const {return fexpname;};
      string             get_legend() const {return flegend;};
      string             get_location() const {return flocation;};
      string             get_comment() const {return fcomment;};
      string             get_xheader() const {return fxheader;};
      string             get_yheader() const {return fyheader;};
      string             get_reaction() const {return freaction;};
      
      // setters
      void set_sqrts(double sqrts) {fsqrts = sqrts;};
      void set_raprange(range raprange) {fraprange.min = raprange.min; fraprange.max = raprange.max;};
      void set_graph(TGraphAsymmErrors *gstat, TGraphAsymmErrors *gsyst);
      void set_name(string name) {fname = name;};
      void set_expname(string expname) {fexpname = expname;};
      void set_legend(string legend) {flegend = legend;};
      void set_location(string location) {flocation = location;};
      void set_comment(string comment) {fcomment = comment;};
      void set_xheader(string xheader) {fxheader = xheader;};
      void set_yheader(string yheader) {fyheader = yheader;};
      void set_reaction(string reaction) {freaction = reaction;};
};

#endif // ifndef dataset_h
