#ifndef ins1409298_tableall
#define ins1409298_tableall

#include "table1.h"
#include "table2.h"

#include "../../include/utils.h"
#include "../../include/dataset.h"

#include <iostream>
#include <sstream>

using namespace std;

dataset dataset_atlas_7tev(double rapidity) {
   using namespace ins1409298table1;

   dataset ans;

   ans.set_sqrts(sqrts);
   ans.set_expname("ATLAS");
   ans.set_location(location);
   ans.set_reaction(reackey);
   ans.set_comment(dscomment);
   ans.set_xheader(xheader);
   ans.set_yheader(yheader);

   double absrap = fabs(rapidity);
   if (absrap<0.25) {
      ans.set_raprange(range(0,0.25));
      ans.set_graph(graph_9066_d1x1y1(true), graph_9066_d1x1y1(false));
   } else if (absrap<0.5) {
      ans.set_raprange(range(0.25,0.5));
      ans.set_graph(graph_9066_d1x1y2(true), graph_9066_d1x1y2(false));
   } else if (absrap<0.75) {
      ans.set_raprange(range(0.5,0.75));
      ans.set_graph(graph_9066_d1x1y3(true), graph_9066_d1x1y3(false));
   } else if (absrap<1.0) {
      ans.set_raprange(range(0.75, 1.0));
      ans.set_graph(graph_9066_d1x1y4(true), graph_9066_d1x1y4(false));
   } else if (absrap<1.25) {
      ans.set_raprange(range(1.0,1.25));
      ans.set_graph(graph_9066_d1x1y5(true), graph_9066_d1x1y5(false));
   } else if (absrap<1.5) {
      ans.set_raprange(range(1.25,1.5));
      ans.set_graph(graph_9066_d1x1y6(true), graph_9066_d1x1y6(false));
   } else if (absrap<1.75) {
      ans.set_raprange(range(1.5,1.75));
      ans.set_graph(graph_9066_d1x1y7(true), graph_9066_d1x1y7(false));
   } else if (absrap<2.0) {
      ans.set_raprange(range(1.75,2.0));
      ans.set_graph(graph_9066_d1x1y8(true), graph_9066_d1x1y8(false));
   } else {
      cout << "ERROR wrong rapidity value " << rapidity << endl;
   }

   ostringstream oss;
   oss << "ATLAS, 7TeV, " << ans.get_raprange().min << "<|y|<" << ans.get_raprange().max;
   ans.set_legend(oss.str());

   return ans;
}

dataset dataset_atlas_8tev(double rapidity) {
   using namespace ins1409298table2;

   dataset ans;

   ans.set_sqrts(sqrts);
   ans.set_expname("ATLAS");
   ans.set_location(location);
   ans.set_reaction(reackey);
   ans.set_comment(dscomment);
   ans.set_xheader(xheader);
   ans.set_yheader(yheader);

   double absrap = fabs(rapidity);
   if (absrap<0.25) {
      ans.set_raprange(range(0,0.25));
      ans.set_graph(graph_9066_d2x1y1(true), graph_9066_d2x1y1(false));
   } else if (absrap<0.5) {
      ans.set_raprange(range(0.25,0.5));
      ans.set_graph(graph_9066_d2x1y2(true), graph_9066_d2x1y2(false));
   } else if (absrap<0.75) {
      ans.set_raprange(range(0.5,0.75));
      ans.set_graph(graph_9066_d2x1y3(true), graph_9066_d2x1y3(false));
   } else if (absrap<1.0) {
      ans.set_raprange(range(0.75, 1.0));
      ans.set_graph(graph_9066_d2x1y4(true), graph_9066_d2x1y4(false));
   } else if (absrap<1.25) {
      ans.set_raprange(range(1.0,1.25));
      ans.set_graph(graph_9066_d2x1y5(true), graph_9066_d2x1y5(false));
   } else if (absrap<1.5) {
      ans.set_raprange(range(1.25,1.5));
      ans.set_graph(graph_9066_d2x1y6(true), graph_9066_d2x1y6(false));
   } else if (absrap<1.75) {
      ans.set_raprange(range(1.5,1.75));
      ans.set_graph(graph_9066_d2x1y7(true), graph_9066_d2x1y7(false));
   } else if (absrap<2.0) {
      ans.set_raprange(range(1.75,2.0));
      ans.set_graph(graph_9066_d2x1y8(true), graph_9066_d2x1y8(false));
   } else {
      cout << "ERROR wrong rapidity value " << rapidity << endl;
   }

   ostringstream oss;
   oss << "ATLAS, 7TeV, " << ans.get_raprange().min << "<|y|<" << ans.get_raprange().max;
   ans.set_legend(oss.str());

   return ans;
}

#endif // #ifndef ins1409298_tableall
