#ifndef ins1205646_tableall
#define ins1205646_tableall

#include "table1.h"

#include "../../include/utils.h"
#include "../../include/dataset.h"

#include <iostream>
#include <sstream>

using namespace std;

dataset dataset_lhcb_276tev(double rapidity) {
   using namespace ins1205646table1;

   dataset ans;

   ans.set_sqrts(sqrts);
   ans.set_expname("LHCb");
   ans.set_location(location);
   ans.set_reaction(reackey);
   ans.set_comment(dscomment);
   ans.set_xheader(xheader);
   ans.set_yheader(yheader);

   double absrap = fabs(rapidity);
   if (absrap<4.5) {
      ans.set_raprange(range(2,4.5));
      ans.set_graph(graph_table1_y1(true), graph_table1_y1(false));
   } else {
      cout << "ERROR wrong rapidity value " << rapidity << endl;
   }

   ostringstream oss;
   oss << "LHCb, 2.76TeV, " << ans.get_raprange().min << "<|y|<" << ans.get_raprange().max;
   ans.set_legend(oss.str());

   return ans;
}


#endif // #ifndef ins1205646_tableall
