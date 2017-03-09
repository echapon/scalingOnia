#ifndef ins891233_tableall
#define ins891233_tableall

#include "table6.h"

#include "../../include/utils.h"
#include "../../include/dataset.h"

#include <iostream>
#include <sstream>

using namespace std;

dataset dataset_lhcb_7tev(double rapidity) {
   using namespace ins891233table6;

   dataset ans;

   ans.set_sqrts(sqrts);
   ans.set_expname("LHCb");
   ans.set_location(location);
   ans.set_reaction(reackey);
   ans.set_comment(dscomment);
   ans.set_xheader(xheader);
   ans.set_yheader(yheader);

   double absrap = fabs(rapidity);
   if (absrap<2.5) {
      ans.set_raprange(range(2,2.5));
      ans.set_graph(graph_table6_y1(true), graph_table6_y1(false));
   } else if (absrap<3) {
      ans.set_raprange(range(2.5,3));
      ans.set_graph(graph_table6_y2(true), graph_table6_y2(false));
   } else if (absrap<3.5) {
      ans.set_raprange(range(3,3.5));
      ans.set_graph(graph_table6_y3(true), graph_table6_y3(false));
   } else if (absrap<4) {
      ans.set_raprange(range(3.5,4));
      ans.set_graph(graph_table6_y4(true), graph_table6_y4(false));
   } else if (absrap<4.5) {
      ans.set_raprange(range(4,4.5));
      ans.set_graph(graph_table6_y5(true), graph_table6_y5(false));
   } else {
      cout << "ERROR wrong rapidity value " << rapidity << endl;
   }

   ostringstream oss;
   oss << "LHCb, 7TeV, " << ans.get_raprange().min << "<|y|<" << ans.get_raprange().max;
   ans.set_legend(oss.str());

   return ans;
}


#endif // #ifndef ins891233_tableall
