#ifndef ins1345023_tableall
#define ins1345023_tableall

#include "table1.h"

#include "../../include/utils.h"
#include "../../include/dataset.h"

#include <iostream>
#include <sstream>

using namespace std;

dataset dataset_cms_7tev(double rapidity) {
   using namespace ins1345023table1;

   dataset ans;

   ans.set_sqrts(sqrts);
   ans.set_expname("CMS");
   ans.set_location(location);
   ans.set_reaction(reackey);
   ans.set_comment(dscomment);
   ans.set_xheader(xheader);
   ans.set_yheader(yheader);

   double absrap = fabs(rapidity);
   if (absrap<0.3) {
      ans.set_raprange(range(0,0.3));
      ans.set_graph(graph_table1_y1(1,true), graph_table1_y(1,false));
   } else if (absrap<0.6) {
      ans.set_raprange(range(0.3,0.6));
      ans.set_graph(graph_table1_y1(2,true), graph_table1_y(2,false));
   } else if (absrap<0.9) {
      ans.set_raprange(range(0.6,0.9));
      ans.set_graph(graph_table1_y1(3,true), graph_table1_y(3,false));
   } else if (absrap<1.2) {
      ans.set_raprange(range(0.9,1.2));
      ans.set_graph(graph_table1_y1(4,true), graph_table1_y(4,false));
   } else {
      cout << "ERROR wrong rapidity value " << rapidity << endl;
   }

   ostringstream oss;
   oss << "CMS, 7TeV, " << ans.get_raprange().min << "<|y|<" << ans.get_raprange().max;
   ans.set_legend(oss.str());

   return ans;
}

#endif // #ifndef ins1345023_tableall
