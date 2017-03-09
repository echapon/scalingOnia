#ifndef ins1230344_tableall
#define ins1230344_tableall

#include "table3.h"
#include "table4.h"
#include "table5.h"
#include "table6.h"
#include "table7.h"

#include "../../include/utils.h"
#include "../../include/dataset.h"

#include <iostream>
#include <sstream>

using namespace std;

dataset dataset_lhcb_8tev(double rapidity) {

   dataset ans;

   ans.set_sqrts(ins1230344table3::sqrts);
   ans.set_expname("LHCb");
   ans.set_location(ins1230344table3::location);
   ans.set_reaction(ins1230344table3::reackey);
   ans.set_comment(ins1230344table3::dscomment);
   ans.set_xheader(ins1230344table3::xheader);
   ans.set_yheader(ins1230344table3::yheader);

   double absrap = fabs(rapidity);
   if (absrap<2.5) {
      using namespace ins1230344table3;
      ans.set_raprange(range(2,2.5));
      ans.set_graph(graph_table3_y1(true), graph_table3_y1(false));
   } else if (absrap<3) {
      using namespace ins1230344table4;
      ans.set_raprange(range(2.5,3));
      ans.set_graph(graph_table4_y1(true), graph_table4_y1(false));
   } else if (absrap<3.5) {
      using namespace ins1230344table5;
      ans.set_raprange(range(3,3.5));
      ans.set_graph(graph_table5_y1(true), graph_table5_y1(false));
   } else if (absrap<4) {
      using namespace ins1230344table6;
      ans.set_raprange(range(3.5,4));
      ans.set_graph(graph_table6_y1(true), graph_table6_y1(false));
   } else if (absrap<4.5) {
      using namespace ins1230344table7;
      ans.set_raprange(range(4,4.5));
      ans.set_graph(graph_table7_y1(true), graph_table7_y1(false));
   } else {
      cout << "ERROR wrong rapidity value " << rapidity << endl;
   }

   ostringstream oss;
   oss << "LHCb, 8TeV, " << ans.get_raprange().min << "<|y|<" << ans.get_raprange().max;
   ans.set_legend(oss.str());

   return ans;
}


#endif // #ifndef ins1230344_tableall
