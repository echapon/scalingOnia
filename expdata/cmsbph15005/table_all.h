#ifndef cmsbph15005_tableall
#define cmsbph15005_tableall

#include "tableA1.h"

#include "../../include/utils.h"
#include "../../include/dataset.h"

#include <iostream>
#include <sstream>

using namespace std;

dataset dataset_cms_13tev(double rapidity) {
   using namespace cmsbph15005tableA1;

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
      ans.set_graph(graph_tableA1(0,true), graph_tableA1(0,false));
   } else if (absrap<0.6) {
      ans.set_raprange(range(0.3,0.6));
      ans.set_graph(graph_tableA1(1,true), graph_tableA1(1,false));
   } else if (absrap<0.9) {
      ans.set_raprange(range(0.6,0.9));
      ans.set_graph(graph_tableA1(2,true), graph_tableA1(2,false));
   } else if (absrap<1.2) {
      ans.set_raprange(range(0.9,1.2));
      ans.set_graph(graph_tableA1(3,true), graph_tableA1(3,false));
   } else {
      ans.set_raprange(range(0,1.2));
      ans.set_graph(graph_tableA1(4,true), graph_tableA1(4,false));
   }

   ostringstream oss;
   oss << "CMS, 13TeV, " << ans.get_raprange().min << "<|y|<" << ans.get_raprange().max;
   ans.set_legend(oss.str());

   return ans;
}


#endif // #ifndef cmsbph15005_tableall
