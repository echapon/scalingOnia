#include "dataset.h"

#include <fstream>

void dataset::set_graph(TGraphAsymmErrors *gstat, TGraphAsymmErrors *gsyst) {
   if (gstat) fstat = gstat;
   if (gsyst) fsyst = gsyst;
   if (gstat && gsyst) ftot = combgraph(gstat,gsyst);
}
