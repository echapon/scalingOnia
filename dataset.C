#ifndef dataset_C
#define dataset_C

#include "include/dataset.h"

#include <fstream>

void dataset::set_graph(TGraphAsymmErrors *gstat, TGraphAsymmErrors *gsyst) {
   if (gstat) fstat = gstat;
   if (gsyst) fsyst = gsyst;
   if (gstat && gsyst) ftot = combgraph(gstat,gsyst);
   else if (gsyst) ftot = gsyst;
   else ftot = gstat;
}

#endif // ifndef dataset_C
