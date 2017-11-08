#include "all.C"

void plotATLASCMS() {
   vector<dataset> data;
   dataset data7;
   data7 = dataset_atlas_7tev(0.1);
   data.push_back(data7); 
   dataset data72;
   data72 = dataset_cms_7tev(0.1);
   data.push_back(data72); 

   vector<dataset> theory;

   plot(data,theory);
}
