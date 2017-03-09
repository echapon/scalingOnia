#include "all.C"

void plot1() {
   vector<dataset> data;
   dataset data1;
   data1 = dataset_atlas_7tev(0);
   data.push_back(data1); 

   dataset data_interpol_lin = dataset_atlas_8tev(0);
   data_interpol_lin.interpolate(20,lin,true);
   data_interpol_lin.set_legend("linear interpolation");
   data.push_back(data_interpol_lin);

   dataset data_interpol_cspline = dataset_atlas_8tev(0);
   data_interpol_cspline.interpolate(20,cspline,true);
   data_interpol_cspline.set_legend("cspline interpolation");
   data.push_back(data_interpol_cspline);

   dataset data_interpol_loglin = dataset_atlas_8tev(0);
   data_interpol_loglin.interpolate(20,loglin,true);
   data_interpol_loglin.set_legend("loglin interpolation");
   data.push_back(data_interpol_loglin);

   dataset data_interpol_logcspline = dataset_atlas_8tev(0);
   data_interpol_logcspline.interpolate(20,logcspline,true);
   data_interpol_logcspline.set_legend("logcspline interpolation");
   data.push_back(data_interpol_logcspline);
   
   vector<dataset> theory;
   plot(data,theory);
}
