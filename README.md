# Some code for quarkonium scaling studies.

## How to run the code (example): in ROOT:
- .L all.C+
- dataset d7tev = dataset_atlas_7tev(0.25); 
- dataset d8tev = dataset_atlas_8tev(0.25); 
- vector<dataset> datas; 
- datas.push_back(d7tev); 
- datas.push_back(d8tev); 
- plot(datas,vector<dataset>());
