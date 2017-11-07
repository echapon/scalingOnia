#include "all.C"

void plotTest() {
   vector<dataset> data;
   dataset data8;
   data8 = dummy_dataset(10,10,100,8000,-4);
   data.push_back(data8); 
   dataset data13;
   data13 = dummy_dataset(10,10,100,13000,-4);
   data.push_back(data13); 

   vector<dataset> theory;
   dataset theory8;
   theory8 = dummy_dataset(10,10,100,8000,-5);
   theory.push_back(theory8); 
   dataset theory13;
   theory13 = dummy_dataset(10,10,100,13000,-5);
   theory.push_back(theory13); 

   plot(data,theory);
}
