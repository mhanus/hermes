#include "problem_data.h"
#include <fstream>
#include <algorithm>

ProblemData::ProblemData(const std::string& datafile)
{
  std::ifstream df(datafile.c_str(), std::ifstream::in);
  
  double sigma;
  double data[G];
  std::string mat_name;
  
  int sz = 0;
  while (df >> ws >> mat_name)
  {
    std::string tmp;
    df >> tmp;
    
    tmp.erase(std::remove( tmp.begin(), tmp.end(), '\'' ), tmp.end());
    
    if (tmp.length() > 0)
      mat_name += "_" + tmp;
    
    std::cout << "Loading data for " << mat_name << std::endl;
    
    mat_list.insert(mat_name);
    assert(mat_list.size() == ++sz);
        
    df >> data[0] >> sigma >> data[1] >> sigma;
    D[mat_name] = rank1(data, data+G);
    df >> data[0] >> sigma >> data[1] >> sigma;
    nu[mat_name] = rank1(data, data+G);
    df >> data[0] >> sigma >> data[1] >> sigma;    
    Sf[mat_name] = rank1(data, data+G);
    df >> data[0] >> sigma >> data[1] >> sigma;    
    Sa[mat_name] = rank1(data, data+G);
    
    Ss[mat_name] = rank2(G);
    for (int g = 0; g < G; g++)
      Ss[mat_name][g] = rank1(G);
    
    df >> Ss[mat_name][0][0] >> sigma >> Ss[mat_name][1][0] >> sigma 
       >> Ss[mat_name][0][1] >> sigma >> Ss[mat_name][1][1] >> sigma;
  }
  
  df.close();
}
