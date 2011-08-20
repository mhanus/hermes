// Utility functions simplifying some repeating tasks in main.cpp.
//
#ifndef HERMES_REPORT_ALL
  #define HERMES_REPORT_ALL
#endif

#include "hermes2d.h"
using namespace Hermes::Hermes2D; 

void report_num_dof(const std::string& msg, const Hermes::vector< Space<double>* >& spaces);
std::string itos(int t);