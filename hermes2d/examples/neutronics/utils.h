#include "hermes2d.h"
#include <istream>
using namespace Hermes::Hermes2D;
using namespace Hermes::Mixins;

class Line : public std::string 
{ 
  friend std::istream & operator>>(std::istream & is, Line & line)
  {   
    return std::getline(is, line);
  }
};

inline double todbl(const std::string& str) { return atof(str.c_str()); }

template<class OutIt>
void read_solution_from_file(std::istream& is, OutIt dest)
{
  typedef std::istream_iterator<Line> InIt;
  std::transform(InIt(is), InIt(), dest, todbl);
}

void report_num_dof(const std::string& msg, const Hermes::vector<SpaceSharedPtr<double> > spaces);
void report_errors(const std::string& msg, const Hermes::vector< double > errors);
std::string itos(int t);
