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
void report_errors(const std::string& msg, const ErrorCalculator<double>& error_calculator);

template<typename dtype> std::string tostr(dtype t)
{
  std::stringstream ss; ss << t;
  return ss.str();
}

void save_algebraic_representation(WeakForm<double>* wf, const Hermes::vector<SpaceSharedPtr<double> >& spaces, const std::string& varname, bool assign_dofs = true);

enum VisualizationOptions 
{
    HERMES_SCALAR_VISUALIZATION   = 0x01,
    HERMES_ANGULAR_VISUALIZATION   = 0x02,
    VTK_SCALAR_VISUALIZATION   = 0x04,
    VTK_ANGULAR_VISUALIZATION   = 0x08
};

void load_solution(const std::string& sln_file, 
                   const Hermes::vector<SpaceSharedPtr<double> >& spaces, 
                   const Neutronics::Common::MaterialProperties::MaterialPropertyMaps* matprop,
                   VisualizationOptions visualization, bool mode_3D = false);

void load_solution(const std::string& sln_file, 
                   const Hermes::vector<SpaceSharedPtr<double> >& spaces, 
                   VisualizationOptions visualization, bool mode_3D = false);
                   
void load_solution(const std::string& sln_file, 
                   const SpaceSharedPtr<double>& space, 
                   VisualizationOptions visualization, bool mode_3D = false);                   