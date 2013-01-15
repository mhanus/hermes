#include "hermes2d.h"
#include <istream>
using namespace Hermes::Hermes2D;
using namespace Hermes::Mixins;

class ConstantableSpacesVector
{
  public:
    ConstantableSpacesVector(Hermes::vector<Space<double> *>* spaces_)
    {
      set(spaces_);
    }
    
    ~ConstantableSpacesVector()
    {
      delete constant;
    }
    
    void set(Hermes::vector<Space<double> *>* spaces_)
    {
      non_constant = spaces_;
      constant = new Hermes::vector<const Space<double> *>(spaces_->size());
      for (Hermes::vector<Space<double> *>::const_iterator it = spaces_->begin(); it != spaces_->end(); ++it)
        constant->push_back(*it);
    }
    
    ConstantableSpacesVector& operator=(const ConstantableSpacesVector& other)
    {
      if (this != &other)
      {
        delete constant;
        
        for (int i = 0; i < non_constant->size(); i++)
          delete non_constant->at(i);
        
        this->set(&other.get());
      }
      
      return *this;
    }
    
    Hermes::vector<const Space<double> *>& get_const() const
    {
      return *constant;
    }
    
    Hermes::vector<Space<double> *>& get() const
    {
      return *non_constant;
    }
    
  private:
    Hermes::vector<const Space<double> *>* constant;
    Hermes::vector<Space<double> *>* non_constant;
};

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

void report_num_dof(const std::string& msg, const Hermes::vector<Space<double> *> spaces);
void report_errors(const std::string& msg, const Hermes::vector< double > errors);
std::string itos(int t);
