#include "hermes2d.h"
using namespace Hermes::Hermes2D;

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
    
    Hermes::vector<const Space<double> *>& get_const()
    {
      return *constant;
    }
    
    Hermes::vector<Space<double> *>& get()
    {
      return *non_constant;
    }
    
  private:
    Hermes::vector<const Space<double> *>* constant;
    Hermes::vector<Space<double> *>* non_constant;
};

void report_num_dof(const std::string& msg, const Hermes::vector<Space<double> *> spaces);
void report_errors(const std::string& msg, const Hermes::vector< double > errors);
std::string itos(int t);