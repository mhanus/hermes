#include "hermes2d.h"
using namespace Hermes::Hermes2D;
using namespace Hermes::Mixins;

class ConstantableSpacesVector
{
  public:
    ConstantableSpacesVector(Hermes::vector<SpaceSharedPtr<double> >* spaces_)
    {
      set(spaces_);
    }
    
    ~ConstantableSpacesVector()
    {
      delete constant;
    }
    
    void set(Hermes::vector<SpaceSharedPtr<double> >* spaces_)
    {
      non_constant = spaces_;
      constant = new Hermes::vector<SpaceSharedPtr<double> >(spaces_->size());
      for (Hermes::vector<SpaceSharedPtr<double> >::const_iterator it = spaces_->begin(); it != spaces_->end(); ++it)
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
    
    Hermes::vector<SpaceSharedPtr<double> >& get_const() const
    {
      return *constant;
    }
    
    Hermes::vector<SpaceSharedPtr<double> >& get() const
    {
      return *non_constant;
    }
    
  private:
    Hermes::vector<SpaceSharedPtr<double> >* constant;
    Hermes::vector<SpaceSharedPtr<double> >* non_constant;
};

void report_num_dof(const std::string& msg, const Hermes::vector<SpaceSharedPtr<double> > spaces);
void report_errors(const std::string& msg, const Hermes::vector< double > errors);
std::string itos(int t);
