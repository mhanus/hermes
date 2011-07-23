#ifndef ___H2D_NEUTRONICS_WEAK_FORM_PARTS_H
#define ___H2D_NEUTRONICS_WEAK_FORM_PARTS_H

#include "common_definitions.h"
#include "material_properties.h"
#include "support_classes.h"

namespace Hermes { namespace Hermes2D { namespace Neutronics { namespace WeakFormParts
{             
  namespace Diffusion
  { 
    using namespace MaterialProperties::Diffusion;
    using DataStructures::rank1;
    using DataStructures::rank2;
    
    class GenericForm
    {
      protected:
        const MaterialPropertyMaps& matprop;
        GeomType geom_type;
        
        GenericForm(const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR)
          : matprop(matprop), geom_type(geom_type) 
        {};
    };
    
    struct VacuumBoundaryCondition
    {
      // TODO: General albedo boundary condition.
      class Jacobian : public WeakForm::MatrixFormSurf
      {
        public:
          Jacobian(unsigned int g, GeomType geom_type = HERMES_PLANAR) 
            : WeakForm::MatrixFormSurf(g,g,HERMES_ANY), 
            g(g), geom_type(geom_type)
          {};
          
          Jacobian(unsigned int g, const std::string& area, GeomType geom_type = HERMES_PLANAR) 
            : WeakForm::MatrixFormSurf(g,g,area),
            g(g), geom_type(geom_type)
          {};
          
          template<typename Real, typename Scalar>
          Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
          
          virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                              Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
            return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::MatrixFormSurf* clone() {
            return new Jacobian(*this);
          }
                        
        private:
          unsigned int g;
          GeomType geom_type;
      };
      
      class Residual : public WeakForm::VectorFormSurf
      {
        public:
          Residual(unsigned int g, GeomType geom_type = HERMES_PLANAR) 
            : WeakForm::VectorFormSurf(g,HERMES_ANY), 
            g(g), geom_type(geom_type)
          {};
          
          Residual(unsigned int g, const std::string& area, GeomType geom_type = HERMES_PLANAR) 
            : WeakForm::VectorFormSurf(g,area),
            g(g), geom_type(geom_type)
          {};
          
          template<typename Real, typename Scalar>
          Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
          
          virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
                              Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
            return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::VectorFormSurf* clone() {
            return new Residual(*this);
          }
                        
        private:
          unsigned int g;
          GeomType geom_type;
      };
    };
    
    struct DiffusionReaction
    {   
      class Jacobian : public WeakForm::MatrixFormVol, protected GenericForm
      {
        public:            
          Jacobian(unsigned int g, 
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
            : WeakForm::MatrixFormVol(g, g, HERMES_ANY, HERMES_SYM),
              GenericForm(matprop, geom_type),
              g(g)
          {};
              
          Jacobian(unsigned int g, const std::string& area,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
            : WeakForm::MatrixFormVol(g, g, area, HERMES_SYM),
              GenericForm(matprop, geom_type),
              g(g)
          {};
          
          template<typename Real, typename Scalar>
          Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const;

          virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const { 
            return  matrix_form<double, scalar> (n, wt, u_ext, u, v, e, ext);
          }
          
          virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
            return  matrix_form<Ord, Ord> (n, wt, u_ext, u, v, e, ext);
          }

          // This is to make the form usable in rk_time_step().
          virtual WeakForm::MatrixFormVol* clone() {
            return new Jacobian(*this);
          }

        private:
          
          unsigned int g;
      };
      
      class Residual : public WeakForm::VectorFormVol, protected GenericForm
      {
        public:
          
          Residual(unsigned int g, 
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
            : WeakForm::VectorFormVol(g, HERMES_ANY),
              GenericForm(matprop, geom_type),
              g(g)
          {};
              
          Residual(unsigned int g, const std::string& area,
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
            : WeakForm::VectorFormVol(g, area),
              GenericForm(matprop, geom_type), 
              g(g)
          {};
          
          template<typename Real, typename Scalar>
          Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
          
          virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<scalar> *ext) const  {
            return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const  {
            return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::VectorFormVol* clone() {
            return new Residual(*this);
          }
          
        private:
          
          unsigned int g;
      };
    };
  
    struct FissionYield
    {
      class Jacobian : public WeakForm::MatrixFormVol, protected GenericForm
      {
        public:
          
          Jacobian( unsigned int gto, unsigned int gfrom, 
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
            : WeakForm::MatrixFormVol(gto, gfrom), 
              GenericForm(matprop, geom_type),
              gto(gto), gfrom(gfrom)
          {};
          
          Jacobian( unsigned int gto, unsigned int gfrom, const std::string& area,
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
            : WeakForm::MatrixFormVol(gto, gfrom, area), 
              GenericForm(matprop, geom_type),
              gto(gto), gfrom(gfrom)
          {};
          
          template<typename Real, typename Scalar>
          Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const;
          
          virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const { 
            return  -1.0 * matrix_form<double, scalar> (n, wt, u_ext, u, v, e, ext);
          }
          
          virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
            return  matrix_form<Ord, Ord> (n, wt, u_ext, u, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::MatrixFormVol* clone() {
            return new Jacobian(*this);
          }
          
        private:
          
          unsigned int gto, gfrom;
      };
  
      class OuterIterationForm : public WeakForm::VectorFormVol, protected GenericForm
      {
        public:
          
          OuterIterationForm( unsigned int g, 
                              const MaterialPropertyMaps& matprop,
                              const Hermes::vector<MeshFunction*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : WeakForm::VectorFormVol(g, HERMES_ANY, iterates),
              GenericForm(matprop, geom_type),
              g(g), keff(keff)
          {
            if (g >= iterates.size())
              error(Messages::E_INVALID_GROUP_INDEX);
          }
          
          OuterIterationForm( unsigned int g, const std::string& area,
                              const MaterialPropertyMaps& matprop,
                              const Hermes::vector<MeshFunction*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : WeakForm::VectorFormVol(g, area, iterates),
              GenericForm(matprop, geom_type),
              g(g), keff(keff)
          {
            if (g >= iterates.size())
              error(Messages::E_INVALID_GROUP_INDEX);
          }
          
          OuterIterationForm( unsigned int g, const Hermes::vector<std::string>& areas,
                              const MaterialPropertyMaps& matprop,
                              const Hermes::vector<MeshFunction*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : WeakForm::VectorFormVol(g, areas, iterates),
              GenericForm(matprop, geom_type),
              g(g), keff(keff)
          {
            if (g >= iterates.size())
              error(Messages::E_INVALID_GROUP_INDEX);
          }
          
          template<typename Real, typename Scalar>
          Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

          virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                              Geom<double> *e, ExtData<scalar> *ext) const {
            return -1.0 * vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const  {
            return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
          }

          // This is to make the form usable in rk_time_step().
          virtual WeakForm::VectorFormVol* clone() {
            return new OuterIterationForm(*this);
          }
          
          void update_keff(double new_keff) { keff = new_keff; }
          
        private:
          
          unsigned int g;
          double keff;
      };
    
      class Residual : public WeakForm::VectorFormVol, protected GenericForm
      {
        public:
          Residual( unsigned int gto, unsigned int gfrom, 
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
            : WeakForm::VectorFormVol(gto), 
              GenericForm(matprop, geom_type),
              gto(gto), gfrom(gfrom)
          {};
          
          Residual( unsigned int gto, unsigned int gfrom, const std::string& area,
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
            : WeakForm::VectorFormVol(gto, area), 
              GenericForm(matprop, geom_type),
              gto(gto), gfrom(gfrom)
          {};
          
          template<typename Real, typename Scalar>
          Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
          
          virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<scalar> *ext) const {
            return -1.0 * vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::VectorFormVol* clone() {
            return new Residual(*this);
          }
          
        private:
          
          unsigned int gto, gfrom;
      };
    };

    struct Scattering
    {      
      class Jacobian : public WeakForm::MatrixFormVol, protected GenericForm
      {
        public:
          
          Jacobian( unsigned int gto, unsigned int gfrom, 
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
            : WeakForm::MatrixFormVol(gto, gfrom), 
              GenericForm(matprop, geom_type),
              gto(gto), gfrom(gfrom)
          {
            this->scaling_factor = (gto != gfrom) ? -1 : 0;
          };
          
          Jacobian( unsigned int gto, unsigned int gfrom, const std::string& area,
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
            : WeakForm::MatrixFormVol(gto, gfrom, area), 
              GenericForm(matprop, geom_type),
              gto(gto), gfrom(gfrom)
          {
            this->scaling_factor = (gto != gfrom) ? -1 : 0;
          };
          
          template<typename Real, typename Scalar>
          Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const;
          
          virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const { 
            return  matrix_form<double, scalar> (n, wt, u_ext, u, v, e, ext);
          }
          
          virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
            return  matrix_form<Ord, Ord> (n, wt, u_ext, u, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::MatrixFormVol* clone() {
            return new Jacobian(*this);
          }
          
        private:
          
          unsigned int gto, gfrom;
      };
    
      class Residual : public WeakForm::VectorFormVol, protected GenericForm
      {
        public:
          Residual( unsigned int gto, unsigned int gfrom, 
                    const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR )
            : WeakForm::VectorFormVol(gto), 
              GenericForm(matprop, geom_type),
              gto(gto), gfrom(gfrom)
          {
            this->scaling_factor = (gto != gfrom) ? -1 : 0;
          };
          
          Residual( unsigned int gto, unsigned int gfrom, const std::string& area,
                    const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR )
            : WeakForm::VectorFormVol(gto, area), 
              GenericForm(matprop, geom_type),
              gto(gto), gfrom(gfrom)
          {
            this->scaling_factor = (gto != gfrom) ? -1 : 0;
          };
          
          template<typename Real, typename Scalar>
          Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
          
          virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<scalar> *ext) const {
            return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::VectorFormVol* clone() {
            return new Residual(*this);
          }
          
        private:
          
          unsigned int gto, gfrom;
      };
    };
    
    struct ExternalSources
    {
      class LinearForm : public WeakForm::VectorFormVol, protected GenericForm
      {
        public:
          
          LinearForm( unsigned int g, 
                      const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
            : WeakForm::VectorFormVol(g), 
              GenericForm(matprop, geom_type),
              g(g)
          {};
          
          LinearForm( unsigned int g, const std::string& area,
                      const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
            : WeakForm::VectorFormVol(g, area), 
              GenericForm(matprop, geom_type),
              g(g)
          {};
          
          template<typename Real, typename Scalar>
          Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
          
          virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<scalar> *ext) const {
            return -1.0 * vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
          }
                          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::VectorFormVol* clone() {
            return new LinearForm(*this);
          }
        
        private:
          
          unsigned int g;      
      }; 
    };
        
  }                

  namespace SPN
  {
    using namespace MaterialProperties::SPN;
    using SupportClasses::SPN::Coeffs;
    using SupportClasses::SPN::MomentGroupFlattener;
    using DataStructures::rank1;
    using DataStructures::rank2;
    using DataStructures::rank3;
    
    //TODO: Make Diffusion::GenericForm only a GenericForm, which takes pointer to
    // MaterialProperties::Common::MaterialPropertyMaps (hence all "matprop." will 
    // have to be changed to "matprop->". The following class will then not be needed.
    class GenericForm
    {
      protected:
        const MaterialPropertyMaps& matprop;
        GeomType geom_type;
        MomentGroupFlattener mg;
        
        GenericForm(const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR)
          : matprop(matprop), geom_type(geom_type), mg(matprop.get_G())
        {};
    };
    
    struct VacuumBoundaryCondition
    {
      // TODO: General albedo boundary condition.
      class Jacobian : protected GenericForm, public WeakForm::MatrixFormSurf
      {
        public:
          Jacobian(unsigned int m, unsigned int n, unsigned int g, 
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(matprop, geom_type),
              WeakForm::MatrixFormSurf(mg.pos(m,g),mg.pos(n,g),HERMES_ANY),
              mrow(m), mcol(n), g(g)
          {};
          
          Jacobian(unsigned int m, unsigned int n,  unsigned int g, const std::string& area, 
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(matprop, geom_type),
              WeakForm::MatrixFormSurf(mg.pos(m,g),mg.pos(n,g),area),
              mrow(m), mcol(n), g(g)
          {};
          
          template<typename Real, typename Scalar>
          Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext ) const;
          
          virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                Func<double> *v, Geom<double> *e, ExtData<scalar> *ext ) const {
            return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::MatrixFormSurf* clone() {
            return new Jacobian(*this);
          }
                        
        private:
          
          unsigned int mrow, mcol;
          unsigned int g;
      };
      
      class Residual : protected GenericForm, public WeakForm::VectorFormSurf
      {
        public:
          Residual(unsigned int m, unsigned int N, unsigned int g,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(matprop, geom_type),
              WeakForm::VectorFormSurf(mg.pos(m,g),HERMES_ANY), 
              mrow(m), N_odd((N+1)/2), g(g)
          {};
          
          Residual(unsigned int m, unsigned int N, unsigned int g, const std::string& area, 
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(matprop, geom_type),
              WeakForm::VectorFormSurf(mg.pos(m,g),area), 
              mrow(m), N_odd((N+1)/2), g(g)
          {};
          
          template<typename Real, typename Scalar>
          Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
          
          virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
                              Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
            return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::VectorFormSurf* clone() {
            return new Residual(*this);
          }
                        
        private:
          unsigned int mrow;
          unsigned int N_odd;
          unsigned int g;
      };
    };
    
    struct DiagonalStreamingAndReactions
    {   
      class Jacobian : protected GenericForm, public WeakForm::MatrixFormVol
      {
        public:            
          Jacobian(unsigned int m, unsigned int g, 
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(matprop, geom_type),
              WeakForm::MatrixFormVol(mg.pos(m,g), mg.pos(m,g), HERMES_ANY, HERMES_SYM),
              mrow(m), g(g)
          {};
              
          Jacobian(unsigned int m, unsigned int g, const std::string& area,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
            : GenericForm(matprop, geom_type),
              WeakForm::MatrixFormVol(mg.pos(m,g), mg.pos(m,g), area, HERMES_SYM),
              mrow(m), g(g)
          {};
          
          template<typename Real, typename Scalar>
          Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const;

          virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const { 
            return  matrix_form<double, scalar> (n, wt, u_ext, u, v, e, ext);
          }
          
          virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
            return  matrix_form<Ord, Ord> (n, wt, u_ext, u, v, e, ext);
          }

          // This is to make the form usable in rk_time_step().
          virtual WeakForm::MatrixFormVol* clone() {
            return new Jacobian(*this);
          }

        private:
          
          unsigned int mrow;
          unsigned int g;
      };
      
      class Residual : protected GenericForm, public WeakForm::VectorFormVol
      {
        public:
          
          Residual(unsigned int m, unsigned int g, 
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(matprop, geom_type),
              WeakForm::VectorFormVol(mg.pos(m,g), HERMES_ANY),
              mrow(m), g(g)
          {};
              
          Residual(unsigned int m, unsigned int g, const std::string& area,
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
            : GenericForm(matprop, geom_type), 
              WeakForm::VectorFormVol(mg.pos(m,g), area),
              mrow(m), g(g)
          {};
          
          template<typename Real, typename Scalar>
          Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
          
          virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<scalar> *ext) const  {
            return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const  {
            return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::VectorFormVol* clone() {
            return new Residual(*this);
          }
          
        private:
          
          unsigned int mrow;
          unsigned int g;
      };
    };
  
    struct FissionYield
    {
      class Jacobian : protected GenericForm, public WeakForm::MatrixFormVol
      {
        public:
          
          Jacobian( unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                    const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(matprop, geom_type),
              WeakForm::MatrixFormVol(mg.pos(m,gto), mg.pos(n,gfrom), HERMES_ANY, sym), 
              mrow(m), mcol(n), gto(gto), gfrom(gfrom)
          {};
          
          Jacobian( unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                    const std::string& area, const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(matprop, geom_type),
              WeakForm::MatrixFormVol(mg.pos(m,gto), mg.pos(n,gfrom), area, sym), 
              mrow(m), mcol(n), gto(gto), gfrom(gfrom)
          {};
          
          template<typename Real, typename Scalar>
          Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const;
          
          virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const { 
            return  matrix_form<double, scalar> (n, wt, u_ext, u, v, e, ext);
          }
          
          virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
            return  matrix_form<Ord, Ord> (n, wt, u_ext, u, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::MatrixFormVol* clone() {
            return new Jacobian(*this);
          }
          
        private:
          
          unsigned int mrow, mcol;
          unsigned int gto, gfrom;
      };
  
      class OuterIterationForm : protected GenericForm, public WeakForm::VectorFormVol
      {
        public:
          
          OuterIterationForm( unsigned int m, unsigned int g,
                              const MaterialPropertyMaps& matprop,
                              const Hermes::vector<MeshFunction*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              WeakForm::VectorFormVol(mg.pos(m,g), HERMES_ANY, iterates),
              mrow(m), g(g), keff(keff)
          {};
          
          OuterIterationForm( unsigned int m, unsigned int g, const std::string& area,
                              const MaterialPropertyMaps& matprop,
                              const Hermes::vector<MeshFunction*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              WeakForm::VectorFormVol(mg.pos(m,g), area, iterates),
              mrow(m), g(g), keff(keff)
          {};
          
          OuterIterationForm( unsigned int m, unsigned int g, const Hermes::vector<std::string>& areas,
                              const MaterialPropertyMaps& matprop,
                              const Hermes::vector<MeshFunction*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              WeakForm::VectorFormVol(mg.pos(m,g), areas, iterates),
              mrow(m), g(g), keff(keff)
          {};
          
          template<typename Real, typename Scalar>
          Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

          virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                              Geom<double> *e, ExtData<scalar> *ext) const {
            return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const  {
            return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
          }

          // This is to make the form usable in rk_time_step().
          virtual WeakForm::VectorFormVol* clone() {
            return new OuterIterationForm(*this);
          }
          
          void update_keff(double new_keff) { keff = new_keff; }
          
        private:
          
          unsigned int mrow;
          unsigned int g;
          double keff;
      };
    
      class Residual : protected GenericForm, public WeakForm::VectorFormVol
      {
        public:
          Residual( unsigned int m, unsigned int N_odd, unsigned int gto,
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              WeakForm::VectorFormVol(mg.pos(m,gto)), 
              mrow(m), N_odd(N_odd), gto(gto)
          {};
          
          Residual( unsigned int m, unsigned int N_odd, unsigned int gto,
                    const std::string& area, const MaterialPropertyMaps& matprop, 
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              WeakForm::VectorFormVol(mg.pos(m,gto), area), 
              mrow(m), N_odd(N_odd), gto(gto)
          {};
          
          template<typename Real, typename Scalar>
          Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
          
          virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<scalar> *ext) const {
            return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::VectorFormVol* clone() {
            return new Residual(*this);
          }
          
        private:
          
          unsigned int mrow;
          unsigned int N_odd;
          unsigned int gto;
      };
    };
    
    struct OffDiagonalStreaming
    {      
      class Jacobian : protected GenericForm, public WeakForm::MatrixFormVol
      {
        public:
          
          Jacobian( unsigned int m, unsigned int gto, unsigned int gfrom,
                    const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(matprop, geom_type),
              WeakForm::MatrixFormVol(mg.pos(m,gto), mg.pos(m,gfrom), HERMES_ANY, sym), 
              mrow(m), gto(gto), gfrom(gfrom)
          {};
          
          Jacobian( unsigned int m, unsigned int gto, unsigned int gfrom,
                    const std::string& area, const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(matprop, geom_type),
              WeakForm::MatrixFormVol(mg.pos(m,gto), mg.pos(m,gfrom), area, sym), 
              mrow(m), gto(gto), gfrom(gfrom)
          {};
          
          template<typename Real, typename Scalar>
          Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const;
          
          virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const { 
            return  matrix_form<double, scalar> (n, wt, u_ext, u, v, e, ext);
          }
          
          virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
            return  matrix_form<Ord, Ord> (n, wt, u_ext, u, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::MatrixFormVol* clone() {
            return new Jacobian(*this);
          }
          
        private:
          
          unsigned int mrow;
          unsigned int gto, gfrom;
      };
    
      class Residual : protected GenericForm, public WeakForm::VectorFormVol
      {
        public:
          Residual( unsigned int m, unsigned int gto,
                    const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              WeakForm::VectorFormVol(mg.pos(m,gto)), 
              mrow(m), gto(gto)
          {};
          
          Residual( unsigned int m, unsigned int gto,
                    const std::string& area, const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              WeakForm::VectorFormVol(mg.pos(m,gto), area), 
              mrow(m), gto(gto)
          {};
          
          template<typename Real, typename Scalar>
          Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
          
          virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<scalar> *ext) const {
            return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::VectorFormVol* clone() {
            return new Residual(*this);
          }
          
        private:
          
          unsigned int mrow;
          unsigned int gto;
      };
    };
    
    struct OffDiagonalReactions
    {      
      class Jacobian : protected GenericForm, public WeakForm::MatrixFormVol
      {
        public:
          
          Jacobian( unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                    const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(matprop, geom_type),
              WeakForm::MatrixFormVol(mg.pos(m,gto), mg.pos(n,gfrom), HERMES_ANY, sym), 
              mrow(m), mcol(n), gto(gto), gfrom(gfrom)
          {};
          
          Jacobian( unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                    const std::string& area, const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(matprop, geom_type),
              WeakForm::MatrixFormVol(mg.pos(m,gto), mg.pos(n,gfrom), area, sym), 
              mrow(m), mcol(n), gto(gto), gfrom(gfrom)
          {};
          
          template<typename Real, typename Scalar>
          Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const;
          
          virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const { 
            return  matrix_form<double, scalar> (n, wt, u_ext, u, v, e, ext);
          }
          
          virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
            return  matrix_form<Ord, Ord> (n, wt, u_ext, u, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::MatrixFormVol* clone() {
            return new Jacobian(*this);
          }
          
        private:
          
          unsigned int mrow, mcol;
          unsigned int gto, gfrom;
      };
    
      class Residual : protected GenericForm, public WeakForm::VectorFormVol
      {
        public:
          Residual( unsigned int m, unsigned int N_odd, unsigned int gto,
                    const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              WeakForm::VectorFormVol(mg.pos(m,gto)), 
              mrow(m), N_odd(N_odd), gto(gto)
          {};
          
          Residual( unsigned int m, unsigned int N_odd, unsigned int gto,
                    const std::string& area, const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              WeakForm::VectorFormVol(mg.pos(m,gto), area), 
              mrow(m), N_odd(N_odd), gto(gto)
          {};
          
          template<typename Real, typename Scalar>
          Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
          
          virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<scalar> *ext) const {
            return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::VectorFormVol* clone() {
            return new Residual(*this);
          }
          
        private:
          
          unsigned int mrow;
          unsigned int N_odd;
          unsigned int gto;
      };
    };
    
    struct ExternalSources
    {
      class LinearForm : protected GenericForm, public WeakForm::VectorFormVol
      {
        public:
          
          LinearForm( unsigned int m, unsigned int g,
                      const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
            : GenericForm(matprop, geom_type),
              WeakForm::VectorFormVol(mg.pos(m,g)), 
              mrow(m), g(g)
          {};
          
          LinearForm( unsigned int m, unsigned int g, const std::string& area,
                      const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
            : GenericForm(matprop, geom_type),
              WeakForm::VectorFormVol(mg.pos(m,g), area), 
              mrow(m), g(g)
          {};
          
          template<typename Real, typename Scalar>
          Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
          
          virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<scalar> *ext) const {
            return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
          }
                          
          // This is to make the form usable in rk_time_step().
          virtual WeakForm::VectorFormVol* clone() {
            return new LinearForm(*this);
          }
        
        private:
          
          unsigned int mrow;
          unsigned int g;      
      }; 
    };
  }
        
/* WeakFormParts */
}
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
}

#endif