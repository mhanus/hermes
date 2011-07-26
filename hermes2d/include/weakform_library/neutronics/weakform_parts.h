#ifndef ___H2D_NEUTRONICS_WEAK_FORM_PARTS_H
#define ___H2D_NEUTRONICS_WEAK_FORM_PARTS_H

#include "support_classes.h"

namespace Hermes { namespace Hermes2D { namespace Neutronics 
{             
  namespace Diffusion { namespace WeakFormParts
  { 
    using namespace MaterialProperties;
    
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
      class Jacobian : public MatrixFormSurf<double>
      {
        public:
          Jacobian(unsigned int g, GeomType geom_type = HERMES_PLANAR) 
            : MatrixFormSurf<double>(g,g,HERMES_ANY), 
            g(g), geom_type(geom_type)
          {};
          
          Jacobian(unsigned int g, const std::string& area, GeomType geom_type = HERMES_PLANAR) 
            : MatrixFormSurf<double>(g,g,area),
            g(g), geom_type(geom_type)
          {};
          
          template<typename Real>
          Real matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                            Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const;
          
          virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                              Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
            return matrix_form<double>(n, wt, u_ext, u, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return matrix_form<Ord>(n, wt, u_ext, u, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual MatrixFormSurf<double>* clone() {
            return new Jacobian(*this);
          }
                        
        private:
          unsigned int g;
          GeomType geom_type;
      };
      
      class Residual : public VectorFormSurf<double>
      {
        public:
          Residual(unsigned int g, GeomType geom_type = HERMES_PLANAR) 
            : VectorFormSurf<double>(g,HERMES_ANY), 
            g(g), geom_type(geom_type)
          {};
          
          Residual(unsigned int g, const std::string& area, GeomType geom_type = HERMES_PLANAR) 
            : VectorFormSurf<double>(g,area),
            g(g), geom_type(geom_type)
          {};
          
          template<typename Real>
          Real vector_form(int n, double *wt, Func<Real> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const;
          
          virtual double value(int n, double *wt, Func<double> *u_ext[],
                              Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
            return vector_form<double>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual VectorFormSurf<double>* clone() {
            return new Residual(*this);
          }
                        
        private:
          unsigned int g;
          GeomType geom_type;
      };
    };
    
    struct DiffusionReaction
    {   
      class Jacobian : public MatrixFormVol<double>, protected GenericForm
      {
        public:            
          Jacobian(unsigned int g, 
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
            : MatrixFormVol<double>(g, g, HERMES_ANY, HERMES_SYM),
              GenericForm(matprop, geom_type),
              g(g)
          {};
              
          Jacobian(unsigned int g, const std::string& area,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
            : MatrixFormVol<double>(g, g, area, HERMES_SYM),
              GenericForm(matprop, geom_type),
              g(g)
          {};
          
          template<typename Real>
          Real matrix_form( int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                              Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext  ) const;

          virtual double value( int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                                Func<double> *v, Geom<double> *e, ExtData<double> *ext  ) const { 
            return  matrix_form<double> (n, wt, u_ext, u, v, e, ext);
          }
          
          virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
            return  matrix_form<Ord> (n, wt, u_ext, u, v, e, ext);
          }

          // This is to make the form usable in rk_time_step().
          virtual MatrixFormVol<double>* clone() {
            return new Jacobian(*this);
          }

        private:
          
          unsigned int g;
      };
      
      class Residual : public VectorFormVol<double>, protected GenericForm
      {
        public:
          
          Residual(unsigned int g, 
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
            : VectorFormVol<double>(g, HERMES_ANY),
              GenericForm(matprop, geom_type),
              g(g)
          {};
              
          Residual(unsigned int g, const std::string& area,
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
            : VectorFormVol<double>(g, area),
              GenericForm(matprop, geom_type), 
              g(g)
          {};
          
          template<typename Real>
          Real vector_form(int n, double *wt, Func<Real> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const;
          
          virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<double> *ext) const  {
            return vector_form<double>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const  {
            return vector_form<Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual VectorFormVol<double>* clone() {
            return new Residual(*this);
          }
          
        private:
          
          unsigned int g;
      };
    };
  
    struct FissionYield
    {
      class Jacobian : public MatrixFormVol<double>, protected GenericForm
      {
        public:
          
          Jacobian( unsigned int gto, unsigned int gfrom, 
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
            : MatrixFormVol<double>(gto, gfrom), 
              GenericForm(matprop, geom_type),
              gto(gto), gfrom(gfrom)
          {};
          
          Jacobian( unsigned int gto, unsigned int gfrom, const std::string& area,
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
            : MatrixFormVol<double>(gto, gfrom, area), 
              GenericForm(matprop, geom_type),
              gto(gto), gfrom(gfrom)
          {};
          
          template<typename Real>
          Real matrix_form( int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                              Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext  ) const;
          
          virtual double value( int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                                Func<double> *v, Geom<double> *e, ExtData<double> *ext  ) const { 
            return  -1.0 * matrix_form<double> (n, wt, u_ext, u, v, e, ext);
          }
          
          virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
            return  matrix_form<Ord> (n, wt, u_ext, u, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual MatrixFormVol<double>* clone() {
            return new Jacobian(*this);
          }
          
        private:
          
          unsigned int gto, gfrom;
      };
  
      class OuterIterationForm : public VectorFormVol<double>, protected GenericForm
      {
        public:
          
          OuterIterationForm( unsigned int g, 
                              const MaterialPropertyMaps& matprop,
                              const Hermes::vector<MeshFunction<double>*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : VectorFormVol<double>(g, HERMES_ANY, iterates),
              GenericForm(matprop, geom_type),
              g(g), keff(keff)
          {
            if (g >= iterates.size())
              error_function(Messages::E_INVALID_GROUP_INDEX);
          }
          
          OuterIterationForm( unsigned int g, const std::string& area,
                              const MaterialPropertyMaps& matprop,
                              const Hermes::vector<MeshFunction<double>*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : VectorFormVol<double>(g, area, iterates),
              GenericForm(matprop, geom_type),
              g(g), keff(keff)
          {
            if (g >= iterates.size())
              error_function(Messages::E_INVALID_GROUP_INDEX);
          }
          
          OuterIterationForm( unsigned int g, const Hermes::vector<std::string>& areas,
                              const MaterialPropertyMaps& matprop,
                              const Hermes::vector<MeshFunction<double>*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : VectorFormVol<double>(g, areas, iterates),
              GenericForm(matprop, geom_type),
              g(g), keff(keff)
          {
            if (g >= iterates.size())
              error_function(Messages::E_INVALID_GROUP_INDEX);
          }
          
          template<typename Real>
          Real vector_form(int n, double *wt, Func<Real> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const;

          virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                              Geom<double> *e, ExtData<double> *ext) const {
            return -1.0 * vector_form<double>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const  {
            return vector_form<Ord>(n, wt, u_ext, v, e, ext);
          }

          // This is to make the form usable in rk_time_step().
          virtual VectorFormVol<double>* clone() {
            return new OuterIterationForm(*this);
          }
          
          void update_keff(double new_keff) { keff = new_keff; }
          
        private:
          
          unsigned int g;
          double keff;
      };
    
      class Residual : public VectorFormVol<double>, protected GenericForm
      {
        public:
          Residual( unsigned int gto, unsigned int gfrom, 
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
            : VectorFormVol<double>(gto), 
              GenericForm(matprop, geom_type),
              gto(gto), gfrom(gfrom)
          {};
          
          Residual( unsigned int gto, unsigned int gfrom, const std::string& area,
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
            : VectorFormVol<double>(gto, area), 
              GenericForm(matprop, geom_type),
              gto(gto), gfrom(gfrom)
          {};
          
          template<typename Real>
          Real vector_form(int n, double *wt, Func<Real> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const;
          
          virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<double> *ext) const {
            return -1.0 * vector_form<double>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual VectorFormVol<double>* clone() {
            return new Residual(*this);
          }
          
        private:
          
          unsigned int gto, gfrom;
      };
    };

    struct Scattering
    {      
      class Jacobian : public MatrixFormVol<double>, protected GenericForm
      {
        public:
          
          Jacobian( unsigned int gto, unsigned int gfrom, 
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
            : MatrixFormVol<double>(gto, gfrom), 
              GenericForm(matprop, geom_type),
              gto(gto), gfrom(gfrom)
          {
            this->scaling_factor = (gto != gfrom) ? -1 : 0;
          };
          
          Jacobian( unsigned int gto, unsigned int gfrom, const std::string& area,
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
            : MatrixFormVol<double>(gto, gfrom, area), 
              GenericForm(matprop, geom_type),
              gto(gto), gfrom(gfrom)
          {
            this->scaling_factor = (gto != gfrom) ? -1 : 0;
          };
          
          template<typename Real>
          Real matrix_form( int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                              Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext  ) const;
          
          virtual double value( int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                                Func<double> *v, Geom<double> *e, ExtData<double> *ext  ) const { 
            return  matrix_form<double> (n, wt, u_ext, u, v, e, ext);
          }
          
          virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
            return  matrix_form<Ord> (n, wt, u_ext, u, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual MatrixFormVol<double>* clone() {
            return new Jacobian(*this);
          }
          
        private:
          
          unsigned int gto, gfrom;
      };
    
      class Residual : public VectorFormVol<double>, protected GenericForm
      {
        public:
          Residual( unsigned int gto, unsigned int gfrom, 
                    const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR )
            : VectorFormVol<double>(gto), 
              GenericForm(matprop, geom_type),
              gto(gto), gfrom(gfrom)
          {
            this->scaling_factor = (gto != gfrom) ? -1 : 0;
          };
          
          Residual( unsigned int gto, unsigned int gfrom, const std::string& area,
                    const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR )
            : VectorFormVol<double>(gto, area), 
              GenericForm(matprop, geom_type),
              gto(gto), gfrom(gfrom)
          {
            this->scaling_factor = (gto != gfrom) ? -1 : 0;
          };
          
          template<typename Real>
          Real vector_form(int n, double *wt, Func<Real> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const;
          
          virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<double> *ext) const {
            return vector_form<double>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual VectorFormVol<double>* clone() {
            return new Residual(*this);
          }
          
        private:
          
          unsigned int gto, gfrom;
      };
    };
    
    struct ExternalSources
    {
      class LinearForm : public VectorFormVol<double>, protected GenericForm
      {
        public:
          
          LinearForm( unsigned int g, 
                      const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
            : VectorFormVol<double>(g), 
              GenericForm(matprop, geom_type),
              g(g)
          {};
          
          LinearForm( unsigned int g, const std::string& area,
                      const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
            : VectorFormVol<double>(g, area), 
              GenericForm(matprop, geom_type),
              g(g)
          {};
          
          template<typename Real>
          Real vector_form(int n, double *wt, Func<Real> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const;
          
          virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<double> *ext) const {
            return -1.0 * vector_form<double>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord>(n, wt, u_ext, v, e, ext);
          }
                          
          // This is to make the form usable in rk_time_step().
          virtual VectorFormVol<double>* clone() {
            return new LinearForm(*this);
          }
        
        private:
          
          unsigned int g;      
      }; 
    };
    
  /* WeakFormParts */
  }
  /* Diffusion */
  }

  namespace SPN { namespace WeakFormParts
  {
    using namespace MaterialProperties;
    using SupportClasses::Coeffs;
    using SupportClasses::MomentGroupFlattener;
    
    //TODO: Make Diffusion::GenericForm only a GenericForm, which takes pointer to
    // Common::MaterialProperties::MaterialPropertyMaps (hence all "matprop." will 
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
      class Jacobian : protected GenericForm, public MatrixFormSurf<double>
      {
        public:
          Jacobian(unsigned int m, unsigned int n, unsigned int g, 
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(matprop, geom_type),
              MatrixFormSurf<double>(mg.pos(m,g),mg.pos(n,g),HERMES_ANY),
              mrow(m), mcol(n), g(g)
          {};
          
          Jacobian(unsigned int m, unsigned int n,  unsigned int g, const std::string& area, 
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(matprop, geom_type),
              MatrixFormSurf<double>(mg.pos(m,g),mg.pos(n,g),area),
              mrow(m), mcol(n), g(g)
          {};
          
          template<typename Real>
          Real matrix_form( int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                              Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext ) const;
          
          virtual double value( int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                                Func<double> *v, Geom<double> *e, ExtData<double> *ext ) const {
            return matrix_form<double>(n, wt, u_ext, u, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return matrix_form<Ord>(n, wt, u_ext, u, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual MatrixFormSurf<double>* clone() {
            return new Jacobian(*this);
          }
                        
        private:
          
          unsigned int mrow, mcol;
          unsigned int g;
      };
      
      class Residual : protected GenericForm, public VectorFormSurf<double>
      {
        public:
          Residual(unsigned int m, unsigned int N, unsigned int g,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(matprop, geom_type),
              VectorFormSurf<double>(mg.pos(m,g),HERMES_ANY), 
              mrow(m), N_odd((N+1)/2), g(g)
          {};
          
          Residual(unsigned int m, unsigned int N, unsigned int g, const std::string& area, 
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(matprop, geom_type),
              VectorFormSurf<double>(mg.pos(m,g),area), 
              mrow(m), N_odd((N+1)/2), g(g)
          {};
          
          template<typename Real>
          Real vector_form(int n, double *wt, Func<Real> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const;
          
          virtual double value(int n, double *wt, Func<double> *u_ext[],
                              Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
            return vector_form<double>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual VectorFormSurf<double>* clone() {
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
      class Jacobian : protected GenericForm, public MatrixFormVol<double>
      {
        public:            
          Jacobian(unsigned int m, unsigned int g, 
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(matprop, geom_type),
              MatrixFormVol<double>(mg.pos(m,g), mg.pos(m,g), HERMES_ANY, HERMES_SYM),
              mrow(m), g(g)
          {};
              
          Jacobian(unsigned int m, unsigned int g, const std::string& area,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
            : GenericForm(matprop, geom_type),
              MatrixFormVol<double>(mg.pos(m,g), mg.pos(m,g), area, HERMES_SYM),
              mrow(m), g(g)
          {};
          
          template<typename Real>
          Real matrix_form( int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                              Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext  ) const;

          virtual double value( int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                                Func<double> *v, Geom<double> *e, ExtData<double> *ext  ) const { 
            return  matrix_form<double> (n, wt, u_ext, u, v, e, ext);
          }
          
          virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
            return  matrix_form<Ord> (n, wt, u_ext, u, v, e, ext);
          }

          // This is to make the form usable in rk_time_step().
          virtual MatrixFormVol<double>* clone() {
            return new Jacobian(*this);
          }

        private:
          
          unsigned int mrow;
          unsigned int g;
      };
      
      class Residual : protected GenericForm, public VectorFormVol<double>
      {
        public:
          
          Residual(unsigned int m, unsigned int g, 
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(matprop, geom_type),
              VectorFormVol<double>(mg.pos(m,g), HERMES_ANY),
              mrow(m), g(g)
          {};
              
          Residual(unsigned int m, unsigned int g, const std::string& area,
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
            : GenericForm(matprop, geom_type), 
              VectorFormVol<double>(mg.pos(m,g), area),
              mrow(m), g(g)
          {};
          
          template<typename Real>
          Real vector_form(int n, double *wt, Func<Real> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const;
          
          virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<double> *ext) const  {
            return vector_form<double>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const  {
            return vector_form<Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual VectorFormVol<double>* clone() {
            return new Residual(*this);
          }
          
        private:
          
          unsigned int mrow;
          unsigned int g;
      };
    };
  
    struct FissionYield
    {
      class Jacobian : protected GenericForm, public MatrixFormVol<double>
      {
        public:
          
          Jacobian( unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                    const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(matprop, geom_type),
              MatrixFormVol<double>(mg.pos(m,gto), mg.pos(n,gfrom), HERMES_ANY, sym), 
              mrow(m), mcol(n), gto(gto), gfrom(gfrom)
          {};
          
          Jacobian( unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                    const std::string& area, const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(matprop, geom_type),
              MatrixFormVol<double>(mg.pos(m,gto), mg.pos(n,gfrom), area, sym), 
              mrow(m), mcol(n), gto(gto), gfrom(gfrom)
          {};
          
          template<typename Real>
          Real matrix_form( int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                              Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext  ) const;
          
          virtual double value( int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                                Func<double> *v, Geom<double> *e, ExtData<double> *ext  ) const { 
            return  matrix_form<double> (n, wt, u_ext, u, v, e, ext);
          }
          
          virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
            return  matrix_form<Ord> (n, wt, u_ext, u, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual MatrixFormVol<double>* clone() {
            return new Jacobian(*this);
          }
          
        private:
          
          unsigned int mrow, mcol;
          unsigned int gto, gfrom;
      };
  
      class OuterIterationForm : protected GenericForm, public VectorFormVol<double>
      {
        public:
          
          OuterIterationForm( unsigned int m, unsigned int g,
                              const MaterialPropertyMaps& matprop,
                              const Hermes::vector<MeshFunction<double>*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              VectorFormVol<double>(mg.pos(m,g), HERMES_ANY, iterates),
              mrow(m), g(g), keff(keff)
          {};
          
          OuterIterationForm( unsigned int m, unsigned int g, const std::string& area,
                              const MaterialPropertyMaps& matprop,
                              const Hermes::vector<MeshFunction<double>*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              VectorFormVol<double>(mg.pos(m,g), area, iterates),
              mrow(m), g(g), keff(keff)
          {};
          
          OuterIterationForm( unsigned int m, unsigned int g, const Hermes::vector<std::string>& areas,
                              const MaterialPropertyMaps& matprop,
                              const Hermes::vector<MeshFunction<double>*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              VectorFormVol<double>(mg.pos(m,g), areas, iterates),
              mrow(m), g(g), keff(keff)
          {};
          
          template<typename Real>
          Real vector_form(int n, double *wt, Func<Real> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const;

          virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                              Geom<double> *e, ExtData<double> *ext) const {
            return vector_form<double>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const  {
            return vector_form<Ord>(n, wt, u_ext, v, e, ext);
          }

          // This is to make the form usable in rk_time_step().
          virtual VectorFormVol<double>* clone() {
            return new OuterIterationForm(*this);
          }
          
          void update_keff(double new_keff) { keff = new_keff; }
          
        private:
          
          unsigned int mrow;
          unsigned int g;
          double keff;
      };
    
      class Residual : protected GenericForm, public VectorFormVol<double>
      {
        public:
          Residual( unsigned int m, unsigned int N_odd, unsigned int gto,
                    const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              VectorFormVol<double>(mg.pos(m,gto)), 
              mrow(m), N_odd(N_odd), gto(gto)
          {};
          
          Residual( unsigned int m, unsigned int N_odd, unsigned int gto,
                    const std::string& area, const MaterialPropertyMaps& matprop, 
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              VectorFormVol<double>(mg.pos(m,gto), area), 
              mrow(m), N_odd(N_odd), gto(gto)
          {};
          
          template<typename Real>
          Real vector_form(int n, double *wt, Func<Real> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const;
          
          virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<double> *ext) const {
            return vector_form<double>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual VectorFormVol<double>* clone() {
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
      class Jacobian : protected GenericForm, public MatrixFormVol<double>
      {
        public:
          
          Jacobian( unsigned int m, unsigned int gto, unsigned int gfrom,
                    const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(matprop, geom_type),
              MatrixFormVol<double>(mg.pos(m,gto), mg.pos(m,gfrom), HERMES_ANY, sym), 
              mrow(m), gto(gto), gfrom(gfrom)
          {};
          
          Jacobian( unsigned int m, unsigned int gto, unsigned int gfrom,
                    const std::string& area, const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(matprop, geom_type),
              MatrixFormVol<double>(mg.pos(m,gto), mg.pos(m,gfrom), area, sym), 
              mrow(m), gto(gto), gfrom(gfrom)
          {};
          
          template<typename Real>
          Real matrix_form( int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                              Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext  ) const;
          
          virtual double value( int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                                Func<double> *v, Geom<double> *e, ExtData<double> *ext  ) const { 
            return  matrix_form<double> (n, wt, u_ext, u, v, e, ext);
          }
          
          virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
            return  matrix_form<Ord> (n, wt, u_ext, u, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual MatrixFormVol<double>* clone() {
            return new Jacobian(*this);
          }
          
        private:
          
          unsigned int mrow;
          unsigned int gto, gfrom;
      };
    
      class Residual : protected GenericForm, public VectorFormVol<double>
      {
        public:
          Residual( unsigned int m, unsigned int gto,
                    const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              VectorFormVol<double>(mg.pos(m,gto)), 
              mrow(m), gto(gto)
          {};
          
          Residual( unsigned int m, unsigned int gto,
                    const std::string& area, const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              VectorFormVol<double>(mg.pos(m,gto), area), 
              mrow(m), gto(gto)
          {};
          
          template<typename Real>
          Real vector_form(int n, double *wt, Func<Real> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const;
          
          virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<double> *ext) const {
            return vector_form<double>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual VectorFormVol<double>* clone() {
            return new Residual(*this);
          }
          
        private:
          
          unsigned int mrow;
          unsigned int gto;
      };
    };
    
    struct OffDiagonalReactions
    {      
      class Jacobian : protected GenericForm, public MatrixFormVol<double>
      {
        public:
          
          Jacobian( unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                    const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(matprop, geom_type),
              MatrixFormVol<double>(mg.pos(m,gto), mg.pos(n,gfrom), HERMES_ANY, sym), 
              mrow(m), mcol(n), gto(gto), gfrom(gfrom)
          {};
          
          Jacobian( unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                    const std::string& area, const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(matprop, geom_type),
              MatrixFormVol<double>(mg.pos(m,gto), mg.pos(n,gfrom), area, sym), 
              mrow(m), mcol(n), gto(gto), gfrom(gfrom)
          {};
          
          template<typename Real>
          Real matrix_form( int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                              Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext  ) const;
          
          virtual double value( int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                                Func<double> *v, Geom<double> *e, ExtData<double> *ext  ) const { 
            return  matrix_form<double> (n, wt, u_ext, u, v, e, ext);
          }
          
          virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
            return  matrix_form<Ord> (n, wt, u_ext, u, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual MatrixFormVol<double>* clone() {
            return new Jacobian(*this);
          }
          
        private:
          
          unsigned int mrow, mcol;
          unsigned int gto, gfrom;
      };
    
      class Residual : protected GenericForm, public VectorFormVol<double>
      {
        public:
          Residual( unsigned int m, unsigned int N_odd, unsigned int gto,
                    const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              VectorFormVol<double>(mg.pos(m,gto)), 
              mrow(m), N_odd(N_odd), gto(gto)
          {};
          
          Residual( unsigned int m, unsigned int N_odd, unsigned int gto,
                    const std::string& area, const MaterialPropertyMaps& matprop,
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(matprop, geom_type),
              VectorFormVol<double>(mg.pos(m,gto), area), 
              mrow(m), N_odd(N_odd), gto(gto)
          {};
          
          template<typename Real>
          Real vector_form(int n, double *wt, Func<Real> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const;
          
          virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<double> *ext) const {
            return vector_form<double>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord>(n, wt, u_ext, v, e, ext);
          }
          
          // This is to make the form usable in rk_time_step().
          virtual VectorFormVol<double>* clone() {
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
      class LinearForm : protected GenericForm, public VectorFormVol<double>
      {
        public:
          
          LinearForm( unsigned int m, unsigned int g,
                      const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
            : GenericForm(matprop, geom_type),
              VectorFormVol<double>(mg.pos(m,g)), 
              mrow(m), g(g)
          {};
          
          LinearForm( unsigned int m, unsigned int g, const std::string& area,
                      const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
            : GenericForm(matprop, geom_type),
              VectorFormVol<double>(mg.pos(m,g), area), 
              mrow(m), g(g)
          {};
          
          template<typename Real>
          Real vector_form(int n, double *wt, Func<Real> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const;
          
          virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                Geom<double> *e, ExtData<double> *ext) const {
            return vector_form<double>(n, wt, u_ext, v, e, ext);
          }

          virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const {
            return vector_form<Ord>(n, wt, u_ext, v, e, ext);
          }
                          
          // This is to make the form usable in rk_time_step().
          virtual VectorFormVol<double>* clone() {
            return new LinearForm(*this);
          }
        
        private:
          
          unsigned int mrow;
          unsigned int g;      
      }; 
    };
  
  /* WeakFormParts */
  }     
  /* SPN */
  }
  
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
}

#endif