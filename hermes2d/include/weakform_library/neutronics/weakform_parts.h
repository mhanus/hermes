#ifndef ___H2D_NEUTRONICS_WEAK_FORM_PARTS_H
#define ___H2D_NEUTRONICS_WEAK_FORM_PARTS_H

#include "support_classes.h"

namespace Hermes { namespace Hermes2D { namespace Neutronics 
{             
  namespace Diffusion { namespace WeakFormParts
  { 
    using namespace MaterialProperties;
        
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
          
          Jacobian(const std::string& area, unsigned int g, GeomType geom_type = HERMES_PLANAR) 
            : MatrixFormSurf<double>(g,g,area),
            g(g), geom_type(geom_type)
          {};
          
          Jacobian(const Hermes::vector<std::string>& areas, unsigned int g, GeomType geom_type = HERMES_PLANAR) 
            : MatrixFormSurf<double>(g,g,areas),
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
          
          Residual(const std::string& area, 
                   unsigned int g, GeomType geom_type = HERMES_PLANAR) 
            : VectorFormSurf<double>(g,area),
            g(g), geom_type(geom_type)
          {};
          
          Residual(const Hermes::vector<std::string>& areas, 
                   unsigned int g, GeomType geom_type = HERMES_PLANAR) 
            : VectorFormSurf<double>(g,areas),
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
      class Jacobian : public MatrixFormVol<double>
      {
        public:            
          Jacobian(unsigned int g, double D, double Sigma_r,
                  GeomType geom_type = HERMES_PLANAR) 
            : MatrixFormVol<double>(g, g, HERMES_ANY, HERMES_SYM),
              D(D), Sigma_r(Sigma_r), geom_type(geom_type)
          {};
              
          Jacobian(const std::string& area, 
                   unsigned int g, double D, double Sigma_r,
                   GeomType geom_type = HERMES_PLANAR)
            : MatrixFormVol<double>(g, g, area, HERMES_SYM),
              D(D), Sigma_r(Sigma_r), geom_type(geom_type)
          {};
          
          Jacobian(const Hermes::vector<std::string>& areas, 
                   unsigned int g, double D, double Sigma_r,
                   GeomType geom_type = HERMES_PLANAR)
            : MatrixFormVol<double>(g, g, areas, HERMES_SYM),
              D(D), Sigma_r(Sigma_r), geom_type(geom_type)
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
          
          double D, Sigma_r;
          GeomType geom_type;
      };
      
      class Residual : public VectorFormVol<double>
      {
        public:
          
          Residual(unsigned int g, double D, double Sigma_r,
                   GeomType geom_type = HERMES_PLANAR) 
            : VectorFormVol<double>(g, HERMES_ANY),
              g(g), D(D), Sigma_r(Sigma_r), geom_type(geom_type)
          {};
              
          Residual(const std::string& area, 
                   unsigned int g, double D, double Sigma_r, 
                   GeomType geom_type = HERMES_PLANAR)
            : VectorFormVol<double>(g, area),
              g(g), D(D), Sigma_r(Sigma_r), geom_type(geom_type)
          {};
          
          Residual(const Hermes::vector<std::string>& areas,
                   unsigned int g, double D, double Sigma_r, 
                   GeomType geom_type = HERMES_PLANAR)
            : VectorFormVol<double>(g, areas),
              g(g), D(D), Sigma_r(Sigma_r), geom_type(geom_type)
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
          double D, Sigma_r;
          GeomType geom_type;
      };
    };
  
    struct FissionYield
    {
      class Jacobian : public MatrixFormVol<double>
      {
        public:
          
          Jacobian( unsigned int gto, unsigned int gfrom,
                    double chi_to, double nu_from, double Sigma_f_from,
                    GeomType geom_type = HERMES_PLANAR )
            : MatrixFormVol<double>(gto, gfrom), 
              gto(gto), gfrom(gfrom), chi_to(chi_to), nu_from(nu_from), Sigma_f_from(Sigma_f_from),
              geom_type(geom_type)
          {};
          
          Jacobian( const std::string& area, 
                    unsigned int gto, unsigned int gfrom, 
                    double chi_to, double nu_from, double Sigma_f_from,
                    GeomType geom_type = HERMES_PLANAR )
            : MatrixFormVol<double>(gto, gfrom, area), 
              gto(gto), gfrom(gfrom), chi_to(chi_to), nu_from(nu_from), Sigma_f_from(Sigma_f_from),
              geom_type(geom_type)
          {};
          
          Jacobian( const Hermes::vector<std::string>& areas,
                    unsigned int gto, unsigned int gfrom, 
                    double chi_to, double nu_from, double Sigma_f_from,
                    GeomType geom_type = HERMES_PLANAR )
            : MatrixFormVol<double>(gto, gfrom, areas), 
              gto(gto), gfrom(gfrom), chi_to(chi_to), nu_from(nu_from), Sigma_f_from(Sigma_f_from),
              geom_type(geom_type)
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
          double chi_to, nu_from, Sigma_f_from;
          GeomType geom_type;
      };
  
      class OuterIterationForm : public VectorFormVol<double>
      {
        public:
          
          OuterIterationForm( unsigned int g, 
                              double chi_to, const rank1& nu, const rank1& Sigma_f, 
                              const Hermes::vector<MeshFunction<double>*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : VectorFormVol<double>(g, HERMES_ANY, iterates),
              g(g), chi_to(chi_to), nu(nu), Sigma_f(Sigma_f),
              keff(keff), geom_type(geom_type)
          {
            if (g >= iterates.size())
              error_function(Messages::E_INVALID_GROUP_INDEX);
          }
          
          OuterIterationForm( const std::string& area, 
                              unsigned int g, 
                              double chi_to, const rank1& nu, const rank1& Sigma_f, 
                              const Hermes::vector<MeshFunction<double>*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : VectorFormVol<double>(g, area, iterates),
              g(g), chi_to(chi_to), nu(nu), Sigma_f(Sigma_f),
              keff(keff), geom_type(geom_type)
          {
            if (g >= iterates.size())
              error_function(Messages::E_INVALID_GROUP_INDEX);
          }
          
          OuterIterationForm( const Hermes::vector<std::string>& areas,
                              unsigned int g, 
                              double chi_to, const rank1& nu, const rank1& Sigma_f, 
                              const Hermes::vector<MeshFunction<double>*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : VectorFormVol<double>(g, areas, iterates),
              g(g), chi_to(chi_to), nu(nu), Sigma_f(Sigma_f),
              keff(keff), geom_type(geom_type)
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
          double chi_to;
          rank1 nu, Sigma_f;
          double keff;
          GeomType geom_type;
      };
    
      class Residual : public VectorFormVol<double>
      {
        public:
          Residual( unsigned int gto, unsigned int gfrom,
                    double chi_to, double nu_from, double Sigma_f_from,
                    GeomType geom_type = HERMES_PLANAR )
            : VectorFormVol<double>(gto), 
              gto(gto), gfrom(gfrom), chi_to(chi_to), nu_from(nu_from), Sigma_f_from(Sigma_f_from),
              geom_type(geom_type)
          {};
          
          Residual( const std::string& area, 
                    unsigned int gto, unsigned int gfrom,
                    double chi_to, double nu_from, double Sigma_f_from,
                    GeomType geom_type = HERMES_PLANAR )
            : VectorFormVol<double>(gto, area), 
              gto(gto), gfrom(gfrom), chi_to(chi_to), nu_from(nu_from), Sigma_f_from(Sigma_f_from),
              geom_type(geom_type)
          {};
          
          Residual( const Hermes::vector<std::string>& areas,
                    unsigned int gto, unsigned int gfrom,
                    double chi_to, double nu_from, double Sigma_f_from,
                    GeomType geom_type = HERMES_PLANAR )
            : VectorFormVol<double>(gto, areas), 
              gto(gto), gfrom(gfrom), chi_to(chi_to), nu_from(nu_from), Sigma_f_from(Sigma_f_from),
              geom_type(geom_type)
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
          double chi_to, nu_from, Sigma_f_from;
          GeomType geom_type;
      };
    };

    struct Scattering
    {      
      class Jacobian : public MatrixFormVol<double>
      {
        public:
          
          Jacobian( unsigned int gto, unsigned int gfrom, double Sigma_s_to_from,
                    GeomType geom_type = HERMES_PLANAR )
            : MatrixFormVol<double>(gto, gfrom), 
              Sigma_s_to_from(Sigma_s_to_from), geom_type(geom_type)
          {
            this->scaling_factor = (gto != gfrom) ? -1 : 0;
          };
          
          Jacobian( const std::string& area,
                    unsigned int gto, unsigned int gfrom, 
                    double Sigma_s_to_from,
                    GeomType geom_type = HERMES_PLANAR )
            : MatrixFormVol<double>(gto, gfrom, area), 
              Sigma_s_to_from(Sigma_s_to_from), geom_type(geom_type)
          {
            this->scaling_factor = (gto != gfrom) ? -1 : 0;
          };
          
          Jacobian( const Hermes::vector<std::string>& areas,
                    unsigned int gto, unsigned int gfrom, 
                    double Sigma_s_to_from,
                    GeomType geom_type = HERMES_PLANAR )
            : MatrixFormVol<double>(gto, gfrom, areas), 
              Sigma_s_to_from(Sigma_s_to_from), geom_type(geom_type)
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
          
          double Sigma_s_to_from;
          GeomType geom_type;
      };
    
      class Residual : public VectorFormVol<double>
      {
        public:
          Residual( unsigned int gto, unsigned int gfrom, 
                    double Sigma_s_to_from,
                    GeomType geom_type = HERMES_PLANAR )
            : VectorFormVol<double>(gto), 
              gfrom(gfrom), Sigma_s_to_from(Sigma_s_to_from), geom_type(geom_type)
          {
            this->scaling_factor = (gto != gfrom) ? -1 : 0;
          };
          
          Residual( const std::string& area,
                    unsigned int gto, unsigned int gfrom, 
                    double Sigma_s_to_from,
                    GeomType geom_type = HERMES_PLANAR )
            : VectorFormVol<double>(gto, area),
              gfrom(gfrom), Sigma_s_to_from(Sigma_s_to_from), geom_type(geom_type)
          {
            this->scaling_factor = (gto != gfrom) ? -1 : 0;
          };
          
          Residual( const Hermes::vector<std::string>& areas,
                    unsigned int gto, unsigned int gfrom, 
                    double Sigma_s_to_from,
                    GeomType geom_type = HERMES_PLANAR )
            : VectorFormVol<double>(gto, areas),
              gfrom(gfrom), Sigma_s_to_from(Sigma_s_to_from), geom_type(geom_type)
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
          
          unsigned int gfrom;
          double Sigma_s_to_from;
          GeomType geom_type;
      };
    };
    
    struct ExternalSources
    {
      class LinearForm : public VectorFormVol<double>
      {
        public:
          
          LinearForm( unsigned int g, double src,
                      GeomType geom_type = HERMES_PLANAR)
            : VectorFormVol<double>(g), 
              src(src), geom_type(geom_type)
          {};
          
          LinearForm( const std::string& area,
                      unsigned int g, double src,
                      GeomType geom_type = HERMES_PLANAR)
            : VectorFormVol<double>(g, area), 
              src(src), geom_type(geom_type)
          {};
          
          LinearForm( const Hermes::vector<std::string>& areas,
                      unsigned int g, double src,
                      GeomType geom_type = HERMES_PLANAR)
            : VectorFormVol<double>(g, areas), 
              src(src), geom_type(geom_type)
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
          
          double src;
          GeomType geom_type;
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
    
    class GenericForm
    {
      protected:
        GeomType geom_type;
        unsigned int G;
        MomentGroupFlattener mg;
        
        GenericForm(double G, GeomType geom_type = HERMES_PLANAR)
          : geom_type(geom_type), G(G), mg(G)
        {};
    };
    
    struct VacuumBoundaryCondition
    {
      // TODO: General albedo boundary condition.
      class Jacobian : protected GenericForm, public MatrixFormSurf<double>
      {
        public:
          Jacobian(unsigned int m, unsigned int n, unsigned int g, 
                  unsigned int G, GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(G, geom_type),
              MatrixFormSurf<double>(mg.pos(m,g),mg.pos(n,g),HERMES_ANY),
              mrow(m), mcol(n), g(g)
          {};
          
          Jacobian(const std::string& area,
                   unsigned int m, unsigned int n, unsigned int g, unsigned int G, 
                   GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(G, geom_type),
              MatrixFormSurf<double>(mg.pos(m,g),mg.pos(n,g),area),
              mrow(m), mcol(n), g(g)
          {};
          
          Jacobian(const Hermes::vector<std::string>& areas,
                   unsigned int m, unsigned int n, unsigned int g, unsigned int G, 
                   GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(G, geom_type),
              MatrixFormSurf<double>(mg.pos(m,g),mg.pos(n,g),areas),
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
                   unsigned int G, GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(G, geom_type),
              VectorFormSurf<double>(mg.pos(m,g),HERMES_ANY), 
              mrow(m), N_odd((N+1)/2), g(g)
          {};
          
          Residual(const std::string& area,
                   unsigned int m, unsigned int N, unsigned int g, unsigned int G, 
                   GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(G, geom_type),
              VectorFormSurf<double>(mg.pos(m,g),area), 
              mrow(m), N_odd((N+1)/2), g(g)
          {};
          
          Residual(const Hermes::vector<std::string>& areas,
                   unsigned int m, unsigned int N, unsigned int g, unsigned int G, 
                   GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(G, geom_type),
              VectorFormSurf<double>(mg.pos(m,g),areas), 
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
          Jacobian(unsigned int m, unsigned int g, unsigned int G,
                   double D, double Sigma_r,
                   GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(G, geom_type),
              MatrixFormVol<double>(mg.pos(m,g), mg.pos(m,g), HERMES_ANY, HERMES_SYM),
              mrow(m), g(g), D(D), Sigma_r(Sigma_r)
          {};
              
          Jacobian(const std::string& area, 
                   unsigned int m, unsigned int g, unsigned int G, 
                   double D, double Sigma_r,
                   GeomType geom_type = HERMES_PLANAR)
            : GenericForm(G, geom_type),
              MatrixFormVol<double>(mg.pos(m,g), mg.pos(m,g), area, HERMES_SYM),
              mrow(m), g(g), D(D), Sigma_r(Sigma_r)
          {};
          
          Jacobian(const Hermes::vector<std::string>& areas,
                   unsigned int m, unsigned int g, unsigned int G, 
                   double D, double Sigma_r,
                   GeomType geom_type = HERMES_PLANAR)
            : GenericForm(G, geom_type),
              MatrixFormVol<double>(mg.pos(m,g), mg.pos(m,g), areas, HERMES_SYM),
              mrow(m), g(g), D(D), Sigma_r(Sigma_r)
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
          double D, Sigma_r;
      };
      
      class Residual : protected GenericForm, public VectorFormVol<double>
      {
        public:
          
          Residual(unsigned int m, unsigned int g, unsigned int G,
                   double D, double Sigma_r,
                   GeomType geom_type = HERMES_PLANAR) 
            : GenericForm(G, geom_type),
              VectorFormVol<double>(mg.pos(m,g), HERMES_ANY),
              mrow(m), g(g), D(D), Sigma_r(Sigma_r)
          {};
              
          Residual(const std::string& area,
                   unsigned int m, unsigned int g, unsigned int G,
                   double D, double Sigma_r, 
                   GeomType geom_type = HERMES_PLANAR)
            : GenericForm(G, geom_type), 
              VectorFormVol<double>(mg.pos(m,g), area),
              mrow(m), g(g), D(D), Sigma_r(Sigma_r)
          {};
          
          Residual(const Hermes::vector<std::string>& areas,
                   unsigned int m, unsigned int g, unsigned int G,
                   double D, double Sigma_r, 
                   GeomType geom_type = HERMES_PLANAR)
            : GenericForm(G, geom_type), 
              VectorFormVol<double>(mg.pos(m,g), areas),
              mrow(m), g(g), D(D), Sigma_r(Sigma_r)
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
          double D, Sigma_r;
      };
    };
  
    struct FissionYield
    {
      class Jacobian : protected GenericForm, public MatrixFormVol<double>
      {
        public:
          
          Jacobian( unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                    unsigned int G, double chi_to, double nu_from, double Sigma_f_from,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(G, geom_type),
              MatrixFormVol<double>(mg.pos(m,gto), mg.pos(n,gfrom), HERMES_ANY, sym), 
              mrow(m), mcol(n), gto(gto), gfrom(gfrom), 
              chi_to(chi_to), nu_from(nu_from), Sigma_f_from(Sigma_f_from)
          {};
          
          Jacobian( const std::string& area, 
                    unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                    unsigned int G, double chi_to, double nu_from, double Sigma_f_from,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(G, geom_type),
              MatrixFormVol<double>(mg.pos(m,gto), mg.pos(n,gfrom), area, sym), 
              mrow(m), mcol(n), gto(gto), gfrom(gfrom),
              chi_to(chi_to), nu_from(nu_from), Sigma_f_from(Sigma_f_from)
          {};
          
          Jacobian( const Hermes::vector<std::string>& areas,
                    unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                    unsigned int G, double chi_to, double nu_from, double Sigma_f_from,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(G, geom_type),
              MatrixFormVol<double>(mg.pos(m,gto), mg.pos(n,gfrom), areas, sym), 
              mrow(m), mcol(n), gto(gto), gfrom(gfrom),
              chi_to(chi_to), nu_from(nu_from), Sigma_f_from(Sigma_f_from)
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
          double chi_to, nu_from, Sigma_f_from;
      };
  
      class OuterIterationForm : protected GenericForm, public VectorFormVol<double>
      {
        public:
          
          OuterIterationForm( unsigned int m, unsigned int g,
                              unsigned int G, double chi_to, const rank1& nu, const rank1& Sigma_f,
                              const Hermes::vector<MeshFunction<double>*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : GenericForm(G, geom_type),
              VectorFormVol<double>(mg.pos(m,g), HERMES_ANY, iterates),
              mrow(m), g(g), keff(keff), chi_to(chi_to), nu(nu), Sigma_f(Sigma_f)
          {};
          
          OuterIterationForm( const std::string& area, 
                              unsigned int m, unsigned int g,  unsigned int G,
                              double chi_to, const rank1& nu, const rank1& Sigma_f,
                              const Hermes::vector<MeshFunction<double>*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : GenericForm(G, geom_type),
              VectorFormVol<double>(mg.pos(m,g), area, iterates),
              mrow(m), g(g), keff(keff), chi_to(chi_to), nu(nu), Sigma_f(Sigma_f)
          {};
          
          OuterIterationForm( const Hermes::vector<std::string>& areas, 
                              unsigned int m, unsigned int g, unsigned int G, 
                              double chi_to, const rank1& nu, const rank1& Sigma_f,
                              const Hermes::vector<MeshFunction<double>*>& iterates,
                              double keff = 1.0,
                              GeomType geom_type = HERMES_PLANAR )
            : GenericForm(G, geom_type),
              VectorFormVol<double>(mg.pos(m,g), areas, iterates),
              mrow(m), g(g), keff(keff), chi_to(chi_to), nu(nu), Sigma_f(Sigma_f)
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
          double chi_to;
          rank1 nu, Sigma_f;
      };
    
      class Residual : protected GenericForm, public VectorFormVol<double>
      {
        public:
          Residual( unsigned int m, unsigned int N_odd, unsigned int gto, unsigned int G, 
                    double chi_to, const rank1& nu, const rank1& Sigma_f,
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(G, geom_type),
              VectorFormVol<double>(mg.pos(m,gto)), 
              mrow(m), N_odd(N_odd), gto(gto), chi_to(chi_to), nu(nu), Sigma_f(Sigma_f)
          {};
          
          Residual( const std::string& area, 
                    unsigned int m, unsigned int N_odd, unsigned int gto, unsigned int G,
                    double chi_to, const rank1& nu, const rank1& Sigma_f,
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(G, geom_type),
              VectorFormVol<double>(mg.pos(m,gto), area), 
              mrow(m), N_odd(N_odd), gto(gto), chi_to(chi_to), nu(nu), Sigma_f(Sigma_f)
          {};
          
          Residual( const Hermes::vector<std::string>& areas,
                    unsigned int m, unsigned int N_odd, unsigned int gto, unsigned int G,
                    double chi_to, const rank1& nu, const rank1& Sigma_f,
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(G, geom_type),
              VectorFormVol<double>(mg.pos(m,gto), areas), 
              mrow(m), N_odd(N_odd), gto(gto), chi_to(chi_to), nu(nu), Sigma_f(Sigma_f)
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
          double chi_to;
          rank1 nu, Sigma_f;
      };
    };
    
    struct OffDiagonalStreaming
    {      
      class Jacobian : protected GenericForm, public MatrixFormVol<double>
      {
        public:
          
          Jacobian( unsigned int m, unsigned int gto, unsigned int gfrom, unsigned int G, 
                    double D,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(G, geom_type),
              MatrixFormVol<double>(mg.pos(m,gto), mg.pos(m,gfrom), HERMES_ANY, sym), 
              mrow(m), gto(gto), gfrom(gfrom), D(D)
          {};
          
          Jacobian( const std::string& area, 
                    unsigned int m, unsigned int gto, unsigned int gfrom, unsigned int G, 
                    double D,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(G, geom_type),
              MatrixFormVol<double>(mg.pos(m,gto), mg.pos(m,gfrom), area, sym), 
              mrow(m), gto(gto), gfrom(gfrom), D(D)
          {};
          
          Jacobian( const Hermes::vector<std::string>& areas,
                    unsigned int m, unsigned int gto, unsigned int gfrom, unsigned int G, 
                    double D,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(G, geom_type),
              MatrixFormVol<double>(mg.pos(m,gto), mg.pos(m,gfrom), areas, sym), 
              mrow(m), gto(gto), gfrom(gfrom), D(D)
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
          double D;
      };
    
      class Residual : protected GenericForm, public VectorFormVol<double>
      {
        public:
          Residual( unsigned int m, unsigned int gto, unsigned int G, 
                    const rank1& D,
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(G, geom_type),
              VectorFormVol<double>(mg.pos(m,gto)), 
              mrow(m), gto(gto), D(D)
          {};
          
          Residual( const std::string& area,
                    unsigned int m, unsigned int gto, unsigned int G, 
                    const rank1& D,
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(G, geom_type),
              VectorFormVol<double>(mg.pos(m,gto), area), 
              mrow(m), gto(gto), D(D)
          {};
          
          Residual( const Hermes::vector<std::string>& areas,
                    unsigned int m, unsigned int gto, unsigned int G, 
                    const rank1& D,
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(G, geom_type),
              VectorFormVol<double>(mg.pos(m,gto), areas), 
              mrow(m), gto(gto), D(D)
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
          rank1 D;
      };
    };
    
    struct OffDiagonalReactions
    {      
      class Jacobian : protected GenericForm, public MatrixFormVol<double>
      {
        public:
          
          Jacobian( unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                    unsigned int G, double Sigma_rn,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(G, geom_type),
              MatrixFormVol<double>(mg.pos(m,gto), mg.pos(n,gfrom), HERMES_ANY, sym), 
              mrow(m), mcol(n), gto(gto), gfrom(gfrom), Sigma_rn(Sigma_rn)
          {};
          
          Jacobian( const std::string& area, 
                    unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                    unsigned int G, double Sigma_rn,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(G, geom_type),
              MatrixFormVol<double>(mg.pos(m,gto), mg.pos(n,gfrom), area, sym), 
              mrow(m), mcol(n), gto(gto), gfrom(gfrom), Sigma_rn(Sigma_rn)
          {};
          
          Jacobian( const Hermes::vector<std::string>& areas,
                    unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                    unsigned int G, double Sigma_rn,
                    GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
            : GenericForm(G, geom_type),
              MatrixFormVol<double>(mg.pos(m,gto), mg.pos(n,gfrom), areas, sym), 
              mrow(m), mcol(n), gto(gto), gfrom(gfrom), Sigma_rn(Sigma_rn)
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
          double Sigma_rn;
      };
    
      class Residual : protected GenericForm, public VectorFormVol<double>
      {
        public:
          Residual( unsigned int m, unsigned int N_odd, unsigned int gto,
                    unsigned int G, const rank3& Sigma_rn,
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(G, geom_type),
              VectorFormVol<double>(mg.pos(m,gto)), 
              mrow(m), N_odd(N_odd), gto(gto), Sigma_rn(Sigma_rn)
          {};
          
          Residual( const std::string& area, 
                    unsigned int m, unsigned int N_odd, unsigned int gto,
                    unsigned int G, const rank3& Sigma_rn,
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(G, geom_type),
              VectorFormVol<double>(mg.pos(m,gto), area), 
              mrow(m), N_odd(N_odd), gto(gto), Sigma_rn(Sigma_rn)
          {};
          
          Residual( const Hermes::vector<std::string>& areas,
                    unsigned int m, unsigned int N_odd, unsigned int gto,
                    unsigned int G, const rank3& Sigma_rn,
                    GeomType geom_type = HERMES_PLANAR )
            : GenericForm(G, geom_type),
              VectorFormVol<double>(mg.pos(m,gto), areas), 
              mrow(m), N_odd(N_odd), gto(gto), Sigma_rn(Sigma_rn)
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
          rank3 Sigma_rn;
      };
    };
    
    struct ExternalSources
    {
      class LinearForm : protected GenericForm, public VectorFormVol<double>
      {
        public:
          
          LinearForm( unsigned int m, unsigned int g, unsigned int G, 
                      double src, 
                      GeomType geom_type = HERMES_PLANAR)
            : GenericForm(G, geom_type),
              VectorFormVol<double>(mg.pos(m,g)), 
              mrow(m), g(g), src(src)
          {};
          
          LinearForm( const std::string& area,
                      unsigned int m, unsigned int g, unsigned int G, 
                      double src, 
                      GeomType geom_type = HERMES_PLANAR)
            : GenericForm(G, geom_type),
              VectorFormVol<double>(mg.pos(m,g), area), 
              mrow(m), g(g), src(src)
          {};
          
          LinearForm( const Hermes::vector<std::string>& areas,
                      unsigned int m, unsigned int g, unsigned int G, 
                      double src, 
                      GeomType geom_type = HERMES_PLANAR)
            : GenericForm(G, geom_type),
              VectorFormVol<double>(mg.pos(m,g), areas), 
              mrow(m), g(g), src(src)
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
          double src;
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