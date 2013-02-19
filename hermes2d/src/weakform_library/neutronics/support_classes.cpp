#include "neutronics/support_classes.h"

namespace Hermes { namespace Hermes2D { namespace Neutronics
{  
  namespace Common { namespace SupportClasses
  {
    void SourceFilter::filter_fn(int n, Hermes::vector<double*> values, double* result)
    {
      int marker = this->get_active_element()->marker;
      std::string material = matprop.get_material(marker, mesh);
      
      memset(result, 0, n*sizeof(double));
      
      if (markers.empty() || markers.find(marker) != markers.end())
      {           
        rank1 Sigma_f = matprop.get_Sigma_f(material);
        rank1 nu = matprop.get_nu(material);
      
        for (int i = 0; i < n; i++) 
          for (unsigned int j = 0; j < values.size(); j++)
            result[i] += nu[j] * Sigma_f[j] * values.at(j)[i];
      }
    }
    
    void SourceFilter::pre_init()
    {
      num = matprop.get_G();
      if(num > 10)
        ErrorHandling::error_function("Unable to create an instance of SourceFilter: Hermes is currently able to handle"
              "only 10 functions in filters.");
      
      for (int i = 0; i < 10; i++)
      {
        item[i] = H2D_FN_VAL & H2D_FN_COMPONENT_0;
        sln[i] = NULL;
      }
      
      num_components = 1;
      have_solutions = false;
    }

    void SourceFilter::assign_solutions(const Hermes::vector< Solution<double>* >& solutions)
    {
      if (solutions.size() != (unsigned) num)
        ErrorHandling::error_function("SourceFilter: Number of solutions does not match the size of data.");
      
      free();
      for (int i = 0; i < num; i++)
        sln[i] = solutions[i];
      init();
      post_init();
    }

    void SourceFilter::assign_solutions(const Hermes::vector< MeshFunction<double>* >& solutions)
    {
      if (solutions.size() != (unsigned) num)
        ErrorHandling::error_function("SourceFilter: Number of solutions does not match the size of data.");
      
      free();
      for (int i = 0; i < num; i++)
        sln[i] = solutions[i];
      init();
      post_init();
    }
    
    void SourceFilter::post_init()
    {
      have_solutions = true;
      
      std::set<std::string>::const_iterator it = source_regions.begin();
      for ( ; it != source_regions.end(); ++it)
        markers.insert(mesh->get_element_markers_conversion().get_internal_marker(*it).marker);
    }
    
    void SourceFilter::set_active_element(Element* e)
    {
      SimpleFilter<double>::set_active_element(e);
      
      order = sln[0]->get_fn_order();
      for (int i = 1; i < num; i++)
        if (sln[i]->get_fn_order() > order)
          order = sln[i]->get_fn_order();
    }

    double SourceFilter::integrate()
    {
      if (!have_solutions)
        return 0.0;
                          
      Quad2D* quad = get_quad_2d(); // Needed for h1_integrate_expression.
      double integral = 0.0;
      Element* e;
                
      for_all_active_elements(e, mesh)
      {
        if (markers.empty() || markers.find(e->marker) != markers.end())
        {
          update_limit_table(e->get_mode());
          this->set_active_element(e);
          RefMap* ru = this->get_refmap();
          int o = this->get_fn_order() + ru->get_inv_ref_order();
          if (geom_type == HERMES_AXISYM_X || geom_type == HERMES_AXISYM_Y)
            o++;
          limit_order(o, e->get_mode());
          this->set_quad_order(o, H2D_FN_VAL);
          double *uval = this->get_fn_values();
          double result = 0.0;
          
          if (geom_type == HERMES_PLANAR)
          {
            h1_integrate_expression(uval[i]);
          }
          else if (geom_type == HERMES_AXISYM_X)
          {
            double* y = ru->get_phys_y(o);
            h1_integrate_expression(y[i] * uval[i]);
          }
          else
          {
            double* x = ru->get_phys_x(o);
            h1_integrate_expression(x[i] * uval[i]);
          }
          
          integral += result;
        }
      }
      
      if (geom_type == HERMES_AXISYM_X || geom_type == HERMES_AXISYM_Y)
        integral *= 2*M_PI;
      
      return integral;
    }
    
    const std::string Visualization::base_title_flux = "Neutron flux: group ";
    const std::string Visualization::base_title_order = "Polynomial orders: group ";
    const std::string Visualization::base_title_mesh = "Core mesh for group ";
    
    void Visualization::init(unsigned int ne, unsigned int ng) 
    {
      n_equations = ne; n_groups = ng;
      
      if (ne > 0 && ng > 0)
      {
        sviews = new Views::ScalarView* [n_equations];
        oviews = new Views::OrderView* [n_equations];
        if (display_meshes)
          mviews = new Views::MeshView* [n_equations];
        else
          mviews = NULL;
      }
      else
      {
        sviews = NULL;
        oviews = NULL;
        mviews = NULL;
      }
    }
    
    Visualization::~Visualization()
    {
      if (sviews != NULL)
      {
#ifndef NOGLUT
        for (unsigned int i = 0; i < n_equations; i++)
          delete sviews[i];
#endif
        delete [] sviews;
      }
      
      if (oviews != NULL)
      {
#ifndef NOGLUT
        for (unsigned int i = 0; i < n_equations; i++)
          delete oviews[i];
#endif
        delete [] oviews;
      }
      
      if (mviews != NULL)
      {
#ifndef NOGLUT
        for (unsigned int i = 0; i < n_equations; i++)
          delete mviews[i];
#endif
        delete [] mviews;
      }
    }

#ifndef NOGLUT
    void Visualization::inspect_meshes(Hermes::vector< Mesh* > meshes)
    {
      if (display_meshes)
      {
        show_meshes(meshes);
        Views::View::wait();
        
        for (unsigned int i = 0; i < n_equations; i++)
          delete mviews[i];
        delete [] mviews;
        
        mviews = NULL;
      }
    }
    
    void Visualization::inspect_solutions(Hermes::vector< Solution<double>* > solutions)
    {
      show_solutions(solutions);
      Views::View::wait();
      
      for (unsigned int i = 0; i < n_equations; i++)
        delete sviews[i];
      delete [] sviews;
      
      sviews = NULL;
    }
    
    void Visualization::inspect_orders(Hermes::vector< Space<double>* > spaces)
    {
      show_orders(spaces);
      Views::View::wait();
      
      for (unsigned int i = 0; i < n_equations; i++)
        delete oviews[i];
      delete [] oviews;
      
      oviews = NULL;
    }
#endif

  /* SupportClasses */
  }
  /* Common */
  }
  
  namespace Diffusion { namespace SupportClasses
  {
#ifndef NOGLUT
    void Visualization::init(unsigned int G, unsigned int width, unsigned int height, bool display_meshes)
    {
      for (unsigned int g = 0; g < n_groups; g++)
      {
        std::string title_flux = base_title_flux + itos(g);
        std::string title_order = base_title_order + itos(g);
        
        sviews[g] = new Views::ScalarView(title_flux.c_str(), new Views::WinGeom(g*(width+2), 0, width, height));
        sviews[g]->show_mesh(false);
        sviews[g]->set_3d_mode(true);
        oviews[g] = new Views::OrderView(title_order.c_str(), new Views::WinGeom(g*(width+2), height+2, width, height));
      }
      
      if (display_meshes)
        for (unsigned int g = 0; g < n_groups; g++)
        {
          std::string title = base_title_mesh + itos(g);
          mviews[g] = new Views::MeshView(title.c_str(), new Views::WinGeom(g*(width+2), 2*(height+2), width, height));
        }
    }
    
    void Visualization::show_meshes(Hermes::vector< Mesh* > meshes)
    {
      if (display_meshes)
        for (unsigned int g = 0; g < n_groups; g++)
          mviews[g]->show(meshes[g]);
    }
    
    void Visualization::show_solutions(Hermes::vector< Solution<double>* > solutions)
    {
      for (unsigned int g = 0; g < n_groups; g++)
        sviews[g]->show(solutions[g]);
    }
    
    void Visualization::show_orders(Hermes::vector< Space<double>* > spaces)
    {
      for (unsigned int g = 0; g < n_groups; g++)
        oviews[g]->show(spaces[g]);
    }
#endif

    void Visualization::save_solutions_vtk(const std::string& base_filename, const std::string& base_varname, 
                                    Hermes::vector< Solution<double>* > solutions, bool mode_3D)
    {
      Views::Linearizer lin;
      for (unsigned int g = 0; g < n_groups; g++)
      {
        std::string appendix = std::string("_group_") + itos(g);
        std::string file = base_filename + appendix + std::string(".vtk");
        std::string var = base_varname + appendix;
        lin.save_solution_vtk(solutions[g], file.c_str(), var.c_str(), mode_3D);
        info("Scalar flux in group %d saved in VTK format to file %s.", g, file.c_str());
      }
    }
    
    void Visualization::save_orders_vtk(const std::string& base_filename, Hermes::vector< Space<double>* > spaces)
    {
      Views::Orderizer ord;
      for (unsigned int g = 0; g < n_groups; g++)
      {
        std::string file = base_filename + std::string("_group_") + itos(g) + std::string(".vtk");
        ord.save_orders_vtk(spaces[g], file.c_str());
        info("Information about approximation space for group %d saved in VTK format to file %s.", g, file.c_str());
      }
    }
    
  /* SupportClasses */
  }
  /* Diffusion */
  }
  
  namespace SPN { namespace SupportClasses
  {
    const double Coeffs::SYSTEM_MATRIX[N_MAX][N_MAX][N_MAX] =            
    {
      {
        {-1., 0, 0, 0, 0}, 
        {2./3., 0, 0, 0, 0}, 
        {-8./15., 0, 0, 0, 0}, 
        {16./35., 0, 0, 0, 0}, 
        {-128./315., 0, 0, 0, 0}
      },
      {
        {-4./9., -5./9., 0, 0, 0}, 
        {16./45., 4./9., 0, 0, 0}, 
        {-32./105., -8./21., 0, 0, 0}, 
        {256./945., 64./189., 0, 0, 0},
        {0, 0, 0, 0, 0}
      },
      {
        {-64./225., -16./45., -9./25., 0, 0}, 
        {128./525., 32./105., 54./175., 0, 0}, 
        {-1024./4725., -256./945., -48./175., 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0}
      },
      {
        {-256./1225., -64./245., -324./1225., -13./49., 0}, 
        {2048./11025., 512./2205., 288./1225., 104./441., 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0}
      },
      {
        {-16384./99225., -4096./19845., -256./1225., -832./3969., -17./81.},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0}
      }
    };
      
    const double Coeffs::D_GRAD_F[N_MAX][N_MAX] = 
    {
      {-1./2., 1./8., -1./16., 5./128., -7./256.},
      {-7./24, 41./384., -1./16., 131./3072., 0},
      {-407./1920., 233./2560., -1777./30720., 0, 0},
      {3023./17920., 90989./1146880., 0, 0, 0},
      {-1456787./10321920., 0, 0, 0, 0}
    };
    
    const double Coeffs::EVEN_MOMENTS[N_MAX][N_MAX] = 
    {
      {1., -2./3., 8./15., -16./35., 128./315.},
      {1./3., -4./15., 8./35., -64./315., 0},
      {1./5., -6./35., 16./105., 0, 0},
      {1./7., -8./63., 0, 0, 0},
      {1./9., 0, 0, 0, 0}
    };
        
    double Coeffs::system_matrix(unsigned int m, unsigned int n, unsigned int k_of_Sigma_t2k)
    {
      if (m >= N_MAX || n >= N_MAX)
        ErrorHandling::error_function("For the maximum implemented SPN order (N = %d), both m, n must lie in the range [0,%d]."
              "Entered (m,n) = (%d,%d).", 2*N_MAX-1, N_MAX-1, m, n);
      
      if (m > n) std::swap(m,n);
      
      if (k_of_Sigma_t2k > n)
        ErrorHandling::error_function("In the m-th SPN equation, m = %d, the coefficients at the unknown generalized fluxes involve"
              "Sigma_tk with k up to %d. Entered k = %d.", n, 2*n, 2*k_of_Sigma_t2k);
      
      return SYSTEM_MATRIX[m][n-m][k_of_Sigma_t2k];
    }
    
    double Coeffs::D_grad_F(unsigned int m, unsigned int n)
    {
      if (m >= N_MAX || n >= N_MAX)
        ErrorHandling::error_function("For the maximum implemented SPN order (N = %d), both m, n must lie in the range [0,%d]."
              "Entered (m,n) = (%d,%d).", 2*N_MAX-1, N_MAX-1, m, n);
      
      if (m > n) std::swap(m,n); 
      return D_GRAD_F[m][n-m];
    }
    
    double Coeffs::even_moment(unsigned int m, unsigned int n)
    {
      if (m >= N_MAX)
        ErrorHandling::error_function("For the maximum implemented SPN order (N = %d), there are %d even moment(s)."
              "Tried to access moment #%d.", 2*N_MAX-1, N_MAX, m+1);
      if (n > N_MAX || n < m)
        ErrorHandling::error_function("The even moment #%d must be expressed in terms of %d odd moment(s), starting with %d."
              "Tried to use odd moment #%d.", m+1, N_MAX-m, m+1, n+1);
              
      return EVEN_MOMENTS[m][n-m];
    }
    
    void MomentFilter::EvenMomentVal::filter_fn(int n, Hermes::vector< double* > values, double* result)
    {         
      for (int i = 0; i < n; i++) 
      {
        result[i] = 0;
        unsigned int exp_mom_idx = req_mom_idx;
        for (unsigned int sol_idx = mg.pos(exp_mom_idx,g); sol_idx < values.size(); sol_idx = mg.pos(++exp_mom_idx,g))
          result[i] += Coeffs::even_moment(req_mom_idx, exp_mom_idx) * values.at(sol_idx)[i];
      }
    }
    
    void MomentFilter::EvenMomentVal::set_active_element(Element* e)
    {
      SimpleFilter<double>::set_active_element(e);

      order = -1;
      
      unsigned int exp_mom_idx = req_mom_idx;
      for (int sol_idx = mg.pos(exp_mom_idx,g); sol_idx < num; sol_idx = mg.pos(++exp_mom_idx,g))
        if (sln[sol_idx]->get_fn_order() > order)
          order = sln[sol_idx]->get_fn_order();
    }
    
    MeshFunction<double>* MomentFilter::EvenMomentVal::clone() const
    {
      Hermes::vector<MeshFunction<double>*> slns;

      for(int i = 0; i < this->num; i++)
        slns.push_back(this->sln[i]->clone());

      return new EvenMomentVal(angular_moment, group, G, slns);
    }
    
    void MomentFilter::EvenMomentValDxDy::filter_fn(int n, Hermes::vector<double *> values, 
                                          Hermes::vector<double *> dx, Hermes::vector<double *> dy, 
                                          double* rslt, double* rslt_dx, double* rslt_dy)
    {
      for (int i = 0; i < n; i++) 
      {
        rslt[i] = rslt_dx[i] = rslt_dy[i] = 0;
        unsigned int exp_mom_idx = req_mom_idx;
        for (unsigned int sol_idx = mg.pos(exp_mom_idx,g); sol_idx < values.size(); sol_idx = mg.pos(++exp_mom_idx,g))
        {
          rslt[i] += Coeffs::even_moment(req_mom_idx, exp_mom_idx) * values.at(sol_idx)[i];
          rslt_dx[i] += Coeffs::even_moment(req_mom_idx, exp_mom_idx) * dx.at(sol_idx)[i];
          rslt_dy[i] += Coeffs::even_moment(req_mom_idx, exp_mom_idx) * dy.at(sol_idx)[i];
        }
      }
    }
    
    void MomentFilter::EvenMomentValDxDy::set_active_element(Element* e)
    {
      DXDYFilter<double>::set_active_element(e);
      
      order = -1;
      
      unsigned int exp_mom_idx = req_mom_idx;
      for (int sol_idx = mg.pos(exp_mom_idx,g); sol_idx < num; sol_idx = mg.pos(++exp_mom_idx,g))
        if (sln[sol_idx]->get_fn_order() > order)
          order = sln[sol_idx]->get_fn_order();
    }
    
    MeshFunction<double>* MomentFilter::EvenMomentValDxDy::clone() const
    {
      Hermes::vector<MeshFunction<double>*> slns;

      for(int i = 0; i < this->num; i++)
        slns.push_back(this->sln[i]->clone());

      return new EvenMomentValDxDy(angular_moment, group, G, slns);
    }
    
    void MomentFilter::OddMomentVal::set_active_element(Element* e)
    {
      Filter<double>::set_active_element(e);

      order = -1;
      
      for (unsigned int gfrom = 0; gfrom < matprop->get_G(); gfrom++)
        if (sln[mg.pos(req_mom_idx,gfrom)]->get_fn_order() > order)
          order = sln[mg.pos(req_mom_idx,gfrom)]->get_fn_order();
    }
    
    MeshFunction<double>* MomentFilter::OddMomentVal::clone() const
    {
      Hermes::vector<MeshFunction<double>*> slns;

      for(int i = 0; i < this->num; i++)
        slns.push_back(this->sln[i]->clone());

      return new OddMomentVal(component, angular_moment, group, G, slns, matprop);
    }
    
    void MomentFilter::OddMomentVal::precalculate(int order, int mask)
    {
      if (mask & (H2D_FN_DX | H2D_FN_DY | H2D_FN_DXX | H2D_FN_DYY | H2D_FN_DXY))
        ErrorHandling::error_function("MomentFilter::OddMomentVal not defined for derivatives.");
      
      Quad2D* quad = this->quads[this->cur_quad];
      int np = quad->get_num_points(order, this->mode);
      Filter<double>::Node* node = new_node(H2D_FN_VAL_0, np);
      
      double **dx = new double* [this->num], **dy = new double* [this->num];
      for (unsigned int gfrom = 0; gfrom < matprop->get_G(); gfrom++)
      {
        unsigned int i = mg.pos(req_mom_idx,gfrom);
        this->sln[i]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);
        this->sln[i]->get_dx_dy_values(dx[i], dy[i]);
      }
   
      std::string material = matprop->get_material(this->element->marker, this->mesh);
      rank2 D = matprop->get_odd_Sigma_rn_inv(material)[req_mom_idx];
      
      for (int i = 0; i < np; i++)
      {
        node->values[0][0][i] = 0.0;
        for (unsigned int gfrom = 0; gfrom < matprop->get_G(); gfrom++)
        {
          if (component == 0)
            node->values[0][0][i] += -D[g][gfrom] * dx[mg.pos(req_mom_idx,gfrom)][i];
          else
            node->values[0][0][i] += -D[g][gfrom] * dy[mg.pos(req_mom_idx,gfrom)][i];
        }
        node->values[0][0][i] *= Coeffs::D(odd_req_mom);
      }
      
      delete [] dx;
      delete [] dy;
      
      if(this->nodes->present(order)) 
      {
        assert(this->nodes->get(order) == this->cur_node);
        ::free(this->nodes->get(order));
      }
      this->nodes->add(node, order);
      this->cur_node = node;
    }
    
    // TODO: Templatize.
    void MomentFilter::get_scalar_fluxes(const Hermes::vector< Solution<double>* >& angular_fluxes, 
                                          Hermes::vector< MeshFunction<double>* >* scalar_fluxes,
                                          unsigned int G)
    {          
      scalar_fluxes->reserve(G);
      for (unsigned int g = 0; g < G; g++)
        scalar_fluxes->push_back(new MomentFilter::EvenMomentVal(0, g, G, angular_fluxes));
    }
    
    void MomentFilter::get_scalar_fluxes(const Hermes::vector< Solution<double>* >& angular_fluxes, 
                                          Hermes::vector< Filter<double>* >* scalar_fluxes,
                                          unsigned int G)
    {          
      scalar_fluxes->reserve(G);
      for (unsigned int g = 0; g < G; g++)
        scalar_fluxes->push_back(new MomentFilter::EvenMomentVal(0, g, G, angular_fluxes));
    }
    
    void MomentFilter::get_scalar_fluxes_with_derivatives(const Hermes::vector< Solution<double>* >& angular_fluxes, 
                                                          Hermes::vector< MeshFunction<double>* >* scalar_fluxes,
                                                          unsigned int G)
    {          
      scalar_fluxes->reserve(G);
      for (unsigned int g = 0; g < G; g++)
        scalar_fluxes->push_back(new MomentFilter::EvenMomentValDxDy(0, g, G, angular_fluxes));
    }
    
    void MomentFilter::get_scalar_fluxes_with_derivatives(const Hermes::vector< Solution<double>* >& angular_fluxes, 
                                                          Hermes::vector< Filter<double>* >* scalar_fluxes,
                                                          unsigned int G)
    {          
      scalar_fluxes->reserve(G);
      for (unsigned int g = 0; g < G; g++)
        scalar_fluxes->push_back(new MomentFilter::EvenMomentValDxDy(0, g, G, angular_fluxes));
    }
    
    void MomentFilter::clear_scalar_fluxes(Hermes::vector< MeshFunction<double>* >* scalar_fluxes)
    {
      Hermes::vector< MeshFunction<double>* >::const_iterator it = scalar_fluxes->begin();
      for( ; it != scalar_fluxes->end(); ++it)
        delete *it;
      scalar_fluxes->clear();
    }
    
    void MomentFilter::clear_scalar_fluxes(Hermes::vector< Filter<double>* >* scalar_fluxes)
    {
      Hermes::vector< Filter<double>* >::const_iterator it = scalar_fluxes->begin();
      for( ; it != scalar_fluxes->end(); ++it)
        delete *it;
      scalar_fluxes->clear();
    }
    
    void SourceFilter::filter_fn(int n, Hermes::vector< double* > values, double* result)
    {
      int marker = this->get_active_element()->marker;
      std::string material = matprop.get_material(marker, mesh);
      
      memset(result, 0, n*sizeof(double));
      
      if (markers.empty() || markers.find(marker) != markers.end())
      {           
        rank1 Sigma_f = matprop.get_Sigma_f(material);
        rank1 nu = matprop.get_nu(material);
        
        for (int i = 0; i < n; i++) 
        {
          for (unsigned int g = 0; g < G; g++)
          {
            double group_scalar_flux = 0;
            
            unsigned int exp_mom_idx = 0;
            for (unsigned int sol_idx = mg.pos(exp_mom_idx,g); sol_idx < values.size(); sol_idx = mg.pos(++exp_mom_idx,g))
              group_scalar_flux += Coeffs::even_moment(0, exp_mom_idx) * values.at(sol_idx)[i];
            
            result[i] += nu[g] * Sigma_f[g] * group_scalar_flux;
          }
        }
      }
    }
    
    void Visualization::init(unsigned int spn_order, unsigned int G, unsigned int width, unsigned int height, bool display_meshes) 
    {
      n_moments = spn_order+1;
      n_odd_moments = n_moments/2;
      Common::SupportClasses::Visualization::init(G * n_odd_moments, G);
      
#ifndef NOGLUT
      for (unsigned int g = 0; g < n_groups; g++)
      {
        std::string title_flux = base_title_flux + itos(g) + std::string(", pseudo-flux #");
        std::string title_order = base_title_order + itos(g) + std::string(", pseudo-flux #");
        for (unsigned int m = 0; m < n_odd_moments; m++)
        {
          unsigned int i = mg.pos(m,g);
          
          sviews[i] = new Views::ScalarView((title_flux + itos(m)).c_str(), new Views::WinGeom(m*(width+2), g*(height+2), width, height));
          sviews[i]->show_mesh(false);
          sviews[i]->set_3d_mode(true);
          oviews[i] = new Views::OrderView((title_order + itos(m)).c_str(), new Views::WinGeom(m*(width+2), n_groups*(height+2) + g*(height+2), width, height));
        }
      }
      
      if (display_meshes)
      {
        for (unsigned int g = 0; g < n_groups; g++)
        {
          std::string title = base_title_mesh + itos(g) + std::string(", moment ");
          for (unsigned int m = 0; m < n_odd_moments; m++)
            mviews[mg.pos(m,g)] = new Views::MeshView((title + itos(m)).c_str(), new Views::WinGeom(m*(width+2), g*(height+2), width, height));
        }
      }
#endif
    }
    
    Visualization::~Visualization()
    {
      assert((sviews_app == NULL && vviews == NULL) || (sviews_app != NULL && vviews != NULL));
      if (sviews_app != NULL)
      {
#ifndef NOGLUT
        for (unsigned int i = 0; i < n_equations; i++)
        {
          delete sviews_app[i];
          delete vviews[i];
        }
#endif
        delete [] sviews_app;
        delete [] vviews;
      }
    }

#ifndef NOGLUT
    void Visualization::show_meshes(Hermes::vector< Mesh* > meshes)
    {
      if (display_meshes)
        for (unsigned int g = 0; g < n_groups; g++)
          for (unsigned int m = 0; m < n_odd_moments; m++)
            mviews[mg.pos(m,g)]->show(meshes[mg.pos(m,g)]);
    }
    
    void Visualization::show_solutions(Hermes::vector< Solution<double>* > solutions)
    {
      for (unsigned int g = 0; g < n_groups; g++)
        for (unsigned int m = 0; m < n_odd_moments; m++)
          sviews[mg.pos(m,g)]->show(solutions[mg.pos(m,g)]);
    }
    
    void Visualization::show_even_flux_moment(unsigned int moment, unsigned int group, Views::ScalarView* sview,
                                              Hermes::vector< Solution<double>* > solutions)
    {
      if ((moment % 2) != 0 || moment >= n_moments)
      {
        ErrorHandling::warning("Invalid moment specified for Visualization::show_even_flux_moment.");
        return;
      }
      if (group >= n_groups)
      {
        ErrorHandling::warning("Invalid group specified for Visualization::show_even_flux_moment.");
        return;
      }
      MomentFilter::EvenMomentVal mf(moment, group, n_groups, solutions);
      sview->show(&mf);
    }
    
    void Visualization::show_odd_flux_moment(unsigned int moment, unsigned int group, Views::VectorView* vview,
                                             Hermes::vector< Solution<double>* > solutions, const MaterialProperties::MaterialPropertyMaps& matprop)
    {
      if ((moment % 2) != 1 || moment >= n_moments)
      {
        ErrorHandling::warning("Invalid moment specified for Visualization::show_odd_flux_moment.");
        return;
      }
      if (group >= n_groups)
      {
        ErrorHandling::warning("Invalid group specified for Visualization::show_odd_flux_moment.");
        return;
      }
      MomentFilter::OddMomentVal mfx(0, moment, group, n_groups, solutions, &matprop);
      MomentFilter::OddMomentVal mfy(1, moment, group, n_groups, solutions, &matprop);
      vview->show(&mfx, &mfy);
    }
    
    void Visualization::show_all_flux_moments(Hermes::vector< Solution<double>* > solutions, const MaterialProperties::MaterialPropertyMaps& matprop)
    {
      assert((sviews_app == NULL && vviews == NULL) || (sviews_app != NULL && vviews != NULL));
      
      if (sviews_app == NULL)
      {
        sviews_app = new Views::ScalarView* [n_equations];
        vviews = new Views::VectorView* [n_equations];
        
        for (unsigned int g = 0; g < n_groups; g++)
        {
          std::string title_flux = base_title_flux + itos(g) + std::string(", flux moment #");
          for (unsigned int m = 0; m < n_odd_moments; m++)
          {
            unsigned int i = mg.pos(m,g);
            
            sviews_app[i] = new Views::ScalarView((title_flux + itos(2*m)).c_str(), new Views::WinGeom(2*m*(width+2), g*(height+2), width, height));
            sviews_app[i]->show_mesh(false);
            sviews_app[i]->set_3d_mode(true);
            vviews[i] = new Views::VectorView((title_flux + itos(2*m+1)).c_str(), new Views::WinGeom((2*m+1)*(width+2), g*(height+2), width, height));
          }
        }
      }
      
      for (unsigned int g = 0; g < n_groups; g++)
      {
        for (unsigned int m = 0; m < n_odd_moments; m++)
        {
          unsigned int i = mg.pos(m,g);
          show_even_flux_moment(2*m, g, sviews_app[i], solutions);
          show_odd_flux_moment(2*m+1, g, vviews[i], solutions, matprop);
        }
      }
    }
    
    void Visualization::inspect_solutions(Hermes::vector< Solution<double>* > solutions)
    {
      show_solutions(solutions);
      Views::View::wait();
      
      for (unsigned int i = 0; i < n_equations; i++)
        delete sviews[i];
      delete [] sviews;
      
      sviews = NULL;
    }

    void Visualization::show_orders(Hermes::vector< Space<double>* > spaces)
    {
      for (unsigned int g = 0; g < n_groups; g++)
        for (unsigned int m = 0; m < n_odd_moments; m++)
          oviews[mg.pos(m,g)]->show(spaces[mg.pos(m,g)]);
    }
#endif
    
    void Visualization::save_solutions_vtk(const std::string& base_filename, const std::string& base_varname, 
                                    Hermes::vector< Solution<double>* > solutions, bool mode_3D)
    { 
      Views::Linearizer lin;
      for (unsigned int g = 0; g < n_groups; g++)
      {
        std::string appendix = std::string("_group_") + itos(g);
        for (unsigned int m = 0; m < n_odd_moments; m++)
        {
          MomentFilter::EvenMomentVal mf(2*m, g, n_groups, solutions);
          
          std::string file = base_filename + std::string("_moment_") + itos(2*m) + appendix + std::string(".vtk");
          std::string var = base_varname + std::string("_moment_") + itos(2*m) + appendix;
          lin.save_solution_vtk(&mf, file.c_str(), var.c_str(), mode_3D);
          info("SP%d moment #%d of solution in group %d saved in VTK format to file %s.", n_moments-1, 2*m, g, file.c_str());
          
          file = base_filename + std::string("_moment_") + itos(2*m+1) + appendix + std::string(".vtk");
          var = base_varname + std::string("_moment_") + itos(2*m+1) + appendix;
          lin.save_solution_vtk(solutions[mg.pos(m,g)], file.c_str(), var.c_str(), mode_3D);
          info("SP%d moment #%d of solution in group %d saved in VTK format to file %s.", n_moments-1, 2*m+1, g, file.c_str());
        }
      }
    }
    
    void Visualization::save_orders_vtk(const std::string& base_filename, Hermes::vector< Space<double>* > spaces)
    {
      Views::Orderizer ord;
      for (unsigned int g = 0; g < n_groups; g++)
        for (unsigned int m = 0; m < n_odd_moments; m++)
        {
          std::string file = base_filename + std::string("_moment_") + itos(2*m+1) + std::string("_group_") + itos(g) + std::string(".vtk");
          ord.save_orders_vtk(spaces[mg.pos(m,g)], file.c_str());
          info("Information about approximation space for moment %d, group %d saved in VTK format to file %s.", 2*m+1, g, file.c_str());
        }
    }
  
  /* SupportClasses */
  }
  /* SPN */
  }
  
  namespace SN { namespace SupportClasses
  {
    double SphericalHarmonic::plgndr(double x) const
    {
      if (abs(x) > 1)
        ErrorHandling::error_function("Invalid arguments for plgndr. |x| > 1");
      
      double pmm = 1;
      
      if (m > 0)
      {
        double prod = 1;
        for (int i = 1; i <= 2*m; i+=2)
          prod *= i;
        pmm = prod * std::pow( sqrt((1-x)*(1+x)), m );
      }
      
      if (l == m)
        return pmm;
      else
      {
        double pmmp1 = (2*m+1)*x*pmm;
        
        if (l == m+1)
          return pmmp1;
        else
        {
          double pll;
          for (int ll = m+2; ll <= l; ll++)
          {
            pll = ( (2*ll-1)*x*pmmp1-(ll+m-1)*pmm ) / ( ll-m );
            pmm = pmmp1;
            pmmp1 = pll;
          }
          return pll;
        }
      }
    }

    double SphericalHarmonic::operator()(double xi, double eta, double mu) const
    {
      if(std::abs(mu*mu+eta*eta+xi*xi-1.0) > 1.0e-5)
        ErrorHandling::error_function("Invalid direction cosines in SphericalHarmonic::operator()");

      if((l == 0) && (m == 0))
        return 1.0;
      else if((l == 1) && (m == -1))
        return eta;
      else if((l == 1) && (m == 0))
        return mu;
      else if((l == 1) && (m == 1))
        return xi;
      else if((l == 2) && (m == -2))
        return sqrt(3.0)*xi*eta;
      else if((l == 2) && (m == -1))
        return sqrt(3.0)*mu*eta;
      else if((l == 2) && (m == 0))
        return 0.5*(3.0*mu*mu-1.0);
      else if((l == 2) && (m == 1))
        return sqrt(3.0)*mu*xi;
      else if((l == 2) && (m == 2))
        return 0.5*sqrt(3.0)*(xi*xi-eta*eta);
      else if((l == 3) && (m == -3))
        return sqrt(5./8.)*eta*(3.0*xi*xi-eta*eta);
      else if((l == 3) && (m == -2))
        return sqrt(15.0)*xi*eta*mu;
      else if((l == 3) && (m == -1))
        return sqrt(3./8.)*eta*(5.0*mu*mu-1.0);
      else if((l == 3) && (m == 0))
        return 0.5*mu*(5.0*mu*mu-3.0);
      else if((l == 3) && (m == 1))
        return sqrt(3./8.)*xi*(5.0*mu*mu-1.0);
      else if((l == 3) && (m == 2))
        return sqrt(15.0/4.0)*mu*(xi*xi-eta*eta);
      else if((l == 3) && (m == 3))
        return sqrt(5./8.)*xi*(xi*xi-3.0*eta*eta);
      else if((l == 4) && (m == -4))
        return 0.5*sqrt(35.)*xi*eta*(xi*xi-eta*eta);
      else if((l == 4) && (m == -3))
        return 0.5*sqrt(0.5*35.)*mu*eta*(3.*xi*xi-eta*eta);
      else if((l == 4) && (m == -2))
        return sqrt(5.)*(21.*mu*mu-3.)*xi*eta/6.;
      else if((l == 4) && (m == -1))
        return 0.5*sqrt(2.5)*mu*eta*(7.*mu*mu-3.);
      else if((l == 4) && (m == 0))
        return (35.*Hermes::sqr(mu)*Hermes::sqr(mu)-30.*mu*mu+3.)/8.;
      else if((l == 4) && (m == 1))
        return 0.5*sqrt(2.5)*mu*xi*(7.*mu*mu-3.);
      else if((l == 4) && (m == 2))
        return sqrt(5.)*(21.*mu*mu-3.)*(xi*xi-eta*eta)/12.;
      else if((l == 4) && (m == 3))
        return 0.5*sqrt(0.5*35.)*mu*xi*(xi*xi-3.*eta*eta);
      else if((l == 4) && (m == 4))
        return sqrt(35.)*(Hermes::sqr(xi)*Hermes::sqr(xi)-6.*Hermes::sqr(xi*eta)+Hermes::sqr(eta)*Hermes::sqr(eta))/8.;
      else if((l == 5) && (m == -5))
        return 21.*eta*(5.*Hermes::sqr(xi)*Hermes::sqr(xi)-10.*Hermes::sqr(xi*eta)+Hermes::sqr(eta)*Hermes::sqr(eta))/(8.*sqrt(14.));
      else if((l == 5) && (m == -4))
        return 0.5*105.*mu*xi*eta*(xi*xi-eta*eta)/sqrt(35.);
      else if((l == 5) && (m == -3))
        return 35.*(9*mu*mu-1.)*eta*(3.*xi*xi-eta*eta)/(8.*sqrt(70.));
      else if((l == 5) && (m == -2))
        return 0.5*sqrt(105.)*mu*(3.*mu*mu-1.)*xi*eta;
      else if((l == 5) && (m == -1))
        return sqrt(15.)*eta*(21.*Hermes::sqr(mu)*Hermes::sqr(mu)-14.*mu*mu+1.)/8.;
      else if((l == 5) && (m == 0))
        return mu*(63.*Hermes::sqr(mu)*Hermes::sqr(mu)-70.*mu*mu+15.)/8.;
      else if((l == 5) && (m == 1))
        return sqrt(15.)*xi*(21.*Hermes::sqr(mu)*Hermes::sqr(mu)-14.*mu*mu+1.)/8.;
      else if((l == 5) && (m == 2))
        return 0.25*sqrt(105.)*mu*(3.*mu*mu-1.)*(xi*xi-eta*eta);
      else if((l == 5) && (m == 3))
        return 35.*(9*mu*mu-1.)*xi*(xi*xi-3.*eta*eta)/(8.*sqrt(70.));
      else if((l == 5) && (m == 4))
        return 105.*mu*(Hermes::sqr(xi)*Hermes::sqr(xi)-6.*Hermes::sqr(xi*eta)+Hermes::sqr(eta)*Hermes::sqr(eta))/(8.*sqrt(35.));
      else if((l == 5) && (m == 5))
        return 21.*xi*(Hermes::sqr(xi)*Hermes::sqr(xi)-10.*Hermes::sqr(xi*eta)+5.*Hermes::sqr(eta)*Hermes::sqr(eta))/(8.*sqrt(14.));
      else
      { 
        double phi = 1./sqrt(1 - Hermes::sqr(mu)) * acos(xi);
        if (eta < 0)
          phi = 2*M_PI - phi;
        
        if (m > 0)
          return sqrt(2*factorial(l - m) / factorial(l + m) ) * plgndr(mu) * cos(m * phi);
        else if (m == 0)
          return plgndr(mu);
        else
          return sqrt(2*factorial(l + m) / factorial(l - m) ) * plgndr(mu) * sin(-m * phi);
      }
    }

    /* Ordinates in the four octants of the upper hemisphere are ordered as
    *
    *      13         12
    *     17  5     4  16
    *    21  9 1   0 8  20
    *            o
    *    22 10 2   3 11 23
    *     18  6     7  19
    *       14        15
    **/
    OrdinatesData::OrdinatesData(unsigned int N, const std::string& filename) : N(N), M(N*(N+2)/2)
    {
      std::ifstream ifs(filename.c_str());
      
      if (ifs.fail())
        Neutronics::ErrorHandling::error_function("Discrete ordinates could not be loaded from %s.", filename.c_str());
      
      std::string tmp;
      std::getline(ifs, tmp);
      std::getline(ifs, tmp);
      
      int read_n;
      do
      {
        ifs >> read_n;
        if (read_n == N)
          break;
        
        std::getline(ifs, tmp);    
        for (int n = 0; n < read_n; n++)
          std::getline(ifs, tmp);
        std::getline(ifs, tmp);
      } 
      while (!ifs.eof());
      
      if (ifs.eof())
        Neutronics::ErrorHandling::error_function("Required set of discrete ordinates could not be found in %s", filename.c_str());
      
      double *mu_base = new double [N/2];
      
      std::getline(ifs, tmp);    
      for (int n = 0; n < N/2; n++)
        ifs >> mu_base[n];
      
      std::getline(ifs, tmp, '"');
      std::getline(ifs, tmp);
      std::getline(ifs, tmp);
      
      do 
      {
        ifs >> read_n;
        if (read_n == N)
          break;
        
        std::getline(ifs, tmp);    
        for (int n = 0; n < read_n; n++)
          std::getline(ifs, tmp);
        std::getline(ifs, tmp);
      }
      while (!ifs.eof());
      
      if (ifs.eof())
        Neutronics::ErrorHandling::error_function("Required set of weights could not be found in %s.", filename.c_str());
      
      double *wt = new double [N/2];
      
      std::getline(ifs, tmp);    
      for (int n = 0; n < N/2; n++)
        ifs >> wt[n];
    
      ifs.close();
      
      xi.reserve(M);
      eta.reserve(M);
      mu.reserve(M);
      pw.reserve(M);
      reflections_about_x.reserve(M);
      reflections_about_y.reserve(M);
      
      int dir = 0;
      
      for (int n = 1; n <= N/2; n++) // for each polar level
      { 
        for (int i = 1; i <= n; i++) // for each azimuthal point in the first octant
        {  
          double omega = (2*n - 2*i + 1)/(2.*n) * M_PI/2.;
          double xi1q = std::sqrt(1-sqr(mu_base[n-1])) * std::cos(omega);
          double eta1q = std::sqrt(1-sqr(mu_base[n-1])) * std::sin(omega);
          
          // octant 1, ordinate index dir
          xi.push_back( xi1q );
          eta.push_back( eta1q );   
          reflections_about_x.push_back( dir + 3 );
          reflections_about_y.push_back( dir + 1 );
          
          // octant 2, ordinate index dir+1
          xi.push_back(-xi1q );
          eta.push_back( eta1q );   
          reflections_about_x.push_back( dir + 2 );
          reflections_about_y.push_back( dir );
          
          // octant 3, ordinate index dir+2
          xi.push_back(-xi1q );
          eta.push_back(-eta1q );   
          reflections_about_x.push_back( dir + 1 );
          reflections_about_y.push_back( dir + 3 );
          
          // octant 4, ordinate index dir+3
          xi.push_back( xi1q );
          eta.push_back(-eta1q );   
          reflections_about_x.push_back( dir );
          reflections_about_y.push_back( dir + 2 );
          
          for (int j = 0; j < 4; j++)
          {
            mu.push_back(mu_base[n-1]);
            pw.push_back( wt[n-1] / (8*n)); // equal weights in each polar level, summing up to 1 over the whole sphere.
          }
          
          dir+=4;
        }
      }
      
      delete [] mu_base;
      delete [] wt;
    }

    std::ostream& operator<<(std::ostream& os, const OrdinatesData& odata)
    {
      using namespace std;
      
      os << "Discrete ordinates" << endl;
      os << "  N = " << odata.N << endl;
      os << "____________________________________" << endl;

      int m = 0;
      std::vector<double>::const_iterator xi = odata.xi.begin();
      std::vector<double>::const_iterator eta = odata.eta.begin();
      std::vector<double>::const_iterator mu = odata.mu.begin();
      std::vector<double>::const_iterator pw = odata.pw.begin(); 
      for ( ; xi != odata.xi.end(); ++xi, ++eta, ++mu, ++pw)
      { 
        os << endl << *xi << ", " << *eta << ", " << *mu << ", " << *pw << endl;
        os << " --- " << odata.xi[odata.reflections_about_x[m]] << ", " << odata.eta[odata.reflections_about_x[m]] << ", " << odata.mu[odata.reflections_about_x[m]] << endl;
        os << "  |  " << odata.xi[odata.reflections_about_y[m]] << ", " << odata.eta[odata.reflections_about_y[m]] << ", " << odata.mu[odata.reflections_about_y[m]] << endl; 
        m++;
      }
      
      os << endl << m << " ordinates loaded (M = " << odata.M << ")." << endl;
      
      double sum = 0.0;
      for (int n = 0; n < odata.M; n++)
        sum += odata.pw[n];
      
      os << "sum of weights over the whole sphere: " << 2*sum << endl;
      
      const char* pwfile = "pw.m";
      FILE* fp;
      fp = fopen(pwfile, "wt");
      fprintf(fp, "pw = [ \n");
      for (int n = 0; n < odata.M; n++)
        fprintf(fp, "\t%1.15f\n", odata.pw[n]);
      fprintf(fp, "];");
      fclose(fp);
      
      os << "weights written to: " << pwfile << endl << endl;
    }


    template<typename Real>
    void OrdinatesData::ordinates_to_moment(unsigned int l, int m, unsigned int g, unsigned int G, const ExtData< Real >* solution_fns, int num_quad_pts, Real* moment_values_at_quad_pts) const
    {
      SphericalHarmonic Rlm(l, m);
      AngleGroupFlattener ag(G);
      
      for (int quad_pt = 0; quad_pt < num_quad_pts; quad_pt++)
      {
        moment_values_at_quad_pts[quad_pt] = Real(0);
        
        for (int n = 0; n < M; n++)
          moment_values_at_quad_pts[quad_pt] += pw[n] * solution_fns->fn[ag.pos(n,g)]->val[quad_pt] * Rlm(xi[n], eta[n], mu[n]);
      
        moment_values_at_quad_pts[quad_pt] *= 2;       // add contribution from the lower hemisphere
        moment_values_at_quad_pts[quad_pt] *= 4*M_PI;  // make integral of the unity function over the whole sphere equal to 4 PI
      }
    }
    
    template void OrdinatesData::ordinates_to_moment<double>(unsigned int l, int m, unsigned int g, unsigned int G, const ExtData< double >* solution_fns, int num_quad_pts, double* moment_values_at_quad_pts) const;
    template void OrdinatesData::ordinates_to_moment<Hermes::Ord>(unsigned int l, int m, unsigned int g, unsigned int G, const ExtData< Hermes::Ord >* solution_fns, int num_quad_pts, Hermes::Ord* moment_values_at_quad_pts) const;

    template<typename Real>
    void OrdinatesData::ordinates_to_moment(unsigned int l, int m, unsigned int g, unsigned int G, const Hermes::vector< Real* >& solution_values_at_quad_pts, int num_quad_pts, Real *moment_values_at_quad_pts) const
    {
      SphericalHarmonic Rlm = SphericalHarmonic(l, m);
      AngleGroupFlattener ag(G);
      
      for (int quad_pt = 0; quad_pt < num_quad_pts; quad_pt++)
      {
        moment_values_at_quad_pts[quad_pt] = Real(0);
         
        for (int n = 0; n < M; n++)
          moment_values_at_quad_pts[quad_pt] += pw[n] * solution_values_at_quad_pts.at(ag.pos(n,g))[quad_pt] * Rlm(xi[n], eta[n], mu[n]);
      
        moment_values_at_quad_pts[quad_pt] *= 2;  // add contribution from the lower hemisphere
        moment_values_at_quad_pts[quad_pt] *= 4*M_PI;  // make integral of the unity function over the whole sphere equal to 4 PI
      }
    }
    
    template void OrdinatesData::ordinates_to_moment<double>(unsigned int l, int m, unsigned int g, unsigned int G, const Hermes::vector< double* >& solution_values_at_quad_pts, int num_quad_pts, double *moment_values_at_quad_pts) const;
    template void OrdinatesData::ordinates_to_moment<Hermes::Ord>(unsigned int l, int m, unsigned int g, unsigned int G, const Hermes::vector< Hermes::Ord* >& solution_values_at_quad_pts, int num_quad_pts, Hermes::Ord *moment_values_at_quad_pts) const;
    
    void MomentFilter::Val::set_active_element(Element* e)
    {
      SimpleFilter<double>::set_active_element(e);

      order = -1;
      
      unsigned int dir = 0;
      for (int sol_idx = ag.pos(dir,g); sol_idx < num; sol_idx = ag.pos(++dir,g))
        if (sln[sol_idx]->get_fn_order() > order)
          order = sln[sol_idx]->get_fn_order();
    }
    
    MeshFunction<double>* MomentFilter::Val::clone()
    {
      Hermes::vector<MeshFunction<double>*> slns;

      for(int i = 0; i < this->num; i++)
        slns.push_back(this->sln[i]->clone());

      return new Val(l, m, g, G, slns, odata);
    }
    
    void MomentFilter::ValDxDy::set_active_element(Element* e)
    {
      DXDYFilter<double>::set_active_element(e);
      
      order = -1;
      
      unsigned int dir = 0;
      for (int sol_idx = ag.pos(dir,g); sol_idx < num; sol_idx = ag.pos(++dir,g))
        if (sln[sol_idx]->get_fn_order() > order)
          order = sln[sol_idx]->get_fn_order();
    }
    
    MeshFunction<double>* MomentFilter::ValDxDy::clone()
    {
      Hermes::vector<MeshFunction<double>*> slns;

      for(int i = 0; i < this->num; i++)
        slns.push_back(this->sln[i]->clone());

      return new ValDxDy(l, m, g, G, slns, odata);
    }
    
    // TODO: Templatize.
    void MomentFilter::get_scalar_fluxes(const Hermes::vector< Solution<double>* >& angular_fluxes, 
                                          Hermes::vector< MeshFunction<double>* >* scalar_fluxes,
                                          unsigned int G,
                                          const OrdinatesData& odata)
    {          
      scalar_fluxes->reserve(G);
      for (unsigned int g = 0; g < G; g++)
        scalar_fluxes->push_back(new MomentFilter::Val(0, 0, g, G, angular_fluxes, odata));
    }
    
    void MomentFilter::get_scalar_fluxes(const Hermes::vector< Solution<double>* >& angular_fluxes, 
                                          Hermes::vector< Filter<double>* >* scalar_fluxes,
                                          unsigned int G,
                                          const OrdinatesData& odata)
    {          
      scalar_fluxes->reserve(G);
      for (unsigned int g = 0; g < G; g++)
        scalar_fluxes->push_back(new MomentFilter::Val(0, 0, g, G, angular_fluxes, odata));
    }
    
    void MomentFilter::get_scalar_fluxes_with_derivatives(const Hermes::vector< Solution<double>* >& angular_fluxes, 
                                                          Hermes::vector< MeshFunction<double>* >* scalar_fluxes,
                                                          unsigned int G,
                                                          const OrdinatesData& odata)
    {          
      scalar_fluxes->reserve(G);
      for (unsigned int g = 0; g < G; g++)
        scalar_fluxes->push_back(new MomentFilter::ValDxDy(0, 0, g, G, angular_fluxes, odata));
    }
    
    void MomentFilter::get_scalar_fluxes_with_derivatives(const Hermes::vector< Solution<double>* >& angular_fluxes, 
                                                          Hermes::vector< Filter<double>* >* scalar_fluxes,
                                                          unsigned int G,
                                                          const OrdinatesData& odata)
    {          
      scalar_fluxes->reserve(G);
      for (unsigned int g = 0; g < G; g++)
        scalar_fluxes->push_back(new MomentFilter::ValDxDy(0, 0, g, G, angular_fluxes, odata));
    }
    
    void MomentFilter::clear_scalar_fluxes(Hermes::vector< MeshFunction<double>* >* scalar_fluxes)
    {
      Hermes::vector< MeshFunction<double>* >::const_iterator it = scalar_fluxes->begin();
      for( ; it != scalar_fluxes->end(); ++it)
        delete *it;
      scalar_fluxes->clear();
    }
    
    void MomentFilter::clear_scalar_fluxes(Hermes::vector< Filter<double>* >* scalar_fluxes)
    {
      Hermes::vector< Filter<double>* >::const_iterator it = scalar_fluxes->begin();
      for( ; it != scalar_fluxes->end(); ++it)
        delete *it;
      scalar_fluxes->clear();
    }
    
  /* SupportClasses */
  }
  /* SN */
  }
  
  double PostProcessor::integrate(MeshFunction<double>* solution, const Hermes::vector<std::string>& areas) const
  {
    Quad2D* quad = &g_quad_2d_std;
    solution->set_quad_2d(quad);
    const Mesh* mesh = solution->get_mesh();
    
    std::set<int> markers;
    Hermes::vector<std::string>::const_iterator it = areas.begin();
    for ( ; it != areas.end(); ++it)
      markers.insert(mesh->get_element_markers_conversion().get_internal_marker(*it).marker);
    
    double integral = 0.0;
    Element* e;
    for_all_active_elements(e, mesh)
    {
      if (markers.empty() || markers.find(e->marker) != markers.end())
      {
        update_limit_table(e->get_mode());
        solution->set_active_element(e);
        RefMap* ru = solution->get_refmap();
        int o = solution->get_fn_order() + ru->get_inv_ref_order();
        limit_order(o, e->get_mode());
        solution->set_quad_order(o, H2D_FN_VAL);
        double *uval = solution->get_fn_values();
        double result = 0.0;
        
        if (geom_type == HERMES_PLANAR)
        {
          h1_integrate_expression(uval[i]);
        }
        else if (geom_type == HERMES_AXISYM_X)
        {
          double* y = ru->get_phys_y(o);
          h1_integrate_expression(y[i] * uval[i]);
        }
        else
        {
          double* x = ru->get_phys_x(o);
          h1_integrate_expression(x[i] * uval[i]);
        }
        
        integral += result;
      }
    }
    
    if (geom_type == HERMES_AXISYM_X || geom_type == HERMES_AXISYM_Y)
      integral *= 2.0*M_PI;
    
    return integral;
  }

  void PostProcessor::normalize_to_unit_fission_source(Hermes::vector< Solution<double>* >* solutions, 
                                                       double integrated_fission_source) const
  {
    if (integrated_fission_source < 1e-12)
      ErrorHandling::error_function("PostProcessor::normalize_to_unit_fission_source : Invalid fission source.");
     
    Hermes::vector< Solution<double>* >::iterator sln = solutions->begin();
    for ( ; sln != solutions->end(); ++sln)
      (*sln)->multiply(1./integrated_fission_source);
  }

  void PostProcessor::normalize_to_unit_fission_source(Hermes::vector< Solution<double>* >* solutions, 
                                                        const Common::MaterialProperties::MaterialPropertyMaps& matprop) const
  {
    Common::SupportClasses::SourceFilter *sf;
    
    if (method == NEUTRONICS_DIFFUSION)
      sf = new Common::SupportClasses::SourceFilter(*solutions, matprop, geom_type);
    else if (method == NEUTRONICS_SPN)
      sf = new SPN::SupportClasses::SourceFilter(*solutions, matprop, geom_type);
    else if (method == NEUTRONICS_SN)
      // TODO
      ErrorHandling::error_function("SN SourceFilter not yet implemented.");
    
    normalize_to_unit_fission_source(solutions, sf->integrate());
    
    delete sf;
  }

  void PostProcessor::normalize_to_unit_power(Hermes::vector< Solution<double>* >* solutions, 
                                              const Common::MaterialProperties::MaterialPropertyMaps& matprop, 
                                              double power_per_fission) const
  {
    // TODO
  }
  
  double PostProcessor::get_integrated_group_reaction_rates_internal( ReactionType reaction, MeshFunction<double>* solution, 
                                                                      const Common::MaterialProperties::MaterialPropertyMaps& matprop, 
                                                                      const Hermes::vector< std::string >& regions,
                                                                      unsigned int this_group, int other_group) const
  {    
    if (this_group > matprop.get_G())
      ErrorHandling::error_function(Messages::E_INVALID_GROUP_INDEX);
    
    Quad2D* quad = &g_quad_2d_std;
    solution->set_quad_2d(quad);
    const Mesh* mesh = solution->get_mesh();
    
    std::set<int> markers;
    Hermes::vector<std::string>::const_iterator it = regions.begin();
    for ( ; it != regions.end(); ++it)
      markers.insert(mesh->get_element_markers_conversion().get_internal_marker(*it).marker);
            
    double integral = 0.0;
    Element* e;
    for_all_active_elements(e, mesh)
    {
      if (markers.empty() || markers.find(e->marker) != markers.end())
      {
        update_limit_table(e->get_mode());
        solution->set_active_element(e);
        RefMap* ru = solution->get_refmap();
        int o = solution->get_fn_order() + ru->get_inv_ref_order();
        limit_order(o, e->get_mode());
        solution->set_quad_order(o, H2D_FN_VAL);
        double *uval = solution->get_fn_values();
        double result = 0.0;
        
        if (geom_type == HERMES_PLANAR)
        {
          h1_integrate_expression(uval[i]);
        }
        else if (geom_type == HERMES_AXISYM_X)
        {
          double* y = ru->get_phys_y(o);
          h1_integrate_expression(y[i] * uval[i]);
        }
        else
        {
          double* x = ru->get_phys_x(o);
          h1_integrate_expression(x[i] * uval[i]);
        }
        
        std::string mat = matprop.get_material(e->marker, mesh);
        double xsec;
        
        switch (reaction)
        {
          case ABSORPTION:
            xsec = matprop.compute_Sigma_a(mat)[this_group];
            break;
          case TOTAL:
            xsec = matprop.compute_Sigma_t(mat)[this_group];
            break;
          case IN_SCATTERING:
            if (other_group > (int) matprop.get_G() || other_group < 0)
              ErrorHandling::error_function(Messages::E_INVALID_GROUP_INDEX);
            
            xsec = matprop.compute_Sigma_s(mat)[this_group][other_group];
            break;
          case SELF_SCATTERING:
            xsec = matprop.compute_Sigma_s(mat)[this_group][this_group];
            break;
          case OUT_SCATTERING:
            if (other_group > (int) matprop.get_G() || other_group < 0)
              ErrorHandling::error_function(Messages::E_INVALID_GROUP_INDEX);
            
            xsec = matprop.compute_Sigma_s(mat)[other_group][this_group];
            break;
          case FISSION:
            xsec = matprop.get_Sigma_f(mat)[this_group];
            break;
          case NU_FISSION:
            xsec = matprop.get_Sigma_f(mat)[this_group] * matprop.get_nu(mat)[this_group];
            break;
        };
        
        integral += xsec * result;
      }
    }
    
    if (geom_type == HERMES_AXISYM_X || geom_type == HERMES_AXISYM_Y)
      integral *= 2.0*M_PI;
    
    return integral;
  }

  
  void PostProcessor::get_integrated_group_reaction_rates(ReactionType reaction, 
                                                          const Hermes::vector< Solution<double>* >& solutions, Hermes::vector< double >* results, 
                                                          const Common::MaterialProperties::MaterialPropertyMaps& matprop, 
                                                          unsigned int group, const Hermes::vector< std::string >& regions) const
  {
    Hermes::vector<MeshFunction<double>*> scalar_fluxes;
    
    if (method == NEUTRONICS_SPN)
      SPN::SupportClasses::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, matprop.get_G());
    else if (method == NEUTRONICS_SN)
      SN::SupportClasses::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, matprop.get_G(), odata);
    else if (method == NEUTRONICS_DIFFUSION)
      for (Hermes::vector<Solution<double>*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
        scalar_fluxes.push_back(*it);
    
    Hermes::vector<std::string>::const_iterator region = regions.begin(); 

    for ( ; region != regions.end(); ++region)
    {
      double result = 0.0;
      
      if (reaction == IN_SCATTERING || reaction == OUT_SCATTERING)
      {
        for (unsigned int g_other = 0; g_other < matprop.get_G(); g_other++)
        {
          if (reaction == IN_SCATTERING && group != g_other)
            result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[g_other], matprop, *region, group, g_other);
          else if (reaction == OUT_SCATTERING && group != g_other)
            result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[group], matprop, *region, group, g_other);
        }
      }
      else
        result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[group], matprop, *region, group);
      
      results->push_back(result);
    }
    
    if (method == NEUTRONICS_SPN)
      SPN::SupportClasses::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
    else if (method == NEUTRONICS_SN)
      SN::SupportClasses::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
  }
  
  double PostProcessor::get_integrated_group_reaction_rates(ReactionType reaction, const Hermes::vector< Solution<double>* >& solutions, 
                                                            const Common::MaterialProperties::MaterialPropertyMaps& matprop, 
                                                            unsigned int group, const Hermes::vector< std::string >& regions) const
  {
    Hermes::vector<MeshFunction<double>*> scalar_fluxes;
    
    if (method == NEUTRONICS_SPN)
      SPN::SupportClasses::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, matprop.get_G());
    else if (method == NEUTRONICS_SN)
      SN::SupportClasses::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, matprop.get_G(), odata);
    else if (method == NEUTRONICS_DIFFUSION)
      for (Hermes::vector<Solution<double>*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
        scalar_fluxes.push_back(*it);
    
    double result = 0.0;
      
    if (reaction == IN_SCATTERING || reaction == OUT_SCATTERING)
    {
      for (unsigned int g_other = 0; g_other < matprop.get_G(); g_other++)
      {
        if (reaction == IN_SCATTERING && group != g_other)
          result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[g_other], matprop, regions, group, g_other);
        else if (reaction == OUT_SCATTERING && group != g_other)
          result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[group], matprop, regions, group, g_other);
      }
    }
    else
      result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[group], matprop, regions, group);
      
    if (method == NEUTRONICS_SPN)
      SPN::SupportClasses::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
    else if (method == NEUTRONICS_SN)
      SN::SupportClasses::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
    
    return result;
  }
  
  void PostProcessor::get_integrated_reaction_rates(ReactionType reaction, const Hermes::vector< Solution<double>* >& solutions, 
                                                    Hermes::vector< double >* results, 
                                                    const Common::MaterialProperties::MaterialPropertyMaps& matprop, 
                                                    const Hermes::vector< std::string >& regions) const
  {
    Hermes::vector<MeshFunction<double>*> scalar_fluxes;
    
    if (method == NEUTRONICS_SPN)
      SPN::SupportClasses::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, matprop.get_G());
    else if (method == NEUTRONICS_SN)
      SN::SupportClasses::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, matprop.get_G(), odata);
    else if (method == NEUTRONICS_DIFFUSION)
      for (Hermes::vector<Solution<double>*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
        scalar_fluxes.push_back(*it);
    
    Hermes::vector<std::string>::const_iterator region = regions.begin(); 
    
    for ( ; region != regions.end(); ++region)
    {
      double result = 0.0;
      for (unsigned int group = 0; group < matprop.get_G(); group++)
      {
        if (reaction == IN_SCATTERING || reaction == OUT_SCATTERING)
        {
          for (unsigned int g_other = 0; g_other < matprop.get_G(); g_other++)
          {
            if (reaction == IN_SCATTERING && group != g_other)
              result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[g_other], matprop, *region, group, g_other);
            else if (reaction == OUT_SCATTERING && group != g_other)
              result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[group], matprop, *region, group, g_other);
          }
        }
        else
          result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[group], matprop, *region, group);
      }
        
      results->push_back(result);
    }

    if (method == NEUTRONICS_SPN)
      SPN::SupportClasses::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
    else if (method == NEUTRONICS_SN)
      SN::SupportClasses::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
  }
  
  double PostProcessor::get_integrated_reaction_rates(ReactionType reaction, const Hermes::vector< Solution<double>* >& solutions, 
                                                      const Common::MaterialProperties::MaterialPropertyMaps& matprop, 
                                                      const Hermes::vector< std::string >& regions) const
  {
    double result = 0.0;
    for (unsigned int group = 0; group < matprop.get_G(); group++)
      result += get_integrated_group_reaction_rates(reaction, solutions, matprop, group, regions);
    return result;
  }

  void PostProcessor::get_integrated_group_scalar_fluxes( const Hermes::vector< Solution<double>* >& solutions, 
                                                          Hermes::vector< double >* results,
                                                          unsigned int group, unsigned int G,
                                                          const Hermes::vector< std::string >& regions) const
  {
    Hermes::vector<MeshFunction<double>*> scalar_fluxes;
    
    if (method == NEUTRONICS_SPN)
      SPN::SupportClasses::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, G);
    else if (method == NEUTRONICS_SN)
      SN::SupportClasses::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, G, odata);
    else if (method == NEUTRONICS_DIFFUSION)
      for (Hermes::vector<Solution<double>*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
        scalar_fluxes.push_back(*it);
    
    Hermes::vector<std::string>::const_iterator region = regions.begin(); 
    for ( ; region != regions.end(); ++region)
      results->push_back(integrate(scalar_fluxes[group], *region));
    
    if (method == NEUTRONICS_SPN)
      SPN::SupportClasses::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
    else if (method == NEUTRONICS_SN)
      SN::SupportClasses::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
  }
  
  double PostProcessor::get_integrated_group_scalar_fluxes(const Hermes::vector< Solution<double>* >& solutions, 
                                                            unsigned int group, unsigned int G,
                                                            const Hermes::vector< std::string >& regions) const
  {
    Hermes::vector<MeshFunction<double>*> scalar_fluxes;
    
    if (method == NEUTRONICS_SPN)
      SPN::SupportClasses::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, G);
    else if (method == NEUTRONICS_SN)
      SN::SupportClasses::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, G, odata);
    else if (method == NEUTRONICS_DIFFUSION)
      for (Hermes::vector<Solution<double>*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
        scalar_fluxes.push_back(*it);
    
    double result = integrate(scalar_fluxes[group], regions);
    
    if (method == NEUTRONICS_SPN)
      SPN::SupportClasses::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
    else if (method == NEUTRONICS_SN)
      SN::SupportClasses::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
    
    return result;
  }
  
  void PostProcessor::get_integrated_scalar_fluxes(const Hermes::vector< Solution<double>* >& solutions, Hermes::vector< double >* results, 
                                                    unsigned int G, const Hermes::vector< std::string >& regions) const
  {
    Hermes::vector<MeshFunction<double>*> scalar_fluxes;
    
    if (method == NEUTRONICS_SPN)
      SPN::SupportClasses::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, G);
    else if (method == NEUTRONICS_SN)
      SN::SupportClasses::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, G, odata);
    else if (method == NEUTRONICS_DIFFUSION)
      for (Hermes::vector<Solution<double>*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
        scalar_fluxes.push_back(*it);
    
    Hermes::vector<std::string>::const_iterator region = regions.begin(); 
    for ( ; region != regions.end(); ++region)
    {
      double result = 0.0;
      for (unsigned int group = 0; group < G; group++)
        result += integrate(scalar_fluxes[group], *region);
      results->push_back(result);
    }
    
    if (method == NEUTRONICS_SPN)
      SPN::SupportClasses::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
    else if (method == NEUTRONICS_SN)
      SN::SupportClasses::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
  }

  double PostProcessor::get_integrated_scalar_fluxes(const Hermes::vector< Solution<double>* >& solutions, unsigned int G,
                                                      const Hermes::vector< std::string >& regions) const
  {
    double result = 0.0;
    for (unsigned int group = 0; group < G; group++)
      result += get_integrated_group_scalar_fluxes(solutions, group, G, regions);
    return result;
  }
  
  void PostProcessor::get_areas(Mesh *mesh, const Hermes::vector<std::string>& regions, Hermes::vector<double>* results) const
  {
    ConstantSolution<double> unity(mesh, 1.0);
    
    Hermes::vector<std::string>::const_iterator region = regions.begin(); 
    for ( ; region != regions.end(); ++region)
      results->push_back(integrate(&unity, *region));
  }

  double PostProcessor::get_area(Mesh* mesh, const Hermes::vector< std::string >& regions) const
  {
    ConstantSolution<double> unity(mesh, 1);
    return integrate(&unity, regions);
  }    

/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
} 
