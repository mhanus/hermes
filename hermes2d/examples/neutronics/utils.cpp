#define HERMES_REPORT_ALL
#include "utils.h"

#include <iterator>

// Utility functions that simplify repeated reporting of number of DOF during adaptivity.
void report_num_dof(const std::string& msg, const Hermes::vector< SpaceSharedPtr<double> > spaces)
{
  std::stringstream ss;
  
  ss << msg << spaces[0]->get_num_dofs();
  
  for (unsigned int i = 1; i < spaces.size(); i++)
    ss << " + " << spaces[i]->get_num_dofs();
  
  if (spaces.size() > 1)
    ss << " = " << Space<double>::get_num_dofs(spaces);
  
  Loggable::Static::info(ss.str().c_str());
}

void report_errors(const std::string& msg, const Hermes::Hermes2D::ErrorCalculator< double >& error_calculator)
{
  std::stringstream ss;
  ss << msg;
  
  for (unsigned int i = 0; i < error_calculator.get_component_count()-1; i++)
    ss << sqrt(error_calculator.get_error_squared(i))*100 << "%%, ";
  
  ss << sqrt(error_calculator.get_error_squared(error_calculator.get_component_count()-1))*100 << "%%";
  
  Loggable::Static::info(ss.str().c_str());
}

void save_algebraic_representation(Hermes::Hermes2D::WeakForm< double >* wf, const Hermes::vector< SpaceSharedPtr<double> >& spaces, const std::string& varname, bool assign_dofs)
{
  DiscreteProblem<double> dp(wf, spaces);
  dp.set_linear();
  dp.set_do_not_use_cache();
  
  if (assign_dofs)
    Space<double>::assign_dofs(spaces);
  
  bool has_matrices = (!wf->get_mfvol().empty() || !wf->get_mfsurf().empty() || !wf->get_mfDG().empty());
  bool has_vectors = (!wf->get_vfvol().empty() || !wf->get_vfsurf().empty() || !wf->get_vfDG().empty());
  
  SparseMatrix<double>* matrix = NULL;
  Vector<double>* vector = NULL;
  
  Mixins::MatrixRhsOutput<double> dumper;
  
  if (has_matrices)
  {
    matrix = create_matrix<double>();
    
    dumper.output_matrix();
    dumper.set_matrix_E_matrix_dump_format(DF_HERMES_BIN);
    dumper.set_matrix_number_format("%1.15f");
    dumper.set_matrix_filename(varname+".dat");
    dumper.set_matrix_varname(varname);
  }
  if (has_vectors)
  {
    vector = create_vector<double>();
    
    dumper.output_rhs();
    dumper.set_rhs_E_matrix_dump_format(DF_HERMES_BIN);;
    dumper.set_rhs_number_format("%1.15f");
    dumper.set_rhs_filename(varname+".dat");
    dumper.set_rhs_varname(varname);
  }
  
  dp.assemble(matrix, vector);
  dumper.process_matrix_output(matrix);
  dumper.process_vector_output(vector);
  
  if (matrix)
    delete matrix;
  if (vector)
    delete vector;
}

void load_solution(const std::string& sln_file, 
                   const Hermes::vector<SpaceSharedPtr<double> >& spaces, 
                   const Neutronics::Common::MaterialProperties::MaterialPropertyMaps* matprop,
                   VisualizationOptions visualization, bool mode_3D)
{  
  int ndof =  Space<double>::get_num_dofs(spaces);
  int G = matprop->get_G();
  Neutronics::NeutronicsMethod method;
  
  if (dynamic_cast<const Neutronics::Diffusion::MaterialProperties::MaterialPropertyMaps*>(matprop))
    method = Neutronics::NEUTRONICS_DIFFUSION;
  else if (dynamic_cast<const Neutronics::SPN::MaterialProperties::MaterialPropertyMaps*>(matprop))
    method = Neutronics::NEUTRONICS_SPN;
  else if (dynamic_cast<const Neutronics::SN::MaterialProperties::MaterialPropertyMaps*>(matprop))
    method = Neutronics::NEUTRONICS_SN;
  
  std::vector<double> x_ext;
  x_ext.reserve(ndof);
  std::ifstream ifs( sln_file.c_str() , std::ifstream::in );
  read_solution_from_file(ifs, std::back_inserter(x_ext));
  ifs.close();
  
  Hermes::vector<MeshFunctionSharedPtr<double> > sol_ext;
  for (unsigned int i = 0; i < spaces.size(); i++) 
    sol_ext.push_back(new Solution<double>()); 
  Solution<double>::vector_to_solutions(x_ext.data(), spaces, sol_ext);
  
  Neutronics::Common::SupportClasses::Visualization *vis;
  
  switch (method)
  {
    case Neutronics::NEUTRONICS_DIFFUSION:
    {
      vis = new Neutronics::Diffusion::SupportClasses::Visualization(G);
      
      if (visualization & HERMES_SCALAR_VISUALIZATION)
        vis->show_scalar_fluxes(sol_ext);
      
      break;
    }
    case Neutronics::NEUTRONICS_SPN:
    {
      int N = 2*spaces.size()/G - 1; 
      Neutronics::SPN::SupportClasses::Visualization* spn_vis = new Neutronics::SPN::SupportClasses::Visualization(N, G);
      vis = spn_vis;
      
      if (visualization & HERMES_SCALAR_VISUALIZATION)
        spn_vis->show_scalar_fluxes(sol_ext);
      
      if (visualization & HERMES_ANGULAR_VISUALIZATION)
        spn_vis->show_all_flux_moments(sol_ext, *dynamic_cast<const Neutronics::SPN::MaterialProperties::MaterialPropertyMaps*>(matprop));
      
      break;
    }
    case Neutronics::NEUTRONICS_SN:
    {
      int M = spaces.size()/G;
      int N = -1 + sqrt(2*M+1);
      vis = new Neutronics::SN::SupportClasses::Visualization(G, N, Neutronics::SN::SupportClasses::OrdinatesData(N, "lgvalues.txt"));
      
      if (visualization & HERMES_SCALAR_VISUALIZATION)
        vis->show_scalar_fluxes(sol_ext);
      
      if (visualization & HERMES_ANGULAR_VISUALIZATION)
        vis->show_solutions(sol_ext);
      
      break;
    }  
    default:
      return;
  }
  
  if(visualization & (VTK_ANGULAR_VISUALIZATION))
    vis->save_solutions_vtk("psi", "psi", sol_ext, mode_3D);
  if(visualization & (VTK_SCALAR_VISUALIZATION))  
    vis->save_scalar_fluxes_vtk("phi", "phi", sol_ext, mode_3D);
  
  if (visualization & (HERMES_SCALAR_VISUALIZATION | HERMES_ANGULAR_VISUALIZATION))
    Views::View::wait();
  
  delete vis;
}

