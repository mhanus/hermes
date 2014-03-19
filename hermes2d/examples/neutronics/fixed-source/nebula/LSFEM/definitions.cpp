#include "definitions.h"
#include <cstdlib>
#include <iomanip>

#define VER1
//#define VER3
#define CUTOFF	1.e-8

double ExtinctionFunction::value(double x, double y) const
{
  double sx = x*x;
  double sy = y*y;
  
  double r;
  if (sx + sy < 1)
  	r = chi0 * std::exp(5. - 5./(1.-(sx+sy)));
  else
  	r = 0;

  return r;//std::max(CUTOFF, r);
}

Ord ExtinctionFunction::value(Ord x, Ord y) const
{
  return Ord::get_max_order();
}

IsotropicScatteringMatrixForm::IsotropicScatteringMatrixForm(double extinction, double thermalization) : WeakForm<double>(1)
{
  add_matrix_form(new IsotropicScatteringMatrixForm::ScatteringMF(ScatteringFunction(thermalization, ExtinctionFunction(extinction))));
}

template<typename Real>
Real IsotropicScatteringMatrixForm::ScatteringMF::matrix_form( int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                                               Func<Real> *v, Geom<Real> *e, Func<Real> **ext  ) const
{
  Real result(0.0);
  for (int i = 0; i < n; i++)
    result += wt[i] * Sigma_s.value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
  return result;
}

SNWeakForm::SNWeakForm(unsigned int N, double extinction, double thermalization,
                       const Hermes::vector<std::string>& reflective_boundaries, const Hermes::vector<std::string>& inflow_boundaries,
                       const char* out_tensor) 
  : WeakForm<double>(N*(N+2)/2), N(N), M(N*(N+2)/2), odata(SupportClasses::OrdinatesData(N, "../lgvalues.txt"))
{  
    std::cout << odata;
    
    bool assemble_all = strlen(out_tensor) == 0;    
    bool assemble_A = !strcmp(out_tensor, "A");
    bool assemble_Q = !strcmp(out_tensor, "Q");
    bool assemble_S = !strcmp(out_tensor, "S");

    if (!assemble_all && !strcmp(out_tensor, "AQ"))
    {
      assemble_A = true;
      assemble_Q = true;
    }
    
    ExtinctionFunction Sigma_t(extinction);
    ScatteringFunction Sigma_s(thermalization, Sigma_t);
    SourceFunction Q(thermalization, Sigma_t);
    
    for (unsigned int n = 0; n < M; n++)
    {
      if (assemble_all || assemble_A)
      {
        if (!reflective_boundaries.empty())
        {
          add_matrix_form_surf(new SpecularReflectionMF_X(odata, n, reflective_boundaries));
          add_matrix_form_surf(new SpecularReflectionMF_Y(odata, n, reflective_boundaries));
        }
        add_matrix_form_surf(new BoundaryStreamingMF(odata, n, 1.0));
        
        add_matrix_form(new VolumetricStreamingAndReactionsMF(odata, n, Sigma_t));
      }
      
      if (assemble_all || assemble_Q) 
      {
        add_vector_form(new VolumetricExternalSourceVF(odata, n, Q, Sigma_s, Sigma_t));
        if (!inflow_boundaries.empty())
        	add_vector_form_surf(new BoundaryStreamingVF(n, inflow_boundaries, 1.0));
      }
      
      if (assemble_all)
        add_vector_form(new VolumetricScatteringSourceVF(odata, n, Sigma_s, Sigma_t));
      
      if (assemble_S)
					for (unsigned int m = n; m < M; m++)
						add_matrix_form(new VolumetricScatteringMF(odata, m, n, Sigma_s, Sigma_t));

    }
}

double SNWeakForm::VolumetricStreamingAndReactionsMF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                  Geom<double> *e, Func<double> **ext) const
{
  double result = 0.0;

  for (int quad_pt = 0; quad_pt < n; quad_pt++)
  {
  	double St = Sigma_t.value(e->x[quad_pt], e->y[quad_pt]);
#if defined(VER1)
  	result += wt[quad_pt] * ( static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, u->dx[quad_pt], u->dy[quad_pt]) +
  	                        	St * u->val[quad_pt] )
													* ( static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, v->dx[quad_pt], v->dy[quad_pt]) +
															St * v->val[quad_pt] );

#else if defined(VER3)
  	if (fabs(St) > 1e-16)
  	{
			result += wt[quad_pt] * 1/St
															 * static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, u->dx[quad_pt], u->dy[quad_pt])
															 * static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, v->dx[quad_pt], v->dy[quad_pt]);
			result += wt[quad_pt] * St * u->val[quad_pt] * v->val[quad_pt];
  	}
  	else
  	{
  		result += wt[quad_pt] * static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, u->dx[quad_pt], u->dy[quad_pt])
  													* static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, v->dx[quad_pt], v->dy[quad_pt]);
  	}
#endif
  }

  //std::cout << "VolumetricStreamingAndReactionsMF :: " << result << std::endl;
  
  result = result * 1/(odata.pw[direction] * 2 * 4*M_PI);
  return result;
}

Ord SNWeakForm::VolumetricStreamingAndReactionsMF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                  Geom<Ord> *e, Func<Ord> **ext) const
{ 
	return Ord::get_max_order();
}
/*
template<typename Real>
Real SNWeakForm::VolumetricStreamingAndReactionsMF::b(Real x, Real y) const
{
  return Real(1);
}
*/

#if defined(VER1)

template<typename Real>
Real SNWeakForm::VolumetricScatteringSourceVF::vector_form(int n, double* wt, Func< Real >* u_ext[], Func< Real >* v, Geom< Real >* e, Func< Real >** ext) const
{
  Real result(0.0);
  Real *moment_values_at_quad_pts = new Real [n];
  SNWeakForm* _wf = static_cast<SNWeakForm*>(wf);

  odata.ordinates_to_moment(0, 0, 0, 1, u_ext, n, moment_values_at_quad_pts);
  SupportClasses::SphericalHarmonic Rlm(0, 0);
  double rlm_in_dir = Rlm(odata.xi[direction], odata.eta[direction], odata.mu[direction]);
    
  for (int quad_pt = 0; quad_pt < n; quad_pt++)
  {
    Real source = 1./(Sigma_t.value(e->x[quad_pt], e->y[quad_pt])) * ( _wf->calculate_a_dot_v(direction, v->dx[quad_pt], v->dy[quad_pt]) + Sigma_t.value(e->x[quad_pt], e->y[quad_pt]) * v->val[quad_pt] )
    								* Sigma_s.value(e->x[quad_pt], e->y[quad_pt]) * moment_values_at_quad_pts[quad_pt] * rlm_in_dir;
    
    // Contribution to spatial integral.
    if (geom_type == HERMES_AXISYM_X)
      result += source * wt[quad_pt] * e->y[quad_pt];
    else if (geom_type == HERMES_AXISYM_Y)
      result += source * wt[quad_pt]  * e->x[quad_pt];
    else
      result += source * wt[quad_pt];
  }
  
  result *= 1 / (4*M_PI);
    
  delete [] moment_values_at_quad_pts;
  
  //std::cout << "VolumetricScatteringSourceVF :: " << result << std::endl;
  
  return result;
}

template<typename Real>
Real SNWeakForm::VolumetricScatteringMF::matrix_form(int n, double* wt, Func< Real >* u_ext[], Func< Real >* u, Func< Real >* v, Geom< Real >* e, Func< Real >** ext) const
{
  Real result(0.0);
  //Real *moment_values_at_quad_pts = new Real [n];
  SNWeakForm* _wf = static_cast<SNWeakForm*>(wf);

  //odata.ordinate_to_moment(direction_from, 0, 0, 0, 1, u, n, moment_values_at_quad_pts);

  //SupportClasses::SphericalHarmonic Rlm(0, 0);
  //double rlm_in_dir = Rlm(odata.xi[direction], odata.eta[direction], odata.mu[direction]);

  for (int quad_pt = 0; quad_pt < n; quad_pt++)
  {
  	Real St = Sigma_t.value(e->x[quad_pt], e->y[quad_pt]);
  	Real Ss = Sigma_s.value(e->x[quad_pt], e->y[quad_pt]);

    Real tmp = 		( _wf->calculate_a_dot_v(direction, v->dx[quad_pt], v->dy[quad_pt]) + St * v->val[quad_pt] )
        						* Ss * u->val[quad_pt]
    						+ ( _wf->calculate_a_dot_v(direction_from, u->dx[quad_pt], u->dy[quad_pt]) + St * u->val[quad_pt] )
    							  * Ss * v->val[quad_pt]
                - Ss * u->val[quad_pt] * Ss * v->val[quad_pt];

    // Contribution to spatial integral.
    if (geom_type == HERMES_AXISYM_X)
      result += tmp * wt[quad_pt] * e->y[quad_pt];
    else if (geom_type == HERMES_AXISYM_Y)
      result += tmp * wt[quad_pt] * e->x[quad_pt];
    else
      result += tmp * wt[quad_pt];
  }

  result *= 1 / (4*M_PI);

  //delete [] moment_values_at_quad_pts;

  //std::cout << "VolumetricScatteringSourceVF :: " << result << std::endl;

  return result;
}

template<typename Real>
Real SNWeakForm::VolumetricExternalSourceVF::vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v,
                                                  Geom<Real> *e, Func<Real> **ext) const
{
  Real result = Real(0.0);
  SNWeakForm* _wf = static_cast<SNWeakForm*>(wf);
  for (int i = 0; i < n; i++)
    result += wt[i] * Q.value(e->x[i], e->y[i])
    								* ( _wf->calculate_a_dot_v(direction, v->dx[i], v->dy[i])
    										+ Sigma_t.value(e->x[i], e->y[i]) * v->val[i]
    										- Sigma_s.value(e->x[i], e->y[i]) * v->val[i] );
  
  //std::cout << "VolumetricExternalSourceVF :: " << result << std::endl;
  
  return result;
}

double SNWeakForm::BoundaryStreamingMF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                      Geom<double> *e, Func<double> **ext) const
{
  double result = 0.0;
  for (int i = 0; i < n; i++)
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[i], e->ny[i]);

    if (a_dot_n < 0)
      result += wt[i] * u->val[i] * (-a_dot_n) * v->val[i];
  }

  //std::cout << "BoundaryStreamingMF :: " << result << std::endl;

  return alpha * result * 1/(odata.pw[direction] * 2 * 4*M_PI);
}

Ord SNWeakForm::BoundaryStreamingMF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                      Geom<Ord> *e, Func<Ord> **ext) const
{
  return u->val[0]*v->val[0];
}

double SNWeakForm::BoundaryStreamingVF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0;

  for (int quad_pt = 0; quad_pt < n; quad_pt++)
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);

    if (a_dot_n < 0)
    {
      double boundary_data = influx<double>(e->x[quad_pt],e->y[quad_pt]);
      result += wt[quad_pt] * boundary_data * (-a_dot_n) * v->val[quad_pt];
    }
  }
  return alpha*result;
}

Ord SNWeakForm::BoundaryStreamingVF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  //std::cout << "BoundaryStreamingVF :: " << v->val[0] * influx<Ord>(e->x[0], e->y[0]) << std::endl;
  return v->val[0] * influx<Ord>(e->x[0], e->y[0]);
}

#else if defined(VER3)

template<typename Real>
Real SNWeakForm::VolumetricScatteringMF::matrix_form(int n, double* wt, Func< Real >* u_ext[], Func< Real >* u, Func< Real >* v, Geom< Real >* e, Func< Real >** ext) const
{
  Real result(0.0);
  //Real *trialfn_moment_values_at_quad_pts = new Real [n];
  //Real *testfn_moment_values_at_quad_pts = new Real [n];
  SNWeakForm* _wf = static_cast<SNWeakForm*>(wf);

        	//odata.ordinate_to_moment(direction_from, 0, 0, 0, 1, u, n, trialfn_moment_values_at_quad_pts);
        	//odata.ordinate_to_moment(direction, 0, 0, 0, 1, v, n, testfn_moment_values_at_quad_pts);


          //SupportClasses::SphericalHarmonic Rlm(0, 0);
          //double rlm_in_dir = Rlm(odata.xi[direction], odata.eta[direction], odata.mu[direction]);
          //double rlm_in_dir_from = Rlm(odata.xi[direction_from], odata.eta[direction_from], odata.mu[direction_from]);

          Real C(1/(4*M_PI));
          //Real C(1.0);

          //std::cout << std::setprecision(15);

          for (int quad_pt = 0; quad_pt < n; quad_pt++)
          {
          	Real Ss = Sigma_s.value(e->x[quad_pt], e->y[quad_pt]);
          	Real St = max(CUTOFF, Sigma_t.value(e->x[quad_pt], e->y[quad_pt]));
          	Real D = Ss / (St*St - St*Ss);

        	  //Real trialfn_group_source_moment = trialfn_moment_values_at_quad_pts[quad_pt] * rlm_in_dir_from * rlm_in_dir;
        	  //Real testfn_group_source_moment = testfn_moment_values_at_quad_pts[quad_pt] * rlm_in_dir_from * rlm_in_dir ;

          	/*if (St > 1E-2) {
          	std::cout << "Ss = " << Ss << std::endl << std::flush;;
          	std::cout << "St = " << St << std::endl << std::flush;;
          	std::cout << "St-Ss = " << St-Ss << std::endl << std::flush;;
          	std::cout << "D = " <<  D << std::endl << std::flush;
          	std::cout << "________________________________________________"  << std::endl << std::endl << std::flush;
          	}*/
        	  Real tmp = -/* odata.pw[direction_from] * 2 * */ C * D
        	          	                                          * _wf->calculate_a_dot_v(direction_from, u->dx[quad_pt], u->dy[quad_pt])
        	          	                                          * _wf->calculate_a_dot_v(direction, v->dx[quad_pt], v->dy[quad_pt]);
						tmp += C * Ss * u->val[quad_pt] * v->val[quad_pt];


            // Contribution to spatial integral.
            if (geom_type == HERMES_AXISYM_X)
              result += tmp * wt[quad_pt] * e->y[quad_pt];
            else if (geom_type == HERMES_AXISYM_Y)
              result += tmp * wt[quad_pt] * e->x[quad_pt];
            else
              result += tmp * wt[quad_pt];
          }

    //delete [] trialfn_moment_values_at_quad_pts;
    //delete [] testfn_moment_values_at_quad_pts;


  return result;
}


template<typename Real>
Real SNWeakForm::VolumetricScatteringSourceVF::vector_form(int n, double* wt, Func< Real >* u_ext[], Func< Real >* v, Geom< Real >* e, Func< Real >** ext) const
{
  Real result(0.0);
  Real *trialfn_moment_values_at_quad_pts = new Real [n];
  Real *testfn_moment_values_at_quad_pts = new Real [n];
  SNWeakForm* _wf = static_cast<SNWeakForm*>(wf);

        	odata.ordinates_to_moment(0, 0, 0, 1, u_ext, n, trialfn_moment_values_at_quad_pts);
        	odata.ordinate_to_moment(direction, 0, 0, 0, 1, v, n, testfn_moment_values_at_quad_pts);


          SupportClasses::SphericalHarmonic Rlm(0, 0);
          double rlm_in_dir = Rlm(odata.xi[direction], odata.eta[direction], odata.mu[direction]);

          Real C(1/(4*M_PI));

          for (int quad_pt = 0; quad_pt < n; quad_pt++)
          {
          	Real Ss = Sigma_s.value(e->x[quad_pt], e->y[quad_pt]);
          	Real St = Sigma_t.value(e->x[quad_pt], e->y[quad_pt]);

          	Func<Real>* u = u_ext[direction];

        	  Real trialfn_group_source_moment = trialfn_moment_values_at_quad_pts[quad_pt] * rlm_in_dir;
        	  Real testfn_group_source_moment = testfn_moment_values_at_quad_pts[quad_pt] * rlm_in_dir;

        	  Real tmp = trialfn_group_source_moment * C * Ss
        	                                          * ( _wf->calculate_a_dot_v(direction, v->dx[quad_pt], v->dy[quad_pt]) +
        	                                          		St * v->val[quad_pt] );
								tmp = testfn_group_source_moment * C * Ss
																* ( _wf->calculate_a_dot_v(direction, u->dx[quad_pt], u->dy[quad_pt]) +
																	St * u->val[quad_pt] );
								tmp -= testfn_group_source_moment * C * Ss * trialfn_group_source_moment * C * Ss;

						tmp = tmp * 1./St;

            // Contribution to spatial integral.
            if (geom_type == HERMES_AXISYM_X)
              result += tmp * wt[quad_pt] * e->y[quad_pt];
            else if (geom_type == HERMES_AXISYM_Y)
              result += tmp * wt[quad_pt] * e->x[quad_pt];
            else
              result += tmp * wt[quad_pt];
          }

    delete [] trialfn_moment_values_at_quad_pts;
    delete [] testfn_moment_values_at_quad_pts;


  return /*1/Sigma_t**/result;
}


template<typename Real>
Real SNWeakForm::VolumetricExternalSourceVF::vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v,
                                                  Geom<Real> *e, Func<Real> **ext) const
{
  Real result = Real(0.0);
  //Real *testfn_moment_values_at_quad_pts = new Real [n];
  SNWeakForm* _wf = static_cast<SNWeakForm*>(wf);


          //odata.ordinate_to_moment(direction, 0, 0, 0, 1, v, n, testfn_moment_values_at_quad_pts);

          //SupportClasses::SphericalHarmonic Rlm(0, 0);
          //double rlm_in_dir = Rlm(odata.xi[direction], odata.eta[direction], odata.mu[direction]);
          Real C(1./(4*M_PI));

          for (int quad_pt = 0; quad_pt < n; quad_pt++)
          {
            //Real testfn_group_source_moment = testfn_moment_values_at_quad_pts[quad_pt] * rlm_in_dir;

          	Real Ss = Sigma_s.value(e->x[quad_pt], e->y[quad_pt]);
          	Real St = max(CUTOFF, Sigma_t.value(e->x[quad_pt], e->y[quad_pt]));
          	Real D = Ss / (St*St - St*Ss);
          	Real Qv = Q.value(e->x[quad_pt], e->y[quad_pt]);

          	Real tmp = Qv * (1/St + D
          	             * _wf->calculate_a_dot_v(direction, v->dx[quad_pt], v->dy[quad_pt])) + Qv * v->val[quad_pt];

            // Contribution to spatial integral.
            if (geom_type == HERMES_AXISYM_X)
              result += tmp * wt[quad_pt] * e->y[quad_pt];
            else if (geom_type == HERMES_AXISYM_Y)
              result += tmp * wt[quad_pt] * e->x[quad_pt];
            else
              result += tmp * wt[quad_pt];

          }

  //delete [] testfn_moment_values_at_quad_pts;

  //std::cout << "VolumetricExternalSourceVF :: " << result << std::endl;

  return result;
}

double SNWeakForm::BoundaryStreamingMF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                      Geom<double> *e, Func<double> **ext) const
{
  double result = 0.0;
  for (int i = 0; i < n; i++)
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[i], e->ny[i]);

    //if (a_dot_n > 0)
      result += wt[i] * u->val[i] * Hermes::abs(a_dot_n) * v->val[i];
  }

  //std::cout << "BoundaryStreamingMF :: " << result << std::endl;

  result = result * 1/(odata.pw[direction] * 2 * 4*M_PI);
  return result;
}

Ord SNWeakForm::BoundaryStreamingMF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                      Geom<Ord> *e, Func<Ord> **ext) const
{
  return u->val[0]*v->val[0];
}

double SNWeakForm::BoundaryStreamingVF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++)
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);

    if (a_dot_n < 0)
    {
      double boundary_data = influx<double>(e->x[quad_pt],e->y[quad_pt]);
      result += wt[quad_pt] * boundary_data * (-a_dot_n) * v->val[quad_pt];
    }
  }
  return 2*result;
}

Ord SNWeakForm::BoundaryStreamingVF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  //std::cout << "BoundaryStreamingVF :: " << v->val[0] * influx<Ord>(e->x[0], e->y[0]) << std::endl;
  return v->val[0] * influx<Ord>(e->x[0], e->y[0]);
}

#endif


double SNWeakForm::SpecularReflectionMF_X::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                      Geom<double> *e, Func<double> **ext) const
{
  double eps = 1e-14;
  double result = 0.0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++)
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);
  
    if (a_dot_n < 0)
      if (fabs(e->ny[quad_pt] - 1.0) < eps || fabs(e->ny[quad_pt] + 1.0) < eps && fabs(e->nx[quad_pt]) < eps)
        result += -wt[quad_pt] * u->val[quad_pt] * (-a_dot_n) * v->val[quad_pt];
  }
  
  result = result * 1/(odata.pw[direction] * 2 * 4*M_PI);
  return result;
}

Ord SNWeakForm::SpecularReflectionMF_X::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                      Geom<Ord> *e, Func<Ord> **ext) const
{ 
  return u->val[0]*v->val[0];
}

double SNWeakForm::SpecularReflectionMF_Y::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                      Geom<double> *e, Func<double> **ext) const
{
  double eps = 1e-14;
  double result = 0.0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++)
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);
    
    if (a_dot_n < 0)
      if (fabs(e->nx[quad_pt] - 1.0) < eps || fabs(e->nx[quad_pt] + 1.0) < eps && fabs(e->ny[quad_pt]) < eps)
        result += -wt[quad_pt] * u->val[quad_pt] * (-a_dot_n) * v->val[quad_pt];
  }
  
  result = result * 1/(odata.pw[direction] * 2 * 4*M_PI);
  return result;
}

Ord SNWeakForm::SpecularReflectionMF_Y::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                      Geom<Ord> *e, Func<Ord> **ext) const
{ 
  return u->val[0]*v->val[0];
}




double SNWeakForm::calculate_a_dot_v(int n, double vx, double vy) const
{
  //return std::cos((n+1)*M_PI/(2.1*N))*vx + std::sin((n+1)*M_PI/(2.1*N))*vy;
  return odata.xi[n]*vx + odata.eta[n]*vy;
}

Ord SNWeakForm::calculate_a_dot_v(int n, Ord vx, Ord vy) const
{
  //return std::cos((n+1)*M_PI/(2.1*N))*vx + std::sin((n+1)*M_PI/(2.1*N))*vy;
  return vx + vy;
}

