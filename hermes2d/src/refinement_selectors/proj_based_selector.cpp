// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "proj_based_selector.h"
#include "hcurl_proj_based_selector.h"
#include <algorithm>
#include "order_permutator.h"
#include "algebra/dense_matrix_operations.h"

using namespace Hermes::Algebra::DenseMatrixOperations;

namespace Hermes
{
  namespace Hermes2D
  {
    namespace RefinementSelectors
    {
      template<typename Scalar>
      ProjBasedSelector<Scalar>::ProjBasedSelector(CandList cand_list, int
        max_order, Shapeset* shapeset, const Range& vertex_order, const
        Range& edge_bubble_order) :
        OptimumSelector<Scalar>(cand_list, max_order, shapeset, vertex_order, edge_bubble_order),
        warn_uniform_orders(false),
        error_weight_h(H2DRS_DEFAULT_ERR_WEIGHT_H),
        error_weight_p(H2DRS_DEFAULT_ERR_WEIGHT_P),
        error_weight_aniso(H2DRS_DEFAULT_ERR_WEIGHT_ANISO)
      {
          cached_shape_vals_valid = new bool[2];
          cached_shape_ortho_vals = new TrfShape[2];
          cached_shape_vals = new TrfShape[2];

          //clean svals initialization state
          std::fill(cached_shape_vals_valid, cached_shape_vals_valid + H2D_NUM_MODES, false);

          //clear matrix cache
          for (int m = 0; m < H2D_NUM_MODES; m++)
          for (int i = 0; i < H2DRS_MAX_ORDER + 2; i++)
          for (int k = 0; k < H2DRS_MAX_ORDER + 2; k++)
            proj_matrix_cache[m][i][k] = nullptr;
        }

      template<typename Scalar>
      ProjBasedSelector<Scalar>::~ProjBasedSelector()
      {
        //delete matrix cache
        for (int m = 0; m < H2D_NUM_MODES; m++)
        {
          for (int i = 0; i < H2DRS_MAX_ORDER + 2; i++)
          for (int k = 0; k < H2DRS_MAX_ORDER + 2; k++)
          {
            if (proj_matrix_cache[m][i][k] != nullptr)
              delete[] proj_matrix_cache[m][i][k];
          }
        }

        delete[] cached_shape_vals_valid;
        delete[] cached_shape_ortho_vals;
        delete[] cached_shape_vals;
      }

      template<typename Scalar>
      void ProjBasedSelector<Scalar>::set_error_weights(double weight_h, double weight_p, double weight_aniso)
      {
        error_weight_h = weight_h;
        error_weight_p = weight_p;
        error_weight_aniso = weight_aniso;
      }

      template<typename Scalar>
      double ProjBasedSelector<Scalar>::get_error_weight_h() const
      {
        return error_weight_h;
      }

      template<typename Scalar>
      double ProjBasedSelector<Scalar>::get_error_weight_p() const
      {
        return error_weight_p;
      }

      template<typename Scalar>
      double ProjBasedSelector<Scalar>::get_error_weight_aniso() const
      {
        return error_weight_aniso;
      }

      template<typename Scalar>
      ProjBasedSelector<Scalar>::TrfShapeExp::TrfShapeExp() : num_gip(0), num_expansion(0), values(nullptr) {};

      template<typename Scalar>
      ProjBasedSelector<Scalar>::TrfShapeExp::~TrfShapeExp()
      {
        free_with_check(values, true);
      }

      template<typename Scalar>
      void ProjBasedSelector<Scalar>::TrfShapeExp::allocate(int num_expansion, int num_gip)
      {
        free_with_check(values, true);
        values = new_matrix<double>(num_expansion, num_gip);
        this->num_expansion = num_expansion;
        this->num_gip = num_gip;
      }

      template<typename Scalar>
      double* ProjBasedSelector<Scalar>::TrfShapeExp::operator[](int inx_expansion)
      {
        return values[inx_expansion];
      }

      template<typename Scalar>
      bool ProjBasedSelector<Scalar>::TrfShapeExp::empty() const
      {
        return values == nullptr;
      }

      template<typename Scalar>
      void ProjBasedSelector<Scalar>::evaluate_cands_error(Hermes::vector<Cand>& candidates, Element* e, MeshFunction<Scalar>* rsln)
      {
        bool tri = e->is_triangle();

        // find range of orders
        typename OptimumSelector<Scalar>::CandsInfo info_h, info_p, info_aniso;
        this->update_cands_info(candidates, info_h, info_p, info_aniso);

        // calculate squared projection errors of elements of candidates
        CandElemProjError herr[4], anisoerr[4], perr;
        calc_projection_errors(e, info_h, info_p, info_aniso, rsln, herr, perr, anisoerr);

        //evaluate errors and dofs
        for (unsigned i = 0; i < candidates.size(); i++)
        {
          Cand& c = candidates[i];
          double error_squared = 0.0;
          if (tri)
          {
            switch (c.split)
            {
            case H2D_REFINEMENT_H:
              error_squared = 0.0;
              for (int j = 0; j < H2D_MAX_ELEMENT_SONS; j++)
              {
                int order = H2D_GET_H_ORDER(c.p[j]);
                c.errors[j] = herr[j][order][order];
                error_squared += c.errors[j];
              }
              error_squared *= 0.25; //element of a candidate occupies 1/4 of the reference domain defined over a candidate
              break;

            case H2D_REFINEMENT_P:
            {
                                   int order = H2D_GET_H_ORDER(c.p[0]);
                                   c.errors[0] = perr[order][order];
                                   error_squared = perr[order][order];
            }
              break;

            default:
              throw Hermes::Exceptions::Exception("Unknown split type \"%d\" at candidate %d", c.split, i);
            }
          }
          else
          {
            //quad
            switch (c.split)
            {
            case H2D_REFINEMENT_H:
              error_squared = 0.0;
              for (int j = 0; j < H2D_MAX_ELEMENT_SONS; j++)
              {
                int order_h = H2D_GET_H_ORDER(c.p[j]), order_v = H2D_GET_V_ORDER(c.p[j]);
                c.errors[j] = herr[j][order_h][order_v];
                error_squared += c.errors[j];
              }
              error_squared *= 0.25; //element of a candidate occupies 1/4 of the reference domain defined over a candidate
              break;

            case H2D_REFINEMENT_ANISO_H:
            case H2D_REFINEMENT_ANISO_V:
            {
                                         error_squared = 0.0;
                                         for (int j = 0; j < 2; j++)
                                         {
                                           c.errors[j] = anisoerr[(c.split == H2D_REFINEMENT_ANISO_H) ? j : j + 2][H2D_GET_H_ORDER(c.p[j])][H2D_GET_V_ORDER(c.p[j])];
                                           error_squared += c.errors[j];
                                         }
                                         error_squared *= 0.5;  //element of a candidate occupies 1/2 of the reference domain defined over a candidate
            }
              break;

            case H2D_REFINEMENT_P:
            {
                                   int order_h = H2D_GET_H_ORDER(c.p[0]), order_v = H2D_GET_V_ORDER(c.p[0]);
                                   c.errors[0] = perr[order_h][order_v];
                                   error_squared = c.errors[0];
            }
              break;

            default:
              throw Hermes::Exceptions::Exception("Unknown split type \"%d\" at candidate %d", c.split, i);
            }
          }

          //calculate error from squared error
          c.error = sqrt(error_squared);

          //apply weights
          switch (c.split)
          {
          case H2D_REFINEMENT_H: c.error *= error_weight_h; break;
          case H2D_REFINEMENT_ANISO_H:
          case H2D_REFINEMENT_ANISO_V: c.error *= error_weight_aniso; break;
          case H2D_REFINEMENT_P: c.error *= error_weight_p; break;
          default: throw Hermes::Exceptions::Exception("Unknown split type \"%d\" at candidate %d", c.split, i);
          }
        }
      }

      template<typename Scalar>
      void ProjBasedSelector<Scalar>::calc_projection_errors(Element* e, const typename OptimumSelector<Scalar>::CandsInfo& info_h, const typename OptimumSelector<Scalar>::CandsInfo& info_p, const  typename OptimumSelector<Scalar>::CandsInfo& info_aniso, MeshFunction<Scalar>* rsln, CandElemProjError herr[H2D_MAX_ELEMENT_SONS], CandElemProjError perr, CandElemProjError anisoerr[H2D_MAX_ELEMENT_SONS])
      {
        ElementMode2D mode = e->get_mode();

        // select quadrature, obtain integration points and weights
        Quad2D* quad = &g_quad_2d_std;
        rsln->set_quad_2d(quad);
        double3* gip_points = quad->get_points(H2DRS_INTR_GIP_ORDER, mode);
        int num_gip_points = quad->get_num_points(H2DRS_INTR_GIP_ORDER, mode);

        // everything is done on the reference domain
        Scalar* rval[H2D_MAX_ELEMENT_SONS][MAX_NUMBER_FUNCTION_VALUES_FOR_SELECTORS];

        Solution<Scalar>* rsln_sln = dynamic_cast<Solution<Scalar>*>(rsln);
        if (rsln_sln != nullptr)
          rsln_sln->enable_transform(false);

        // obtain reference solution values on all four refined sons
        Element* base_element = rsln->get_mesh()->get_element(e->id);

        // value on base element.
        if (base_element->active)
        {
          for (int son = 0; son < H2D_MAX_ELEMENT_SONS; son++)
          {
            rsln->set_active_element(base_element);
            rsln->push_transform(son);
            rsln->set_quad_order(H2DRS_INTR_GIP_ORDER);

            //obtain precalculated values
            precalc_ref_solution(son, rsln, e, H2DRS_INTR_GIP_ORDER, rval);
          }
        }
        else
        {
          for (int son = 0; son < H2D_MAX_ELEMENT_SONS; son++)
          {
            rsln->set_active_element(base_element->sons[son]);
            rsln->set_quad_order(H2DRS_INTR_GIP_ORDER);

            //obtain precalculated values
            precalc_ref_solution(son, rsln, e, H2DRS_INTR_GIP_ORDER, rval);
          }
        }

        //retrieve transformations
        Trf* trfs = nullptr;
        int num_noni_trfs = 0;
        if (mode == HERMES_MODE_TRIANGLE)
        {
          trfs = tri_trf;
          num_noni_trfs = H2D_TRF_TRI_NUM;
        }
        else
        {
          trfs = quad_trf;
          num_noni_trfs = H2D_TRF_QUAD_NUM;
        }

        // precalculate values of shape functions
        TrfShape empty_shape_vals;

        if (!cached_shape_vals_valid[mode])
        {
#pragma omp critical (cached_shape_vals_valid)
          if (!cached_shape_vals_valid[mode])
          {
            precalc_ortho_shapes(gip_points, num_gip_points, trfs, num_noni_trfs, this->shape_indices[mode], this->max_shape_inx[mode], cached_shape_ortho_vals[mode], mode);
            precalc_shapes(gip_points, num_gip_points, trfs, num_noni_trfs, this->shape_indices[mode], this->max_shape_inx[mode], cached_shape_vals[mode], mode);
            cached_shape_vals_valid[mode] = true;
          }
        }

        //issue a warning if ortho values are defined and the selected cand_list might benefit from that but it cannot because elements do not have uniform orders
        if (!warn_uniform_orders && mode == HERMES_MODE_QUAD && !cached_shape_ortho_vals[mode][H2D_TRF_IDENTITY].empty())
        {
          warn_uniform_orders = true;
          if (this->cand_list == H2D_H_ISO || this->cand_list == H2D_H_ANISO || this->cand_list == H2D_P_ISO || this->cand_list == H2D_HP_ISO || this->cand_list == H2D_HP_ANISO_H)
          {
            this->warn_if(!info_h.uniform_orders || !info_aniso.uniform_orders || !info_p.uniform_orders, "Possible inefficiency: %s might be more efficient if the input mesh contains elements with uniform orders strictly.", get_cand_list_str(this->cand_list));
          }
        }

        TrfShape& svals = cached_shape_vals[mode];
        TrfShape& ortho_svals = cached_shape_ortho_vals[mode];

#pragma region candidatesEvaluation
        //H-candidates
        if (!info_h.is_empty())
        {
          if (base_element->active)
          {
            Trf* sub_trfs[4] = { &trfs[0], &trfs[1], &trfs[2], &trfs[3] };
            Hermes::vector<TrfShapeExp>* p_trf_svals[4] = { &svals[0], &svals[1], &svals[2], &svals[3] };
            Hermes::vector<TrfShapeExp>* p_trf_ortho_svals[4] = { &ortho_svals[0], &ortho_svals[1], &ortho_svals[2], &ortho_svals[3] };
            for (int son = 0; son < H2D_MAX_ELEMENT_SONS; son++)
            {
              int sub_rval[1] = { son };
              calc_error_cand_element(mode, gip_points, num_gip_points
                , 1, &base_element, &sub_trfs[son], sub_rval
                , &p_trf_svals[son], &p_trf_ortho_svals[son]
                , info_h, herr[son], rval);
            }
          }
          else
          {
            Trf* p_trf_identity[1] = { &trfs[H2D_TRF_IDENTITY] };
            Hermes::vector<TrfShapeExp>* p_trf_svals[1] = { &svals[H2D_TRF_IDENTITY] };
            Hermes::vector<TrfShapeExp>* p_trf_ortho_svals[1] = { &ortho_svals[H2D_TRF_IDENTITY] };
            for (int son = 0; son < H2D_MAX_ELEMENT_SONS; son++)
            {
              int sub_rval[1] = { son };
              calc_error_cand_element(mode, gip_points, num_gip_points
                , 1, &base_element->sons[son], p_trf_identity, sub_rval
                , p_trf_svals, p_trf_ortho_svals
                , info_h, herr[son], rval);
            }
          }
        }

        //ANISO-candidates
        if (!info_aniso.is_empty())
        {
          if (base_element->active)
          {
            const int tr[4][2] = { { 0, 1 }, { 3, 2 }, { 0, 3 }, { 1, 2 } };
            for (int version = 0; version < 4; version++)
            {
              Trf* sub_trfs[2] = { &trfs[tr[version][0]], &trfs[tr[version][1]] };
              Element* sub_domains[2] = { base_element, base_element };
              int sub_rval[2] = { tr[version][0], tr[version][1] };
              Hermes::vector<TrfShapeExp>* sub_svals[2] = { &svals[tr[version][0]], &svals[tr[version][1]] };
              Hermes::vector<TrfShapeExp>* sub_ortho_svals[2] = { &ortho_svals[tr[version][0]], &ortho_svals[tr[version][1]] };
              calc_error_cand_element(mode, gip_points, num_gip_points
                , 2, sub_domains, sub_trfs, sub_rval
                , sub_svals, sub_ortho_svals
                , info_aniso, anisoerr[version], rval);
            }
          }
          else
          {
            const int sons[4][2] = { { 0, 1 }, { 3, 2 }, { 0, 3 }, { 1, 2 } }; //indices of sons for sub-areas
            const int tr[4][2] = { { 6, 7 }, { 6, 7 }, { 4, 5 }, { 4, 5 } }; //indices of ref. domain transformations for sub-areas
            for (int version = 0; version < 4; version++)
            { // 2 elements for vertical split, 2 elements for horizontal split
              Trf* sub_trfs[2] = { &trfs[tr[version][0]], &trfs[tr[version][1]] };
              Element* sub_domains[2] = { base_element->sons[sons[version][0]], base_element->sons[sons[version][1]] };
              int sub_rval[2] = { sons[version][0], sons[version][1] };
              Hermes::vector<TrfShapeExp>* sub_svals[2] = { &svals[tr[version][0]], &svals[tr[version][1]] };
              Hermes::vector<TrfShapeExp>* sub_ortho_svals[2] = { &ortho_svals[tr[version][0]], &ortho_svals[tr[version][1]] };
              calc_error_cand_element(mode, gip_points, num_gip_points
                , 2, sub_domains, sub_trfs, sub_rval
                , sub_svals, sub_ortho_svals
                , info_aniso, anisoerr[version], rval);
            }
          }
        }

        //P-candidates
        if (!info_p.is_empty())
        {
          if (base_element->active)
          {
            Trf* sub_trfs[4] = { &trfs[0], &trfs[1], &trfs[2], &trfs[3] };
            int sub_rval[4] = { 0, 1, 2, 3 };
            Hermes::vector<TrfShapeExp>* sub_svals[4] = { &svals[0], &svals[1], &svals[2], &svals[3] };
            Hermes::vector<TrfShapeExp>* sub_ortho_svals[4] = { &ortho_svals[0], &ortho_svals[1], &ortho_svals[2], &ortho_svals[3] };
            Element* sub_domains[4] = { base_element, base_element, base_element, base_element };

            calc_error_cand_element(mode, gip_points, num_gip_points
              , 4, sub_domains, sub_trfs, sub_rval
              , sub_svals, sub_ortho_svals
              , info_p, perr, rval);
          }
          else
          {
            Trf* sub_trfs[4] = { &trfs[0], &trfs[1], &trfs[2], &trfs[3] };
            int sub_rval[4] = { 0, 1, 2, 3 };
            Hermes::vector<TrfShapeExp>* sub_svals[4] = { &svals[0], &svals[1], &svals[2], &svals[3] };
            Hermes::vector<TrfShapeExp>* sub_ortho_svals[4] = { &ortho_svals[0], &ortho_svals[1], &ortho_svals[2], &ortho_svals[3] };

            calc_error_cand_element(mode, gip_points, num_gip_points
              , 4, base_element->sons, sub_trfs, sub_rval
              , sub_svals, sub_ortho_svals
              , info_p, perr, rval);
          }
        }
#pragma endregion

        for (int son = 0; son < H2D_MAX_ELEMENT_SONS; son++)
          this->free_ref_solution_data(son, rval);
      }

      template<typename Scalar>
      void ProjBasedSelector<Scalar>::calc_error_cand_element(const ElementMode2D mode
        , double3* gip_points, int num_gip_points
        , const int num_sub, Element** sub_domains, Trf** sub_trfs, int* sons
        , Hermes::vector<TrfShapeExp>** sub_nonortho_svals, Hermes::vector<TrfShapeExp>** sub_ortho_svals
        , const typename OptimumSelector<Scalar>::CandsInfo& info
        , CandElemProjError errors_squared, Scalar* rval[H2D_MAX_ELEMENT_SONS][MAX_NUMBER_FUNCTION_VALUES_FOR_SELECTORS]
        )
      {
        //allocate space
        int max_num_shapes = this->next_order_shape[mode][this->max_order == H2DRS_DEFAULT_ORDER ? H2DRS_MAX_ORDER : this->max_order];
        Scalar* right_side = new Scalar[max_num_shapes];
        int* shape_inxs = new int[max_num_shapes];
        int* indx = new int[max_num_shapes]; //solver data
        double* d = new double[max_num_shapes]; //solver data
        double** proj_matrix = new_matrix<double>(max_num_shapes, max_num_shapes);
        ProjMatrixCache& proj_matrices = proj_matrix_cache[mode];
        Hermes::vector<typename OptimumSelector<Scalar>::ShapeInx>& full_shape_indices = this->shape_indices[mode];

        //check whether ortho-svals are available
        bool ortho_svals_available = true;
        for (int i = 0; i < num_sub && ortho_svals_available; i++)
          ortho_svals_available &= !sub_ortho_svals[i]->empty();

        /// An array of cached right-hand side values.
        Hermes::vector< ValueCacheItem<Scalar> > nonortho_rhs_cache;
        Hermes::vector< ValueCacheItem<Scalar> > ortho_rhs_cache;
        for (int i = 0; i <= this->max_shape_inx[mode]; i++)
        {
          nonortho_rhs_cache.push_back(ValueCacheItem<Scalar>());
          ortho_rhs_cache.push_back(ValueCacheItem<Scalar>());
        }

        //calculate for all orders
        double sub_area_corr_coef = 1.0 / num_sub;
        OrderPermutator order_perm(info.min_quad_order, info.max_quad_order, mode == HERMES_MODE_TRIANGLE || info.uniform_orders);
        do
        {
          int quad_order = order_perm.get_quad_order();
          int order_h = H2D_GET_H_ORDER(quad_order), order_v = H2D_GET_V_ORDER(quad_order);

          //build a list of shape indices from the full list
          int num_shapes = 0;
          unsigned int inx_shape = 0;
          while (inx_shape < full_shape_indices.size())
          {
            typename OptimumSelector<Scalar>::ShapeInx& shape = full_shape_indices[inx_shape];
            if (order_h >= shape.order_h && order_v >= shape.order_v)
            {
              if (num_shapes >= max_num_shapes)
                throw Exceptions::Exception("more shapes than predicted, possible incosistency");
              shape_inxs[num_shapes] = shape.inx;
              num_shapes++;
            }
            inx_shape++;
          }

          //continue only if there are shapes to process
          if (num_shapes == 0)
            continue;

          bool use_ortho = ortho_svals_available && order_perm.get_order_h() == order_perm.get_order_v();

          //select a cache
          Hermes::vector< ValueCacheItem<Scalar> >& rhs_cache = use_ortho ? ortho_rhs_cache : nonortho_rhs_cache;
          Hermes::vector<TrfShapeExp>** sub_svals = use_ortho ? sub_ortho_svals : sub_nonortho_svals;

          //calculate projection matrix iff no ortho is used
          if (!use_ortho)
          {
            if (!proj_matrices[order_h][order_v])
            {
#pragma omp critical
              {
                if (!proj_matrices[order_h][order_v])
                {
                  proj_matrices[order_h][order_v] = build_projection_matrix(gip_points, num_gip_points, shape_inxs, num_shapes, mode);
                }
              }
            }

            copy_matrix(proj_matrix, proj_matrices[order_h][order_v], num_shapes, num_shapes); //copy projection matrix because original matrix will be modified
          }

          //build right side (fill cache values that are missing)
          for (int inx_sub = 0; inx_sub < num_sub; inx_sub++)
          {
            Element* this_sub_domain = sub_domains[inx_sub];
            ElemSubTrf this_sub_trf = { sub_trfs[inx_sub], 1 / sub_trfs[inx_sub]->m[0], 1 / sub_trfs[inx_sub]->m[1] };
            ElemGIP this_sub_gip = { gip_points, num_gip_points };
            Hermes::vector<TrfShapeExp>& this_sub_svals = *(sub_svals[inx_sub]);

            for (int k = 0; k < num_shapes; k++)
            {
              int shape_inx = shape_inxs[k];
              ValueCacheItem<Scalar>& shape_rhs_cache = rhs_cache[shape_inx];
              if (!shape_rhs_cache.is_valid())
              {
                TrfShapeExp empty_sub_vals;
                ElemSubShapeFunc this_sub_shape = { shape_inx, this_sub_svals.empty() ? empty_sub_vals : this_sub_svals[shape_inx] };
                shape_rhs_cache.set(shape_rhs_cache.get() + evaluate_rhs_subdomain(this_sub_domain, this_sub_gip, sons[inx_sub], this_sub_trf, this_sub_shape, rval));
              }
            }
          }

          //copy values from cache and apply area correction coefficient
          for (int k = 0; k < num_shapes; k++)
          {
            ValueCacheItem<Scalar>& rhs_cache_value = rhs_cache[shape_inxs[k]];
            right_side[k] = sub_area_corr_coef * rhs_cache_value.get();
            rhs_cache_value.mark();
          }

          //solve iff no ortho is used
          if (!use_ortho)
          {
            ludcmp(proj_matrix, num_shapes, indx, d);
            lubksb<double, Scalar>(proj_matrix, num_shapes, indx, right_side);
          }

          //calculate error
          double error_squared = 0;
          for (int inx_sub = 0; inx_sub < num_sub; inx_sub++)
          {
            Element* this_sub_domain = sub_domains[inx_sub];
            ElemSubTrf this_sub_trf = { sub_trfs[inx_sub], 1 / sub_trfs[inx_sub]->m[0], 1 / sub_trfs[inx_sub]->m[1] };
            ElemGIP this_sub_gip = { gip_points, num_gip_points };
            ElemProj elem_proj = { shape_inxs, num_shapes, *(sub_svals[inx_sub]), right_side, quad_order };

            error_squared += evaluate_error_squared_subdomain(this_sub_domain, this_sub_gip, sons[inx_sub], this_sub_trf, elem_proj, rval);
          }
          errors_squared[order_h][order_v] = error_squared * sub_area_corr_coef; //apply area correction coefficient
        } while (order_perm.next());

        free_with_check(proj_matrix, true);
        delete[] right_side;
        delete[] shape_inxs;
        delete[] indx;
        delete[] d;
      }

      template class HERMES_API ProjBasedSelector<double>;
      template class HERMES_API ProjBasedSelector<std::complex<double> >;
    }
  }
}
