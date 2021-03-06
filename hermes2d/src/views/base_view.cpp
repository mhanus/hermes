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
// $Id: view4.cpp 1086 2008-10-21 09:05:44Z jakub $

#include "base_view.h"

#ifndef NOGLUT

#include <GL/freeglut.h>
#include "global.h"
#include "filter.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      template<typename Scalar>
      BaseView<Scalar>::BaseView(const char* title, WinGeom* wg)
        : ScalarView((char*)title, wg), sln(MeshFunctionSharedPtr<Scalar>(new Solution<Scalar>())), complex_filter(MeshFunctionSharedPtr<double>())
      {
        pss = nullptr;
        this->show_edges = true;
        basic_title.assign(title);
      }

      template<typename Scalar>
      BaseView<Scalar>::BaseView(char* title, WinGeom* wg)
        : ScalarView(title, wg), sln(MeshFunctionSharedPtr<Scalar>(new Solution<Scalar>())), complex_filter(MeshFunctionSharedPtr<double>())
      {
        pss = nullptr;
        this->show_edges = true;
        basic_title.assign(title);
      }

      template<typename Scalar>
      void BaseView<Scalar>::show(SpaceSharedPtr<Scalar> space, int item)
      {
        this->space = space;
        free();
        pss = new PrecalcShapeset(space->shapeset);
        ndof = this->space->get_num_dofs();
        base_index = 0;
        this->eps = eps;
        this->item = item;
        update_solution();
      }

      template<typename Scalar>
      void BaseView<Scalar>::free()
      {
        if (pss != nullptr) { delete pss; pss = nullptr; }
      }

      template<>
      void BaseView<double>::update_solution()
      {
        double* coeffs = malloc_with_check<double>(ndof);
        memset(coeffs, 0, sizeof(double)* ndof);
        if (base_index >= 0)
        {
          if (base_index < ndof)
            coeffs[base_index] = 1.0;
          Solution<double>::vector_to_solution(coeffs, space, sln, pss, false);
        }
        else
        {
          Solution<double>::vector_to_solution(coeffs, space, sln, pss, true);
        }

        ScalarView::show(sln, item);
        update_title();

        free_with_check(coeffs);
      }
      template<>
      void BaseView<std::complex<double> >::update_solution()
      {
        std::complex<double>* coeffs = malloc_with_check<std::complex<double> >(ndof);
        memset(coeffs, 0, sizeof(std::complex<double>) * ndof);
        if (base_index >= 0)
        {
          if (base_index < ndof) coeffs[base_index] = 1.0;
          Solution<std::complex<double> >::vector_to_solution(coeffs, space, sln, pss, false);
        }
        else
        {
          Solution<std::complex<double> >::vector_to_solution(coeffs, space, sln, pss, true);
        }


        complex_filter = MeshFunctionSharedPtr<double>(new Hermes::Hermes2D::RealFilter(sln.get()));
        ScalarView::show(complex_filter, item);
        update_title();

        free_with_check(coeffs);
      }

      template<typename Scalar>
      void BaseView<Scalar>::update_title()
      {
        std::stringstream str;
        str << basic_title << " - dof = " << base_index;
        if (base_index < 0)
          str << " (Dirichlet lift)";
        View::set_title(str.str().c_str());
      }

      template<typename Scalar>
      void BaseView<Scalar>::on_special_key(int key, int x, int y)
      {
        switch (key)
        {
        case GLUT_KEY_LEFT:
          if (base_index > -1) base_index--;
          update_solution();
          break;

        case GLUT_KEY_RIGHT:
          if (base_index < ndof - 1) base_index++;
          update_solution();
          break;

        default:
          ScalarView::on_special_key(key, x, y);
        }
      }

      template<typename Scalar>
      const char* BaseView<Scalar>::get_help_text() const
      {
        return
          "BaseView\n\n"
          "Controls:\n"
          "  Left mouse - pan\n"
          "  Right mouse - zoom\n"
          "  Left arrow - previous basis function\n"
          "  Right arrow - next basis function\n"
          "  3 - toggle 3D mode\n"
          "  C - center image\n"
          "  F - toggle smooth palette\n"
          "  H - render high-quality frame\n"
          "  M - toggle mesh\n"
          "  P - cycle palettes\n"
          "  S - save screenshot\n"
          "  F1 - this help\n"
          "  Esc, Q - quit\n\n"
          "3D mode:\n"
          "  Left mouse - rotate\n"
          "  Right mouse - zoom\n"
          "  Middle mouse - pan\n"
          "  * - increase Z scale\n"
          "  / - decrease Z scale";
      }
    }
  }
}
#endif

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      template class HERMES_API BaseView<double>;
      template class HERMES_API BaseView<std::complex<double> >;
    }
  }
}
