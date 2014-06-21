// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file epetra.cpp
\brief EpetraMatrix and EpetraVector storage classes for Amesos, AztecOO, ... .
*/
#include "config.h"
#ifdef HAVE_EPETRA
#include "epetra.h"
#include "callstack.h"

namespace Hermes
{
  namespace Algebra
  {
    /// \brief A communicator for Epetra objects (serial version).
    static Epetra_SerialComm seq_comm;

    template<typename Scalar>
    EpetraMatrix<Scalar>::EpetraMatrix()
    {
      this->mat = nullptr;
      this->mat_im = nullptr;
      this->grph = nullptr;
      this->std_map = nullptr;
      this->owner = true;

      this->row_storage = true;
      this->col_storage = false;
    }

    template<typename Scalar>
    EpetraMatrix<Scalar>::EpetraMatrix(const EpetraMatrix<Scalar> &op) : SparseMatrix<Scalar>(op)
    {
      this->mat = new Epetra_CrsMatrix(*op.mat);
      if (op.mat_im)
        this->mat_im = new Epetra_CrsMatrix(*op.mat_im);
      else
        this->mat_im = nullptr;

      assert(mat != nullptr);
      this->grph = new Epetra_CrsGraph(this->mat->Graph());
      this->std_map = new Epetra_BlockMap(this->grph->Map());
      this->owner = true;

      this->row_storage = op.row_storage;
      this->col_storage = op.col_storage;
    }

    template<typename Scalar>
    EpetraMatrix<Scalar>::EpetraMatrix(Epetra_RowMatrix &op)
    {
      this->mat = dynamic_cast<Epetra_CrsMatrix*>(&op);
      assert(mat != nullptr);

      this->size = this->mat->NumGlobalRows();
      this->grph = const_cast<Epetra_CrsGraph*>(&this->mat->Graph());
      this->std_map = const_cast<Epetra_BlockMap*>(&this->grph->Map());
      this->owner = false;

      this->row_storage = true;
      this->col_storage = false;
    }

    template<typename Scalar>
    EpetraMatrix<Scalar>::~EpetraMatrix()
    {
      free();
    }

    template<typename Scalar>
    void EpetraMatrix<Scalar>::prealloc(unsigned int n)
    {
      this->size = n;
      // alloc trilinos structs
      std_map = new Epetra_Map(n, 0, seq_comm);
      grph = new Epetra_CrsGraph(Copy, *std_map, 0);
    }

    template<typename Scalar>
    void EpetraMatrix<Scalar>::pre_add_ij(unsigned int row, unsigned int col)
    {
      int col_to_pass = col;
      grph->InsertGlobalIndices(row, 1, &col_to_pass);
    }

    template<>
    void EpetraMatrix<double>::finish()
    {
      mat->FillComplete();
    }

    template<>
    void EpetraMatrix<std::complex<double> >::finish()
    {
      mat->FillComplete();
      mat_im->FillComplete();
    }

    template<>
    void EpetraMatrix<double>::alloc()
    {
      // create the matrix
      grph->FillComplete();
      mat = new Epetra_CrsMatrix(Copy, *grph);
    }

    template<>
    void EpetraMatrix<std::complex<double> >::alloc()
    {
      grph->FillComplete();
      // create the matrix
      mat = new Epetra_CrsMatrix(Copy, *grph);
      mat_im = new Epetra_CrsMatrix(Copy, *grph);
    }

    template<>
    void EpetraMatrix<double>::free()
    {
      if (owner)
      {
        delete mat; mat = nullptr;
        delete grph; grph = nullptr;
        delete std_map; std_map = nullptr;
      }
    }

    template<>
    void EpetraMatrix<std::complex<double> >::free()
    {
      if (owner)
      {
        delete mat; mat = nullptr;
        delete mat_im; mat_im = nullptr;
        delete grph; grph = nullptr;
        delete std_map; std_map = nullptr;
      }
    }

    template<typename Scalar>
    Scalar EpetraMatrix<Scalar>::get(unsigned int m, unsigned int n) const
    {
      int n_entries = mat->NumGlobalEntries(m);
      Hermes::vector<double> vals(n_entries);
      Hermes::vector<int> idxs(n_entries);
      mat->ExtractGlobalRowCopy(m, n_entries, n_entries, &vals[0], &idxs[0]);
      for (int i = 0; i < n_entries; i++)
      if (idxs[i] == (int)n)
        return vals[i];
      return 0.0;
    }

    template<typename Scalar>
    int EpetraMatrix<Scalar>::get_num_row_entries(unsigned int row)
    {
      return mat->NumGlobalEntries(row);
    }

    template<typename Scalar>
    void EpetraMatrix<Scalar>::extract_row_copy(unsigned int row, unsigned int len, unsigned int &n_entries, double *vals, unsigned int *idxs)
    {
      int* idxs_to_pass = malloc_with_check<EpetraMatrix<Scalar>, int>(len, this);
      for (unsigned int i = 0; i < len; i++)
        idxs_to_pass[i] = idxs[i];
      int n_entries_to_pass = n_entries;
      mat->ExtractGlobalRowCopy(row, len, n_entries_to_pass, vals, idxs_to_pass);
      free_with_check(idxs_to_pass);
    }

    template<>
    void EpetraMatrix<double>::zero()
    {
      mat->PutScalar(0.0);
    }

    template<>
    void EpetraMatrix<std::complex<double> >::zero()
    {
      mat->PutScalar(0.0);
      mat_im->PutScalar(0.0);
    }

    template<>
    void EpetraMatrix<double>::add(unsigned int m, unsigned int n, double v)
    {
      if (v != 0.0)
      {    // ignore zero values
        int n_to_pass = n;
#pragma omp critical (EpetraMatrixAdd)
        {
          int ierr = mat->SumIntoMyValues(m, 1, &v, &n_to_pass);
          if (ierr != 0)
            throw Hermes::Exceptions::Exception("Failed to insert into Epetra matrix");
        }
      }
    }

    template<>
    void EpetraMatrix<std::complex<double> >::add(unsigned int m, unsigned int n, std::complex<double> v)
    {
      if (v != 0.0)
      {    // ignore zero values
        double v_r = std::real<double>(v);
        int n_to_pass = n;
        int ierr = mat->SumIntoMyValues(m, 1, &v_r, &n_to_pass);
        if (ierr != 0) throw Hermes::Exceptions::Exception("Failed to insert into Epetra matrix");
        assert(ierr == 0);
        double v_i = std::imag<double>(v);
        ierr = mat_im->SumIntoMyValues(m, 1, &v_i, &n_to_pass);
        assert(ierr == 0);
      }
    }

    template<>
    void EpetraMatrix<std::complex<double> >::multiply_with_vector(std::complex<double>* vector_in, std::complex<double>*& vector_out, bool vector_out_initialized) const
    {
      SparseMatrix<std::complex<double> >::multiply_with_vector(vector_in, vector_out, vector_out_initialized);
    }

    template<>
    void EpetraMatrix<double>::multiply_with_vector(double* vector_in, double*& vector_out, bool vector_out_initialized) const
    {
      Epetra_Vector x(View, mat->OperatorDomainMap(), vector_in);
      Epetra_Vector y(mat->OperatorRangeMap());
      mat->Multiply(false, x, y);
      y.ExtractCopy(vector_out);
    }

    template<typename Scalar>
    void EpetraMatrix<Scalar>::add(unsigned int m, unsigned int n, Scalar *mat, int *rows, int *cols, const int size)
    {
      for (unsigned int i = 0; i < m; i++)        // rows
      for (unsigned int j = 0; j < n; j++)      // cols
      if (rows[i] >= 0 && cols[j] >= 0) // not Dir. dofs.
        add(rows[i], cols[j], mat[i * size + j]);
    }

    template<>
    void EpetraMatrix<double>::multiply_with_Scalar(double value)
    {
      mat->Scale(value);
    }

    template<>
    void EpetraMatrix<std::complex<double> >::multiply_with_Scalar(std::complex<double> value)
    {
      // TODO
    }

    template<typename Scalar>
    void EpetraMatrix<Scalar>::add_sparse_matrix(SparseMatrix<Scalar>* mat)
    {
      EpetraMatrix<Scalar>* ep_mat = dynamic_cast<EpetraMatrix<Scalar>*>(mat);
      assert(ep_mat);

      EpetraExt::MatrixMatrix::Add(*ep_mat->mat, false, 1.0, *this->mat, 1.0);
      if (this->mat_im && ep_mat->mat_im)
        EpetraExt::MatrixMatrix::Add(*ep_mat->mat_im, false, 1.0, *this->mat_im, 1.0);
    }

    template<>
    void EpetraMatrix<double>::export_to_file(const char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format)
    {
      throw Exceptions::MethodNotImplementedException("EpetraMatrix<double>::export_to_file");
      /*
      switch (fmt)
      {
      case DF_MATLAB_SPARSE:
      case EXPORT_FORMAT_PLAIN_ASCII:
      EpetraExt::RowMatrixToHandle(file, *this->mat);
      }
      */
    }

    template<>
    void EpetraMatrix<std::complex<double> >::export_to_file(const char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format)
    {
      throw Exceptions::MethodNotImplementedException("EpetraMatrix<double>::export_to_file");
    }

    template<typename Scalar>
    unsigned int EpetraMatrix<Scalar>::get_matrix_size() const
    {
      return this->size;
    }

    template<typename Scalar>
    double EpetraMatrix<Scalar>::get_fill_in() const
    {
      return mat->NumGlobalNonzeros() / ((double)this->size*this->size);
    }

    template<typename Scalar>
    unsigned int EpetraMatrix<Scalar>::get_nnz() const
    {
      return mat->NumGlobalNonzeros();
    }

    template<typename Scalar>
    EpetraVector<Scalar>::EpetraVector()
    {
      this->std_map = nullptr;
      this->vec = nullptr;
      this->vec_im = nullptr;
      this->size = 0;
      this->owner = true;
    }

    template<typename Scalar>
    EpetraVector<Scalar>::EpetraVector(const Epetra_Vector &v)
    {
      this->vec = (Epetra_Vector *)&v;
      this->vec_im = nullptr;
      this->std_map = (Epetra_BlockMap *)&v.Map();
      this->size = v.MyLength();
      this->owner = false;
    }

    template<typename Scalar>
    EpetraVector<Scalar>::~EpetraVector()
    {
      if (owner) free();
    }

    template<typename Scalar>
    void EpetraVector<Scalar>::alloc(unsigned int n)
    {
      free();
      this->size = n;
      std_map = new Epetra_Map(this->size, 0, seq_comm);
      vec = new Epetra_Vector(*std_map);
      vec_im = new Epetra_Vector(*std_map);
      zero();
    }

    template<typename Scalar>
    void EpetraVector<Scalar>::zero()
    {
      for (unsigned int i = 0; i < this->size; i++) (*vec)[i] = 0.0;
      if (vec_im)
      for (unsigned int i = 0; i < this->size; i++) (*vec_im)[i] = 0.0;
    }

    template<typename Scalar>
    Vector<Scalar>* EpetraVector<Scalar>::change_sign()
    {
      for (unsigned int i = 0; i < this->size; i++)
        (*vec)[i] *= -1.;
      for (unsigned int i = 0; i < this->size; i++)
        (*vec_im)[i] *= -1.;
      return this;
    }

    template<typename Scalar>
    void EpetraVector<Scalar>::free()
    {
      if (this->owner)
      {
        delete std_map; std_map = nullptr;
        delete vec; vec = nullptr;
        delete vec_im; vec_im = nullptr;
        this->size = 0;
      }
    }

    template<>
    void EpetraVector<double>::set(unsigned int idx, double y)
    {
      (*vec)[idx] = y;
    }

    template<>
    void EpetraVector<std::complex<double> >::set(unsigned int idx, std::complex<double> y)
    {
      (*vec)[idx] = std::real(y);
      (*vec_im)[idx] = std::imag(y);
    }

    template<>
    void EpetraVector<double>::add(unsigned int idx, double y)
    {
      (*vec)[idx] += y;
    }

    template<>
    void EpetraVector<std::complex<double> >::add(unsigned int idx, std::complex<double> y)
    {
      (*vec)[idx] += std::real(y);
      (*vec_im)[idx] += std::imag(y);
    }

    template<typename Scalar>
    void EpetraVector<Scalar>::add(unsigned int n, unsigned int *idx, Scalar *y)
    {
      for (unsigned int i = 0; i < n; i++)
        add(idx[i], y[i]);
    }

    template<typename Scalar>
    Scalar EpetraVector<Scalar>::get(unsigned int idx) const
    {
      return (*vec)[idx];
    }

    template<typename Scalar>
    void EpetraVector<Scalar>::extract(Scalar *v) const
    {
      vec->ExtractCopy((double *)v); ///< \todo this can't be used with complex numbers
    }

    template<typename Scalar>
    void EpetraVector<Scalar>::export_to_file(const char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format)
    {
      switch (fmt)
      {
      case EXPORT_FORMAT_MATRIX_MARKET:
        EpetraExt::VectorToMatrixMarketFile(filename, *this->vec);
        return;
      case EXPORT_FORMAT_MATLAB_MATIO:
      case EXPORT_FORMAT_PLAIN_ASCII:
        EpetraExt::VectorToMatlabFile(filename, *this->vec);
        return;
      }
    }

    template class HERMES_API EpetraMatrix<double>;
    template class HERMES_API EpetraMatrix<std::complex<double> >;
    template class HERMES_API EpetraVector<double>;
    template class HERMES_API EpetraVector<std::complex<double> >;
  }
}
#endif