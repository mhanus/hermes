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
/*! \file cs_matrix.cpp
\brief Basic cs (Compressed sparse) matrix classes and operations.
*/
#include "cs_matrix.h"

namespace Hermes
{
  namespace Algebra
  {
    static int find_position(int *Ai, int Alen, int idx)
    {
      assert(Ai != NULL);
      assert(Alen > 0);
      assert(idx >= 0);

      register int lo = 0, hi = Alen - 1, mid;

      while (true)
      {
        mid = (lo + hi) >> 1;

        if(idx < Ai[mid]) hi = mid - 1;
        else if(idx > Ai[mid]) lo = mid + 1;
        else break;

        // Sparse matrix entry not found (raise an error when trying to add
        // value to this position, return 0 when obtaining value there).
        if(lo > hi)
        {
          mid = -1;
          break;
        }
      }
      return mid;
    }

    template<typename Scalar>
    CSMatrix<Scalar>::CSMatrix()
    {
      this->size = 0; nnz = 0;
      Ap = NULL;
      Ai = NULL;
      Ax = NULL;
    }

    template<typename Scalar>
    CSMatrix<Scalar>::CSMatrix(unsigned int size)
    {
      this->size = size;
      this->alloc();
    }

    template<typename Scalar>
    CSMatrix<Scalar>::~CSMatrix()
    {
      free();
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::multiply_with_Scalar(Scalar value)
    {
      for (unsigned int i = 0; i < this->nnz; i++)
        Ax[i] *= value;
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::alloc()
    {
      assert(this->pages != NULL);

      // initialize the arrays Ap and Ai
      Ap = new int[this->size + 1];
      int aisize = this->get_num_indices();
      Ai = new int[aisize];

      // sort the indices and remove duplicities, insert into Ai
      unsigned int i;
      int pos = 0;
      for (i = 0; i < this->size; i++)
      {
        Ap[i] = pos;
        pos += this->sort_and_store_indices(this->pages[i], Ai + pos, Ai + aisize);
      }
      Ap[i] = pos;

      delete [] this->pages;
      this->pages = NULL;

      nnz = Ap[this->size];

      Ax = new Scalar[nnz];
      memset(Ax, 0, sizeof(Scalar) * nnz);
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::free()
    {
      nnz = 0;
      if(Ap)
      {
        delete [] Ap;
        Ap = NULL;
      }
      if(Ai)
      {
        delete [] Ai;
        Ai = NULL;
      }
      if(Ax)
      {
        delete [] Ax;
        Ax = NULL;
      }
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::set_row_zero(unsigned int n)
    {
      for(int i = 0; i < Ap[n + 1] - Ap[n]; i++)
        Ax[Ap[n] + i] = Scalar(0);
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::zero()
    {
      memset(Ax, 0, sizeof(Scalar) * nnz);
    }

    template<typename Scalar>
    unsigned int CSMatrix<Scalar>::get_nnz() const
    {
      return this->nnz;
    }

    double inline real(double x)
    {
      return x;
    }

    double inline imag(double x)
    {
      return 0;
    }

    double inline real(std::complex<double> x)
    {
      return x.real();
    }

    double inline imag(std::complex<double> x)
    {
      return x.imag();;
    }

    template<typename Scalar>
    double CSMatrix<Scalar>::get_fill_in() const
    {
      return nnz / (double) (this->size * this->size);
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::create(unsigned int size, unsigned int nnz, int* ap, int* ai, Scalar* ax)
    {
      this->nnz = nnz;
      this->size = size;
      this->Ap = new int[this->size + 1];
      this->Ai = new int[nnz];
      this->Ax = new Scalar[nnz];
      for (unsigned int i = 0; i < this->size + 1; i++)
        this->Ap[i] = ap[i];
      for (unsigned int i = 0; i < nnz; i++)
      {
        this->Ax[i] = ax[i];
        this->Ai[i] = ai[i];
      }
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::switch_orientation()
    {
      // The variable names are so to reflect CSC -> CSR direction.
      // From the "Ap indexed by columns" to "Ap indexed by rows".
      int* tempAp = new int[this->size + 1];
      int* tempAi = new int[nnz];
      Scalar* tempAx = new Scalar[nnz];

      int run_i = 0;
      for(int target_row = 0; target_row < this->size; target_row++)
      {
        tempAp[target_row] = run_i;
        for(int src_column = 0; src_column < this->size; src_column++)
        {
          for(int src_row = this->Ap[src_column]; src_row < this->Ap[src_column + 1]; src_row++)
          {
            if(this->Ai[src_row] == target_row)
            {
              tempAi[run_i] = src_column;
              tempAx[run_i++] = this->Ax[src_row];
            }
          }
        }
      }

      tempAp[this->size] = this->nnz;
      memcpy(this->Ai, tempAi, sizeof(int) * nnz);
      memcpy(this->Ap, tempAp, sizeof(int) * (this->size + 1));
      memcpy(this->Ax, tempAx, sizeof(Scalar) * nnz);
      delete [] tempAi;
      delete [] tempAx;
      delete [] tempAp;
    }

    template<typename Scalar>
    int *CSMatrix<Scalar>::get_Ap() const
    {
      return this->Ap;
    }

    template<typename Scalar>
    int *CSMatrix<Scalar>::get_Ai() const
    {
      return this->Ai;
    }

    template<typename Scalar>
    Scalar *CSMatrix<Scalar>::get_Ax() const
    {
      return this->Ax;
    }

    template<>
    void CSMatrix<double>::add(unsigned int m, unsigned int n, double v)
    {
      if(v != 0.0)   // ignore zero values.
      {
        // Find m-th row in the n-th column.
        int pos = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
        // Make sure we are adding to an existing non-zero entry.
        if(pos < 0)
        {
          this->info("CSMatrix<Scalar>::add(): i = %d, j = %d.", m, n);
          throw Hermes::Exceptions::Exception("Sparse matrix entry not found: [%i, %i]", m, n);
        }

#pragma omp atomic
        Ax[Ap[n] + pos] += v;
      }
    }

    template<>
    void CSMatrix<std::complex<double> >::add(unsigned int m, unsigned int n, std::complex<double> v)
    {
      if(v != 0.0)   // ignore zero values.
      {
        // Find m-th row in the n-th column.
        int pos = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
        // Make sure we are adding to an existing non-zero entry.
        if(pos < 0)
        {
          this->info("CSMatrix<Scalar>::add(): i = %d, j = %d.", m, n);
          throw Hermes::Exceptions::Exception("Sparse matrix entry not found: [%i, %i]", m, n);
        }

#pragma omp critical
        Ax[Ap[n] + pos] += v;
      }
    }

    template<typename Scalar>
    Scalar CSMatrix<Scalar>::get(unsigned int m, unsigned int n) const
    {
      // Find m-th row in the n-th column.
      int mid = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);

      if(mid < 0) // if the entry has not been found
        return 0.0;
      else
        return Ax[Ap[n] + mid];
    }

    template<typename Scalar>
    CSCMatrix<Scalar>::CSCMatrix() : CSMatrix<Scalar>()
    {
    }

    template<typename Scalar>
    CSCMatrix<Scalar>::CSCMatrix(unsigned int size) : CSMatrix<Scalar>(size)
    {
    }

    template<typename Scalar>
    CSCMatrix<Scalar>::~CSCMatrix()
    {
    }

    template<>
    void CSCMatrix<double>::add(unsigned int m, unsigned int n, double v)
    {
      CSMatrix<double>::add(m, n, v);
    }

    template<>
    void CSCMatrix<std::complex<double> >::add(unsigned int m, unsigned int n, std::complex<double> v)
    {
      CSMatrix<std::complex<double> >::add(m, n, v);
    }
    
    template<typename Scalar>
    void CSCMatrix<Scalar>::add_as_block(unsigned int offset_i, unsigned int offset_j, SparseMatrix<Scalar>* mat)
    {
      CSCMatrix<Scalar>* mat_cast = dynamic_cast<CSCMatrix<Scalar>*>(mat);
      if(!mat_cast)
        throw Hermes::Exceptions::Exception("Wrong matrix type detected in CSCMatrix<Scalar>::add_as_block().");

      CSCMatrix<Scalar>::Iterator mat_it(mat_cast);
      CSCMatrix<Scalar>::Iterator this_it(this);

      // Sanity check.
      bool this_not_empty = this_it.init();
      if(!this_not_empty)
        throw Hermes::Exceptions::Exception("Empty matrix detected in CSCMatrix<Scalar>::add_as_block_csc().");

      // Iterate through the small matrix column by column and add all nonzeros
      // to the large one.
      bool mat_not_finished = mat_it.init();
      if(!mat_not_finished)
        throw Hermes::Exceptions::Exception("Empty matrix detected in CSCMatrix<Scalar>::add_as_block_csc().");

      int mat_i, mat_j;
      Scalar mat_val;
      while(mat_not_finished)
      {
        mat_it.get_current_position(mat_i, mat_j, mat_val);
        bool found = this_it.move_to_position(mat_i + offset_i, mat_j + offset_j);
        if(!found)
          throw Hermes::Exceptions::Exception("Nonzero matrix entry at %d, %d not found in CSCMatrix<Scalar>::add_as_block().",
          mat_i + offset_i, mat_j + offset_j);
        this_it.add_to_current_position(mat_val);
        mat_not_finished = mat_it.move_ptr();
      }
    }

    template<typename Scalar>
    void CSCMatrix<Scalar>::add_sparse_matrix(SparseMatrix<Scalar>* mat)
    {
      assert(this->get_size() == mat->get_size());
      
      CSCMatrix<Scalar>* mat_cast = dynamic_cast<CSCMatrix<Scalar>*>(mat);
      if(!mat_cast)
        throw Hermes::Exceptions::Exception("Wrong matrix type detected in CSCMatrix<Scalar>::add_sparse_matrix().");
      
      // Create iterators for both matrices.
      CSCMatrix<Scalar>::Iterator mat_it(mat_cast);
      CSCMatrix<Scalar>::Iterator this_it(this);
      
      int mat_i, mat_j;
      Scalar mat_val;
      int this_i, this_j;
      Scalar this_val;

      bool mat_not_finished = mat_it.init();
      bool this_not_finished = this_it.init();
      while(mat_not_finished && this_not_finished)
      {
        mat_it.get_current_position(mat_i, mat_j, mat_val);
        //printf("mat: current position %d %d %g\n", mat_i, mat_j, mat_val);
        this_it.get_current_position(this_i, this_j, this_val);
        //printf("this: current position %d %d %g\n", this_i, this_j, this_val);
        while(mat_i != this_i || mat_j != this_j)
        {
          //printf("SHOULD NOT BE HERE\n");
          this_not_finished = this_it.move_ptr();
          if(!this_not_finished)
          {
            printf("Entry %d %d does not exist in the matrix to which it is contributed.\n", mat_i, mat_j);
            throw Hermes::Exceptions::Exception("Incompatible matrices in add_csc_matrix().");
          }
          this_it.get_current_position(this_i, this_j, this_val);
        }
        this_it.add_to_current_position(mat_val);
        mat_not_finished = mat_it.move_ptr();
        this_not_finished = this_it.move_ptr();
        if(mat_not_finished && !this_not_finished)
          throw Hermes::Exceptions::Exception("Incompatible matrices in add_csc_matrix().");
      }
    }

    template<typename Scalar>
    Scalar CSCMatrix<Scalar>::get(unsigned int m, unsigned int n) const
    {
      return CSMatrix<Scalar>::get(m, n);
    }

    template<typename Scalar>
    void CSCMatrix<Scalar>::multiply_with_vector(Scalar* vector_in, Scalar*& vector_out, bool vector_out_initialized) const
    {
      if(!vector_out_initialized)
        vector_out = new Scalar[this->size];
      memset(vector_out, 0, sizeof(Scalar) * this->size);
      {
        for(int i = 0; i < this->size; i++)
        {
          for(int j = 0; j < this->Ap[i + 1] - this->Ap[i]; j++)
           vector_out[this->Ai[this->Ap[i] + j]] += this->Ax[this->Ap[i] + j] * vector_in[i];
        }
      }
    }

    static int i_coordinate(int i, int j, bool invert)
    {
      if(invert)
        return i;
      else
        return j;
    }

    static int j_coordinate(int i, int j, bool invert)
    {
      if(invert)
        return j;
      else
        return i;
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::export_to_file(const char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format, bool invert_storage)
    {
      switch (fmt)
      {
      case EXPORT_FORMAT_MATRIX_MARKET:
        {
          FILE* file = fopen(filename, "w");
          if(!file)
            throw Exceptions::IOException(Exceptions::IOException::Write, filename);
          if(Hermes::Helpers::TypeIsReal<Scalar>::value)
            fprintf(file, "%%%%Matrix<Scalar>Market matrix coordinate real\n");
          else
            fprintf(file, "%%%%Matrix<Scalar>Market matrix coordinate complex\n");

          fprintf(file, "%d %d %d\n", this->size, this->size, this->nnz);

          if(invert_storage)
            this->switch_orientation();
          for (unsigned int j = 0; j < this->size; j++)
          {
            for (int i = Ap[j]; i < Ap[j + 1]; i++)
            {
              Hermes::Helpers::fprint_coordinate_num(file, i_coordinate(Ai[i] + 1, j + 1, invert_storage), j_coordinate(Ai[i] + 1, j + 1, invert_storage), Ax[i], number_format);
              fprintf(file, "\n");
            }
          }
          if(invert_storage)
            this->switch_orientation();

          fclose(file);
        }
        break;

      case EXPORT_FORMAT_MATLAB_MATIO:
        {
#ifdef WITH_MATIO
          mat_sparse_t sparse;
          sparse.nzmax = this->nnz;
          if(invert_storage)
            this->switch_orientation();

          sparse.nir = this->nnz;
          sparse.ir = Ai;
          sparse.njc = this->size + 1;
          sparse.jc = (int *) Ap;
          sparse.ndata = this->nnz;

          size_t dims[2];
          dims[0] = this->size;
          dims[1] = this->size;

          mat_t *mat = Mat_CreateVer(filename, "", MAT_FT_MAT5);

          matvar_t *matvar;

          // For complex. No allocation here.
          double* Ax_re = NULL;
          double* Ax_im = NULL;

          // For real.
          if(Hermes::Helpers::TypeIsReal<Scalar>::value)
          {
            sparse.data = Ax;
            matvar = Mat_VarCreate(var_name, MAT_C_SPARSE, MAT_T_DOUBLE, 2, dims, &sparse, MAT_F_DONT_COPY_DATA);
          }
          else
          {
            // For complex.
            Ax_re = new double[this->nnz];
            Ax_im = new double[this->nnz];
            struct mat_complex_split_t z = {Ax_re, Ax_im};

            for(int i = 0; i < this->nnz; i++)
            {
              Ax_re[i] = ((std::complex<double>)(this->Ax[i])).real();
              Ax_im[i] = ((std::complex<double>)(this->Ax[i])).imag();
              sparse.data = &z;
            }
            matvar = Mat_VarCreate(var_name, MAT_C_SPARSE, MAT_T_DOUBLE, 2, dims, &sparse, MAT_F_DONT_COPY_DATA | MAT_F_COMPLEX);
          }

          if (matvar)
          {
            Mat_VarWrite(mat, matvar, MAT_COMPRESSION_ZLIB);
            Mat_VarFree(matvar);
          }
          if(invert_storage)
            this->switch_orientation();
          if(Ax_re)
            delete [] Ax_re;
          if(Ax_im)
            delete [] Ax_im;
          Mat_Close(mat);

          if(!matvar)
            throw Exceptions::IOException(Exceptions::IOException::Write, filename);
#else
          int ssize = sizeof(double);
          int njc = this->size+1;
     
          hermes_fwrite(&this->size, sizeof(int), 1, file);
          hermes_fwrite(&nnz, sizeof(int), 1, file);
          hermes_fwrite(&njc, sizeof(int), 1, file); 
          hermes_fwrite(&nnz, sizeof(int), 1, file);
          hermes_fwrite(&nnz, sizeof(int), 1, file);
          hermes_fwrite(Ap, sizeof(int), this->size + 1, file);
          hermes_fwrite(Ai, sizeof(int), nnz, file);
          hermes_fwrite(Ax, sizeof(double), nnz, file);	
#endif
        }
        break;

      case EXPORT_FORMAT_PLAIN_ASCII:
        {
          FILE* file = fopen(filename, "w");
          if(!file)
            throw Exceptions::IOException(Exceptions::IOException::Write, filename);

          if(invert_storage)
            this->switch_orientation();
          for (unsigned int j = 0; j < this->size; j++)
          {
            for (int i = Ap[j]; i < Ap[j + 1]; i++)
            {
              Helpers::fprint_coordinate_num(file, i_coordinate(Ai[i], j, invert_storage), j_coordinate(Ai[i], j, invert_storage), Ax[i], number_format);
              fprintf(file, "\n");
            }
          }
          if(invert_storage)
            this->switch_orientation();

          fclose(file);
        }
      }
    }

    template<typename Scalar>
    void CSCMatrix<Scalar>::export_to_file(const char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format)
    {
      CSMatrix<Scalar>::export_to_file(filename, var_name, fmt, number_format, false);
    }

    template<typename Scalar>
    void CSRMatrix<Scalar>::export_to_file(const char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format)
    {
      CSMatrix<Scalar>::export_to_file(filename, var_name, fmt, number_format, true);
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::import_from_file(const char *filename, const char *var_name, MatrixExportFormat fmt, bool invert_storage)
    {
      switch (fmt)
      {
      case EXPORT_FORMAT_MATRIX_MARKET:
        throw Exceptions::MethodNotOverridenException("CSMatrix<Scalar>::import_from_file - Matrix Market");
        break;
      case EXPORT_FORMAT_MATLAB_MATIO:
        {
#ifdef WITH_MATIO
          mat_t    *matfp;
          matvar_t *matvar;

          matfp = Mat_Open(filename,MAT_ACC_RDONLY);

          if (!matfp )
          {
            throw Exceptions::IOException(Exceptions::IOException::Read, filename);
            return;
          }

          matvar = Mat_VarRead(matfp, var_name);

          if (matvar)
          {
            mat_sparse_t *sparse = (mat_sparse_t *)matvar->data;

            this->nnz = sparse->nir;
            this->Ax = new Scalar[this->nnz];
            this->Ai = new int[this->nnz];
            this->size = sparse->njc - 1;
            this->Ap = new int[this->size + 1];

            void* data = NULL;
            if(Hermes::Helpers::TypeIsReal<Scalar>::value)
              data = sparse->data;
            else
            {
              std::complex<double>* complex_data = new std::complex<double>[this->nnz];
              double* real_array = (double*)((mat_complex_split_t*)sparse->data)->Re;
              double* imag_array = (double*)((mat_complex_split_t*)sparse->data)->Im;
              for(int i = 0; i < this->nnz; i++)
                complex_data[i] = std::complex<double>(real_array[i], imag_array[i]);
              data = (void*)complex_data;
            }
            memcpy(this->Ax, data, this->nnz * sizeof(Scalar));
            if(!Hermes::Helpers::TypeIsReal<Scalar>::value)
              delete [] data;
            memcpy(this->Ap, sparse->jc, (this->size + 1) * sizeof(Scalar));
            memcpy(this->Ai, sparse->ir, this->nnz * sizeof(int));

            if(invert_storage)
              this->switch_orientation();
          }

          Mat_Close(matfp);

          if(!matvar)
            throw Exceptions::IOException(Exceptions::IOException::Read, filename);
#else
          throw Exceptions::Exception("MATIO not included.");
#endif
        }
        break;

      case EXPORT_FORMAT_PLAIN_ASCII:
        throw Exceptions::MethodNotOverridenException("CSMatrix<Scalar>::import_from_file - Simple format");
        break;
      }
    }

    template<typename Scalar>
    void CSCMatrix<Scalar>::import_from_file(const char *filename, const char *var_name, MatrixExportFormat fmt)
    {
      CSMatrix<Scalar>::import_from_file(filename, var_name, fmt, false);
    }

    template<typename Scalar>
    void CSRMatrix<Scalar>::import_from_file(const char *filename, const char *var_name, MatrixExportFormat fmt)
    {
      CSMatrix<Scalar>::import_from_file(filename, var_name, fmt, true);
    }

    template<typename Scalar>
    CSMatrix<Scalar>* CSCMatrix<Scalar>::duplicate() const
    {
      CSCMatrix<Scalar>* new_matrix = new CSCMatrix<Scalar>();
      new_matrix->create(this->get_size(), this->get_nnz(), this->get_Ap(),  this->get_Ai(),  this->get_Ax());
      return new_matrix;
    }

    template<typename Scalar>
    CSRMatrix<Scalar>::CSRMatrix() : CSMatrix<Scalar>()
    {
    }

    template<typename Scalar>
    CSRMatrix<Scalar>::CSRMatrix(unsigned int size) : CSMatrix<Scalar>(size)
    {
    }

    template<typename Scalar>
    CSRMatrix<Scalar>::~CSRMatrix()
    {
    }

    template<>
    void CSRMatrix<double>::add(unsigned int m, unsigned int n, double v)
    {
      CSMatrix<double>::add(n, m, v);
    }

    template<>
    void CSRMatrix<std::complex<double> >::add(unsigned int m, unsigned int n, std::complex<double> v)
    {
      CSMatrix<std::complex<double> >::add(n, m, v);
    }

    template<typename Scalar>
    Scalar CSRMatrix<Scalar>::get(unsigned int m, unsigned int n) const
    {
      return CSMatrix<Scalar>::get(n, m);
    }

    template<typename Scalar>
    void CSRMatrix<Scalar>::pre_add_ij(unsigned int row, unsigned int col)
    {
      if(this->pages[row] == NULL || this->pages[row]->count >= SparseMatrix<Scalar>::PAGE_SIZE)
      {
        typename SparseMatrix<Scalar>::Page *new_page = new typename SparseMatrix<Scalar>::Page;
        new_page->count = 0;
        new_page->next = this->pages[row];
        this->pages[row] = new_page;
      }
      this->pages[row]->idx[this->pages[row]->count++] = col;
    }

    template<typename Scalar>
    SparseMatrix<Scalar>* CSRMatrix<Scalar>::duplicate() const
    {
      CSRMatrix<Scalar>* new_matrix = new CSRMatrix<Scalar>();
      new_matrix->create(this->get_size(), this->get_nnz(), this->get_Ap(),  this->get_Ai(),  this->get_Ax());
      return new_matrix;
    }
   
    template<typename Scalar>
    bool CSMatrix<Scalar>::Iterator::init()
    {
      if(this->size == 0 || this->nnz == 0) return false;
      this->Ap_pos = 0;
      this->Ai_pos = 0;
      return true;
    }

    template<typename Scalar>
    CSMatrix<Scalar>::Iterator::Iterator(CSMatrix<Scalar>* mat)
    {
      this->size = mat->get_size();
      this->nnz = mat->get_nnz();
      this->Ai = mat->get_Ai();
      this->Ap = mat->get_Ap();
      this->Ax = mat->get_Ax();
      this->Ai_pos = 0;
      this->Ap_pos = 0;
    }

    template<typename Scalar>
    bool CSMatrix<Scalar>::Iterator::move_ptr()
    {
      if(Ai_pos >= nnz - 1) return false; // It is no longer possible to find next element.
      if(Ai_pos + 1 >= Ap[Ap_pos + 1])
      {
        Ap_pos++;
      }
      Ai_pos++;
      return true;
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::Iterator::add_to_current_position(Scalar val)
    {
      this->Ax[this->Ai_pos] += val;
    }

    template<typename Scalar>
    CSCMatrix<Scalar>::Iterator::Iterator(CSCMatrix<Scalar>* mat) : CSMatrix<Scalar>::Iterator(mat)
    {
    }

    template<typename Scalar>
    void CSCMatrix<Scalar>::Iterator::get_current_position(int& i, int& j, Scalar& val)
    {
      i = this->Ai[this->Ai_pos];
      j = this->Ap_pos;
      val = this->Ax[this->Ai_pos];
    }

    template<typename Scalar>
    bool CSCMatrix<Scalar>::Iterator::move_to_position(int i, int j)
    {
      int ii, jj;
      Scalar val;
      get_current_position(ii, jj, val);
      while (!(ii == i && jj == j))
      {
        if(!this->move_ptr()) return false;
        get_current_position(ii, jj, val);
      }
      return true;
    }

    template<typename Scalar>
    CSRMatrix<Scalar>::Iterator::Iterator(CSRMatrix<Scalar>* mat) : CSMatrix<Scalar>::Iterator(mat)
    {
    }

    template<typename Scalar>
    void CSRMatrix<Scalar>::Iterator::get_current_position(int& i, int& j, Scalar& val)
    {
      i = this->Ai[this->Ai_pos];
      j = this->Ap_pos;
      val = this->Ax[this->Ai_pos];
    }

    template<typename Scalar>
    bool CSRMatrix<Scalar>::Iterator::move_to_position(int i, int j)
    {
      int ii, jj;
      Scalar val;
      get_current_position(ii, jj, val);
      while (!(ii == i && jj == j))
      {
        if(!this->move_ptr()) return false;
        get_current_position(ii, jj, val);
      }
      return true;
    }
  }
}

template class HERMES_API Hermes::Algebra::CSMatrix<double>;
template class HERMES_API Hermes::Algebra::CSMatrix<std::complex<double> >;

template class HERMES_API Hermes::Algebra::CSCMatrix<double>;
template class HERMES_API Hermes::Algebra::CSCMatrix<std::complex<double> >;

template class HERMES_API Hermes::Algebra::CSRMatrix<double>;
template class HERMES_API Hermes::Algebra::CSRMatrix<std::complex<double> >;