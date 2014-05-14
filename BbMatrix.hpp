#ifndef GUARD_BBMATRIX_HPP
#define GUARD_BBMATRIX_HPP

#include "BbVector.hpp"
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <iostream>
#include <algorithm>

namespace BbMath
{
	// Forward declarations
	template<typename T>
	class BbVector;

	//========================================================================================//
	//==================================== CLASS BbMatrix ====================================//
	//========================================================================================//

	template<typename T>
	class BbMatrix
	{
	public:
		// Typedefs for the user
		typedef T				value_type;
		typedef T&				reference;
		typedef const T&		const_reference;
		typedef T*				RowIterator;
		typedef const T*		Const_RowIterator;
		typedef std::size_t		size_type;
		typedef std::ptrdiff_t	difference_type;

	private:
		// Pointer to the matrix
		T** myMatrix;
		size_type nRows;
		size_type nColumns;
		bool hasSize;

		std::allocator<T>  auxAlloc;
		std::allocator<T*> matrixAlloc;

	private:
		void create()
		{
			myMatrix = nullptr;
			nRows    = 0;
			nColumns = 0;
			hasSize  = false;
		}

		void create(size_type r, size_type c, const T& val)
		{
			if (r == 0 || c == 0)
			{
				create();
				return;
			}

			myMatrix = matrixAlloc.allocate(r + 1);

			for (size_type i = 1; i <= r; ++i)
			{
				myMatrix[i] = auxAlloc.allocate(c + 1);
				std::uninitialized_fill(myMatrix[i] + 1, myMatrix[i] + 1 + c, val);
			}

			nRows    = r;
			nColumns = c;
			hasSize  = true;
		}

		void create(const BbMatrix<T>& other)
		{
			nRows    = other.nRows;
			nColumns = other.nColumns;

			myMatrix = matrixAlloc.allocate(nRows + 1);

			for (size_type i = 1; i <= nRows; ++i)
			{
				myMatrix[i] = auxAlloc.allocate(nColumns + 1);
				std::uninitialized_copy(other.myMatrix[i] + 1, other.myMatrix[i] + 1 + nColumns, myMatrix[i] + 1);
			}

			hasSize = true;
		}

		void uncreate()
		{
			if (myMatrix)
			{
				delete_matrix();
			}

			myMatrix = nullptr;
			nRows    = 0;
			nColumns = 0;
			hasSize  = false;
		}

		void delete_matrix()
		{
			T** it = myMatrix + nRows + 1;
			size_type i = nRows + 1;

			while (it != myMatrix + 1)
			{
				delete_array(--i);
				matrixAlloc.destroy(--it);
			}

			matrixAlloc.deallocate(myMatrix, nRows + 1);
		}

		void delete_array(size_type i)
		{
			RowIterator it = myMatrix[i] + nColumns + 1;

			while (it != myMatrix[i] + 1)
			{
				auxAlloc.destroy(--it);
			}

			auxAlloc.deallocate(myMatrix[i], nColumns + 1);
		}

		void shrink(size_type newRows, size_type newColumns)
		{
			T** nm = matrixAlloc.allocate(newRows + 1);

			for (size_type i = 1; i <= newRows; ++i)
			{
				nm[i] = auxAlloc.allocate(newColumns + 1);
				std::uninitialized_copy(myMatrix[i] + 1, myMatrix[i] + 1 + newColumns, nm[i] + 1);
			}

			uncreate();

			myMatrix = nm;
			nRows    = newRows;
			nColumns = newColumns;
			hasSize  = true;
		}

		bool range_check(size_type i, size_type j, const std::string& msg = "Error") const
		{
			if (i < 1 || i > nRows)
				throw std::out_of_range(msg + " - Index i is out of range!!\n");

			if (j < 1 || j > nColumns)
				throw std::out_of_range(msg + " - Index j is out of range!!\n");

			return true;
		}

	public:
		// ****************** Constructors
		BbMatrix()
		{
			create();
		}

		BbMatrix(size_type r, size_type c, const T& val = static_cast<T>(0))
		{
			create(r, c, val);
		}

		BbMatrix(const BbMatrix& other)
		{
			create(other);
		}

		BbMatrix(BbMatrix&& other)
			: myMatrix(other.myMatrix), nRows(other.nRows), nColumns(other.nColumns), hasSize(true)
		{
			other.myMatrix = nullptr; other.hasSize = false;
		}

		BbMatrix(const std::initializer_list<std::initializer_list<T>>& lists)
		{
			// Read number of rows
			nRows = lists.size();

			if (nRows == 0)
			{
				create();
				return;
			}

			// Read number of columns
			bool allEqual = true;
			nColumns = lists.begin()->size();
			std::for_each(lists.begin(), lists.end(), [&allEqual, this](const std::initializer_list<T>& l){ if (l.size() != nColumns) allEqual = false; });

			if (!allEqual)
			{
				throw std::length_error("Error in BbMatrix(const initializer_list<initializer_list<T>>&) - The length of all the rows of the matrix must be equal!!\n");
			}

			if (nColumns == 0)
			{
				create();
				return;
			}

			// Construct the matrix
			myMatrix = matrixAlloc.allocate(nRows + 1);

			int i = 1;
			for (auto it = lists.begin(); it != lists.end(); ++it, ++i)
			{
				myMatrix[i] = auxAlloc.allocate(nColumns + 1);
				std::uninitialized_copy(it->begin(), it->end(), myMatrix[i] + 1);
			}

			hasSize = true;
		}

		explicit BbMatrix(const std::vector<std::vector<T>>& vec)
		{
			// Read number of rows
			nRows = vec.size();

			if (nRows == 0)
			{
				create();
				return;
			}

			// Read number of columns
			bool allEqual = true;
			nColumns = vec[0].size();
			std::for_each(vec.begin(), vec.end(), [&allEqual, this](const std::vector<T>& v){ if (v.size() != nColumns) allEqual = false; });

			if (!allEqual)
			{
				throw std::length_error("Error in BbMatrix(const vector<vector<T>>&) - The length of all the rows of the matrix must be equal!");
			}

			if (nColumns == 0)
			{
				create();
				return;
			}

			// Construct the matrix
			myMatrix = matrixAlloc.allocate(nRows + 1);

			for (size_type i = 1; i <= nRows; ++i)
			{
				myMatrix[i] = auxAlloc.allocate(nColumns + 1);
				std::uninitialized_copy(vec[i - 1].begin(), vec[i - 1].end(), myMatrix[i] + 1);
			}

			hasSize = true;
		}

		// ****************** Destructor
		~BbMatrix()
		{
			uncreate();
		}

		// ****************** Indexing - with and without range check
			  T& operator()(size_type i, size_type j)		{ range_check(i, j); return myMatrix[i][j]; }
		const T& operator()(size_type i, size_type j) const { range_check(i, j); return myMatrix[i][j]; }
			  T* operator[](size_type k)					{ return myMatrix[k]; }
		const T* operator[](size_type k)			  const { return myMatrix[k]; }

		// ****************** Dynamic allocation
		void resize(size_type newRows, size_type newColumns)
		{
			// If new dimensions match the old ones, don't do anything
			if (newRows == nRows && newColumns == nColumns)
				return;

			// Nullify the matrix
			if (newRows == 0 || newColumns == 0)
			{
				uncreate();
				create();
				return;
			}

			// Give dimensions to the matrix
			if (!hasSize)
			{
				create(newRows, newColumns, static_cast<T>(0));
				return;
			}

			// Both new dimensions are <= than the old ones: shrink
			if ((newRows < nRows && newColumns <= nColumns) || (newColumns < nColumns && newRows <= nRows))
			{
				shrink(newRows, newColumns);
			}
			// Dimensions are different: create a new matrix initialized to static_cast<T>(0)
			else
			{
				uncreate();
				create(newRows, newColumns, static_cast<T>(0));
			}
		}

		bool has_size() const
		{
			return hasSize;
		}

		size_type n_rows() const
		{
			return nRows;
		}

		size_type n_columns() const
		{
			return nColumns;
		}

		bool is_square() const
		{
			return nRows == nColumns;
		}

		// ****************** Assignment
		BbMatrix& operator=(T val)
		{
			for (size_type i = 1; i <= nRows; ++i)
				std::fill_n(myMatrix[i] + 1, nColumns, val);

			return *this;
		}

		BbMatrix& operator=(const BbMatrix<T>& rhs)
		{
			// Check for self assignment
			if (myMatrix != rhs.myMatrix)
			{
				// If dimensions match, just copy without dealing with memory allocation
				if (nRows == rhs.nRows && nColumns == rhs.nColumns)
				{
					for (size_type i = 1; i <= nRows; ++i)
						std::copy(rhs.myMatrix[i] + 1, rhs.myMatrix[i] + 1 + nColumns, myMatrix[i] + 1);
				}
				// Otherwise, destroy *this and copy rhs
				else
				{
					uncreate();
					create(rhs);
				}
			}

			return *this;
		}

		BbMatrix& operator=(BbMatrix<T>&& rhs)
		{
			// Check for self assignment
			if (myMatrix != rhs.myMatrix)
			{
				uncreate();

				// Take ownership of the rhs
				myMatrix = rhs.myMatrix;
				nRows    = rhs.nRows;
				nColumns = rhs.nColumns;
				hasSize  = true;

				// Nullify the rhs to prevent destruction
				rhs.myMatrix = nullptr;
				rhs.hasSize  = false;
			}

			return *this;
		}

		BbVector<T> get_row(size_type i) const
		{
			static const std::string msg = "Error in get_row(i)";
			range_check(i, 1, msg);

			return BbVector<T>(row_cbegin(i), row_cend(i));
		}

		void set_row(size_type i, const T& val)
		{
			static const std::string msg = "Error in set_row(i, const T& val)";
			range_check(i, 1, msg);

			std::fill(row_begin(i), row_end(i), val);
		}

		void set_row(size_type i, const BbVector<T>& row)
		{
			static const std::string msg = "Error in set_row(i, const BbVector& val)";
			range_check(i, 1, msg);

			if (row.size() != nColumns)
				throw std::length_error(msg + " - Sizes don't match!!\n");

			std::copy(row.cbegin(), row.cend(), row_begin(i));
		}

		BbVector<T> get_column(size_type j) const
		{
			static const std::string msg = "Error in get_column(j)";
			range_check(1, j, msg);

			return BbVector<T>(col_begin(j), col_end(j));
		}

		void set_column(size_type j, const T& val)
		{
			static const std::string msg = "Error in set_column(j, const T& val)";
			range_check(1, j, msg);

			std::fill(col_begin(j), col_end(j), val);
		}

		void set_column(size_type j, const BbVector<T>& col)
		{
			static const std::string msg = "Error in set_column(j, const BbVector<T>& col)";
			range_check(1, j, msg);

			if (col.size() != nRows)
				throw std::length_error(msg + " - Sizes don't match!!\n");

			std::copy(col.cbegin(), col.cend(), col_begin(j));
		}

		BbVector<T> get_diagonal() const
		{
			if (nRows != nColumns)
				throw std::length_error("Error in get_diagonal() - Matrix isn't square!!\n");

			BbVector<T> ret(nRows);

			for (size_type k = 1; k <= nRows; ++k)
				ret[k] = myMatrix[k][k];

			return ret;
		}

		void set_diagonal(const T& val)
		{
			if (nRows != nColumns)
				throw std::length_error("Error in set_diagonal(const T& val) - Matrix isn't square!!\n");

			for (size_type k = 1; k <= nRows; ++k)
				myMatrix[k][k] = val;
		}

		void set_diagonal(const BbVector<T>& diag)
		{
			if (nRows != nColumns)
				throw std::length_error("Error in set_diagonal(const BbVector& diag) - Matrix isn't square!!\n");

			if (nRows != diag.size())
				throw std::length_error("Error in set_diagonal(const BbVector& diag) - Sizes don't match!!\n");

			for (size_type k = 1; k <= nRows; ++k)
				myMatrix[k][k] = diag[k];
		}

		void swap_rows(size_type swapped, size_type with)
		{
			if (swapped == with)
				return;

			static const std::string msg = "Error in swap_rows()";
			range_check(swapped, 1, msg);
			range_check(with,    1, msg);

			std::swap(myMatrix[swapped], myMatrix[with]);
		}

		void swap_columns(size_type swapped, size_type with)
		{
			if (swapped == with)
				return;

			static const std::string msg = "Error in swap_rows()";
			range_check(1, swapped, msg);
			range_check(1,    with, msg);

			for (size_type i = 1; i <= nRows; ++i)
				std::swap(myMatrix[i][swapped], myMatrix[i][with]);
		}

		void seek_pivot_in_column(size_type column)
		{
			const size_type j = column;
			size_type index = j;

			T max = myMatrix[j][j];

			for (size_type i = j + 1; i <= nRows; ++i)
			{
				if (myMatrix[i][j] > max)
				{
					max = myMatrix[i][j];
					index = i;
				}
			}

			if (index != j)
				this->swap_rows(j, index);
		}

		void insert_row(size_type i, const BbVector<T>& row)
		{
			if (i < 1 || i > nRows + 1)
				throw std::out_of_range("Error in insert_row(i) - Index i is out of range!!\n");

			if (row.size() != nColumns)
				throw std::length_error("Error in insert_row(i) - Sizes don't match!!\n");

			T** nm = matrixAlloc.allocate(nRows + 1 + 1);

			for (size_type k = 1; k < i; ++k)
			{
				nm[k] = auxAlloc.allocate(nColumns + 1);
				std::uninitialized_copy(myMatrix[k] + 1, myMatrix[k] + 1 + nColumns, nm[k] + 1);
			}

			nm[i] = auxAlloc.allocate(row.size() + 1);
			std::uninitialized_copy(row.start, row.avail, nm[i] + 1);

			for (size_type k = i; k <= nRows; ++k)
			{
				nm[k + 1] = auxAlloc.allocate(nColumns + 1);
				std::uninitialized_copy(myMatrix[k] + 1, myMatrix[k] + 1 + nColumns, nm[k + 1] + 1);
			}

			delete_matrix();

			myMatrix = nm;
			++nRows;
		}

		void insert_column(size_type j, const BbVector<T>& col)
		{
			if (j < 1 || j > nRows + 1)
				throw std::out_of_range("Error in insert_column(j) - Index i is out of range!!\n");

			if (col.size() != nRows)
				throw std::length_error("Error in insert_row(i) - Sizes don't match!!\n");

			T** nm = matrixAlloc.allocate(nRows + 1);

			for (size_type k = 1; k <= nRows; ++k)
			{
				nm[k] = auxAlloc.allocate(nColumns + 1 + 1);
				std::uninitialized_copy(myMatrix[k] + 1, myMatrix[k] + 1 + j - 1, nm[k] + 1);
				auxAlloc.construct(nm[k] + 1 + j - 1, col[k]);
				std::uninitialized_copy(myMatrix[k] + 1 + j - 1 , myMatrix[k] + 1 + nColumns, nm[k] + 1 + j);
			}

			delete_matrix();

			myMatrix = nm;
			++nColumns;
		}

		void push_back_column(const BbVector<T>& col)
		{
			insert_column(nColumns + 1, col);
		}

		void push_back_row(const BbVector<T>& row)
		{
			insert_row(nRows + 1, row);
		}

		void delete_row(size_type i)
		{
			static const std::string msg = "Error in delete_row(i)";
			range_check(i, 1, msg);

			size_type newRows = nRows - 1;
			T** nm = matrixAlloc.allocate(newRows + 1);

			// Copy the 1st part of the matrix into the new one
			for (size_type k = 1; k < i; ++k)
			{
				nm[k] = auxAlloc.allocate(nColumns + 1);
				std::uninitialized_copy(myMatrix[k] + 1, myMatrix[k] + 1 + nColumns, nm[k] + 1);
			}
			// Copy the 2nd part of the matrix into the new one
			for (size_type k = i + 1; k <= nRows; ++k)
			{
				nm[k - 1] = auxAlloc.allocate(nColumns + 1);
				std::uninitialized_copy(myMatrix[k] + 1, myMatrix[k] + 1 + nColumns, nm[k - 1] + 1);
			}

			delete_matrix();

			myMatrix = nm;
			nRows = newRows;
		}

		void delete_column(size_type j)
		{
			static const std::string msg = "Error in delete_column(j)";
			range_check(1, j, msg);

			size_type newColumns = nColumns - 1;

			for (size_type k = 1; k <= nRows; ++k)
			{
				T* temp = auxAlloc.allocate(newColumns + 1);
				std::uninitialized_copy(myMatrix[k] + 1, myMatrix[k] + j, temp + 1);
				std::uninitialized_copy(myMatrix[k] + 1 + j, myMatrix[k] + 1 + nColumns, temp + j);
				delete_array(k);
				myMatrix[k] = temp;
			}

			nColumns = newColumns;
		}

		BbVector<T> extract_row(size_type i)
		{
			BbVector<T> ret{ get_row(i) };
			delete_row(i);

			return ret;
		}

		BbVector<T> extract_column(size_type j)
		{
			BbVector<T> ret{ get_column(j) };
			delete_column(j);

			return ret;
		}

		void chop_row()
		{
			delete_row(nRows);
		}

		void chop_column()
		{
			delete_column(nColumns);
		}

		void normalize_row(size_type i)
		{
			static const std::string msg = "Error in normalize_row(i)";
			range_check(i, 1, msg);

			// TODO: divide by the closest power of 2, to preserve mantissa and just act on the exponent
			T denom = this->get_row(i).get_max_abs();
			if (denom < 1.E-20)
			{
				throw std::runtime_error("Error in normalize_row(i) - Singular matrix!!\n");
			}

			std::for_each(row_begin(i), row_end(i), [&denom](T& val){ val /= denom; });
		}

		void normalize_column(size_type j)
		{
			static const std::string msg = "Error in normalize_column(j)";
			range_check(1, j, msg);

			T denom = this->get_column(j).get_max_abs();
			if (abs(denom) < 1.E-20)
			{
				throw std::runtime_error("Error in normalize_column(j) - Singular matrix!!\n");
			}

			std::for_each(col_begin(j), col_end(j), [&denom](T& val){ val /= denom; });
		}

		ColumnIterator<T> col_begin	(size_type j) const { return ColumnIterator<T>(myMatrix, 1, j); }
		ColumnIterator<T> col_end	(size_type j) const { return ColumnIterator<T>(myMatrix, nRows + 1, j); }
		RowIterator row_begin		(size_type i)		{ if (!hasSize) { return nullptr; } return myMatrix[i] + 1; }
		RowIterator row_end			(size_type i)		{ if (!hasSize) { return nullptr; } return myMatrix[i] + nColumns + 1; }
		Const_RowIterator row_cbegin(size_type i) const { if (!hasSize) { return nullptr; } return myMatrix[i] + 1; }
		Const_RowIterator row_cend  (size_type i) const { if (!hasSize) { return nullptr; } return myMatrix[i] + nColumns + 1; }

		void transpose()
		{
			if (!hasSize)
				return;

			// Case of non square matrix: create a the new transposed one and take ownership
			if (nRows != nColumns)
			{
				BbMatrix<T> temp(nColumns, nRows);

				for (size_type i = 1; i <= nColumns; ++i)
				{
					for (size_type j = 1; j <= nRows; ++j)
					{
						temp[i][j] = myMatrix[j][i];
					}
				}

				uncreate();
				*this = std::move(temp);
			}
			// Case of square matrix: fast transposition
			else
			{
				for (size_type i = 1; i <= nRows; ++i)
				{
					for (size_type j = i + 1; j <= nColumns; ++j)
					{
						std::swap(myMatrix[i][j], myMatrix[j][i]);
					}
				}
			}
		}

		template<typename C>
		friend std::ostream& operator<<(std::ostream& os, const BbMatrix<C>& mat);
	};

	template<typename T>
	std::ostream& operator<<(std::ostream& os, const BbMatrix<T>& mat)
	{
		constexpr size_t maxSize = 40;

		if (!mat.has_size())
		{
			os << std::endl << "Empty BbMatrix!!" << std::endl;
		}
		else if (mat.n_rows() <= maxSize && mat.n_columns() <= maxSize)
		{
			int nRows	 = mat.n_rows();
			int nColumns = mat.n_columns();

			std::for_each(mat.myMatrix + 1, mat.myMatrix + 1 + nRows, [nColumns, &os](const T* row)
			{
				std::for_each(row + 1, row + 1 + nColumns, [&os](const T& val)
				{
					os.precision(3);
					os.width(8);
					os << val << " ";
				});

				os << std::endl;
			});
		}
		else
		{
			os << std::endl << "The matrix is too big to be printed!!" << std::endl;
		}

		return os;
	}

	template<typename T>
	BbMatrix<T> product_xyT(const BbVector<T>& lval, const BbVector<T>& rval)
	{
		if (lval.size() == 0 || rval.size() == 0)
			throw std::length_error("Error in product_xyT() - One of the vectors is empty!!\n");

		BbMatrix<T> ret(lval.size(), rval.size());

		for (size_t i = 1; i <= lval.size(); ++i)
		{
			for (size_t j = 1; j <= rval.size(); ++j)
			{
				ret[i][j] = lval[i] * rval[j];
			}
		}

		return ret;
	}

	template<typename T>
	BbVector<T> product_Ax(const BbMatrix<T>& lval, const BbVector<T>& rval)
	{
		if (lval.n_columns() != rval.size())
			throw std::length_error("Error in product_Ax() - Sizes don't match!!\n");

		BbVector<T> ret(lval.n_rows());
		BbVector<T> aux;

		for (size_t i = 1; i <= lval.n_rows(); ++i)
		{
			aux = lval.get_row(i);
			ret[i] = product_xTy(aux, rval);
		}

		return ret;
	}

	template<typename T>
	BbMatrix<T> product_AB(const BbMatrix<T>& lval, const BbMatrix<T>& rval)
	{
		if (lval.n_columns() != rval.n_rows())
			throw std::length_error("Error in product_AB() - Columns and rows don't match!!\n");

		BbMatrix<T> ret(lval.n_rows(), rval.n_columns());

		for (size_t i = 1; i <= lval.n_rows(); ++i)
		{
			for (size_t j = 1; j <= rval.n_columns(); ++j)
			{
				for (size_t k = 1; k <= lval.n_columns(); ++k)
				{
					ret[i][j] += lval[i][k] * rval[k][j];
				}
			}
		}

		return ret;
	}

	template<typename T>
	BbMatrix<T> transpose(const BbMatrix<T>& mat)
	{
		BbMatrix<T> ret(mat);
		ret.transpose();

		return ret;
	}

	using BbMatrixD = BbMatrix<double>;
	using BbMatrixI = BbMatrix<int>;

	BbVectorD solve_Dx(const BbMatrixD& D, const BbVectorD& b);
	BbVectorD solve_Rx(const BbMatrixD& R, const BbVectorD& b);
	BbVectorD solve_Lx(const BbMatrixD& L, const BbVectorD& b);
	BbVectorD solve_gauss_elimination(const BbMatrixD& A, const BbVectorD& b);
	BbVectorD solve_gauss_factorization(const BbMatrixD& A, const BbVectorD& b);
	void	  solve_gauss_factorization(const BbMatrixD& A, BbVectorD* bx);
	BbMatrixD inverse(const BbMatrixD& other);
}

#endif
