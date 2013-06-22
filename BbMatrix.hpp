#ifndef GUARD_BBMATRIX_HPP
#define GUARD_BBMATRIX_HPP

#include "ColumnIterator.hpp"
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>

namespace BbMath
{
	// Forward declarations
	class BbVector;

	//========================================================================================//
	//==================================== CLASS BbMatrix ====================================//
	//========================================================================================//

	class BbMatrix
	{
	public:
		// Typedefs for the user
		typedef double		   value_type;
		typedef double&		   reference;
		typedef const double&  const_reference;
		typedef double*		   RowIterator;
		typedef const double*  Const_RowIterator;
		typedef std::size_t	   size_type;
		typedef std::ptrdiff_t difference_type;

	private:
		// Pointer to the matrix
		double** myMatrix;
		size_type nRows;
		size_type nColumns;
		bool hasSize;

		std::allocator<double>  auxAlloc;
		std::allocator<double*> matrixAlloc;

	private:
		void create();
		void create(size_type r, size_type c, const double& val);
		void create(const BbMatrix& other);
		void uncreate();
		void delete_matrix();
		void delete_array(size_type i);

		void shrink(size_type newRows, size_type newColumns);

		bool range_check(size_type i, size_type j, const std::string& msg = "") const;

	public:
		// ****************** Constructors
		BbMatrix()												   { create(); }
		BbMatrix(size_type r, size_type c, const double& val = 0.) { create(r, c, val); }
		BbMatrix(const BbMatrix& other)							   { create(other); }
		BbMatrix(BbMatrix&& other)
			: myMatrix(other.myMatrix), nRows(other.nRows), nColumns(other.nColumns), hasSize(true)
		{ other.myMatrix = nullptr; other.hasSize = false;	}
		BbMatrix(const std::initializer_list<std::initializer_list<double>>& lists);
		explicit BbMatrix(const std::vector<std::vector<double>>& vec);

		// ****************** Destructor
		~BbMatrix() { uncreate(); }

		// ****************** Indexing - with and without range check
			  double& operator()(size_type i, size_type j)		 { range_check(i, j); return myMatrix[i][j]; }
		const double& operator()(size_type i, size_type j) const { range_check(i, j); return myMatrix[i][j]; }
			  double* operator[](size_type k)					 { return myMatrix[k]; }
		const double* operator[](size_type k) const				 { return myMatrix[k]; }

		// ****************** Dynamic allocation
		void resize(size_type newRows, size_type newColumns);

		bool has_size() const		{ return hasSize; }
		size_type n_rows() const	{ return nRows; }
		size_type n_columns() const	{ return nColumns; }
		bool is_square() const		{ return nRows == nColumns; }

		// ****************** Assignment
		BbMatrix& operator=(double val);
		BbMatrix& operator=(const BbMatrix& other);
		BbMatrix& operator=(BbMatrix&& other);

		BbVector get_row(size_type i) const;
		void set_row(size_type i, const double& val);
		void set_row(size_type i, const BbVector& row);
		BbVector get_column(size_type j) const;
		void set_column(size_type j, const double& val);
		void set_column(size_type j, const BbVector& col);
		BbVector get_diagonal() const;
		void set_diagonal(const double& val);
		void set_diagonal(const BbVector& diag);

		void swap_rows(size_type swapped, size_type with);
		void swap_columns(size_type swapped, size_type with);
		void seek_pivot_in_column(size_type column);
		void insert_row(size_type i, const BbVector& row);
		void insert_column(size_type j, const BbVector& col);
		void push_back_row(const BbVector& row);
		void push_back_column(const BbVector& col);
		void delete_row(size_type i);
		void delete_column(size_type j);
		BbVector extract_row(size_type i);
		BbVector extract_column(size_type j);
		void chop_row();
		void chop_column();
		void normalize_row(size_type i);
		void normalize_column(size_type j);

		ColumnIterator col_begin	(size_type j) const { return ColumnIterator(myMatrix, 1, j); }
		ColumnIterator col_end		(size_type j) const { return ColumnIterator(myMatrix, nRows + 1, j); }
		RowIterator row_begin		(size_type i)		{ if (!hasSize) { return nullptr; } return myMatrix[i] + 1; }
		RowIterator row_end			(size_type i)		{ if (!hasSize) { return nullptr; } return myMatrix[i] + nColumns + 1; }
		Const_RowIterator row_cbegin(size_type i) const { if (!hasSize) { return nullptr; } return myMatrix[i] + 1; }
		Const_RowIterator row_cend  (size_type i) const { if (!hasSize) { return nullptr; } return myMatrix[i] + nColumns + 1; }

		void transpose();

		friend std::ostream& operator<<(std::ostream& os, const BbMatrix& mat);
	};

	inline void BbMatrix::delete_array(size_type i)
	{
		RowIterator it = myMatrix[i] + nColumns + 1;

		while (it != myMatrix[i] + 1)
		{
			auxAlloc.destroy(--it);
		}

		auxAlloc.deallocate(myMatrix[i], nColumns + 1);
	}

	inline bool BbMatrix::range_check(size_type i, size_type j, const std::string& msg) const
	{		
		if (i < 1 || i > nRows)
			throw std::out_of_range(msg + " - Index i is out of range!!\n");

		if (j < 1 || j > nColumns)
			throw std::out_of_range(msg + " - Index j is out of range!!\n");

		return true;
	}

	std::ostream& operator<<(std::ostream& os, const BbMatrix& mat);
	BbVector product_Ax(const BbMatrix& lval, const BbVector& rval);
	BbMatrix product_AB(const BbMatrix& lval, const BbMatrix& rval);
	BbMatrix transpose(const BbMatrix& mat);

	BbVector solve_Dx(const BbMatrix& D, const BbVector& b);
	BbVector solve_Rx(const BbMatrix& R, const BbVector& b);
	BbVector solve_Lx(const BbMatrix& L, const BbVector& b);
	BbVector solve_gauss_elimination(const BbMatrix& A, const BbVector& b);
	BbVector solve_gauss_factorization(const BbMatrix& A, const BbVector& b);
    void	 solve_gauss_factorization(const BbMatrix& A, BbVector* bx);
}

#endif
