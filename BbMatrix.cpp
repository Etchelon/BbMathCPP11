#include "BbVector.hpp"
#include "BbMatrix.hpp"
#include <iostream>
#include <algorithm>
#include <stdexcept>

using namespace std;

namespace BbMath
{
	// ******************************************************************************************************** //
	// *** PRIVATE CONSTRUCTING FUNCTIONS ********************************************************************* //
	// ******************************************************************************************************** //

	void BbMatrix::create()
	{
		myMatrix = nullptr;
		nRows    = 0;
		nColumns = 0;
		hasSize  = false;
	}

	void BbMatrix::create(size_type r, size_type c, const double& val)
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
			uninitialized_fill(myMatrix[i] + 1, myMatrix[i] + 1 + c, val);
		}

		nRows    = r;
		nColumns = c;
		hasSize  = true;
	}

	void BbMatrix::create(const BbMatrix& other)
	{
		nRows    = other.nRows;
		nColumns = other.nColumns;

		myMatrix = matrixAlloc.allocate(nRows + 1);

		for (size_type i = 1; i <= nRows; ++i)
		{
			myMatrix[i] = auxAlloc.allocate(nColumns + 1);
			uninitialized_copy(other.myMatrix[i] + 1, other.myMatrix[i] + 1 + nColumns, myMatrix[i] + 1);
		}

		hasSize = true;
	}

	void BbMatrix::delete_matrix()
	{
		double** it = myMatrix + nRows + 1;
		size_type i = nRows + 1;

		while (it != myMatrix + 1)
		{
			delete_array(--i);
			matrixAlloc.destroy(--it);
		}

		matrixAlloc.deallocate(myMatrix, nRows + 1);
	}

	void BbMatrix::uncreate()
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

	// ******************************************************************************************************** //
	// ***************************************************************** END PRIVATE CONSTRUCTING FUNCTIONS *** //
	// ******************************************************************************************************** //

	// ******************************************************************************************************** //
	// *** PUBLIC CONSTRUCTORS ******************************************************************************** //
	// ******************************************************************************************************** //

	BbMatrix::BbMatrix(const initializer_list<initializer_list<double>>& lists)
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
		for_each(lists.begin(), lists.end(), [&allEqual, this](const initializer_list<double>& l){ if (l.size() != nColumns) allEqual = false; });

		if (!allEqual)
		{
			throw length_error("Error in BbMatrix(const initializer_list<initializer_list<double>>&) - The length of all the rows of the matrix must be equal!!\n");
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
			uninitialized_copy(it->begin(), it->end(), myMatrix[i] + 1);
		}

		hasSize = true;
	}

	BbMatrix::BbMatrix(const vector<vector<double>>& vec)
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
		for_each(vec.begin(), vec.end(), [&allEqual, this](const vector<double>& v){ if (v.size() != nColumns) allEqual = false; });

		if (!allEqual)
		{
			throw length_error("Error in BbMatrix(const vector<vector<double>>&) - The length of all the rows of the matrix must be equal!");
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
			uninitialized_copy(vec[i - 1].begin(), vec[i - 1].end(), myMatrix[i] + 1);
		}

		hasSize = true;
	}

	// ******************************************************************************************************** //
	// **************************************************************************** END PUBLIC CONSTRUCTORS *** //
	// ******************************************************************************************************** //

	// ******************************************************************************************************** //
	// *** ASSIGNMENT OPERATORS ******************************************************************************* //
	// ******************************************************************************************************** //

	BbMatrix& BbMatrix::operator=(double val)
	{
		for (size_type i = 1; i <= nRows; ++i)
		{
			fill_n(myMatrix[i] + 1, nColumns, val);
		}

		return *this;
	}

	BbMatrix& BbMatrix::operator=(const BbMatrix& rhs)
	{
		// Check for self assignment
		if (myMatrix != rhs.myMatrix)
		{
			// If dimensions match, just copy without dealing with memory allocation
			if (nRows == rhs.nRows && nColumns == rhs.nColumns)
			{
				for (size_type i = 1; i <= nRows; ++i)
				{
					copy(rhs.myMatrix[i] + 1, rhs.myMatrix[i] + 1 + nColumns, myMatrix[i] + 1);
				}
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

	BbMatrix& BbMatrix::operator=(BbMatrix&& rhs)
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

	// ******************************************************************************************************** //
	// *************************************************************************** END ASSIGNMENT OPERATORS *** //
	// ******************************************************************************************************** //

	// ******************************************************************************************************** //
	// *** MISCELLANEOUS FUNCTIONS **************************************************************************** //
	// ******************************************************************************************************** //

	ostream& operator<<(ostream& os, const BbMatrix& mat)
	{
        constexpr size_t maxSize = 40;

		if (!mat.has_size())
		{
			os << endl << "Empty BbMatrix!!" << endl;
		}
		else if (mat.n_rows() <= maxSize && mat.n_columns() <= maxSize)
		{
			int nRows	 = mat.n_rows();
			int nColumns = mat.n_columns();

			for_each(mat.myMatrix + 1, mat.myMatrix + 1 + nRows, [nColumns, &os](const double* row)
			{
				for_each(row + 1, row + 1 + nColumns, [&os](const double& val)
				{
					os.precision(3);
					os.width(8);
					os << val << " ";
				});

				os << endl;
			});
		}
		else
		{
			os << endl << "The matrix is too big to be printed!!" << endl;
		}

		return os;
	}

	// ******************************************************************************************************** //
	// *** ACCESS FUNCTIONS *********************************************************************************** //
	// ******************************************************************************************************** //

	BbVector BbMatrix::get_row(size_type i) const
	{
		static const string msg = "Error in get_row(i)";
		range_check(i, 1, msg);

		return BbVector(row_cbegin(i), row_cend(i));
	}

	void BbMatrix::set_row(size_type i, const double& val)
	{
		static const string msg = "Error in set_row(i, const double& val)";
		range_check(i, 1, msg);

		std::fill(row_begin(i), row_end(i), val);
	}

	void BbMatrix::set_row(size_type i, const BbVector& row)
	{
		static const string msg = "Error in set_row(i, const BbVector& val)";
		range_check(i, 1, msg);

		if (row.size() != nColumns)
		{
			throw length_error(msg + " - Sizes don't match!!\n");
		}

		std::copy(row.cbegin(), row.cend(), row_begin(i));
	}

	BbVector BbMatrix::get_column(size_type j) const
	{
		static const string msg = "Error in get_column(j)";
		range_check(1, j, msg);

		return BbVector(col_begin(j), col_end(j));
	}

	void BbMatrix::set_column(size_type j, const double& val)
	{
		static const string msg = "Error in set_column(j, const double& val)";
		range_check(1, j, msg);

		std::fill(col_begin(j), col_end(j), val);
	}

	void BbMatrix::set_column(size_type j, const BbVector& col)
	{
		static const string msg = "Error in set_column(j, const BbVector& col)";
		range_check(1, j, msg);

		if (col.size() != nRows)
		{
			throw length_error(msg + " - Sizes don't match!!\n");
		}

		std::copy(col.cbegin(), col.cend(), col_begin(j));
	}

	BbVector BbMatrix::get_diagonal() const
	{
		if (nRows != nColumns)
		{
			throw length_error("Error in get_diagonal() - Matrix isn't square!!\n");
		}

		BbVector ret(nRows);

		for (size_type k = 1; k <= nRows; ++k)
		{
			ret[k] = myMatrix[k][k];
		}

		return ret;
	}

	BbVector product_Ax(const BbMatrix& lval, const BbVector& rval)
	{
		if (lval.n_columns() != rval.size())
		{
			throw length_error("Error in product_Ax() - Sizes don't match!!\n");
		}

		BbVector ret(lval.n_rows());
		BbVector aux;

		for (size_t i = 1; i <= lval.n_rows(); ++i)
		{
			aux = lval.get_row(i);
			ret[i] = product_xTy(aux, rval);
		}

		return ret;
	}

	BbMatrix product_xyT(const BbVector& lval, const BbVector& rval)
	{
		if (lval.size() == 0 || rval.size() == 0)
		{
			throw length_error("Error in product_xyT() - One of the vectors is empty!!\n");
		}

		BbMatrix ret(lval.size(), rval.size());

		for (size_t i = 1; i <= lval.size(); ++i)
		{
			for (size_t j = 1; j <= rval.size(); ++j)
			{
				ret[i][j] = lval[i] * rval[j];
			}
		}

		return ret;
	}

	BbMatrix product_AB(const BbMatrix& lval, const BbMatrix& rval)
	{
		if (lval.n_columns() != rval.n_rows())
		{
			throw length_error("Error in product_AB() - Columns and rows don't match!!\n");
		}

		BbMatrix ret(lval.n_rows(), rval.n_columns());

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

	void BbMatrix::resize(size_type newRows, size_type newColumns)
	{
		// If new dimensions match the old ones, don't do anything
		if (newRows == nRows && newColumns == nColumns)
		{
			return;
		}

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
			create(newRows, newColumns, 0.);
			return;
		}

		// Both new dimensions are <= than the old ones: shrink
		if ((newRows < nRows && newColumns <= nColumns) || (newColumns < nColumns && newRows <= nRows))
		{
			shrink(newRows, newColumns);
		}
		// Dimensions are different: create a new matrix initialized to 0.
		else
		{
			uncreate();
			create(newRows, newColumns, 0.);
		}
	}

	void BbMatrix::shrink(size_type newRows, size_type newColumns)
	{
		double** nm = matrixAlloc.allocate(newRows + 1);

		for (size_type i = 1; i <= newRows; ++i)
		{
			nm[i] = auxAlloc.allocate(newColumns + 1);
			uninitialized_copy(myMatrix[i] + 1, myMatrix[i] + 1 + newColumns, nm[i] + 1);
		}

		uncreate();

		myMatrix = nm;
		nRows    = newRows;
		nColumns = newColumns;
		hasSize  = true;
	}

	void BbMatrix::transpose()
	{
		if (!hasSize)
		{
			return;
		}

		// Case of non square matrix: create a the new transposed one and take ownership
		if (nRows != nColumns)
		{
			BbMatrix temp(nColumns, nRows);

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

	void BbMatrix::swap_rows(size_type swapped, size_type with)
	{
		if (swapped == with)
		{
			return;
		}

		static const string msg = "Error in swap_rows()";
		range_check(swapped, 1, msg);
		range_check(with,    1, msg);

		std::swap(myMatrix[swapped], myMatrix[with]);
	}

	void BbMatrix::swap_columns(size_type swapped, size_type with)
	{
		if (swapped == with)
		{
			return;
		}

		static const string msg = "Error in swap_rows()";
		range_check(1, swapped, msg);
		range_check(1,    with, msg);

		for (size_type i = 1; i <= nRows; ++i)
		{
			std::swap(myMatrix[i][swapped], myMatrix[i][with]);
		}
	}

	void BbMatrix::set_diagonal(const double& val)
	{
		if (nRows != nColumns)
		{
			throw length_error("Error in set_diagonal(const double& val) - Matrix isn't square!!\n");
		}

		for (size_type k = 1; k <= nRows; ++k)
		{
			myMatrix[k][k] = val;
		}
	}

	void BbMatrix::set_diagonal(const BbVector& diag)
	{
		if (nRows != nColumns)
		{
			throw length_error("Error in set_diagonal(const BbVector& diag) - Matrix isn't square!!\n");
		}

		if (nRows != diag.size())
		{
			throw length_error("Error in set_diagonal(const BbVector& diag) - Sizes don't match!!\n");
		}

		for (size_type k = 1; k <= nRows; ++k)
		{
			myMatrix[k][k] = diag[k];
		}
	}

	BbMatrix transpose(const BbMatrix& mat)
	{
		BbMatrix ret(mat);
		ret.transpose();

		return ret;
	}

	// Solve linear systems
	BbVector solve_Dx(const BbMatrix& D, const BbVector& b)
	{
		if (!D.is_square())
		{
			throw length_error("Error in solve_Dx() - Matrix D isn't square!!\n");
		}

		if (D.n_rows() != b.size())
		{
			throw length_error("Error in solve_Dx() - Sizes don't match!!\n");
		}

		size_t n = b.size();
		BbVector x(n);

		for (size_t i = 1; i <= n; ++i)
		{
			x[i] = b[i] / D[i][i];
		}

		return x;
	}

	BbVector solve_Rx(const BbMatrix& R, const BbVector& b)
	{
		if (R.n_rows() != R.n_columns())
		{
			throw length_error("Error in solve_Rx() - Matrix R isn't square!!\n");
		}

		if (R.n_rows() != b.size())
		{
			throw length_error("Error in solve_Dx() - Sizes don't match!!\n");
		}

		size_t n = b.size();
		BbVector x(b);

		for (size_t i = n; i >= 1; --i)
		{
			for (size_t k = i + 1; k <= n; ++k)
			{
				x[i] -= R[i][k] * x[k];
			}

			x[i] /= R[i][i];
		}

		return x;
	}

	BbVector solve_Lx(const BbMatrix& L, const BbVector& b)
	{
		if (L.n_rows() != L.n_columns())
		{
			throw length_error("Error in solve_Lx() - Matrix L isn't square!!\n");
		}

		if (L.n_rows() != b.size())
		{
			throw length_error("Error in solve_Lx() - Sizes don't match!!\n");
		}

		size_t n = b.size();
		BbVector x(b);

		for (size_t i = 1; i <= n; ++i)
		{
			for (size_t k = 1; k <= i - 1; ++k)
			{
				x[i] -= L[i][k] * x[k];
			}

			x[i] /= L[i][i];
		}

		return x;
	}

	void BbMatrix::seek_pivot_in_column(size_type column)
	{
		const size_type j = column;
		size_type index = j;

		double max = myMatrix[j][j];

		for (size_type i = j + 1; i <= nRows; ++i)
		{
			if (myMatrix[i][j] > max)
			{
				max = myMatrix[i][j];
				index = i;
			}
		}

		if (index != j)
		{
			this->swap_rows(j, index);
		}
	}

	BbVector solve_gauss_elimination(const BbMatrix& A, const BbVector& b)
	{
		if (A.n_rows() != A.n_columns())
		{
			throw length_error("Error in solve_gauss_elimination() - Matrix A isn't square!!\n");
		}

		if (A.n_rows() != b.size())
		{
			throw length_error("Error in solve_gauss_elimination() - Sizes don't match!!\n");
		}

		size_t r = b.size();
		size_t c = r + 1;

		BbMatrix Ab(A);
		Ab.push_back_column(b);

		for (size_t k = 1; k <= r - 1; ++k)
		{
			// TODO: seek pivot

			for (size_t i = k + 1; i <= r; ++i)
			{
				if (Ab[i][k] == 0.)
				{
					continue;
				}

				double l = Ab[i][k] / Ab[k][k];

				for (size_t j = 1; j <= c; ++j)
				{
					Ab[i][j] -= l * Ab[k][j];
				}
			}
		}

		BbVector bb(Ab.get_column(c));
		Ab.chop_column();

		return solve_Rx(Ab, bb);
	}

	void BbMatrix::push_back_column(const BbVector& col)
	{
		insert_column(nColumns + 1, col);
	}

	void BbMatrix::push_back_row(const BbVector& row)
	{
		insert_row(nRows + 1, row);
	}

	void BbMatrix::delete_column(size_type j)
	{
		static const string msg = "Error in delete_column(j)";
		range_check(1, j, msg);

		size_type newColumns = nColumns - 1;

		for (size_type k = 1; k <= nRows; ++k)
		{
			double* temp = auxAlloc.allocate(newColumns + 1);
			uninitialized_copy(myMatrix[k] + 1, myMatrix[k] + j, temp + 1);
			uninitialized_copy(myMatrix[k] + 1 + j, myMatrix[k] + 1 + nColumns, temp + j);
			delete_array(k);
			myMatrix[k] = temp;
		}

		nColumns = newColumns;
	}

	BbVector BbMatrix::extract_row(BbMatrix::size_type i)
	{
		BbVector ret{ get_row(i) };
		delete_row(i);

		return ret;
	}

	BbVector BbMatrix::extract_column(BbMatrix::size_type j)
	{
		BbVector ret{ get_column(j) };
		delete_column(j);

		return ret;
	}

	void BbMatrix::chop_column()
	{
		delete_column(nColumns);
	}

	void BbMatrix::delete_row(size_type i)
	{
		static const string msg = "Error in delete_row(i)";
		range_check(i, 1, msg);

		size_type newRows = nRows - 1;
		double** nm = matrixAlloc.allocate(newRows + 1);

		// Copy the 1st part of the matrix into the new one
		for (size_type k = 1; k < i; ++k)
		{
			nm[k] = auxAlloc.allocate(nColumns + 1);
			uninitialized_copy(myMatrix[k] + 1, myMatrix[k] + 1 + nColumns, nm[k] + 1);
		}
		// Copy the 2nd part of the matrix into the new one
		for (size_type k = i + 1; k <= nRows; ++k)
		{
			nm[k - 1] = auxAlloc.allocate(nColumns + 1);
			uninitialized_copy(myMatrix[k] + 1, myMatrix[k] + 1 + nColumns, nm[k - 1] + 1);
		}

		delete_matrix();

		myMatrix = nm;
		nRows = newRows;
	}

	void BbMatrix::chop_row()
	{
		delete_row(nRows);
	}

	BbVector solve_gauss_factorization(const BbMatrix& A, const BbVector& b)
	{
		// Check requirements
		if (!A.is_square())
		{
			throw length_error("Error in solve_gauss_factorization() - Matrix A isn't square!!\n");
		}

		if (A.n_rows() != b.size())
		{
			throw length_error("Error in solve_gauss_factorization() - Sizes don't match!!\n");
		}

		size_t n = b.size();
		BbMatrix Ab(A);			// Work matrix that will be factorized
		BbVector coeff(n);		// Vector of the scaling factor for the unknowns

		// Bring the system in its standard form
		// Augment the matrix with the vector of constant terms
		Ab.push_back_column(b);

		// Normalize rows to 1
		for (size_t i = 1; i <= n; ++i)
		{
			Ab.normalize_row(i);
		}

		// Normalize columns and save the scaling factor
		for (size_t j = 1; j <= n; ++j)
		{
			double denom = abs(*(std::max_element(Ab.col_begin(j), Ab.col_end(j), [](const double& a, const double& b){ return abs(a) < abs(b); })));
			if (denom < 1.E-20)
			{
				throw runtime_error("Error in normalize_column(j) - Division by 0 in column normalization!!\n");
			}

			for_each(Ab.col_begin(j), Ab.col_end(j), [&denom](double& val){ val /= denom; });
			coeff[j] = denom;
		}

		// Factorize the matrix
		for (size_t k = 1; k <= n - 1; ++k)
		{
			Ab.seek_pivot_in_column(k);

			// Check for null pivot
			if (abs(Ab[k][k]) < 1.E-20)
			{
				throw runtime_error("Error in normalize_column(j) - Singular matrix!!\n");
			}

			// Factorize the rows
			for (size_t i = k + 1; i <= n; ++i)
			{
				if (Ab[i][k] == 0.)
				{
					continue;
				}

				Ab[i][k] /= Ab[k][k];

				for (size_t j = k + 1; j <= n; ++j)
				{
					Ab[i][j] -= Ab[i][k] * Ab[k][j];
				}
			}
		}

		// Bring the matrix back to NxN and prepare the output vector
		BbVector x = Ab.extract_column(n + 1);

		BbVector backup = Ab.get_diagonal();

		// Solve the left system
		Ab.set_diagonal(1.);
		x = std::move(solve_Lx(Ab, x));

		// Solve the right system
		Ab.set_diagonal(backup);
		x = std::move(solve_Rx(Ab, x));

		// Scale the solution back
        for (size_t i = 1; i <= n; ++i)
		{
			x[i] /= coeff[i];
		}

		return x;
	}

    void solve_gauss_factorization(const BbMatrix& A, BbVector* bx)
    {
        BbVector& b = *bx;		// Alias for the input vector

        // Check requirements
        if (!A.is_square())
        {
            throw length_error("Error in solve_gauss_factorization() - Matrix A isn't square!!\n");
        }

        if (A.n_rows() != b.size())
        {
            throw length_error("Error in solve_gauss_factorization() - Sizes don't match!!\n");
        }

        size_t n = b.size();
        BbMatrix Ab(A);			// Work matrix that will be factorized
        BbVector coeff(n);		// Vector of the scaling factor for the unknowns

        // Bring the system in its standard form
        // Augment the matrix with the vector of constant terms
        Ab.push_back_column(b);

        // Normalize rows to 1
        for (size_t i = 1; i <= n; ++i)
        {
            Ab.normalize_row(i);
        }

        // Normalize columns and save the scaling factor
        for (size_t j = 1; j <= n; ++j)
        {
            double denom = abs(*(std::max_element(Ab.col_begin(j), Ab.col_end(j), [](const double& a, const double& b){ return abs(a) < abs(b); })));
            if (denom < 1.E-20)
            {
                throw runtime_error("Error in normalize_column(j) - Division by 0 in column normalization!!\n");
            }

            for_each(Ab.col_begin(j), Ab.col_end(j), [&denom](double& val){ val /= denom; });
            coeff[j] = denom;
        }

        // Factorize the matrix
        for (size_t k = 1; k <= n - 1; ++k)
        {
            Ab.seek_pivot_in_column(k);

            // Check for null pivot
            if (abs(Ab[k][k]) < 1.E-20)
            {
                throw runtime_error("Error in normalize_column(j) - Singular matrix!!\n");
            }

            // Factorize the rows
            for (size_t i = k + 1; i <= n; ++i)
            {
                if (Ab[i][k] == 0.)
                {
                    continue;
                }

                Ab[i][k] /= Ab[k][k];

                for (size_t j = k + 1; j <= n; ++j)
                {
                    Ab[i][j] -= Ab[i][k] * Ab[k][j];
                }
            }
        }

        // Bring the matrix back to NxN and prepare the output vector
        b = std::move(Ab.extract_column(n + 1));

        BbVector backup = Ab.get_diagonal();

        // Solve the left system
        Ab.set_diagonal(1.);
        b = std::move(solve_Lx(Ab, b));

        // Solve the right system
        Ab.set_diagonal(backup);
        b = std::move(solve_Rx(Ab, b));

        // Scale the solution back
        for (size_t i = 1; i <= n; ++i)
        {
            b[i] /= coeff[i];
        }
    }

    void BbMatrix::insert_row(size_type i, const BbVector& row)
	{
		if (i < 1 || i > nRows + 1)
		{
			throw out_of_range("Error in insert_row(i) - Index i is out of range!!\n");
		}

		if (row.size() != nColumns)
		{
			throw length_error("Error in insert_row(i) - Sizes don't match!!\n");
		}

		double** nm = matrixAlloc.allocate(nRows + 1 + 1);

		for (size_type k = 1; k < i; ++k)
		{
			nm[k] = auxAlloc.allocate(nColumns + 1);
			uninitialized_copy(myMatrix[k] + 1, myMatrix[k] + 1 + nColumns, nm[k] + 1);
		}

		nm[i] = auxAlloc.allocate(row.size() + 1);
		uninitialized_copy(row.start, row.avail, nm[i] + 1);

		for (size_type k = i; k <= nRows; ++k)
		{
			nm[k + 1] = auxAlloc.allocate(nColumns + 1);
			uninitialized_copy(myMatrix[k] + 1, myMatrix[k] + 1 + nColumns, nm[k + 1] + 1);
		}

		delete_matrix();

		myMatrix = nm;
		++nRows;
	}

	void BbMatrix::insert_column(size_type j, const BbVector& col)
	{
		if (j < 1 || j > nRows + 1)
		{
			throw out_of_range("Error in insert_column(j) - Index i is out of range!!\n");
		}

		if (col.size() != nRows)
		{
			throw length_error("Error in insert_row(i) - Sizes don't match!!\n");
		}

		double** nm = matrixAlloc.allocate(nRows + 1);

		for (size_type k = 1; k <= nRows; ++k)
		{
			nm[k] = auxAlloc.allocate(nColumns + 1 + 1);
			uninitialized_copy(myMatrix[k] + 1, myMatrix[k] + 1 + j - 1, nm[k] + 1);
			auxAlloc.construct(nm[k] + 1 + j - 1, col[k]);
			uninitialized_copy(myMatrix[k] + 1 + j - 1 , myMatrix[k] + 1 + nColumns, nm[k] + 1 + j);
		}

		delete_matrix();

		myMatrix = nm;
		++nColumns;
	}

	void BbMatrix::normalize_row(size_type i)
	{
		static const string msg = "Error in normalize_row(i)";
		range_check(i, 1, msg);

		// TODO: divide by the closest power of 2, to preserve mantissa and just act on the exponent
		double denom = this->get_row(i).get_max_abs();
		if (denom < 1.E-20)
		{
			throw runtime_error("Error in normalize_row(i) - Singular matrix!!\n");
		}

		for_each(row_begin(i), row_end(i), [&denom](double& val){ val /= denom; });
	}

	void BbMatrix::normalize_column(size_type j)
	{
		static const string msg = "Error in normalize_column(j)";
		range_check(1, j, msg);

		double denom = this->get_column(j).get_max_abs();
		if (abs(denom) < 1.E-20)
		{
			throw runtime_error("Error in normalize_column(j) - Singular matrix!!\n");
		}

		for_each(col_begin(j), col_end(j), [&denom](double& val){ val /= denom; });
	}
}
