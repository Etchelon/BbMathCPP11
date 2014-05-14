#include "BbVector.hpp"
#include "BbMatrix.hpp"

using namespace std;

namespace BbMath
{
	// Solve linear systems
	BbVectorD solve_Dx(const BbMatrixD& D, const BbVectorD& b)
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
		BbVectorD x(n);

		for (size_t i = 1; i <= n; ++i)
		{
			x[i] = b[i] / D[i][i];
		}

		return x;
	}

	BbVectorD solve_Rx(const BbMatrixD& R, const BbVectorD& b)
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
		BbVectorD x(b);

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

	BbVectorD solve_Lx(const BbMatrixD& L, const BbVectorD& b)
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
		BbVectorD x(b);

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

	BbVectorD solve_gauss_elimination(const BbMatrixD& A, const BbVectorD& b)
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

		BbMatrixD Ab(A);
		Ab.push_back_column(b);

		for (size_t k = 1; k <= r - 1; ++k)
		{
			// TODO: seek pivot

			for (size_t i = k + 1; i <= r; ++i)
			{
				if (abs(Ab[i][k]) < 1.E-20)
					continue;

				double l = Ab[i][k] / Ab[k][k];

				for (size_t j = 1; j <= c; ++j)
				{
					Ab[i][j] -= l * Ab[k][j];
				}
			}
		}

		BbVectorD bb(Ab.get_column(c));
		Ab.chop_column();

		return solve_Rx(Ab, bb);
	}

	BbVectorD solve_gauss_factorization(const BbMatrixD& A, const BbVectorD& b)
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
		BbMatrixD Ab(A);			// Work matrix that will be factorized
		BbVectorD coeff(n);		// Vector of the scaling factor for the unknowns

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
			double denom = abs(*(max_element(Ab.col_begin(j), Ab.col_end(j), [](const double& a, const double& b){ return abs(a) < abs(b); })));
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
				if (abs(Ab[i][k]) < 1.E-20)
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
		BbVectorD x = Ab.extract_column(n + 1);

		BbVectorD backup = Ab.get_diagonal();

		// Solve the left system
		Ab.set_diagonal(1.);
		x = move(solve_Lx(Ab, x));

		// Solve the right system
		Ab.set_diagonal(backup);
		x = move(solve_Rx(Ab, x));

		// Scale the solution back
		for (size_t i = 1; i <= n; ++i)
		{
			x[i] /= coeff[i];
		}

		return x;
	}

	void solve_gauss_factorization(const BbMatrixD& A, BbVectorD* bx)
	{
		BbVectorD& b = *bx;		// Alias for the input vector

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
		BbMatrixD Ab(A);			// Work matrix that will be factorized
		BbVectorD coeff(n);			// Vector of the scaling factor for the unknowns

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
			double denom = abs(*(max_element(Ab.col_begin(j), Ab.col_end(j), [](const double& a, const double& b){ return abs(a) < abs(b); })));
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
				if (abs(Ab[i][k]) < 1.E-20)
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
		b = move(Ab.extract_column(n + 1));

		BbVectorD backup = Ab.get_diagonal();

		// Solve the left system
		Ab.set_diagonal(1.);
		b = move(solve_Lx(Ab, b));

		// Solve the right system
		Ab.set_diagonal(backup);
		b = move(solve_Rx(Ab, b));

		// Scale the solution back
		for (size_t i = 1; i <= n; ++i)
		{
			b[i] /= coeff[i];
		}
	}

	BbMatrixD inverse(const BbMatrixD& other)
	{
		// Temp: only square matrices
		if (other.n_rows() != other.n_columns())
			throw length_error("Only square matrices allowed for the moment in inverse()!\n");

		int n = other.n_rows();
		BbMatrixD ret(n, n);

		for (int i = 1; i <= n; ++i)
		{
			// Create the identity vector of solutions
			BbVectorD b(n, 0.);
			b[i] = 1;

			// Solve the system and get the i-th column of the inverse matrix
			BbVectorD x = solve_gauss_factorization(other, b);
			ret.set_column(i, x);
		}

		return ret;
	}
}
