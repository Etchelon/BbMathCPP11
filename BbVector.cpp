#include "BbVector.hpp"

using namespace std;

namespace BbMath
{
	// Scalar product of 2 BbVectors
	double product_xTy(const BbVectorD& lval, const BbVectorD& rval)
	{
		if (lval.size() == 0 || rval.size() == 0)
		{
			return 0.;
		}

		if (lval.size() != rval.size())
		{
			throw std::length_error("Sizes don't match!!\n");
		}

		return std::inner_product(lval.cbegin(), lval.cend(), rval.cbegin(), 0.);
	}

	// Element by element product of 2 BbVectors
	BbVectorD dot_product(const BbVectorD& lval, const BbVectorD& rval)
	{
		if (lval.size() != rval.size())
		{
			throw std::length_error("Sizes don't match!!\n");
		}

		size_t n = lval.size();
		BbVectorD ret(n);

		for (size_t i = 1; i <= n; ++i)
		{
			ret[i] = lval[i]*rval[i];
		}

		return ret;
	}

	// Element by element division of 2 BbVectors
	BbVectorD dot_division(const BbVectorD& lval, const BbVectorD& rval)
	{
		if (lval.size() != rval.size())
		{
			throw std::length_error("Sizes don't match!!\n");
		}

		size_t n = lval.size();
		BbVectorD ret(n);

		for (size_t i = 1; i <= n; ++i)
		{
			ret[i] = lval[i]/rval[i];
		}

		return ret;
	}
}
