#include "BbVector.hpp"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <stdexcept>

using namespace std;

namespace BbMath
{
	// ******************************************************************************************************** //
	// *** PRIVATE CONSTRUCTING FUNCTIONS ********************************************************************* //
	// ******************************************************************************************************** //

	// Private function for memory management: create and empty BbVector
	void BbVector::create()
	{
		mySize = 0;
		myVector = nullptr;
		limit = avail = start = nullptr;
	}

	// Private function for memory management: create a BbVector by allocating n+1 doubles and setting their value to val
	void BbVector::create(size_type n, const double& val)
	{
		if (n == 0)
		{
			create();
			return;
		}

		mySize = n;
		myVector = alloc.allocate(n + 1);
		start = myVector + 1;
		limit = avail = start + mySize;

		uninitialized_fill(start, avail, val);
	}

	// Private function for memory management: delete the current BbVector by destroying all doubles allocated into myVector and then deallocating such memory
	void BbVector::uncreate()
	{
		if (myVector)
		{
			iterator it = avail;

			while (it != myVector)
			{
				alloc.destroy(--it);
			}

			alloc.deallocate(myVector, limit - myVector);
		}

		mySize = 0;
		myVector = nullptr;
		limit = avail = start = nullptr;
	}

	// ******************************************************************************************************** //
	// ***************************************************************** END PRIVATE CONSTRUCTING FUNCTIONS *** //
	// ******************************************************************************************************** //

	// ******************************************************************************************************** //
	// *** PUBLIC CONSTRUCTORS ******************************************************************************** //
	// ******************************************************************************************************** //
	BbVector::BbVector(const initializer_list<double>& list)
	{
		create(list.begin(), list.end());
	}

	BbVector::BbVector(const vector<double>& vec)
	{
		create(vec.cbegin(), vec.cend());
	}

	BbVector::BbVector(const ColumnIterator &begin, const ColumnIterator &end)
	{
		create(begin, end);
	}

	// ******************************************************************************************************** //
	// **************************************************************************** END PUBLIC CONSTRUCTORS *** //
	// ******************************************************************************************************** //

	// ******************************************************************************************************** //
	// *** ASSIGNMENT OPERATIONS ****************************************************************************** //
	// ******************************************************************************************************** //

	// Assignment operator: copy rhs into this, by deleting the current BbVector using uncreate() and copying the rhs into BbVector using create()
	BbVector& BbVector::operator=(const BbVector& rhs)
	{
		// Check for self assignment
		if (&rhs != this)
		{
			// Free the array in left hand side
			uncreate();
			create(rhs.cbegin(), rhs.cend());
		}

		return *this;
	}

	// Assignment operator: move rhs to this
	BbVector& BbVector::operator=(BbVector&& rhs)
	{
		if (&rhs != this)
		{
			// Free the array in left hand side
			uncreate();

			// Take ownership of the rhs
			myVector = rhs.myVector;
			start	 = rhs.start;
			avail	 = rhs.avail;
			limit	 = rhs.limit;
			mySize	 = rhs.mySize;

			// Nullify the rhs
			rhs.myVector = nullptr;
			rhs.start	 = nullptr;
			rhs.avail	 = nullptr;
			rhs.limit	 = nullptr;
			rhs.mySize	 = 0;
		}

		return *this;
	}

	// Assignment operator: set all elements of the BbVector to value val
	BbVector& BbVector::operator=(const double& val)
	{
		if (mySize)
		{
			fill(start, avail, val);
		}

		return *this;
	}

	// Like std::vector.push_back(): add 1 element to the end of the current BbVector
	void BbVector::push_back(double val)
	{
		if (avail == limit)
		{
			grow();
		}

		unchecked_append(val);
	}

	// ******************************************************************************************************** //
	// ************************************************************************** END ASSIGNMENT OPERATIONS *** //
	// ******************************************************************************************************** //

	// Private function for memory management: shrink the size of BbVector, by copying all elements up to newSize into a temporary double*, deleting and freeing myVector and assigning the temp to myVector
	void BbVector::shrink(size_type newSize)
	{
		double* nv = alloc.allocate(newSize + 1);
		double* ns = nv + 1;
		double* na = uninitialized_copy(start, start + newSize, ns);
		uncreate();

		myVector = nv;
		start    = ns;
		limit    = avail = na;
		mySize   = newSize;
	}

	// Private function for memory management
	// The vector is delimited by 2 pointers: avail points to 1 past the last allocated element of the vector, limit points to 1 past the maximum size of the vector
	// C# equivalent in (some) collection classes are Count and Capacity (if I remember correctly)
	// The point is: when I call push_back to append an element, I don't want the resizing happening every time because dynamic allocation takes resources, so the first time
	// I do it I want to allocate twice as much space than the original size, so for ex if mySize == 5, I can append 5 elements while resizing the vector only once.
	// After I call push_back() the 1st time, mySize == 6 == (avail - myVector), but the space allocated to the vector is 10 == (limit - myVector)
	// (in all cases, add 1 unused element to that BbVector's index starts from 1 instead of 0)
	void BbVector::grow()
	{
		size_type newCapacity = max((ptrdiff_t)(1.5*(limit - myVector - 1)), ptrdiff_t(1));

		double* newVector = alloc.allocate(newCapacity + 1);
		double* newStart  = newVector + 1;
		double* newAvail  = uninitialized_copy(start, avail, newStart);
		uncreate();

		myVector = newVector;
		start    = newStart;
		avail    = newAvail;
		limit    = start + newCapacity;
		mySize   = avail - start;
	}

	// Private function for memory management: called when avail < limit and I have some spare allocated memory where to append a new element
	void BbVector::unchecked_append(const double& val)
	{
		alloc.construct(avail, val);
		++avail;
		++mySize;
	}

	// Print the BbVector to the console
	ostream& operator<<(ostream& os, const BbVector& vec)
	{
		if (vec.size() == 0)
		{
			os << endl << "Empty BbVector!!" << endl;
		}
		else
		{
			for_each(vec.cbegin(), vec.cend(), [&os](const double& val)
			{
				os << val << " ";
			});
			os << endl;
		}

		return os;
	}

	// Infinite norm: the max value of the vector
	double BbVector::normInf() const
	{
		static const string msg = "Error in normInf()";
		is_vector_empty(msg);

		return get_max_abs();
	}

	// Norm-1: the sum of the absolute values of the vector
	double BbVector::norm1() const
	{
		static const string msg = "Error in norm1()";
		is_vector_empty(msg);

		return std::accumulate(start, avail, 0., [](double sum, double elem){ return sum + abs(elem); });
	}

	// Norm-2: the square root of the SSQ of the vector
	double BbVector::norm2() const
	{
		static const string msg = "Error in norm2()";
		is_vector_empty(msg);

		return sqrt(std::accumulate(start, avail, 0., [](double sum, double elem){ return sum + elem*elem; }));
	}

	// Norm-M: Norm-1 divided by the no. of elements of the vector
	double BbVector::normM() const
	{
		static const string msg = "Error in normM()";
		is_vector_empty(msg);

		return mySize ? this->norm1() / mySize : 0.;
	}

	// Return the minimum value of the vector
	double BbVector::get_min() const
	{
		static const string msg = "Error in get_min()";
		is_vector_empty(msg);

		return *(std::min_element(start, avail));
	}

	// Return the minimum absolute value of the vector
	double BbVector::get_min_abs() const
	{
		static const string msg = "Error in get_min_abs()";
		is_vector_empty(msg);

		return abs(*(std::min_element(start, avail, [](const double& a, const double& b){ return abs(a) < abs(b); })));
	}

	// Return the maximum value of the vector
	double BbVector::get_max() const
	{
		static const string msg = "Error in get_max()";
		is_vector_empty(msg);

		return *(std::max_element(start, avail));
	}

	// Return the maximum absolute value of the vector
	double BbVector::get_max_abs() const
	{
		static const string msg = "Error in get_max_abs()";
		is_vector_empty(msg);

		return abs(*(std::max_element(start, avail, [](const double& a, const double& b){ return abs(a) < abs(b); })));
	}

	double BbVector::get_sum() const
	{
		static const string msg = "Error in get_sum()";
		is_vector_empty(msg);

		return accumulate(start, avail, 0.);
	}

	double BbVector::get_sum_abs() const
	{
		static const string msg = "Error in get_sum_abs()";
		is_vector_empty(msg);

		return accumulate(start, avail, 0., [](const double& a, const double& b){ return a + abs(b); });
	}

	// Sort the vector in ascending order
	void BbVector::sort(SortMode mode)
	{
		static const string msg = "Error in sort()";
		is_vector_empty(msg);

		std::sort(start, avail, [&mode](const double& a, const double& b){ return mode == SortMode::Ascending ? a < b : a > b; });
	}

	// Resize the vector to new size
	// If newsize <= 0, delete the vector;
	// If newsize > mySize, create an empty, bigger vector
	void BbVector::resize(size_type newSize)
	{
		if (mySize == newSize)
		{
			return;
		}

		if (newSize <= 0)
		{
			uncreate();
			return;
		}

		if (mySize == 0)
		{
			create(newSize, 0.);
			return;
		}

		// Shrink
		if (newSize < mySize)
		{
			shrink(newSize);
		}
		// Expand
		else if (newSize > mySize)
		{
			uncreate();
			create(newSize, 0.);
		}
	}

	// Scalar product of 2 BbVectors
	double product_xTy(const BbVector& lval, const BbVector& rval)
	{
		if (lval.size() == 0 || rval.size() == 0)
		{
			return 0.;
		}

		if (lval.size() != rval.size())
		{
			throw length_error("Sizes don't match!!\n");
		}

		return std::inner_product(lval.cbegin(), lval.cend(), rval.cbegin(), 0.);
	}

	BbVector dot_product(const BbVector& lval, const BbVector& rval)
	{
		if (lval.size() != rval.size())
		{
			throw length_error("Sizes don't match!!\n");
		}

		size_t n = lval.size();
		BbVector ret(n);

		for (size_t i = 1; i <= n; ++i)
		{
			ret[i] = lval[i]*rval[i];
		}

		return ret;
	}

	BbVector dot_division(const BbVector& lval, const BbVector& rval)
	{
		if (lval.size() != rval.size())
		{
			throw length_error("Sizes don't match!!\n");
		}

		size_t n = lval.size();
		BbVector ret(n);

		for (size_t i = 1; i <= n; ++i)
		{
			ret[i] = lval[i]/rval[i];
		}

		return ret;
	}
}
