#ifndef GUARD_BBVECTOR_HPP
#define GUARD_BBVECTOR_HPP

#include "ColumnIterator.hpp"
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>

namespace BbMath
{
	// Forward declarations
	class BbMatrix;

	//========================================================================================//
	//==================================== CLASS BbVector ====================================//
	//========================================================================================//

	class BbVector
	{
		friend class BbMatrix;

	public:
		// ****************** Enums
		enum class SortMode
		{
			Ascending,
			Descending
		};

	public:
		// ****************** Typedefs for the user
		typedef double		   value_type;
		typedef double&		   reference;
		typedef const double&  const_reference;
		typedef double*		   iterator;
		typedef const double*  const_iterator;
		typedef std::size_t	   size_type;
		typedef std::ptrdiff_t difference_type;

	private:
		// ****************** Data
		double* myVector;	// Points to the 1st chunk of memory allocated to the vector, which remains unused because BbVector's index is 1-based; myVector == BbVector(0), unaccessible
		double* start;		// Points to the 1st accessible element of the vector, start == BbVector(1)
		double* avail;		// Points 1 past the last allocated element of the vector (see comments to grow() function)
		double* limit;		// Points 1 past the last chunk of memory allocated to the vector (see comments to grow() function)
		size_type mySize;	// Returns the no. of elements in the vector; mySize == avail - start

		// Allocator object to dynamically manage memory
		std::allocator<double> alloc;

	private:
		// ****************** Constructor functions
		void create();
		void create(size_type n, const double& val);

		template<class Iter>
		void create(Iter b, Iter e);

		void uncreate();

		// Support functions for dynamic resizing
		void grow();
		void unchecked_append(const double& val);
		void shrink(size_type newSize);

        bool range_check(size_type i, const std::string& msg = "") const;
        bool is_vector_empty(const std::string& msg = "") const;

	public:
		// ****************** Constructors
		BbVector() { create(); }
		explicit BbVector(size_type n, const double& val = 0.) { create(n, val); }
		BbVector(size_type n, const double* arr)			   { create(arr, arr + n); }
		BbVector(const_iterator begin, const_iterator end)	   { create(begin, end); }
		BbVector(const BbVector& other)						   { create(other.cbegin(), other.cend()); }
		BbVector(BbVector&& other)
            : myVector{other.myVector}, start{other.start}, avail{other.avail}, limit{other.limit}, mySize{other.mySize}
		{ other.myVector = other.start = other.avail = other.limit = nullptr; }
		BbVector(const std::initializer_list<double>& list);
		explicit BbVector(const std::vector<double>& vec);
		BbVector(const ColumnIterator& begin, const ColumnIterator& end);

		// ****************** Destructor
		~BbVector() { uncreate(); }

		// ****************** Assignment
		BbVector& operator=(const BbVector& other);
		BbVector& operator=(BbVector&& other);
		BbVector& operator=(const double& val);

		// ****************** Dynamic allocation
		void resize(size_type newSize);
		void push_back(double val);
		size_type size() const	   { return mySize; }
		size_type capacity() const { return limit - start; }

		// ****************** Indexing - with and without range check
			  double& operator()(size_type i)		{ range_check(i); return myVector[i]; }
		const double& operator()(size_type i) const { range_check(i); return myVector[i]; }
			  double& operator[](size_type i)		{ return myVector[i]; }
		const double& operator[](size_type i) const { return myVector[i]; }

		// ****************** Return delimiting pointers
			  iterator begin()		  { return start; }
		const_iterator cbegin() const { return start; }
			  iterator end()		  { return avail; }
		const_iterator cend()   const { return avail; }

		// ****************** Mathematic functions
		double get_min() const;
		double get_max() const;
		double get_min_abs() const;
		double get_max_abs() const;
		double get_sum() const;
		double get_sum_abs() const;
		void sort(SortMode mode);

		// Norms
		double normInf() const;
		double norm1() const;
		double norm2() const;
		double normM() const;
	};

	inline
	bool BbVector::range_check(size_type i, const std::string& msg) const
	{
		if (i < 1 || i > mySize)
			throw std::out_of_range(msg + " - Index i is out of range!!\n");

		return true;
	}

	inline
	bool BbVector::is_vector_empty(const std::string& msg) const
	{
		if (mySize == 0)
			throw std::out_of_range(msg + " - Empty BbVector!!\n");

		return false;
	}

	std::ostream& operator<<(std::ostream& os, const BbVector& vec);

	// Products
	double product_xTy(const BbVector& lval, const BbVector& rval);
	BbMatrix product_xyT(const BbVector& lval, const BbVector& rval);
	BbVector dot_product(const BbVector& lval, const BbVector& rval);
	BbVector dot_division(const BbVector& lval, const BbVector& rval);
}

#endif // GUARD_BBVECTOR_HPP
