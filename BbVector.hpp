#ifndef GUARD_BBVECTOR_HPP
#define GUARD_BBVECTOR_HPP

#include "ColumnIterator.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace BbMath
{
	// Forward declarations
	template<typename T>
	class BbMatrix;

	//========================================================================================//
	//==================================== CLASS BbVector ====================================//
	//========================================================================================//

	template<typename T>
	class BbVector
	{
		friend class BbMatrix<T>;

	public:
		// ******************************************************************************************************** //
		// *** ENUMS ********************************************************************************************** //
		// ******************************************************************************************************** //
		enum class SortMode
		{
			Ascending,
			Descending
		};

	public:
		// ******************************************************************************************************** //
		// *** TYPEDEFS FOR THE USER ****************************************************************************** //
		// ******************************************************************************************************** //
		typedef T				value_type;
		typedef T&				reference;
		typedef const T&		const_reference;
		typedef T*				iterator;
		typedef const T*		const_iterator;
		typedef std::size_t		size_type;
		typedef std::ptrdiff_t	difference_type;

	private:
		// ******************************************************************************************************** //
		// *** DATA *********************************************************************************************** //
		// ******************************************************************************************************** //
		T* myVector;		// Points to the 1st chunk of memory allocated to the vector, which remains unused because BbVector's index is 1-based; myVector == BbVector(0), unaccessible
		T* start;			// Points to the 1st accessible element of the vector, start == BbVector(1)
		T* avail;			// Points 1 past the last allocated element of the vector (see comments to grow() function)
		T* limit;			// Points 1 past the last chunk of memory allocated to the vector (see comments to grow() function)
		size_type mySize;	// Returns the no. of elements in the vector; mySize == avail - start

		// Allocator object to dynamically manage memory
		std::allocator<T> alloc;

	private:
		// ******************************************************************************************************** //
		// *** CONSTRUCTING FUNCTIONS ***************************************************************************** //
		// ******************************************************************************************************** //
		// Private function for memory management: create and empty BbVector
		void create()
		{
			mySize = 0;
			myVector = nullptr;
			limit = avail = start = nullptr;
		}

		// Private function for memory management: create a BbVector by allocating n+1 doubles and setting their value to val
		void create(size_type n, const T& val)
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

			std::uninitialized_fill(start, avail, val);
		}

		// Private function for memory management: create a BbVector by copying another BbVector using iterators (similar to one of std::vector constructors)
		template<typename Iter>
		void create(Iter b, Iter e)
		{
			if (b == e)
			{
				create();
				return;
			}

			mySize = e - b;
			myVector = alloc.allocate(mySize + 1);
			start = myVector + 1;
			limit = avail = start + mySize;

			std::uninitialized_copy(b, e, start);
		}

		// Private function for memory management: delete the current BbVector by destroying all doubles allocated into myVector and then deallocating such memory
		void uncreate()
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
		// *** SUPPORT FUNCTION FOR DYNAMIC RESIZING ************************************************************** //
		// ******************************************************************************************************** //
		// Private function for memory management
		// The vector is delimited by 2 pointers: avail points to 1 past the last allocated element of the vector, limit points to 1 past the maximum size of the vector
		// C# equivalent in (some) collection classes are Count and Capacity (if I remember correctly)
		// The point is: when I call push_back to append an element, I don't want the resizing happening every time because dynamic allocation takes resources, so the first time
		// I do it I want to allocate twice as much space than the original size, so for ex if mySize == 5, I can append 5 elements while resizing the vector only once.
		// After I call push_back() the 1st time, mySize == 6 == (avail - myVector), but the space allocated to the vector is 10 == (limit - myVector)
		// (in all cases, add 1 unused element to that BbVector's index starts from 1 instead of 0)
		void grow()
		{
			size_type newCapacity = std::max(static_cast<std::ptrdiff_t>(2 * (limit - myVector - 1)), std::ptrdiff_t(1));

			T* newVector = alloc.allocate(newCapacity + 1);
			T* newStart  = newVector + 1;
			T* newAvail  = std::uninitialized_copy(start, avail, newStart);
			uncreate();

			myVector = newVector;
			start    = newStart;
			avail    = newAvail;
			limit    = start + newCapacity;
			mySize   = avail - start;
		}

		// Private function for memory management: called when avail < limit and I have some spare allocated memory where to append a new element
		void unchecked_append(const T& val)
		{
			alloc.construct(avail, val);
			++avail;
			++mySize;
		}

		// Private function for memory management: shrink the size of BbVector, by copying all elements up to newSize into a temporary double*, deleting and freeing myVector and assigning the temp to myVector
		void shrink(size_type newSize)
		{
			T* nv = alloc.allocate(newSize + 1);
			T* ns = nv + 1;
			T* na = std::uninitialized_copy(start, start + newSize, ns);
			uncreate();

			myVector = nv;
			start    = ns;
			limit    = avail = na;
			mySize   = newSize;
		}

		bool range_check(size_type i, const std::string& msg = "Error") const
		{
			if (i < 1 || i > mySize)
				throw std::out_of_range(msg + " - Index i is out of range!!\n");

			return true;
		}

		bool is_vector_empty(const std::string& msg) const
		{
			if (mySize == 0)
				throw std::out_of_range(msg + " - Empty BbVector!!\n");

			return false;
		}

	public:
		// ******************************************************************************************************** //
		// *** CONSTRUCTORS *************************************************************************************** //
		// ******************************************************************************************************** //
		BbVector()
		{
			create();
		}

		explicit BbVector(size_type n, const T& val = static_cast<T>(0))
		{
			create(n, val);
		}

		BbVector(size_type n, const T* arr)
		{
			create(arr, arr + n);
		}

		BbVector(const_iterator begin, const_iterator end)
		{
			create(begin, end);
		}

		BbVector(const BbVector& other)
		{
			create(other.cbegin(), other.cend());
		}

		BbVector(BbVector&& other)
			: myVector{other.myVector}, start{other.start}, avail{other.avail}, limit{other.limit}, mySize{other.mySize}
		{
			other.myVector = other.start = other.avail = other.limit = nullptr;
		}

		BbVector(const std::initializer_list<T>& list)
		{
			create(list.begin(), list.end());
		}

		explicit BbVector(const std::vector<T>& vec)
		{
			create(vec.cbegin(), vec.cend());
		}

		BbVector(const ColumnIterator<T> &begin, const ColumnIterator<T> &end)
		{
			create(begin, end);
		}

		~BbVector()
		{
			uncreate();
		}

		// ******************************************************************************************************** //
		// *** ASSIGNMENT OPERATIONS ****************************************************************************** //
		// ******************************************************************************************************** //
		// Assignment operator: copy rhs into this, by deleting the current BbVector using uncreate() and copying the rhs into BbVector using create()
		BbVector& operator=(const BbVector<T>& rhs)
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
		BbVector& operator=(BbVector<T>&& rhs)
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
		BbVector& operator=(const T& val)
		{
			if (mySize)
			{
				std::fill(start, avail, val);
			}

			return *this;
		}

		// ******************************************************************************************************** //
		// *** DYNAMIC ALLOCATION ********************************************************************************* //
		// ******************************************************************************************************** //
		// Resize the vector to new size
		// If newsize <= 0, delete the vector;
		// If newsize > mySize, create an empty, bigger vector
		void resize(size_type newSize)
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
				create(newSize, static_cast<T>(0));
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
				create(newSize, static_cast<T>(0));
			}
		}

		// Like std::vector.push_back(): add 1 element to the end of the current BbVector
		void push_back(T val)
		{
			if (avail == limit)
			{
				grow();
			}

			unchecked_append(val);
		}
		size_type size() const
		{
			return mySize;
		}

		size_type capacity() const
		{
			return limit - start;
		}

		// ******************************************************************************************************** //
		// *** INDEXING AND ACCESS ******************************************************************************** //
		// ******************************************************************************************************** //
			  T& operator()(size_type i)	   { range_check(i); return myVector[i]; }
		const T& operator()(size_type i) const { range_check(i); return myVector[i]; }
			  T& operator[](size_type i)	   { return myVector[i]; }
		const T& operator[](size_type i) const { return myVector[i]; }

		// ****************** Return delimiting pointers
			  iterator begin()		  { return start; }
		const_iterator cbegin() const { return start; }
			  iterator end()		  { return avail; }
		const_iterator cend()   const { return avail; }

		// ******************************************************************************************************** //
		// *** MATHEMATIC FUNCTIONS ******************************************************************************* //
		// ******************************************************************************************************** //
		// Return the minimum value of the vector
		T get_min() const
		{
			static const std::string msg = "Error in get_min()";
			is_vector_empty(msg);

			return *(std::min_element(start, avail));
		}

		// Return the minimum absolute value of the vector
		T get_min_abs() const
		{
			static const std::string msg = "Error in get_min_abs()";
			is_vector_empty(msg);

			return std::abs(*(std::min_element(start, avail, [](const T& a, const T& b){ return std::abs(a) < std::abs(b); })));
		}

		// Return the maximum value of the vector
		T get_max() const
		{
			static const std::string msg = "Error in get_max()";
			is_vector_empty(msg);

			return *(std::max_element(start, avail));
		}

		// Return the maximum absolute value of the vector
		T get_max_abs() const
		{
			static const std::string msg = "Error in get_max_abs()";
			is_vector_empty(msg);

			return std::abs(*(std::max_element(start, avail, [](const T& a, const T& b){ return std::abs(a) < std::abs(b); })));
		}

		T get_sum() const
		{
			static const std::string msg = "Error in get_sum()";
			is_vector_empty(msg);

			return std::accumulate(start, avail, static_cast<T>(0));
		}

		T get_sum_abs() const
		{
			static const std::string msg = "Error in get_sum_abs()";
			is_vector_empty(msg);

			return std::accumulate(start, avail, static_cast<T>(0), [](const T& a, const T& b){ return a + abs(b); });
		}

		// Sort the vector in ascending order
		void sort(SortMode mode)
		{
			static const std::string msg = "Error in sort()";
			is_vector_empty(msg);

			std::sort(start, avail, [&mode](const T& a, const T& b){ return mode == SortMode::Ascending ? a < b : a > b; });
		}

		// Infinite norm: the max value of the vector
		T normInf() const
		{
			static const std::string msg = "Error in normInf()";
			is_vector_empty(msg);

			return get_max_abs();
		}

		// Norm-1: the sum of the absolute values of the vector
		T norm1() const
		{
			static const std::string msg = "Error in norm1()";
			is_vector_empty(msg);

			return std::accumulate(start, avail, static_cast<T>(0), [](T sum, T elem){ return sum + std::abs(elem); });
		}

		// Norm-2: the square root of the SSQ of the vector
		T norm2() const
		{
			static const std::string msg = "Error in norm2()";
			is_vector_empty(msg);

			return std::sqrt(std::accumulate(start, avail, static_cast<T>(0), [](T sum, T elem){ return sum + elem*elem; }));
		}

		// Norm-M: Norm-1 divided by the no. of elements of the vector
		T normM() const
		{
			static const std::string msg = "Error in normM()";
			is_vector_empty(msg);

			return mySize ? this->norm1() / mySize : static_cast<T>(0);
		}
	};

	// Print the BbVector to the output stream
	template<typename T>
	std::ostream& operator<<(std::ostream& os, const BbVector<T>& vec)
	{
		if (vec.size() == 0)
		{
			os << std::endl << "Empty BbVector!!" << std::endl;
		}
		else
		{
			std::for_each(vec.cbegin(), vec.cend(), [&os](const T& val)
			{
				os << val << " ";
			});
			os << std::endl;
		}

		return os;
	}

	using BbVectorD = BbVector<double>;
	using BbVectorI = BbVector<int>;

	// ******************************************************************************************************** //
	// *** VECTOR/MATRIX PRODUCTS ***************************************************************************** //
	// ******************************************************************************************************** //
	double product_xTy(const BbVectorD& lval, const BbVectorD& rval);
	BbVectorD dot_product(const BbVectorD& lval, const BbVectorD& rval);
	BbVectorD dot_division(const BbVectorD& lval, const BbVectorD& rval);
}

#endif // GUARD_BBVECTOR_HPP
