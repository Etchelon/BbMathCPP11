#ifndef GUARD_COLUMNITERATOR_HPP
#define GUARD_COLUMNITERATOR_HPP

#include <memory>

namespace BbMath
{
//========================================================================================//
//================================= CLASS ColumnIterator =================================//
//========================================================================================//
template <typename T>
class ColumnIterator
	: public std::iterator<std::random_access_iterator_tag, T> // Inherit typedefs
{
public:
	typedef std::size_t size_type;
	typedef std::ptrdiff_t difference_type;

private:
	T** mat;		  // Points to the matrix over which the iterator will operate
	size_type row;	// This is the index that is incremented/decremented during iterations
	size_type column; // Should never change as the iterator iterates over the rows of a given column

public:
	// Constructors
	ColumnIterator(T** p, const size_type& i, const size_type& j)
		: mat(p), row(i), column(j) {}
	ColumnIterator(const ColumnIterator& other)
		: mat(other.mat), row(other.row), column(other.column) {}
	ColumnIterator(ColumnIterator&& other)
		: mat(other.mat), row(other.row), column(other.column)
	{
		other.mat = nullptr;
	}
	~ColumnIterator() {}

	ColumnIterator& operator=(const ColumnIterator& rhs)
	{
		if (*this != rhs) {
			mat = rhs.mat;
			column = rhs.column;
			row = rhs.row;
		}
		return *this;
	}
	ColumnIterator& operator=(ColumnIterator&& rhs)
	{
		if (*this != rhs) {
			mat = rhs.mat;
			column = rhs.column;
			row = rhs.row;
			rhs.mat = nullptr;
		}
		return *this;
	}

	// Access operators
	T& operator*()
	{
		return mat[row][column];
	}
	//		const T& operator*() const			   { return mat[row][column]; }
	T& operator[](size_type i)
	{
		return mat[i][column];
	}
	//		const T& operator[](size_type i) const { return mat[i][column]; }

	// Increment/decrement operators
	ColumnIterator& operator++()
	{
		++row;
		return *this;
	}
	ColumnIterator& operator--()
	{
		--row;
		return *this;
	}
	ColumnIterator operator+(const size_type& n) const
	{
		return ColumnIterator(mat, row + n, column);
	}
	ColumnIterator operator-(const size_type& n) const
	{
		return ColumnIterator(mat, row - n, column);
	}
	size_type operator-(const ColumnIterator& rhs) const
	{
		return row - rhs.row;
	}
	template <typename C>
	friend ColumnIterator<C> operator+(size_type n, const ColumnIterator<C>& rhs);

	// Comparison operators
	bool operator==(const ColumnIterator& rhs) const
	{
		return mat == rhs.mat && column == rhs.column && row == rhs.row;
	}
	bool operator!=(const ColumnIterator& rhs) const
	{
		return !(*this == rhs);
	}
	bool operator<(const ColumnIterator& rhs) const;
	bool operator<=(const ColumnIterator& rhs) const;
	bool operator>(const ColumnIterator& rhs) const;
	bool operator>=(const ColumnIterator& rhs) const;
};

template <typename T>
inline ColumnIterator<T> operator+(std::size_t n, const ColumnIterator<T>& rhs)
{
	return rhs + n;
}

template <typename T>
inline bool ColumnIterator<T>::operator<(const ColumnIterator<T>& rhs) const
{
	if (mat != rhs.mat || column != rhs.column) {
		return false;
	}

	return row < rhs.row;
}

template <typename T>
inline bool ColumnIterator<T>::operator<=(const ColumnIterator<T>& rhs) const
{
	if (mat != rhs.mat || column != rhs.column) {
		return false;
	}

	return row <= rhs.row;
}

template <typename T>
inline bool ColumnIterator<T>::operator>(const ColumnIterator<T>& rhs) const
{
	if (mat != rhs.mat || column != rhs.column) {
		return false;
	}

	return row > rhs.row;
}

template <typename T>
inline bool ColumnIterator<T>::operator>=(const ColumnIterator<T>& rhs) const
{
	if (mat != rhs.mat || column != rhs.column) {
		return false;
	}

	return row >= rhs.row;
}
} // namespace BbMath

#endif // GUARD_COLUMNITERATOR_HPP
