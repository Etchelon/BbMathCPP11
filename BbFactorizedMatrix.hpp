#ifndef GUARD_BBFACTORIZEDMATRIX_HPP
#define GUARD_BBFACTORIZEDMATRIX_HPP

#include <vector>
#include <memory>
#include "BbMatrix.hpp"

class BbVector;
class BbMatrix;

class BbFactorizedMatrix
{
private:
	int size;
	BbMatrix original;
	BbMatrix factorized;

public:
	explicit BbFactorizedMatrix(const BbMatrix& m);

};

#endif // GUARD_BBFACTORIZEDMATRIX_HPP
