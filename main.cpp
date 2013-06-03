#include "BbVector.hpp"
#include "BbMatrix.hpp"
#include <iostream>
#include <vector>
#include <chrono>
#include <string>
#include <algorithm>
#include <random>
#include <complex>

using namespace std;
using namespace BbMath;

#include "BbStopWatch.hpp"

BbVector create_BbVector()
{
	return BbVector(6, 2.);
}

int main(int argc, char* argv[])
{
	// Stopwatch used for the tests
	BbStopWatch sw(BbStopWatch::TimeUOM::milliseconds);

/*
#if __cplusplus >= 201103L
	{
		vector<double> vec{ 1., 2., 3., -55. };
		cout << BbVector(vec).get_sum_abs() << endl;

		BbVector bbv{ 1., 5. };
		cout << bbv << endl;

		BbVector a{ 2., 3., 5. };
		cout << "BbVector a, constructed with an initializer_list: " << endl << a << endl;

		BbMatrix A
		{
			{ 1., 3., 4. }, { 2., 7., 9. }, { -11., 45., 67. }
		};
		cout << "BbMatrix A, constructed from initializer_lists: " << endl << A << endl;
	}
#endif

	{
		BbVector a;
		a.push_back(-1.);
		a.push_back(1.);
		a.push_back(0.);
		a.push_back(-666.);
		a.push_back(1.E+4);
		a.push_back(.1111);

		cout << "a = " << a << endl;
		cout << "a.get_min() = " << a.get_min() << endl;
		cout << "a.get_min_abs() = " << a.get_min_abs() << endl;
		cout << "a.get_max() = " << a.get_max() << endl;
		cout << "a.get_max_abs() = " << a.get_max_abs() << endl;

		a.sort(BbVector::SortMode::Descending);
		cout << "Sorted a = " << a << endl;

		cout << "Norm1 = " << a.norm1() << endl;
		cout << "a * a = " << product_xTy(a, a) << endl;

		BbVector aa, bb = create_BbVector();
		cout << "aa = " << aa << endl;
		cout << "bb = " << bb << endl;
		aa = std::move(bb);
		bb.~BbVector();
		cout << "aa = " << aa << endl;
		cout << "bb = " << bb << endl;

		int aaa = 2;
		int bbb = 3;
		int ccc = 5;
		ccc = std::move(aaa);
		cout << endl;
	}

	{
		BbMatrix A(10, 10, 666.);
		int i;

		sw.restart();
		BbVector c;
		for (i = 0; i < 10000000; ++i)
		{
			c = A.get_column(6);
			c[6] = 667.;
		}

		cout << c << endl;
		cout << "get_column() called " << i << " times, took " << sw.split() << " ms" << endl;

		sw.restart();
		BbVector r;
		for (i = 0; i < 10000000; ++i)
		{
			r = A.get_row(6);
			r[6] = 667.;
		}

		cout << r << endl;
		cout << "get_row() called " << i << " times, took " << sw.split() << " ms" << endl;
	}

	{
		double list[5] = { 1., -5., 4., 53., -0.66 };
		vector<double> vec(5000000, 6.);

		cout << "Initialization:" << endl;
		BbVector bbV(5, list);

		BbVector bbA(vec);
		cout << bbV << endl;

		cout << "Norms" << endl;

		cout << "Assignment:" << endl;
		bbV  = 666.;

		bbA  = bbV;

		cout << bbV << endl;
		cout << bbA << endl;

		BbVector aaa(6, 5.);
		aaa = 5.5;
		cout << "aaa = " << aaa << endl;
		aaa.resize(2);
		cout << "aaa = " << aaa << endl;

		BbMatrix B(5, 7), C(3, 4, 2.);
		B = 666.;

		cout << "Matrix B =" << endl << B;
		cout << "Matrix C =" << endl << C;

		C = B;

		cout << "Matrix C =" << endl << C;

		BbMatrix AA(3, 3, 1.);
		AA(1,1) = 1.;
		AA(1,2) = 2.;
		AA(1,3) = 3.;
		BbVector bb(3);
		bb(1) = 4.;
		bb(2) = 3.;
		bb(3) = 2.;

		cout << "AA*bb = " << product_Ax(AA, bb) << endl;
		cout << product_xyT(bb, bb) << endl;
		cout << product_xTy(bb, bb) << endl;
		cout << transpose(product_AB(AA, AA)) << endl;

		AA.resize(2, 2);
		cout << AA << endl;
		cout << transpose(AA) << endl;
	}

	{
		sw.restart();
		BbMatrix uff(10000, 10000, 4.);
		cout << "I have initialized a 10000x10000 BbMatrix in " << sw.split() << " ms!" << endl;

		transpose(uff);
		cout << "I have transposed a 10000x10000 BbMatrix in " << sw.split() << " ms!" << endl;
	}
*/
	{
		BbMatrix R(3, 3);
		BbVector b(3);
		R(1,1) = 1.;
		R(1,2) = 2.;
		R(1,3) = 3.;
		R(2,1) = -1.;
		R(2,2) = 0.;
		R(2,3) = 1.;
		R(3,1) = 3.;
		R(3,2) = 1.;
		R(3,3) = 3.;
		b(1) = -5.;
		b(2) = -3.;
		b(3) = -3.;

//		BbMatrix R
//		{
//			{ 10., 1.E+20 },
//			{ 1., -.9 }
//		};
//		BbVector b{ 1.E+20, .1 };

//		BbMatrix R
//		{
//			{ 10., 1., 0. },
//			{ .5, 1.E-20, 1. },
//			{ 0., 0., 1. }
//		};
//		BbVector b{ 1.E+20, 1., 1. };

//		BbMatrix R
//		{
//			{ -1., 1., -1. },
//			{ 10., -5., 10. },
//			{ 2., 1., 1. }
//		};
//		BbVector b{ 1., -3., 2. };

		cout << "Gauss elimination:" << endl;

		sw.start();
		for (int i = 1; i <= 1000000; ++i)
			solve_gauss_elimination(R, b);

		cout << "Solved a 3x3 linear system by Gauss Elimination for 1 million times. It took " << sw.split() << " ms" << endl;
		cout << solve_gauss_elimination(R, b) << endl;

		cout << "Gauss factorization:" << endl;

		sw.restart();
		for (int i = 1; i <= 1000000; ++i)
			solve_gauss_factorization(R, b);

		cout << "Solved a 3x3 linear system by Gauss Factorization for 1 million times. It took " << sw.split() << " ms" << endl;
		cout << solve_gauss_factorization(R, b) << endl;

		cout << "Gauss factorization with overlapping:" << endl;

		sw.restart();
		for (int i = 1; i <= 1000000; ++i)
		{
			BbVector b{ -5., -3., -3. };
			solve_gauss_factorization(R, &b);

			if (i == 1000000)
			{
				cout << "Solved a 3x3 linear system by Gauss Factorization for 1 million times. It took " << sw.split() << " ms" << endl;
				cout << b << endl;
			}
		}
	}
/*
	{
		cout << "Push back row: " << endl;

		BbMatrix R(3, 3);
		R(1,1) = 1.;
		R(1,2) = 2.;
		R(1,3) = 3.;
		R(2,1) = -1.;
		R(2,2) = 0.;
		R(2,3) = 1.;
		R(3,1) = 3.;
		R(3,2) = 1.;
		R(3,3) = 3.;

		sw.restart();
		for (int i = 1; i <= 5000000; ++i)
		{
			R.push_back_row(R.get_row(1));
			R.chop_row();
		}
		cout << "Took " << sw.split() << " ms" << endl;

		cout << "Push back column: " << endl;

		sw.restart();
		for (int i = 1; i <= 5000000; ++i)
		{
			R.push_back_column(R.get_row(1));
			R.chop_column();
		}
		cout << "Took " << sw.split() << " ms" << endl;
	}

	{
		cout << "Insert row:" << endl;
		BbMatrix AAA(5, 3, 2.);
		BbVector bbb(3, 3.);
		cout << "AAA before row insertion: " << endl << AAA << endl;
		AAA.insert_row(3, bbb);
		AAA.push_back_row(bbb);
		cout << "AAA after row insertion: " << endl << AAA << endl;

		cout << "Insert column:" << endl;
		BbVector bbbb(7, 4.);
		cout << "AAA before column insertion: " << endl << AAA << endl;
		AAA.insert_column(2, bbbb);
		AAA.push_back_column(bbbb);
		cout << "AAA after column insertion: " << endl << AAA << endl;

		cout << "Delete row:" << endl;
		cout << "AAA before row deletion: " << endl << AAA << endl;
		AAA.delete_row(1);
		AAA.chop_row();
		cout << "AAA after row deletion: " << endl << AAA << endl;

		cout << "Delete column:" << endl;
		cout << "AAA before column deletion: " << endl << AAA << endl;
		AAA.delete_column(2);
		AAA.chop_column();
		cout << "AAA after column deletion: " << endl << AAA << endl;

		cout << "Swap rows:" << endl;

		sw.restart();
		for (int i = 0; i < 100000000; ++i)
		{
			AAA.swap_rows(2, 3);
		}
		cout << "AAA after row swap: " << endl << AAA << endl << "Took " << sw.split() << " ms" << endl;
#if __cplusplus >= 201103L
		AAA.push_back_column({ 666., 666., 666., 666., 666. });
		cout << "AAA after addition of a column: " << endl << AAA << endl;

		cout << "Swap columns:" << endl;
		sw.restart();
		for (int i = 0; i < 100000000; ++i)
		{
			AAA.swap_columns(1, 4);
		}
		cout << "AAA after column swap: " << endl << AAA << endl << "Took " << sw.split() << " ms" << endl;
#endif
	}

	{
		BbMatrix A(10, 10, 666.);

		sw.restart();
		for (int i = 1; i <= 1000000; ++i)
		{
			A.transpose();
//			if (i < 5)
//				cout << A << endl;
		}
		cout << "Transposed a 10x10 matrix for 1M times in " << sw.split() << " ms" << endl;
	}
*//*
	{
		sw.restart();
		for (int i = 1; i <= 10000000; ++i)
		{
			BbMatrix A(10, 10, 666.);
			A[6][6] = 667.;
		}

		cout << "Created and deleted a 10x10 matrix for 10M times in " << sw.split() << " ms" << endl;
	}

	{
		BbMatrix* A;

		sw.restart();
		for (int i = 1; i <= 10000000; ++i)
		{
			A = new BbMatrix(10, 10, 666.);
			(*A)[6][6] = 667.;
			delete A;
		}

		cout << "Created and deleted a 10x10 matrix for 10M times in " << sw.split() << " ms" << endl;
	}
*//*
	{
#if __cplusplus >= 201103L
		// Test iterators
		std::default_random_engine generator;
		std::uniform_int_distribution<int> distribution(-666, 666);
		auto rnd = std::bind(distribution, generator);
		BbMatrix A{ { rnd(), rnd() }, { rnd(), rnd() } };

		cout << "A = " << endl << A << endl;
		for (ColumnIterator it = A.cbegin(2); it != A.cend(2); ++it)
		{
			cout << *it << " ";
		}
		cout << endl;

		for_each(A.cbegin(1), A.cend(1), [](double a){ cout << a << " "; }); cout << endl;
		for_each(A.rbegin(1), A.rend(1), [](double a){ cout << a << " "; }); cout << endl;
		cout << "A[2][2] = " << *(A.cbegin(2) + 1) << endl;
#endif
		BbMatrix B(10, 10, 666.);
		cout << "B = " << endl << B << endl;

		BbMatrix::RowIterator it = const_cast<BbMatrix::RowIterator>(B.crbegin(5));
		cout << "B[1][5] = " << *it << endl;
		*it = 15.;
		cout << "B[1][5] = " << *it << endl;

		BbMatrix::RowIterator it2 = it + 4;
		*it2 = 55.;
		cout << "B[5][5] = " << *it2 << endl;

		BbMatrix::RowIterator it3 = it2;
		cout << "B[6][5] = " << *++it3 << endl;

		cout << "B = " << endl << B << endl;

		sw.restart();
		for (int i = 1; i < 1.E+6; ++i)
		{
			B.get_column(5).get_sum();
		}
		cout << "Called get_column(5) for 1M times: " << B.get_column(5) << endl << "Took " << sw.split() << " ms" << endl;

		cout << "B before row normalization: " << endl << B << endl;
		sw.restart();
		for (int i = 1; i < 1E+7; ++i)
			B.normalize_row(5);
		cout << "Took " << sw.split() << " ms" << endl;
		cout << "B after row normalization: " << endl << B << endl;

		cout << "B before column normalization: " << endl << B << endl;
		sw.restart();
		for (int i = 1; i < 1E+7; ++i)
			B.normalize_column(5);
		cout << "Took " << sw.split() << " ms" << endl;
		cout << "B after column normalization: " << endl << B << endl;

		cout << (B.cend(1) < B.cbegin(1)) << " " << (B.cbegin(1) < B.cend(1)) << endl;
		cout << (B.cend(1) > B.cbegin(1)) << " " << (B.cbegin(1) > B.cend(1)) << endl;
		B.resize(0, 0);
		cout << (B.cend(1) < B.cbegin(1)) << " " << (B.cbegin(1) <= B.cend(1)) << endl;
		cout << (B.cend(1) >= B.cbegin(1)) << " " << (B.cbegin(1) > B.cend(1)) << endl;
	}
*//*
	{
		BbMatrix A(5, 5, 666.);
		BbVector v(5, 3.);
		BbVector wrong(4);

		cout << "A = " << endl << A << endl;

		sw.restart();

		for (int i = 1; i < 1E+7; ++i)
		{
			A.set_row(3, v);
		}
		cout << "set_row() with BbVector took " << sw.split() << " ms" << endl;

		for (int i = 1; i < 1E+7; ++i)
		{
			A.set_row(2, 5.);
		}
		cout << "set_row() with double took " << sw.split() << " ms" << endl;

		for (int i = 1; i < 1E+7; ++i)
		{
			A.set_column(3, v);
		}
		cout << "set_column() with BbVector took " << sw.split() << " ms" << endl;

		for (int i = 1; i < 1E+7; ++i)
		{
			A.set_column(4, 7.);
		}
		cout << "set_column() with double took " << sw.split() << " ms" << endl;

		for (int i = 1; i < 1E+7; ++i)
		{
			A.set_diagonal({ 0., 0., 0., 0., 0. });
		}
		cout << "set_diagonal() with initializer_list took " << sw.split() << " ms" << endl;

		for (int i = 1; i < 1E+7; ++i)
		{
			A.set_diagonal(v);
		}
		cout << "set_diagonal() with BbVector took " << sw.split() << " ms" << endl;

		cout << "A = " << endl << A << endl;
	}
*/
	return 0;
}
