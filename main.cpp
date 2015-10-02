#include "KineticAlphaCompex.h"
#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <CGAL/Kinetic/Cartesian.h>
#include "ShortCalculationStructs.h"
#include "ExtendedExactTraits.h"
#include "AlphaComplexKineticKernel.h"

#include <CGAL/determinant.h>

int main()
{
	typedef ExtendedExactTraits Traits;
	Traits tr(0, 1000);

	typedef Traits::Kinetic_kernel KK;
	typedef KK::Point_2 Argument;
	typedef KK::Point_3 Argument3;
	typedef KK::Certificate_function result_type;
	typedef Traits::Kinetic_kernel::Motion_function F;
	std::vector<F::NT> coefsAx;
	coefsAx.push_back(F::NT(3.0));
	coefsAx.push_back(F::NT(7.0));
	F ax(coefsAx.begin(), coefsAx.end());

	std::vector<F::NT> coefsAy;
	coefsAy.push_back(F::NT(1.0));
	coefsAy.push_back(F::NT(2.0));
	F ay(coefsAy.begin(), coefsAy.end());

	std::vector<F::NT> coefsBx;
	coefsBx.push_back(F::NT(2.0));
	coefsBx.push_back(F::NT(3.0));
	F bx(coefsBx.begin(), coefsBx.end());

	std::vector<F::NT> coefsBy;
	coefsBy.push_back(F::NT(6.0));
	coefsBy.push_back(F::NT(8.0));
	F by(coefsBy.begin(), coefsBy.end());
	
	std::vector<F::NT> coefsone;
	coefsone.push_back(F::NT(1.0));
	F one(coefsone.begin(), coefsone.end());

	Argument a(ax,ay);
	Argument b(bx,by);

	result_type M1240 =
			CGAL::determinant(ax, bx, ay, one,
							  by, ax, by, one,
							  ay, bx, ax, one,
							  ay, by, bx, one);
	std::cout<< (a.x() + b.x()) / 2<<"\n";
	std::cout<< M1240 *M1240 <<"\n";
	KineticAlphaComplex x(tr, 2.0);
	std::cout<<"It runs!";

	return 0;
}