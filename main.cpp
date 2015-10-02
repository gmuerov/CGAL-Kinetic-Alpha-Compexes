#include "KineticAlphaCompex.h"
#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <CGAL/Kinetic/Cartesian.h>
#include "ShortCalculationStructs.h"
#include "ExtendedExactTraits.h"
#include "AlphaComplexKineticKernel.h"

int main()
{
	typedef CGAL::Kinetic::Exact_simulation_traits Traits;
	Traits tr(0, 1000);
	typedef AlphaComplexKineticKernel<Traits::Function_kernel>::Side
	KineticAlphaComplex a(tr, 2.0);
	std::cout<<"It runs!";

	return 0;
}