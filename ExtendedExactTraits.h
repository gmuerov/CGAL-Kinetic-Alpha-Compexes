#ifndef EXTENDED_EXACT_TRAITS_H
#define EXTENDED_EXACT_TRAITS_H

#include <CGAL\Kinetic\Exact_simulation_traits.h>
#include "AlphaComplexKineticKernel.h"

class ExtendedExactTraits: public CGAL::Kinetic::Exact_simulation_traits
{
	typedef CGAL::Kinetic::Exact_simulation_traits Base;


public:

	typedef AlphaComplexKineticKernel<Base::Function_kernel> Kinetic_kernel;

	const Kinetic_kernel& kinetic_kernel_object() const {return nkk;}

	ExtendedExactTraits(const Base::Simulator::Time &lb, const Base::Simulator::Time &ub):Base(lb, ub){};
	
protected:

	Kinetic_kernel nkk;
	
};

#endif
