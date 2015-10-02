#ifndef ALPHA_COMPLEX_KINETIC_KERNEL_H
#define ALPHA_COMPLEX_KINETIC_KERNEL_H

#include <CGAL\Kinetic\Certificate_generator.h>
#include <CGAL\Kinetic\Cartesian.h>
#include "ShortCalculationStructs.h"

template<class FunctionKernel>
class AlphaComplexKineticKernel:
	public CGAL::Kinetic::Cartesian<FunctionKernel>
{
protected:

	typedef typename AlphaComplexKineticKernel<FunctionKernel> This;
	typedef typename CGAL::Kinetic::Certificate_generator<This, ShortEdgeCheck2D<This>> SEC; 
	typedef typename CGAL::Kinetic::Certificate_generator<This,ShortTriangleCheck2D<This>> STC;
	
	typedef typename CGAL::Kinetic::Certificate_generator<This, ShortEdgeCheck3D<This>> SEC3;
	typedef typename CGAL::Kinetic::Certificate_generator<This, ShortTriangleCheck3D<This>> STC3;
	typedef typename CGAL::Kinetic::Certificate_generator<This, ShortTetrahedronCheck3D<This>> S4C3;

	//add gabriel certificates generators
	 
public :	
	SEC3 ShortEdgeCheck3DObject()
	{
		return SEC3(CGAL::Kinetic::Cartesian<FunctionKernel>::function_kernel_object());
	}

	STC3 ShortTriangle3DObject()
	{
		return STC3(CGAL::Kinetic::Cartesian<FunctionKernel>::function_kernel_object());
	}

	S4C3 ShortTetrahedronCheck3DObject()
	{
		return S4C3(CGAL::Kinetic::Cartesian<FunctionKernel>::function_kernel_object());
	}
};

#endif