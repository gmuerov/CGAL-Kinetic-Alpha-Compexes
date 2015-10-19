#ifndef ALPHA_COMPLEX_KINETIC_KERNEL_H
#define ALPHA_COMPLEX_KINETIC_KERNEL_H

#include <CGAL\Kinetic\Certificate_generator.h>
#include <CGAL\Kinetic\Cartesian.h>
#include "ShortCalculationStructs.h"

template<class FunctionKernel>
class AlphaComplexKineticKernel:
	public CGAL::Kinetic::Cartesian<FunctionKernel>
{
	typedef typename AlphaComplexKineticKernel<FunctionKernel> This;

public :

	typedef typename CGAL::Kinetic::Certificate_generator<This, ShortEdgeCheck2D<This>> SEC; 
	typedef typename CGAL::Kinetic::Certificate_generator<This, ShortTriangleCheck2D<This>> STC;
	
	typedef typename CGAL::Kinetic::Certificate_generator<This, ShortEdgeCheck3D<This>> ShortEdgeCheck3;
	typedef typename CGAL::Kinetic::Certificate_generator<This, ShortTriangleCheck3D<This>> ShortTriangleCheck3;
	typedef typename CGAL::Kinetic::Certificate_generator<This, ShortTetrahedronCheck3D<This>> ShortTetrahedronCheck3;
	 	
	ShortEdgeCheck3 ShortEdgeCheck3DObject()
	{
		return ShortEdgeCheck3(CGAL::Kinetic::Cartesian<FunctionKernel>::function_kernel_object());
	}

	ShortTriangleCheck3 ShortTriangleCheck3DObject()
	{
		return ShortTriangleCheck3(CGAL::Kinetic::Cartesian<FunctionKernel>::function_kernel_object());
	}

	ShortTetrahedronCheck3 ShortTetrahedronCheck3DObject()
	{
		return ShortTetrahedronCheck3(CGAL::Kinetic::Cartesian<FunctionKernel>::function_kernel_object());
	}

	//GEC3 GabrielEdgeCheck3DObject()
	//{
	//	return GEC3(CGAL::Kinetic::Cartesian<FunctionKernel>::function_kernel_object());
	//}

	//GTC3 GabrielTriangleCheck3DObject()
	//{
	//	return GTC3(CGAL::Kinetic::Cartesian<FunctionKernel>::function_kernel_object());
	//}
};

#endif