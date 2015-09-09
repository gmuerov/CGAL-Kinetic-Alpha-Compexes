#ifndef SHORT_CALCULATIONS_STRUCT_H
#define SHORT_CALCULATIONS_STRUCT_H

#include <CGAL/Kinetic/basic.h>
#include <CGAL/determinant.h>

template <class KK>
struct ShortEdgesCheck2D:
{
	ShortEdgesCheck2D(){}
	typedef typename KK::Certificate_function result_type;
	typedef typename KK::Point_2 Argument;

	result_type operator()(Argument a, Argument b, double alpha)
	{
		typedef typename KK::Certificate_function FT;
		
		FT ax = a.x();
		FT ay = a.y();
		FT bx = b.x();

		return result_type();
	}
}

template<class KK>
struct ShortTriangleCheck2D:
{
	ShortEdgesCheck2D(){}
	typedef typename KK::Certificate_function result_type;
	typedef typename KK::Point_2 Argument;

	result_type operator()(Argument a, Argument b, Argument c, double alpha)
	{
		return result_type();
	}

}

#endif