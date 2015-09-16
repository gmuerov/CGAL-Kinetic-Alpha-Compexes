#ifndef ALPHA_COMPLEX_TRIANGULATION_TRAITS_2_H
#define ALPHA_COMPLEX_TRIANGULATION_TRAITS_2_H

#include<CGAL/Kinetic/basic.h>

template <class Simulation_traits_t, class Triangulation_t>
class AlphaComplexTriangulationDefaultTraits2:
	public CGAL::Kinetic::Delaunay_triangulation_default_traits_2<Simulation_traits_t,Triangulation_t>
{

public:
	
	typedef typename Simulation_traits_t::Kinetic_Kernel KKernel;
	typedef Certificate_generateor<KKernel,ShortEdgeCheck2D<KKernel>> SEC; 
	typedef Certificate_generateor<KKernel,ShortTriangleCheck2D<KKernel>> SEC; 
	 
		
	bool shortEdgeCertificateFailureTime(Edge e, Point_key ks[2], double squarelpha, Time &t, Certificate_data &s) 
	{
		s= sec(
			point(ks[0]), 
			point(ks[1]),
			squarelpha,
			st_.simulator_handle()->current_time(),
			st_.simulator_handle()->end_time());

		return return_certificate_failure_time(e, t, s);
	}

	bool shortTriangleCertificateFailureTime(Edge e, Point_key ks[3], double squarelpha, Time &t, Certificate_data &s) 
	{
		s= stc(
			point(ks[0]), 
			point(ks[1]),
			point(ks[2]),
			squarelpha,
			st_.simulator_handle()->current_time(),
			st_.simulator_handle()->end_time());

		return return_certificate_failure_time(e, t, s);
	}

}

#endif