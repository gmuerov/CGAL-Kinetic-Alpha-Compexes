#ifndef ALPHA_COMPLEX_DEFAULT_TRAITS_3
#define ALPHA_COMPLEX_DEFAULT_TRAITS_3

#include<CGAL/Kinetic/basic.h>
#include<CGAL/Kinetic/Certificate_generator.h>
#include<CGAL/Kinetic/Delaunay_triangulation_default_traits_2.h>
#include "ShortCalculationStructs.h"

template <class Simulation_traits_t, class Triangulation_t>
class AlphaComplexTriangulationDefaultTraits2:
	public CGAL::Kinetic::Delaunay_triangulation_default_traits_2<Simulation_traits_t,Triangulation_t>
{

public:
	
	typedef typename Simulation_traits_t::Kinetic_Kernel KKernel;
	typedef Certificate_generator<KKernel, ShortEdgeCheck3D       <KKernel>> SEC; 
	typedef Certificate_generator<KKernel, ShortTriangleCheck3D   <KKernel>> STC;
	typedef Certificate_generator<KKernel, ShortTetrahedronCheck3D<KKernel>> S4C;
	 
		
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

	bool shortTetrahedronCertificateFailureTime(Edge e, Point_key ks[4], double squarelpha, Time &t, Certificate_data &s) 
	{
		s= stetc(
			point(ks[0]), 
			point(ks[1]),
			point(ks[2]),
			point(ks[3]),
			squarelpha,
			st_.simulator_handle()->current_time(),
			st_.simulator_handle()->end_time());

		return return_certificate_failure_time(e, t, s);
	}

protected:
	SEC sec;
	STC stc;
	S4C stetc; 
}

#endif