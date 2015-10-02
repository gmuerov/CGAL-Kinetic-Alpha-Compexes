#ifndef KINETIC_ALPHA_COMPLEX_TRIANGULATION_BASE_3_H
#define KINETIC_ALPHA_COMPLEX_TRIANGULATION_BASE_3_H

template <class TraitsT, class Visitor>
class KineticAphaComplexTriangulationBase:
	public CGAL::Kinetic::internal::Delaunay_triangulation_base_3<class TraitsT, class Visitor>
{
	typedef CGAL::Kinetic::internal::Delaunay_triangulation_base_3::Face_handle Face;
	
	void hideShowFace(Face f){
		
	}
}
#endif