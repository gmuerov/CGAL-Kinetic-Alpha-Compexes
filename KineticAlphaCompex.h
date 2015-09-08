#ifndef KINETIC_ALPHA_COMPLEX_H
#define KINETIC_ALPHA_COMPLEX_H

#include <CGAL/Kinetic/Delaunay_triangulation_2.h>
#include <CGAL/Kinetic/Event_base.h>
#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <CGAL/triangulation_data_structure_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Point_2.h>

class KineticAlphaComplex: 
	public CGAL::Kinetic::Delaunay_triangulation_2<CGAL::Kinetic::Exact_simulation_traits>
{
	typedef CGAL::Kinetic::Exact_simulation_traits Traits;
	typedef CGAL::Kinetic::Delaunay_triangulation_2<Traits> KDel;
	typedef KDel::Triangulation StaticTriangulation;
	typedef StaticTriangulation::Triangulation_data_structure StaticTDS;
	typedef KDel::Face_handle Face_handle;
	typedef KDel::Vertex_handle Vertex_handle;

	private:
		///The faces that need to be deleted when the alpha complex is requested.
		std::vector<Face_handle> hiddenFaces;

		///The value that determines if an edge is short
		double alphaValue;

		///Keeping the squared alpha to save some computation time
		double squaredAlpha;

	public:
		KineticAlphaComplex(Traits st,
							double alpha,
							StaticTriangulation del);

		KineticAlphaComplex(Traits st,
							double alpha);

		StaticTDS AlphaComplex();
};

#endif