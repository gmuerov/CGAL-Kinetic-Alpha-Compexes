#include "KineticAlphaCompex.h"

KineticAlphaComplex::KineticAlphaComplex(Traits st,
										 double alpha,
										 StaticTriangulation del):
	alphaValue(alpha),
	squaredAlpha(alpha * alpha),
	CGAL::Kinetic::Delaunay_triangulation_2<Traits>(st, del),
	hiddenFaces(del_.number_of_faces())
{
	typedef StaticTriangulation::Finite_faces_iterator FaceIt;

	//setup the initial hidden faces
	for (FaceIt iterator = del.finite_faces_begin();
			iterator != del.finite_faces_end();
			iterator ++)
	{
		std::cout << iterator->dimension();
		switch (iterator->dimension())
		{
		case 0: 
			break;
		case 1:
			{
				StaticTriangulation::Point end1 = iterator->vertex(0)->point();
				StaticTriangulation::Point end2 = iterator->vertex(1)->point();

				typedef StaticTriangulation::Geom_traits GTraits;
			}
			break;
		case 2: 
			{
				StaticTriangulation::Point point1 = iterator->vertex(0)->point();
				StaticTriangulation::Point point2 = iterator->vertex(1)->point();
				StaticTriangulation::Point point3 = iterator->vertex(2)->point();


			}
			break;
		default: break;
		};

	}
}

KineticAlphaComplex::KineticAlphaComplex(Traits st,
										 double alpha):
	alphaValue(alpha), 
	squaredAlpha(alpha * alpha),
	CGAL::Kinetic::Delaunay_triangulation_2<Traits>(st),
	hiddenFaces()
{

}

	
/// Returns the complex as a Triangulation Data Structure
/// The type of the structure can be obtained from AlphaComplex::StaticTDS

KineticAlphaComplex::StaticTDS KineticAlphaComplex::AlphaComplex()
{
	StaticTDS del = KDel::triangulation_data_structure();

	for(std::vector<Face_handle>::iterator iter = hiddenFaces.begin();
		iter != hiddenFaces.end(); iter++)
	{
		del.delete_face(*iter);
	}

	return del;
}