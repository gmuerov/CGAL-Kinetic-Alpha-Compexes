#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <CGAL/Kinetic/Cartesian.h>
#include "ShortCalculationStructs.h"
#include "ExtendedExactTraits.h"
#include "KineticAlphaComplexTriangulation3.h"
#include <CGAL/Kinetic/Delaunay_triangulation_3.h>
#include <CGAL/Kinetic/Delaunay_triangulation_visitor_base_3.h>
#include <CGAL/Random.h>
#include <CGAL/determinant.h>
#include <math.h>
#include "UtilityFunctions.h"
//#include "TrajectoryChangeEvent.h"

#define CGAL_CHECK_EXACTNESS
#define CGAL_KINETIC_CHECK_EXPENSIVE

typedef Helper::Traits Traits;

typedef Helper::Traits::Kinetic_kernel KK;
typedef KK::Certificate_function result_type;
typedef KK::Motion_function F;
typedef Helper::Point Point;
typedef Helper::StaticPoint StaticPoint;
typedef CGAL::Kinetic::Delaunay_triangulation_3<Traits> KineticTriangulation;
typedef Helper::Simulator Simulator;
typedef Traits::Function_kernel::Construct_function CF;
typedef Traits::Active_points_3_table::Key Point_key;

///Creates a random linear trajectory.
///If a point and time are given the new trajectory is guaranteed to pass through the same point at that time.
Point makeRandomMovement(CGAL::Random& rand, Point source = Point(), F::NT time = F::NT(0))
{
    std::vector<F::NT> x;
    std::vector<F::NT> y;
    std::vector<F::NT> z;
    for(int j = 0; j < 2; j++)
    {
        x.push_back(rand.get_double(0, 3) * 5 - 2);
        y.push_back(rand.get_double(0, 3) * 7 - 5);
        z.push_back(rand.get_double(0, 3) * 6 - 3); 
    }
    
    
    if (source.x() != Point().x() &&
        source.y() != Point().y() &&
        source.z() != Point().z())
    {
        
        F::NT xs = x[0] + x[1]*time;
        F::NT ys = y[0] + y[1]*time;
        F::NT zs = z[0] + z[1]*time;

        F::NT sourcex = source.x()(time);
        F::NT sourcey = source.y()(time);
        F::NT sourcez = source.z()(time);

        x[0] = x[0] + sourcex - xs;
        y[0] = y[0] + sourcey - ys;
        z[0] = z[0] + sourcez - zs;
    }
    x[0] = x[0] / 10;
    y[0] = y[0] / 10;
    z[0] = z[0] / 10;
    F x_func(x.begin(), x.end());
    F y_func(y.begin(), y.end());
    F z_func(z.begin(), z.end());

    return Point(x_func, y_func, z_func);
}

int main()
{
	Traits tr(0, 1000);

    typedef KineticAlphaComplexTriangulation3<Traits> AC;
    AC beef(tr, Traits::Simulator::NT(2.0));

	KineticAlphaComplexTriangulation3<Traits>::Cell_handle ddd;
    std::vector<StaticPoint> initialPoints;

    CGAL::Random rand;
    for(int i = 0; i < 6; i++)
    {
        StaticPoint new_point = Helper::sampleCube(rand, 2);
        initialPoints.push_back(new_point);
    }

 //   std::vector<int> VisitedIndexes;

 //   std::vector<Point_key> AlphaComplexKeys;

 //   for(int i = 0; i < 6; i++)
 //   {
 //       /*int f = rand.get_int(0, 5);
 //       while (std::find(VisitedIndexes.begin(), VisitedIndexes.end(), f) 
 //                       != VisitedIndexes.end())
 //           f = rand.get_int(0, 5);

 //       VisitedIndexes.push_back(f);*/
 //       AlphaComplexKeys.push_back(Helper::makePointRotate(initialPoints[i],
 //                   &beef, tr.simulator_handle(), 5));
 //   }

 //   /*for(int i = 0; i < 6; i++)
 //   {
 //       if (std::find(VisitedIndexes.begin(), VisitedIndexes.end(), i) 
 //                       != VisitedIndexes.end())
 //       {
 //           AlphaComplexKeys.push_back(tr.active_points_3_table_handle()->insert(
 //               Point(initialPoints[i])));
 //       }
 //   }*/

 //   beef.set_has_certificates(true);
	//Traits::Simulator::Handle sp= tr.simulator_handle();

	//while (sp->next_event_time() != sp->end_time()) 
 //   {
	//	printf("Current event %d\n",sp->current_event_number());
 //       beef.WriteVerticesAndEdges();
 //       sp->set_current_event_number(sp->current_event_number()+1);
 //   }

 //   printf("Simulator time %d\n",tr.simulator_handle()->current_time());

    for(int i = 0; i< 6; i++)
    {

    }


	return 0;
}