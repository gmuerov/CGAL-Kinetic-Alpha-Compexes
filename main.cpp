#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <CGAL/Kinetic/Cartesian.h>
#include "ShortCalculationStructs.h"
#include "ExtendedExactTraits.h"
#include "KineticAlphaComplexTriangulation3.h"
#include <CGAL/Kinetic/Delaunay_triangulation_3.h>
#include <CGAL/Kinetic/Delaunay_triangulation_visitor_base_3.h>
#include <CGAL/Random.h>
#include <CGAL/determinant.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "UtilityFunctions.h"
#include <boost/chrono.hpp>
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


int main()
{
	boost::chrono::high_resolution_clock::time_point runTimeStart = boost::chrono::high_resolution_clock::now();
	Traits tr(0, 50);

    typedef KineticAlphaComplexTriangulation3<Traits> AC;
    AC beef(tr, Traits::Simulator::NT(2.0));

	KineticAlphaComplexTriangulation3<Traits>::Cell_handle ddd;
    std::vector<StaticPoint> initialPoints;

    int numP = 50;

    CGAL::Random rand;

	//reading
	std::ifstream input( "Points.txt" );

	int nrOfPoints = 100;
	int startOfRotation = 50;
	int nrOfPointsMoving = 25;
	int nrOfCorners = 12;


	int allPoints = 0;
	for( std::string line; getline( input, line ); )
	{
		double x, y, z;
		input >> x >> y >> z;
		StaticPoint new_point = StaticPoint(x,y,z);		
        initialPoints.push_back(new_point);
		allPoints++;
		
	}

	//nrOfPoints = allPoints;

    std::set<int> VisitedIndexes;

    std::vector<Point_key> movingPoints;

    std::vector<Point_key> AlphaComplexKeys;
    
    Simulator::NT angle = Simulator::NT(2) * M_PI/Simulator::NT(nrOfCorners);

    StaticPoint center = initialPoints[startOfRotation-1];
	
	boost::chrono::high_resolution_clock::time_point initialisationTimeStart = boost::chrono::high_resolution_clock::now();

    for(int i = 0; i < nrOfPointsMoving; i++)
    {   
        /*rendam points selection
		int f = rand.get_int(0, nrOfPoints-1);
        while (std::find(VisitedIndexes.begin(), VisitedIndexes.end(), f) 
                        != VisitedIndexes.end())
            f = rand.get_int(0, nrOfPoints-1);
			*/
		int f = startOfRotation + i;			
        VisitedIndexes.insert(f);
        Point_key new_key = Helper::makePointRotate(initialPoints[f], center,
                    &tr, angle);
        movingPoints.push_back(new_key);
        AlphaComplexKeys.push_back(new_key);
    }

	Traits::Simulator::Handle sp= tr.simulator_handle();

    sp->new_event(sp->next_time_representable_as_nt() + Simulator::NT(Helper::dt),
        TrajectoryChangeEvent<AC>(&beef, sp, movingPoints, center, angle));

    for(int i = 0; i <nrOfPoints; i++)
    {   
        if (std::find(VisitedIndexes.begin(), VisitedIndexes.end(), i) 
                        == VisitedIndexes.end())
        {
            AlphaComplexKeys.push_back(tr.active_points_3_table_handle()->insert(
                Point(initialPoints[i])));
        }
    }

	boost::chrono::high_resolution_clock::time_point initialisationTimeEnd = boost::chrono::high_resolution_clock::now();
	std::cout <<"Initialisation Triangulation: ";
	std::cout << boost::chrono::duration_cast<boost::chrono::milliseconds>(initialisationTimeEnd-initialisationTimeStart) << "\n";

    beef.set_has_certificates(true);
	std::ofstream outputFile;
    //outputFile.open ("out.txt");
	while (sp->next_event_time() != sp->end_time()) 
    {
		//boost::chrono::high_resolution_clock::time_point frameTimeStart = boost::chrono::high_resolution_clock::now();
		printf("Frame %i",sp->current_event_number());
		std::cout<<std::endl;
        //beef.WriteVerticesAndEdges(outputFile);
        sp->set_current_event_number(sp->current_event_number()+1);
		/*boost::chrono::high_resolution_clock::time_point frameTimeEnd = boost::chrono::high_resolution_clock::now();
		std::cout << boost::chrono::duration_cast<boost::chrono::milliseconds>(frameTimeEnd-frameTimeStart) << "\n";*/
    }
	
	//outputFile.close();	
	printf("Nr Of Edge Flips: %d\n",beef.getNrOfEdgeFlips());
	printf("Nr Of Facet Flips: %d\n",beef.getNrOfFacetFlips());
	printf("Nr Of Short Edge: %d\n",beef.getNrOfShortEdge());
	printf("Nr Of Short Facet: %d\n",beef.getNrOfShortFacet());
	printf("Nr Of Short Cell: %d\n",beef.getNrOfShortCell());
	boost::chrono::high_resolution_clock::time_point runTimeEnd = boost::chrono::high_resolution_clock::now();
	
	std::cout <<"Initialisation Triangulation: ";
	std::cout << boost::chrono::duration_cast<boost::chrono::milliseconds>(initialisationTimeEnd-initialisationTimeStart) << "\n";
	std::cout <<"Run Time: ";
	std::cout << boost::chrono::duration_cast<boost::chrono::milliseconds>(runTimeEnd-runTimeStart) << "\n";

    //printf("Simulator time: %d\n",tr.simulator_handle()->current_time());

	return 0;
}