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
	// If testWith50RandomPoints is true the program will run a test with random points
	bool testWith50RandomPoints = true;
	// If testWith50ProteinPoints is true the program will run a test with points from the protein 1KT3.pdb database
	bool testWith50ProteinPoints = false;
	// The indexis of the points that will be in motion
	int nrOfPointsMoving[] = {5,6,7,8,9};

	int nrOfPoints,startOfRotation,nrOfCorners;
	double alpha;
	std::string outputTxtFile;
	bool useProteinPoints,useRandomPoints,makeNewRandomPoints;

	if(testWith50RandomPoints)
	{
		 // Nr of points used
		 nrOfPoints = 50;
		 // The rotation point
		 startOfRotation = 25;	
		 // Number of corners is used to simulate a circle	
		 nrOfCorners = 5;
		 // The Alpha
		 alpha = 2.0;
		 // If useProteinPoints true, read points from the protein 1KT3.pdb database
		 useProteinPoints = false;
		 // If useRandomPoints true, read points from the radom points file
		 useRandomPoints = true;
		 // If makeNewRandomPoints true, make and read a file with  random points
		 makeNewRandomPoints = false;
		 //Output file name
		 outputTxtFile = "out.txt";
	}

	if(testWith50ProteinPoints)
	{
		 // Nr of points used
		 nrOfPoints = 50;
		 // The rotation point
		 startOfRotation = 4;
		 // Number of corners is used to simulate a circle
		 nrOfCorners = 5;
		 // The Alpha
		 alpha = 2.0;
		 // If useProteinPoints true, read points from the protein 1KT3.pdb database
		 useProteinPoints = true;
		 // If useRandomPoints true, read points from the radom points file
		 useRandomPoints = false;
		 // If makeNewRandomPoints true, make and read a file with  random points
		 makeNewRandomPoints = false;
		 //Output file name
		 outputTxtFile = "out.txt";
	}

	boost::chrono::high_resolution_clock::time_point runTimeStart = boost::chrono::high_resolution_clock::now();
	Traits tr(0, nrOfCorners*10);

    typedef KineticAlphaComplexTriangulation3<Traits> AC;

	KineticAlphaComplexTriangulation3<Traits>::Cell_handle ddd;
    std::vector<StaticPoint> initialPoints;

    CGAL::Random rand;

    AC beef(tr, Traits::Simulator::NT(alpha));
	
	//Reading Protein Points
	if(useProteinPoints)
	{
		std::ifstream input( "Points.txt" );
		for( std::string line; getline( input, line ); )
		{
			double x, y, z;
			input >> x >> y >> z;
			StaticPoint new_point = StaticPoint(x,y,z);		
			initialPoints.push_back(new_point);
		}
	
	}

	//Reading Random Points
	if(useRandomPoints)
	{
		std::ifstream input( "randPoint.txt" );
		for( std::string line; getline( input, line ); )
		{
			double x, y, z;
			input >> x >> y >> z;
			StaticPoint new_point = StaticPoint(x,y,z);		
			initialPoints.push_back(new_point);
		}
	
	}
	
	//Making and Reading Random Points
	if(makeNewRandomPoints)
	{
		std::ofstream outputPiontFile;
		outputPiontFile.open ("randPoint.txt");
		CGAL::Random rand;
		for(int i=0;i<=nrOfPoints;i++ )
		{	
			StaticPoint new_point = Helper::sampleCube(rand,5);		
			initialPoints.push_back(new_point);
			outputPiontFile <<new_point.x()<<" "<<new_point.y()<<" "<<new_point.z()<<std::endl;
		}
		
		outputPiontFile.close();
	}
	
    std::set<int> VisitedIndexes;

    std::vector<Point_key> movingPoints;

    std::vector<Point_key> AlphaComplexKeys;
    
    Simulator::NT angle = Simulator::NT(2) * M_PI/Simulator::NT(nrOfCorners);

    StaticPoint center = initialPoints[startOfRotation-1];

	boost::chrono::high_resolution_clock::time_point initialisationTimeStart = 
        boost::chrono::high_resolution_clock::now();

    for(int i = 0; i < sizeof(nrOfPointsMoving) / sizeof(int); i++)
    {   

		int f = nrOfPointsMoving[i];			
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
    beef.set_has_certificates(true);


	std::ofstream outputFile;
    outputFile.open (outputTxtFile);
	
	while (sp->next_event_time() != sp->end_time()) 
    {
		printf("Events %i",sp->current_event_number());
		std::cout<<std::endl;
        beef.WriteVerticesAndEdges(outputFile);
        sp->set_current_event_number(sp->current_event_number()+1);
    }
	
	outputFile.close();
	
	std::cout <<"Nr Of Points: ";
	std::cout << nrOfPoints << "\n";
	std::cout <<"Nr Of Rotation: ";
	std::cout << startOfRotation << "\n";
	std::cout <<"Nr Of Point Moving: ";
	for(int i = 0; i < sizeof(nrOfPointsMoving) / sizeof(int); i++)
		std::cout << nrOfPointsMoving[i] << ",";
	std::cout << "\n";
	std::cout <<"Nr Of Corners: ";
	std::cout << nrOfCorners << "\n";

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

	return 0;
}