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
	Traits tr(0, 1000);

    typedef KineticAlphaComplexTriangulation3<Traits> AC;
    AC beef(tr, Traits::Simulator::NT(2.0));

	KineticAlphaComplexTriangulation3<Traits>::Cell_handle ddd;
    std::vector<StaticPoint> initialPoints;

    CGAL::Random rand;
	

	//reading
	std::ifstream input( "Points.txt" );
	std::cout<<"Frame"<<std::endl;
	std::cout<<"Vertices"<<std::endl;
	int nrOfPoints = 0;
	for( std::string line; getline( input, line ); )
	{
		double x, y, z;
		input >> x >> y >> z;
		StaticPoint new_point = StaticPoint(x,y,z);		
        initialPoints.push_back(new_point);
		std::cout<<x<<" ";
		std::cout<<y<<" ";
		std::cout<<z<<std::endl;
		nrOfPoints++;
		
	}

    std::vector<int> VisitedIndexes;

    std::vector<Point_key> movingPoints;

    std::vector<Point_key> AlphaComplexKeys;
    
    Simulator::NT angle = Simulator::NT(2) * M_PI/Simulator::NT(5);

    for(int i = 0; i < 6; i++)
    {
        int f = rand.get_int(0, 99/*nrOfPoints-1*/);
        while (std::find(VisitedIndexes.begin(), VisitedIndexes.end(), f) 
                        != VisitedIndexes.end())
            f = rand.get_int(0, 5);

        VisitedIndexes.push_back(f);
        Point_key new_key = Helper::makePointRotate(initialPoints[f], 
                    &beef, &tr, tr.simulator_handle(), angle);
        movingPoints.push_back(new_key);
        AlphaComplexKeys.push_back(new_key);
    }

	Traits::Simulator::Handle sp= tr.simulator_handle();

    sp->new_event(sp->next_time_representable_as_nt() + Simulator::NT(Helper::dt),
        TrajectoryChangeEvent<AC>(&beef, sp, movingPoints, angle));

    for(int i = 0; i < /*nrOfPoints*/100; i++)
    {   //std::cout<<i<<std::endl;
        if (std::find(VisitedIndexes.begin(), VisitedIndexes.end(), i) 
                        == VisitedIndexes.end())
        {
			//std::cout<<"samting "<<std::endl;
            AlphaComplexKeys.push_back(tr.active_points_3_table_handle()->insert(
                Point(initialPoints[i])));
        }
    }


	AC::Triangulation tri = beef.triangulation();

	
	std::cout<<"Edges"<<std::endl;
        for (AC::Triangulation::All_edges_iterator eit = tri.all_edges_begin();
			eit != tri.all_edges_end(); ++eit) 
		{			
			if(eit->first->vertex(eit->second)->point().is_valid() && eit->first->vertex(eit->third)->point().is_valid())
			{
				std::cout<<eit->first->vertex(eit->second)->point()<<
						 eit->first->vertex(eit->third )->point()<< std::endl;				 		
			}			
        }

	std::cout<<"Facet"<<std::endl;
		for (AC::Triangulation::All_facets_iterator fit = tri.all_facets_begin();
	                 fit != tri.all_facets_end(); ++fit)
		{ 
			bool pointsValid = true;
			
			for(int i=0; i<4; i++)
				if(i != fit->second)
					if(!fit->first->vertex(i)->point().is_valid())
						pointsValid = false;

			if(pointsValid)
			{
				for(int i=0; i<4; i++)
					if(i != fit->second)
						std::cout<<fit->first->vertex(i)->point();
				std::cout<<std::endl;
			}
			
		}

	std::cout<<"Cell"<<std::endl;
		for (AC::Triangulation::All_cells_iterator cit = tri.all_cells_begin();
			cit != tri.all_cells_end(); ++cit)
		{
			
			bool pointsValid = true;
			
			for(int i=0; i<4; i++)
					if(!cit->vertex(i)->point().is_valid())
						pointsValid = false;

			if(pointsValid)
			{
				for(int i=0; i<4; i++)
					std::cout<<cit->vertex(i)->point();
				std::cout<<std::endl;
				
			}
			
		}




    beef.set_has_certificates(true);

	while (sp->next_event_time() != sp->end_time()) 
    {
		printf("Frame",sp->current_event_number());
		std::cout<<std::endl;
        beef.WriteVerticesAndEdges();
        sp->set_current_event_number(sp->current_event_number()+1);
    }

    printf("Simulator time %d\n",tr.simulator_handle()->current_time());

    for(int i = 0; i< 6; i++)
    {

    }


	return 0;
}