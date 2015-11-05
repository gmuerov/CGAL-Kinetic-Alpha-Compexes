#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <CGAL/Kinetic/Cartesian.h>
#include "ShortCalculationStructs.h"
#include "ExtendedExactTraits.h"
#include "KineticAlphaComplexTriangulation3.h"
#include <CGAL/Kinetic/Delaunay_triangulation_3.h>
#include <CGAL/Kinetic/Delaunay_triangulation_visitor_base_3.h>
#include <CGAL/Random.h>
#include <CGAL/determinant.h>

typedef ExtendedExactTraits Traits;

int main()
{
	Traits tr(0, 1000);

	typedef Traits::Kinetic_kernel KK;
	typedef KK::Point_3 Argument3;
	typedef KK::Certificate_function result_type;
	typedef KK::Motion_function F;
    typedef KK::Point_3 Point;
    typedef CGAL::Kinetic::Delaunay_triangulation_3<Traits> KineticTriangulation;

    KineticAlphaComplexTriangulation3<Traits> beef(tr, Traits::Simulator::NT(2.0));

    //KineticTriangulation triang(tr);
	KineticAlphaComplexTriangulation3<Traits>::Cell_handle ddd;
	//std::cout<<ddd;
	/*KineticAlphaComplexTriangulation3<Traits>::Facet d;
	std::cout<<d.first->vertex(1)->point();*/

    CGAL::Random rand;
    for(int i = 0; i < 6; i++)
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

        F x_func(x.begin(), x.end());
        std::cout<<"X function: "<<x_func<<std::endl;
        F y_func(y.begin(), y.end());
        std::cout<<"Y function: "<<y_func<<std::endl;
        F z_func(z.begin(), z.end());
        std::cout<<"Z function: "<<z_func<<std::endl;

        Point new_point(x_func, y_func, z_func);

        tr.active_points_3_table_handle()->insert(new_point);
    }

    beef.set_has_certificates(true);
    //triang.set_has_certificates(true);
	 
	Traits::Simulator::Handle sp= tr.simulator_handle();
	while (sp->next_event_time() != sp->end_time()) {
		printf("Current event %d\n",sp->current_event_number());
		beef.WriteVerticesAndEdges();
        sp->set_current_event_number(sp->current_event_number()+1);
    }
    printf("Simulator time %d\n",tr.simulator_handle()->current_time());

	return 0;
}



/* Motion function tests


    std::vector<F::NT> coefsAx;
	coefsAx.push_back(F::NT(3.0));
	coefsAx.push_back(F::NT(7.0));
	F ax(coefsAx.begin(), coefsAx.end());

	std::vector<F::NT> coefsAy;
	coefsAy.push_back(F::NT(1.0));
	coefsAy.push_back(F::NT(2.0));
	F ay(coefsAy.begin(), coefsAy.end());

    Argument3 additive(ax, ay, ax);
    Argument3 additive2(ax, ay, ay);

	std::vector<F::NT> coefsBx;
	coefsBx.push_back(F::NT(2.0));
	coefsBx.push_back(F::NT(3.0));
	F bx(coefsBx.begin(), coefsBx.end());

	std::vector<F::NT> coefsBy;
	coefsBy.push_back(F::NT(6.0));
	coefsBy.push_back(F::NT(8.0));
	F by(coefsBy.begin(), coefsBy.end());
	
	std::vector<F::NT> coefsone;
	coefsone.push_back(F::NT(1.0));
	F one(coefsone.begin(), coefsone.end());

	Argument a(ax,ay);
	Argument b(bx,by);

	result_type M1240 =
			CGAL::determinant(ax, bx, ay, one,
							  by, ax, by, one,
							  ay, bx, ax, one,
							  ay, by, bx, one);
	std::cout<< (a.x() + b.x()) / 2<<"\n";
	std::cout<< M1240 *M1240 <<"\n";
	KineticAlphaComplex x(tr, 2.0);
	*/