#include "KineticAlphaCompex.h"
#include <CGAL/Kinetic/Exact_simulation_traits.h>

int main()
{
	CGAL::Kinetic::Exact_simulation_traits tr(0, 1000);



	KineticAlphaComplex a(tr, 2.0);
	std::cout<<"It runs!";

	return 0;
}