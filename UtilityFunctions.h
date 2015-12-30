#ifndef UTILITY_FS
#define UTILITY_FS

#include "ExtendedExactTraits.h"
#include "TrajectoryChangeEvent.h"
#include "KineticAlphaComplexTriangulation3.h"
#include <CGAL\Random.h>
#include <math.h>


namespace Helper{

    typedef ExtendedExactTraits Traits;
    typedef Traits::Simulator Simulator;
    typedef Traits::Kinetic_kernel::Point_3 Point;
    typedef Traits::Static_kernel::Point_3 StaticPoint;
    typedef Traits::Kinetic_kernel::Motion_function F;
    typedef Traits::Active_points_3_table::Key Point_key;

    int dt = 10;

    StaticPoint sampleCube(CGAL::Random& rand, double side)
    {
        Simulator::NT x = Simulator::NT(rand.get_double(-side/2, side/2));
        Simulator::NT y = Simulator::NT(rand.get_double(-side/2, side/2));
        Simulator::NT z = Simulator::NT(rand.get_double(-side/2, side/2));
    
        return StaticPoint(x, y, z);
    }

    StaticPoint rotatePoint(StaticPoint stuff, Simulator::NT angle)
    {
        double approx = angle.exact().to_double();

        double sineA = sin(approx);
        double cosineA = cos(approx);

        return StaticPoint(stuff.x() * cosineA - stuff.y() * sineA,
                           stuff.x() * sineA   + stuff.y() * cosineA,
                           stuff.z());
    }

    Point makeMovement(StaticPoint first, StaticPoint second, Simulator::NT dt, Simulator::NT time)
    {
        std::vector<Simulator::NT> x(2);
        std::vector<Simulator::NT> y(2);

        x[1] = (first.x() - second.x())/ dt;
        x[0] = ((dt + time)*first.x() - time * second.x())/dt;

        y[1] = (first.y() - second.y())/ dt;
        y[0] = ((dt + time)*first.y() - time * second.y())/dt;

        F xp(x.begin(), x.end());
        F yp(y.begin(), y.end());
    
        std::vector<Simulator::NT> z(1, first.z());
        return Point(xp, yp, F(z.begin(), z.end()));
    }

    typedef KineticAlphaComplexTriangulation3<Traits> AC;

    Point_key makePointRotate(StaticPoint source, AC* kac,
                          Simulator::Handle sim, int sides)
    {
        Simulator::NT angle = Simulator::NT(2)/Simulator::NT(sides);

        StaticPoint end = rotatePoint(source, angle);
        
        sim->new_event(sim->next_time_representable_as_nt() + Simulator::NT(dt),
            TrajectoryChangeEvent<AC>(kac, end, sides));


        Point moving = makeMovement(source, end, Simulator::NT(10), 
                                    sim->next_time_representable_as_nt());

        return kac->moving_object_table()->insert(moving);
    }

}

#endif