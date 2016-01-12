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
    

    Simulator::NT quickDistance(StaticPoint a, StaticPoint b)
    {
        Simulator::NT dx = a.x() - b.x();
        Simulator::NT dy = a.y() - b.y();

        return dx * dx + dy * dy;
    }

    StaticPoint rotatePoint(StaticPoint source, StaticPoint center, Simulator::NT angle)
    {
        double approx = angle.exact().to_double();

        double sineA   = sin(approx);
        double cosineA = cos(approx);

        Simulator::NT dx = source.x() - center.x();
        Simulator::NT dy = source.y() - center.y();

        StaticPoint newPoint((dx * cosineA) - (dy * sineA  ) + center.x(),
                             (dx * sineA  ) + (dy * cosineA) + center.y(),
                             source.z());

        return newPoint;
    }

    Point makeMovement(StaticPoint first, StaticPoint second, Simulator::NT dt, Simulator::NT time)
    {
        std::vector<Simulator::NT> x(2);
        std::vector<Simulator::NT> y(2);

        x[1] = (second.x() - first.x())/ dt;
        x[0] = ((dt + time)*first.x() - time * second.x())/dt;

        y[1] = (second.y() - first.y())/ dt;
        y[0] = ((dt + time)*first.y() - time * second.y())/dt;

        F xp(x.begin(), x.end());
        F yp(y.begin(), y.end());
    
        std::vector<Simulator::NT> z(1, first.z());

        Point result(xp, yp, F(z.begin(), z.end()));

        CGAL_assertion(CGAL::compare(first.x(), result.x().value_at(time)) == CGAL::EQUAL);
        CGAL_assertion(CGAL::compare(first.y(), result.y().value_at(time)) == CGAL::EQUAL);
        return result;
    }

    typedef KineticAlphaComplexTriangulation3<Traits> AC;

    Point_key makePointRotate(StaticPoint source, StaticPoint center, Traits* tr,
                              Simulator::NT angle)
    {
        StaticPoint end = rotatePoint(source, center, angle);

        Simulator::Handle sim = tr->simulator_handle();

        Point moving = makeMovement(source, end, Simulator::NT(dt),
                                    sim->next_time_representable_as_nt());

        return tr->active_points_3_table_handle()->insert(moving);
    }

    void ModifyPoints(std::vector<Point_key> points, StaticPoint center, 
                    AC* kac, Simulator::Handle sim, Simulator::NT angle)
    {
        Simulator::NT t = sim->next_time_representable_as_nt();

        for(size_t i = 0; i < points.size(); i++)
        {
            Point current = kac->pointEx(points[i]);

            StaticPoint start(current.x().value_at(t),
                              current.y().value_at(t),
                              current.z().value_at(t));
            
            StaticPoint end = rotatePoint(start, center, angle);

            Point moving = makeMovement(start, end, Simulator::NT(dt), t);

            CGAL_assertion(CGAL::compare(current.x().value_at(t), moving.x().value_at(t)) == CGAL::EQUAL);
            CGAL_assertion(CGAL::compare(current.y().value_at(t), moving.y().value_at(t)) == CGAL::EQUAL);
            CGAL_assertion(CGAL::compare(current.z().value_at(t), moving.z().value_at(t)) == CGAL::EQUAL);

            kac->moving_object_table()->set(points[i], moving);
        }

        sim->new_event(sim->next_time_representable_as_nt() + Simulator::NT(dt),
                        TrajectoryChangeEvent<AC>(kac, sim, points, center, angle));
    }
}

#endif