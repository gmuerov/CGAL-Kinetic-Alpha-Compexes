#ifndef TRAJECTORY_CHANGE_H
#define TRAJECTORY_CHANGE_H

#include <CGAL\Kinetic\Event_base.h>
#include "UtilityFunctions.h"

template <class KDS>
class TrajectoryChangeEvent: public CGAL::Kinetic::Event_base<KDS*>
{
    typedef CGAL::Kinetic::Event_base<KDS*> E;

    public:
        typedef typename KDS::Traits Traits;

        TrajectoryChangeEvent(KDS* link, typename Traits::Simulator::Handle init,
                              std::vector<typename Traits::Active_points_3_table::Key> point,
                              typename Traits::Simulator::NT polSides):E(link), sim(init),
                              angle(polSides),toUpdate(point)
        {}

        TrajectoryChangeEvent()
        {}

        void process()
        {
            Helper::ModifyPoints(toUpdate, kds_, sim, angle);
        }
        
    protected:
       const typename Traits::Simulator::Handle sim;
       const std::vector<typename Traits::Active_points_3_table::Key> toUpdate;
       const typename Traits::Simulator::NT angle;
};

#endif