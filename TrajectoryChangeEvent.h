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
        TrajectoryChangeEvent(KDS* link,
                              typename Traits::Static_kernel::Point_3 point,
                              int polSides):E(link), 
                              sides(polSides),toUpdate(point)
        {}

        TrajectoryChangeEvent()
        {}

        void process()
        {
            std::cout<<"Rotating point: "<<toUpdate<<std::endl;
            Helper::makePointRotate(toUpdate, kds_, kds_->simulator(), sides);
        }
        
    protected:
       const typename Traits::Static_kernel::Point_3 toUpdate;
       const int sides;
};

#endif