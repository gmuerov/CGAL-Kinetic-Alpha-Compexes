#ifndef EVENT_SHORT_H
#define EVENT_SHORT_H

#include <CGAL\Kinetic\internal\Delaunay_triangulation_base_3.h>

template <class KD, class RS>
class EventShortEdge:
	public CGAL::Kinetic::internal::Delaunay_event_base_3<KD, RS>
{
	typedef Delaunay_event_base_3<KD, RS> baseEvent;

public:
	//kdel -> it's the Kinetic Alpha Complex
	EventShortEdge(const RS &s, const typename KD::StoredEdge &e, KD *kdel):
		baseEvent(s, kdel), e_(e){}

	void process(){
		kdel()->hideShowFace(e_);
    }

    void audit(typename KD::Event_key k)const{
        KD::Event_key key = kdel()->GetEventKey(e_);
        if (key != k) {
            CGAL_ERROR("Mismatch, for label " << k << " had event " << key);
        }
        CGAL_assertion(key == k);
    }

protected:
  const typename KD::StoredEdge e_;

};

template <class KD, class RS>
class EventShortFacet:
	public CGAL::Kinetic::internal::Delaunay_event_base_3<KD, RS>
{
	typedef Delaunay_event_base_3<KD, RS> baseEvent;
public:
	//kdel -> it's the Kinetic Alpha Complex
	EventShortFacet(const RS &s, const typename KD::Facet &f, KD *kdel):
		baseEvent(s, kdel), f_(f){}

	void process(){
		kdel()->hideShowFace(f_);
	}

     void audit(typename KD::Event_key k)const{
        KD::Event_key key = kdel()->GetEventKey(f_);
        if (key != k) {
            CGAL_ERROR("Mismatch, for label " << k << " had event " << key);
        }
        CGAL_assertion(key == k);
    }

protected:
  const typename KD::Facet f_;

};

template <class KD, class RS>
class EventShortCell:
	public CGAL::Kinetic::internal::Delaunay_event_base_3<KD, RS>
{
	typedef Delaunay_event_base_3<KD, RS> baseEvent;
public:
	//kdel -> it's the Kinetic Alpha Complex
	EventShortCell(const RS &s, const typename KD::Cell_handle &c, KD *kdel):
		baseEvent(s, kdel), c_(c){}

	void process(){
		kdel()->hideShowFace(c_);
	}

     void audit(typename KD::Event_key k)const{
        KD::Event_key key = kdel()->GetEventKey(c_);
        if (key != k) {
            CGAL_ERROR("Mismatch, for label " << k << " had event "<<key);
        }
        CGAL_assertion(key == k);
    }

protected:
  const typename KD::Cell_handle c_;

};
#endif