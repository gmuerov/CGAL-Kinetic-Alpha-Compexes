#ifndef EVENT_SHORT_H
#define EVENT_SHORT_H

template <class KD, class RS>
class EventShortEdge:
	public CGAL::Kinetic::internal::Delaunay_event_base_3<KD, RS>
{
	typedef Delaunay_event_base_3<KD, RS> baseEvent;
	//kdel -> it's the Kinetic Alfha Complex
	EventShortEdge(const RS &s,const typename KD::Edge &e,KD *kdel):baseEvent(s, kdel), e_(e){}
	void process(){
		kdel()->hideShowFace(e_);
	}

protected:
  const typename KD::Edge e_;

};

template <class KD, class RS>
class EventShortFacet:
	public CGAL::Kinetic::internal::Delaunay_event_base_3<KD, RS>
{
	typedef Delaunay_event_base_3<KD, RS> baseEvent;
	//kdel -> it's the Kinetic Alfha Complex
	EventShortFacet(const RS &s,const typename KD::Facet &f,KD *kdel):baseEvent(s, kdel), f_(f){}
	void process(){
		kdel()->hideShowFace(f_);
	}

protected:
  const typename KD::Facet f_;

};

template <class KD, class RS>
class EventShortCell:
	public CGAL::Kinetic::internal::Delaunay_event_base_3<KD, RS>
{
	typedef Delaunay_event_base_3<KD, RS> baseEvent;
	//kdel -> it's the Kinetic Alfha Complex
	EventShortCell(const RS &s,const typename KD::Cell_handle &c,KD *kdel):baseEvent(s, kdel), c_(c){}
	void process(){
		kdel()->hideShowFace(c_);
	}

protected:
  const typename KD::Cell_handle c_;

};
#endif