#ifndef EVANT_GABRIEL_H
#define EVANT_GABRIEL_H

template <class KD, class RS>
class EvantGabriel:
	public CGAL::Kinetic::internal::Delaunay_event_base_3<class KD, class RS>

	typedef Delaunay_event_base_3<KD, RS> baseEvent;
	//kdel -> it's the Kinetic Alfha Complex Kernel
	EvantGabriel(const RS &s,const typename KD::Face &f,KD *kdel):baseEvent(s, kdel), f_(f){}
	void process(){
		kdel()->hideShowFacet(f_);
	}

protected:
  const typename KD::Facet f_;
#endif