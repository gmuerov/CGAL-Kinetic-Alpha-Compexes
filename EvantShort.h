#ifndef EVANT_SHORT_H
#define EVANT_SHORT_H

template <class KD, class RS>
class EvantShort:
	public CGAL::Kinetic::internal::Delaunay_event_base_3<class KD, class RS>
{
	typedef Delaunay_event_base_3<KD, RS> baseEvent;
	//kdel -> it's the Kinetic Alfha Complex Kernel
	EvantShort(const RS &s,const typename KD::Face_handle &f,KD *kdel):baseEvent(s, kdel), f_(f){}
	void process(){
		kdel()->hideShowFace(f_);
	}

protected:
  const typename KD::Face_handle f_;

}


#endif