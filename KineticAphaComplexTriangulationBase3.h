#ifndef KINETIC_ALPHA_COMPLEX_TRIANGULATION_BASE_3_H
#define KINETIC_ALPHA_COMPLEX_TRIANGULATION_BASE_3_H

#include <unordered_map>

template <class TraitsT, class Visitor>
class KineticAphaComplexTriangulationBase:
	public CGAL::Kinetic::internal::Delaunay_triangulation_base_3<class TraitsT, class Visitor>
{
	typedef CGAL::Kinetic::internal::Delaunay_triangulation_base_3<TraitsT, Visitor>::Facet Facet;

	typedef typename Moving_object_table::Key Point_key;

	typedef TraitsT::S4C3 s4C3;
	typedef TraitsT::STC3 sTC3;
	typedef TraitsT::SEC3 sEC3;

	std::set<Cell_handle> hiddenCellList;
	std::set<Facet> hiddenFaceList;
	std::set<Edge> hiddenEdgeList;

	std::unordered_map<Cell_handle, Simulator::Event_key> cellsList;
	std::unordered_map<Facet, Simulator::Event_key> fectsList;
	std::unordered_map<Edge, Simulator::Event_key>> edgesList;
	
	void hideShowFacet(Facet f){
		
	}

  template <class Oit>
  void facetPoint(const Facet &f, Oit out) const
  {
	for (unsigned int i=0; i<4; ++i)
     	if(i != f.second)
		{	
			Point_key k= f.first->vertex(i)->point();
			
			if(k.is_valid())
			{
				*out = k;
				++out;
			};			
	   };
  }

  template <class Oit>
  void edgePoint(const Edge &e, Oit out) const
  {
	for (unsigned int i=0; i<4; ++i)
     	if(i == f.second || i == f.third )
		{	
			Point_key k= e.first->vertex(i)->point();
			
			if(k.is_valid())
			{
				*out = k;
				++out;
			};	
	   };
  }

  template <class Oit>
  void cellPoint(const Cell_handle &c, Oit out) const
  {
	for (unsigned int i=0; i<4; ++i)
	{	
		Point_key k= c.first->vertex(i)->point();
		
		if(k.is_valid())
		{
			*out = k;
			++out;
		};	
	};
  }

Certificate cellRootStack(const Cell_handle &c,
			 const typename Simulator::Time &st) const
  {
    std::vector<Point_key> ids;
    cellPoint(c, ids);
	
	if (ids.size()==4) 
	{
     
      return s4C3(point(ids[0]),
		 point(ids[1]),
		 point(ids[2]),
		 point(ids[3]),
		 st,
		 simulator()->end_time());
    }

    CGAL_postcondition(0);
    return Certificate();
  }

Certificate edgeRootStack(const Edge &e,
			 const typename Simulator::Time &st) const
  {
    std::vector<Point_key> ids;
    edgePoint(e, ids);
	
	if (ids.size()==2) 
	{
     
      return sEC3(point(ids[0]),
		 point(ids[1]),
		 st,
		 simulator()->end_time());
    }

    CGAL_postcondition(0);
    return Certificate();
  }

Certificate facetRootStack(const Facet &f,
			 const typename Simulator::Time &st) const
  {
    std::vector<Point_key> ids;
    facetPoint(f, ids);
	
	if (ids.size()==3) 
	{
     
      return sTC3(point(ids[0]),
		 point(ids[1]),
		 point(ids[2]),
		 st,
		 simulator()->end_time());
    }

    CGAL_postcondition(0);
    return Certificate();
  }

 void makeShortCertificate( const Facet &f,
			 const typename Simulator::Time &st) {
    CGAL_precondition(!has_event(f));
    CGAL_precondition(!has_degree_3_edge(f));
    
    Certificate cert = facetRootStack(f, st);
    if (cert.will_fail()) {
      typename Simulator::Time t= cert.failure_time();
      cert.pop_failure_time();
      typename Simulator::Event_key k = simulator()->new_event(t, facetShortEvent(cert, f, tr_.wrapper_handle()));
	  fectsList.insert(f,k);
    }
  }

   void makeShortCertificate( const Facet &f) {
     makeShortCertificate(f,
		      simulation_traits_object().simulator_handle()->current_time());
   }

 void makeShortCertificate( const Edge &e,
			 const typename Simulator::Time &st) {
    CGAL_precondition(!has_event(e));
    
    Certificate cert = edgeRootStack(e, st);
    if (cert.will_fail()) {
      typename Simulator::Time t= cert.failure_time();
      cert.pop_failure_time();
      typename Simulator::Event_key k = simulator()->new_event(t, edgeShortEvent(cert, e, tr_.wrapper_handle()));
	  edgesList.insert(e,k);
    }
  }

   void makeShortCertificate( const Edge &e) {
     makeShortCertificate(e,
		      simulation_traits_object().simulator_handle()->current_time());
   }

void makeShortCertificate( const Cell_handle &c,
			 const typename Simulator::Time &st) {
    CGAL_precondition(!has_event(c));
    
    Certificate cert = cellRootStack(c, st);
    if (cert.will_fail()) {
      typename Simulator::Time t= cert.failure_time();
      cert.pop_failure_time();
      typename Simulator::Event_key k = simulator()->new_event(t, cellShortEvent(cert, c, tr_.wrapper_handle()));
	  cellsList.insert(c,k);
    }
  }

   void makeShortCertificate( const Cell_handle &c) {
     makeShortCertificate(c,
		      simulation_traits_object().simulator_handle()->current_time());
   }


    void create_all_certificates() {
    CGAL_precondition(!has_certificates_);
 
    for (All_edges_iterator eit = triangulation_.all_edges_begin();
	 eit != triangulation_.all_edges_end(); ++eit) {
      if (is_degree_3(*eit) && !has_degree_4_vertex(*eit)) {
		make_certificate(*eit);
		makeShortCertificate(*eit);
      }
    }
    for (All_facets_iterator eit = triangulation_.all_facets_begin();
	 eit != triangulation_.all_facets_end(); ++eit) {
      if (!has_degree_3_edge(*eit)) {
		make_certificate(*eit);
		makeShortCertificate(*eit);
      }
    }
    for (All_cells_iterator cit= triangulation_.all_cells_begin(); 
	 cit != triangulation_.all_cells_end(); ++cit) {
      v_.create_cell(cit);
	  makeShortCertificate(*eit);
    }
  }

}
#endif