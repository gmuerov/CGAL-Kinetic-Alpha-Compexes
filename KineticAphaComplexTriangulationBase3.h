#ifndef KINETIC_ALPHA_COMPLEX_TRIANGULATION_BASE_3_H
#define KINETIC_ALPHA_COMPLEX_TRIANGULATION_BASE_3_H

#include <unordered_map>
#include <unordered_set>
#include "EventShort.h"
#include <CGAL\Kinetic\internal\Delaunay_triangulation_base_3.h>

template <class TraitsT, class Visitor>
class KineticAphaComplexTriangulationBase:
	public CGAL::Kinetic::internal::Delaunay_triangulation_base_3<TraitsT, Visitor>
{
	typedef typename CGAL::Kinetic::internal::Delaunay_triangulation_base_3<TraitsT, Visitor> Base;

    typedef typename Base::Cell_circulator CellCirculator;

    //typedef typename Base::Simulator::NT NT;

	typedef typename Base::Facet       Facet;
	typedef typename Base::Cell_handle Cell_handle;
	typedef typename Base::Edge		   Edge;

	typedef typename Base::Moving_object_table::Key Point_key;

	typedef typename TraitsT::S4C3 s4C3;
	typedef typename TraitsT::STC3 sTC3;
	typedef typename TraitsT::SEC3 sEC3;

public:
Facet flip(const Edge &e)
{
    CellCirculator ccir = triangulation_.incident_cells(e);
    CellCirculator end = ccir;

    std::vector<cell_handle> cells;
    int degree = 0;

    while(ccir!= end)
    {
        if (*ccir != deletedcell)
            cells.push_back(cell_handle(ccir));

        ccir++;
        degree++;
    }
        
    if (degree > 3)
        return Facet();

    Base::flip(e);

    Cell_handle deletedcell = e.first();

    Event_key deletedkey = cellslist[deletedcell];
    cellslist.remove(deletedcell);

    simulator()->delete_event(deletedkey);

    for(std::vector<cell_handle>::iterator it = cells.begin();
        it != cells.end(); it++)
    {
        removeShortCertificate(it*);

        CheckHiddenAndAddCertificates(it*);
    }
}

protected:

void Initialization()
{
	Base::Initialization();

    for (Base::All_edges_iterator eit = triangulation_.all_edges_begin();
		eit != triangulation_.all_edges_end(); ++eit)
	{
		if(!CheckShortEdge(eit))
		{
			hiddenEdgeList.insert(eit);
		}
	}
			

    for (Base::All_facets_iterator fit = triangulation_.all_edges_begin();
		fit != triangulation_.all_edges_end(); ++fit)
	    {
			if(!CheckShortFacet(fit))
			{
				hiddenFaceList.insert(fit);
			}
		}

    for (Base::All_cells_iterator cit = triangulation_.all_edges_begin();
		fit != triangulation_.all_edges_end(); ++cit)
	    {
			if(!CheckShortCell(cit))
			{
				hiddenCellList.insert(cit);
			}
		}
    
}

 bool CheckShortCell(Cell_handle cell)
{
	std::vector<Point_key> ids;
	cellPoint(cell, ids);

	CGAL::Sign certSign = s4C3.sign_at(
		point(ids[0]),
		point(ids[1]),
		point(ids[2]),
		point(ids[3]),
        simulator->current_time());

        
	return certSign != CGAL::NEGATIVE;
}

 bool CheckShortFacet(Facet facet)
{
	std::vector<Point_key> ids;
	facetPoint(facet, ids);

	CGAL::Sign certSign = sTC3.sign_at(
		point(ids[0]),
		point(ids[1]),
		point(ids[2]),
        simulator->current_time());
        
	return certSign != CGAL::NEGATIVE;
}

 bool CheckShortEdge(Edge edge)
{
	std::vector<Point_key> ids;
	facetPoint(edge, ids);

	CGAL::Sign certSign = sEC3.sign_at(
		point(ids[0]),
		point(ids[1]),
        simulator->current_time());
        
	return certSign != CGAL::NEGATIVE;
}

void CheckHiddenAndAddCertificates(Cell_handle cell)
{
    std::vector<Point_key> ids;
    cellPoint(c, ids);

    CGAL::Sign certSign = s4C3.sign_at(point(ids[0]),
		point(ids[1]),
		point(ids[2]),
		point(ids[3]),
        simulator->current_time());

        
    bool shortFace = certSign == CGAL::NEGATIVE;

    if (!shortFace)
        hiddenCellList.insert(cell);

    makeShortCertificate(cell);

    for(int i = 0; i < 4; i++)
    {
        Facet facet(cell, i);

    }
}

void audit()
{
    Base::audit();

    if (!has_certificates_)
    {
        for (Base::All_edges_iterator eit = triangulation_.all_edges_begin();
	        eit != triangulation_.all_edges_end(); ++eit)
	            CGAL_assertion(edgesList.count(*eit) <= 0);
        for (Base::All_facets_iterator fit = triangulation_.all_edges_begin();
	        fit != triangulation_.all_edges_end(); ++fit)
	            CGAL_assertion(!facesList.count(*fit) <= 0);
        for (Base::All_cells_iterator cit = triangulation_.all_edges_begin();
	        fit != triangulation_.all_edges_end(); ++cit)
	            CGAL_assertion(!cellsList.count(*cit) <= 0);
    }
    else
    {
        for (Base::Finite_edges_iterator eit = triangulation_.finite_edges_begin();
	            eit != triangulation_.finite_edges_end(); ++eit)
                    simulator()->audit_event(edgesList[*eit]);

        for (Base::All_facets_iterator fit = triangulation_.all_edges_begin();
	        fit != triangulation_.all_edges_end(); ++fit)
	            simulator()->audit_event(facesList.count[*fit]);

        for (Base::All_cells_iterator cit = triangulation_.all_edges_begin();
	        fit != triangulation_.all_edges_end(); ++cit)
	            simulator()->audit_event(cellsList.count[*cit]);
    }
}

typename Simulator::Event_key GetEventKey(const Facet f)
{
    return facetsList[f];
}

typename Simulator::Event_key GetEventKey(const Edge f)
{
    return edgesList[f];
}

typename Simulator::Event_key GetEventKey(const Cell_handle f)
{
    return cellsList[f];
}

void removeShortCertificate(Cell_handle cell)
{
    Event_key renewed_key = cellsList[cell];
    simulator()->delete_event(renewed_key);
        
    cellsList.remove(cell);
}

void removeShortCertificate(Edge edge)
{
    Event_key renewed_key = edgesList[edge];
    simulator()->delete_event(renewed_key);
        
    edgesList.remove(edge);
}

void removeShortCertificate(Facet facet)
{
    Event_key renewed_key = facetsList[facet];
    simulator()->delete_event(renewed_key);
        
    facetsList.remove(facet);
}
    
#pragma region Hide/Show functions

	void hideShowFace(Facet f){

		//check if the Facet is contained in the set
		int face = hiddenFaceList.count(f);
		
		//check for the mirror Facet within the set
		Facet mirror = triangulation_.mirror_facet(f);
		int mirroredFace = hiddenFaceList.count(mirror);

		if(face > 0 || mirroredFace > 0)
		{
			if (face > 0)
				hiddenFaceList.remove(f);
			else
				hiddenFaceList.remove(mirror);
		}
		else
			hiddenFaceList.insert(f);

		typename Simulator::Event_key failed = facetsList[f];

		Base::Certificate core = extract_root_stack(failed);

		 if (core.will_fail()) {
			typename Simulator::Time t= core.failure_time();
		    core.pop_failure_time();
			typename Simulator::Event_key k= simulator()->new_event(t, 
					EventShortFacet(core, f, tr_.wrapper_handle()));
      
			facetsList.insert(f, k);
		} else {
			facetsList.insert(f, simulator()->null_event());
		}
	}

	void hideShowFace(Edge e){
		//Check if the edge is hidden
		int edge = hiddenEdgeList.count(e);

		//Look for the mirrored edge as well
		Edge mirror= triangulation_.mirror_edge(e);
		int mirroredEdge = hiddenEdgeList.count(mirror);

		if (edge > 0 || mirroredEdge > 0)
		{
			if (edge > 0)
				hiddenEdgeList.remove(e);
			else
				hiddenEdgeList.remove(mirror);
		}
		else
			hiddenEdgeList.insert(e);
		
		typename Simulator::Event_key failed = edgesList[e];

		Base::Certificate core = extract_root_stack(failed);

		 if (core.will_fail()) {
			typename Simulator::Time t= core.failure_time();
		    core.pop_failure_time();
			typename Simulator::Event_key k= simulator()->new_event(t, 
										  EventShortEdge(core, e, tr_.wrapper_handle()));
      
			edgesList.insert(e, k);
		} else {
			edgesList.insert(e, simulator()->null_event());
		}
	}

	void hideShowFace(Cell_handle c){
		//Check if the cell is hidden
		int cell = hiddenEdgeList.count(c);

		if (cell > 0)
			hiddenCellList.remove(c);
		else
			hiddenCellList.insert(c);
		
		typename Simulator::Event_key failed = facetsList[c];

		Base::Certificate core = extract_root_stack(failed);

		 if (core.will_fail()) {
			typename Simulator::Time t= core.failure_time();
		    core.pop_failure_time();
			typename Simulator::Event_key k= simulator()->new_event(t, 
										  EventShortCell(core, c, tr_.wrapper_handle()));
			cellsList.insert(c, k);
		} else {
			cellsList.insert(c, simulator()->null_event());
		}
	}

#pragma endregion The tree functions for hiding and revealing faces.

#pragma region Point Extraction functions
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
#pragma endregion Functions for extracting the points out of facets, cells and edges

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
	  makeShortCertificate(*cit);
    }
  }

protected:
//    NT squared_alpha;

	std::set<Cell_handle> hiddenCellList;
	std::set<Facet>		  hiddenFaceList;
	std::set<Edge>	      hiddenEdgeList;

	std::unordered_map<Cell_handle, typename Simulator::Event_key> cellsList;
	std::unordered_map<Facet,       typename Simulator::Event_key> facetsList;
	std::unordered_map<Edge,        typename Simulator::Event_key> edgesList;

};
#endif