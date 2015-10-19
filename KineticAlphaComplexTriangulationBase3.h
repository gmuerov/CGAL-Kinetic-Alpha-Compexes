#ifndef KINETIC_ALPHA_COMPLEX_TRIANGULATION_BASE_3_H
#define KINETIC_ALPHA_COMPLEX_TRIANGULATION_BASE_3_H

//#include <unordered_map>
#include <map>
#include <unordered_set>
#include "EventShort.h"
#include <CGAL\Kinetic\internal\Delaunay_triangulation_base_3.h>

template <class TraitsT, class Visitor>
class KineticAlphaComplexTriangulationBase:
	public CGAL::Kinetic::internal::Delaunay_triangulation_base_3<TraitsT, Visitor>
{
	typedef typename CGAL::Kinetic::internal::Delaunay_triangulation_base_3<TraitsT, Visitor> Base;

public:
    typedef typename Base::Cell_circulator CellCirculator;

    typedef typename Base::Simulator::NT NT;

	typedef typename Base::Facet       Facet;
	typedef typename Base::Cell_handle Cell_handle;
	typedef typename Base::Edge		   Edge;

	typedef typename Base::Moving_object_table::Key Point_key;

	typedef typename TraitsT::ShortTetrahedronCheck S4C3;
	typedef typename TraitsT::ShortTriangleCheck STC3;
	typedef typename TraitsT::ShortEdgeCheck SEC3;

    Facet flip(const Edge &e)
{
    Cell_handle deletedcell = e.first();

    CellCirculator ccir = triangulation_.incident_cells(e);
    CellCirculator end  = ccir;

    std::vector<cell_handle> cells;
    int degree = 0;
    do
    {
        if (ccir != deletedcell)
            cells.push_back(Cell_handle(ccir));

        ccir++;
        degree++;
    }
    while(ccir != end)
        
    if (degree > 3)
        return Facet();

    Facet returnedFacet = Base::flip(e);

    Event_key deletedkey = cellslist[deletedcell];
    cellslist.remove(deletedcell);

    simulator()->delete_event(deletedkey);

    for(std::vector<Cell_handle>::iterator it = cells.begin();
        it != cells.end(); it++)
    {
        removeShortCertificate(it*);
        
        bool shortCell = CheckShortCell(it*);
        if (!shortCell)
            hiddenCellList.insert(it*);

        makeShortCertificate(it*);

        for(int i = 0; i < 4; i++)
        {
            Facet facet(it*, i);
            Facet mirror = triangulation_.mirror_facet(facet);

            if (!shortCell && CheckShortFacet(facet))
                hiddenFaceList.insert(facet);

            removeShortCertificate(facet);
            makeShortCertificate(facet);

            for(int j = i + 1; j < 4; j++)
            {
                if (j != i)
                {
                    Edge e(it*, i, j);

                    if (!shortCell && CheckShortEdge(e))
                        hiddenEdgeList.insert(e);

                    removeShortCertificate(e);
                    makeShortCertificate(e);
                }
            }
        }
    }

    return returnedFacet;
}

    Edge flip(const Facet &f)
    {
        Cell_handle oldCell = f.first();

        Edge returnedEdge = Base::flip(f);

        CellCirculator edgeCell = triangulation_.incident_cells(returnedEdge);
        CellCirculator done = edgeCell;
        
        do
        {
            removeShortCertificate(edgeCell);
            makeShortCertificate(edgeCell);

            for(int i = 0; 
        }
        while(edgeCells != done)

        return returnedEdge;

    }

    void audit() const
    {
        Base::audit();

        if (!has_certificates_)
        {
            for (Base::All_edges_iterator eit = triangulation_.all_edges_begin();
	            eit != triangulation_.all_edges_end(); ++eit)
	                CGAL_assertion(edgesList.count(*eit) <= 0);
            for (Base::All_facets_iterator fit = triangulation_.all_facets_begin();
	            fit != triangulation_.all_facets_end(); ++fit)
	                CGAL_assertion(!facetsList.count(*fit) <= 0);
            for (Base::All_cells_iterator cit = triangulation_.all_cells_begin();
	            cit != triangulation_.all_cells_end(); ++cit)
	                CGAL_assertion(!cellsList.count(cit) <= 0);
        }
        else
        {
            for (Base::Finite_edges_iterator eit = triangulation_.finite_edges_begin();
	                eit != triangulation_.finite_edges_end(); ++eit)
                        simulator()->audit_event(edgesList.at(*eit));

            for (Base::All_facets_iterator fit = triangulation_.all_facets_begin();
	            fit != triangulation_.all_facets_end(); ++fit)
	                simulator()->audit_event(facetsList.at(*fit));

            for (Base::All_cells_iterator cit = triangulation_.all_cells_begin();
	            cit != triangulation_.all_cells_end(); ++cit)
	                simulator()->audit_event(cellsList.at(cit));
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
    
    KineticAlphaComplexTriangulationBase(TraitsT tr, Visitor v = Visitor()):
        Base(tr, v)
    {
        Initialization();
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
        //not correct, fix that
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

protected:
    void Initialization()
    {
        for (Base::All_edges_iterator eit = triangulation_.all_edges_begin();
		    eit != triangulation_.all_edges_end(); ++eit)
	    {
		    if(!CheckShortEdge(*eit))
		    {
			    hiddenEdgeList.insert(*eit);
		    }
	    }
			

        for (Base::All_facets_iterator fit = triangulation_.all_facets_begin();
		    fit != triangulation_.all_facets_end(); ++fit)
	        {
			    if(!CheckShortFacet(*fit))
			    {
				    hiddenFaceList.insert(*fit);
			    }
		    }

        for (Base::All_cells_iterator cit = triangulation_.all_cells_begin();
		    cit != triangulation_.all_cells_end(); ++cit)
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
	    cellPoint(cell, std::back_insert_iterator<std::vector<Point_key> >(ids));

	    CGAL::Sign certSign = s4C3.sign_at(
		    point(ids[0]),
		    point(ids[1]),
		    point(ids[2]),
		    point(ids[3]),
            squared_alpha,
            simulator()->current_time());

        
	    return certSign != CGAL::NEGATIVE;
    }

    bool CheckShortFacet(Facet facet)
    {
	    std::vector<Point_key> ids;
	    facetPoint(facet, std::back_insert_iterator<std::vector<Point_key> >(ids));

	    CGAL::Sign certSign = sTC3.sign_at(
		    point(ids[0]),
		    point(ids[1]),
		    point(ids[2]),
            squared_alpha,
            simulator()->current_time());
        
	    return certSign != CGAL::NEGATIVE;
    }

    bool CheckShortEdge(Edge edge)
    {
	    std::vector<Point_key> ids;
        edgePoint(edge, std::back_insert_iterator<std::vector<Point_key> >(ids));

	    CGAL::Sign certSign = sEC3.sign_at(
		    point(ids[0]),
		    point(ids[1]),
            squared_alpha,
            simulator()->current_time());
        
	    return certSign != CGAL::NEGATIVE;
    }

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
		Point_key k= e.first->vertex(e.second)->point();
			
		if(k.is_valid())
		{
			*out = k;
			++out;
		}

        Point_key k1 = e.first->vertex(e.third)->point();
			
		if(k1.is_valid())
		{
			*out = k1;
			++out;
		}
	
    }

    template <class Oit>
    void cellPoint(const Cell_handle &c, Oit out) const
    {
	    for (unsigned int i=0; i<4; ++i)
	    {	
		    Point_key k= c->vertex(i)->point();
		
		    if(k.is_valid())
		    {
			    *out = k;
			    out++;
		    }	
	    }
    }

#pragma endregion Functions for extracting the points out of facets, cells and edges

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

    Certificate cellRootStack(const Cell_handle &c,
			     const typename Simulator::Time &st) const
    {
    std::vector<Point_key> ids(4);
    cellPoint(c, std::back_insert_iterator<std::vector<Point_key> >(ids));
	
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
        edgePoint(e, std::back_insert_iterator<std::vector<Point_key> >(ids));
	
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
        facetPoint(f, std::back_insert_iterator<std::vector<Point_key> >(ids));
	
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
			      const typename Simulator::Time &st) 
    {
   
        CGAL_precondition(!hasShortCertificate(f));
    
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
        CGAL_precondition(!hasShortCertificate(e));
    
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
			     const typename Simulator::Time &st) 
    {
        CGAL_precondition(!hasShortCertificate(c));
    
        Certificate cert = cellRootStack(c, st);
        if (cert.will_fail()) {
          typename Simulator::Time t= cert.failure_time();
          cert.pop_failure_time();
          typename Simulator::Event_key k = simulator()->new_event(t, cellShortEvent(cert, c, tr_.wrapper_handle()));
	      cellsList.insert(c,k);
        }
    }

    void makeShortCertificate( const Cell_handle &c) 
    {
         makeShortCertificate(c,
		          simulation_traits_object().simulator_handle()->current_time());
    }

    void create_all_certificates() 
    {
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

    bool hasShortCertificate(const Edge& e)
    {
        CellCirculator cc = triangulation_.incident_cells(e);
        CellCirculator sent = cc;
    
        int count = 0;

        do
        {
            int i = cc->index(e.first()->vertex(e.first()));
            int j = cc->index(e.first()->vertex(e.second()));
            count += edgesList.count(Edge(cc, i, j));
            if(count > 0)
                return true;
            cc++;
        }
        while (cc != sent);

        return false;
    }

    bool hasShortCertificate(const Facet& f)
    {
        Facet mirrored = triangulation_.mirror_facet(f);
        int count = facetsList.count(f) + facetsList.count(mirrored);

        return count > 0;
    }

    bool hasShortCertificate(const Cell_handle& c)
{
    return cellsList.count(c) > 0;
}

    double squared_alpha;
    
    S4C3 s4C3;
    STC3 sTC3;
    SEC3 sEC3;

    std::set<Cell_handle> hiddenCellList;
    std::set<Facet>		  hiddenFaceList;
    std::set<Edge>	      hiddenEdgeList;

    std::map<Cell_handle, typename Simulator::Event_key> cellsList;
    std::map<Facet,       typename Simulator::Event_key> facetsList;
    std::map<Edge,        typename Simulator::Event_key> edgesList;

};
#endif