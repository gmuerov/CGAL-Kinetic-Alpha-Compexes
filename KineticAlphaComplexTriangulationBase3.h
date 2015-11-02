#ifndef KINETIC_ALPHA_COMPLEX_TRIANGULATION_BASE_3_H
#define KINETIC_ALPHA_COMPLEX_TRIANGULATION_BASE_3_H

//#include <unordered_map>
#include <map>
#include <unordered_set>
#include "EventShort.h"
#include <CGAL\Kinetic\internal\Delaunay_triangulation_base_3.h>
#include <ostream>

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

	typedef typename TraitsT::eventShortEdge edgeShortEvent; 
	typedef typename TraitsT::eventShortFacet facetShortEvent; 
	typedef typename TraitsT::eventShortCell cellShortEvent; 
	typedef typename TraitsT::Simulator Simulator;

    Facet flip(const Edge &e)
	{

        Cell_handle deletedCell = e.first;

		Facet returned = Base::flip(e);
        
        if (returned == Facet())
            return returned;

        std::vector<Cell_handle> cells;
        cells.push_back(returned.first);
        cells.push_back(triangulation_.mirror_facet(returned).first);

		Event_key deletedkey = cellsList[deletedCell];
		cellsList.erase(deletedCell);
		if(deletedkey != Simulator::Event_key())
		{
			simulator()->delete_event(deletedkey);
		}
		
		for(std::vector<Cell_handle>::iterator cit = cells.begin();
			cit != cells.end(); ++cit)
		{
			removeShortCertificate(*cit);
			bool shortCell = CheckShortCell(*cit);

			if (!shortCell && (*cit)->is_valid())
				hiddenCellList.insert(*cit);
				

			makeShortCertificate(*cit);

			for(int i = 0; i < 4; i++)
			{
				Facet facet(*cit, i);
				Facet mirror = triangulation_.mirror_facet(facet);

                if (!shortCell && CheckShortFacet(facet) && hiddenFaceList.count(mirror) == 0)
					hiddenFaceList.insert(facet);

				removeShortCertificate(facet);
				makeShortCertificate(facet);

				for(int j = i + 1; j < 4; j++)
				{
					if (j != i)
					{
						Edge e(*cit, i, j);

                        
						if (!shortCell && CheckShortEdge(e))
							hiddenEdgeList.insert(e);

						removeShortCertificate(e);
						makeShortCertificate(e);
					}
				}
			}
		}
        return returned;
	}

    Edge flip(const Facet &f)
    {
        Cell_handle oldCell = f.first;

        Edge returned = Base::flip(f);

        if (returned == Edge())
            return returned;

        CellCirculator edgeCirc = triangulation_.incident_cells(returned);
        CellCirculator done = edgeCirc;

        do
        {
            removeShortCertificate(edgeCirc);
            makeShortCertificate(edgeCirc);

            bool cellShort = CheckShortCell(edgeCirc);

            if(!cellShort && edgeCirc->is_valid())
                hiddenCellList.insert(Cell_handle(edgeCirc));

            for(int i = 0; i < 4; i++)
            {
                Facet f(edgeCirc, i);

                removeShortCertificate(f);
                makeShortCertificate(f);

                if (!cellShort && CheckShortFacet(f))
                    hiddenFaceList.insert(f);

                for(int j = i + 1; j < 4; j++)
                {
                    if(i != j)
                    {
                        Edge e(edgeCirc, i, j);

                        removeShortCertificate(e);
                        makeShortCertificate(e);

                        if (!cellShort && CheckShortEdge(e))
                            hiddenEdgeList.insert(e);
                    }
                }
            }
			/*for(int i=0;i<4;i++)
				std::cout<<edgeCirc->vertex(i)->point()<<"  ";
			std::cout<<std::endl;*/
            edgeCirc++;
			
        }
        while(edgeCirc != done);
		return returned;
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
	                CGAL_assertion(facetsList.count(*fit) <= 0);
            for (Base::All_cells_iterator cit = triangulation_.all_cells_begin();
	            cit != triangulation_.all_cells_end(); ++cit)
	                CGAL_assertion(cellsList.count(cit) <= 0);
        }
        else
        {
            for (Base::Finite_edges_iterator eit = triangulation_.finite_edges_begin();
	                eit != triangulation_.finite_edges_end(); ++eit)
			{
				if(edgesList.at(*eit) != Simulator::Event_key() )	
					simulator()->audit_event(edgesList.at(*eit));
				
			}     

            for (Base::All_facets_iterator fit = triangulation_.all_facets_begin();
	            fit != triangulation_.all_facets_end(); ++fit)
			{
				if(facetsList.at(*fit) != Simulator::Event_key())
					simulator()->audit_event(facetsList.at(*fit));
			}
	                

            for (Base::All_cells_iterator cit = triangulation_.all_cells_begin();
	            cit != triangulation_.all_cells_end(); ++cit)
			{
				if(cellsList.at(cit) != Simulator::Event_key())
					simulator()->audit_event(cellsList.at(cit));
			}
	                
        }
    }

    typename Simulator::Event_key GetEventKey(const Facet f)
    {
        if(facetsList.count(f) > 0)
            return facetsList[f];
        else
        {
            Facet mir = triangulation_.mirror_facet(f);
            if (facetsList.count(mir) > 0)
                return facetsList[triangulation_.mirror_facet(f)];
            else
                CGAL_ERROR("No key found for facet");
            
            return Simulator::Event_key();
        }
    }

    typename Simulator::Event_key GetEventKey(const Edge e)
    {
        CellCirculator cc = triangulation_.incident_cells(e);
        CellCirculator done = cc;

        do
        {
            Edge current = GetEdgeInCell(cc, e);
            if (edgesList.count(current) > 0)
                return edgesList[current];
            cc++;
        }while(cc!=done);

        CGAL_ERROR("Couldn't find key for edge");
        return Simulator::Event_key();
    }

    typename Simulator::Event_key GetEventKey(const Cell_handle c)
    {
        return cellsList[c];
    }
    
    KineticAlphaComplexTriangulationBase(TraitsT tr, NT alpha, Visitor v = Visitor()):
        Base(tr, v), squared_alpha(alpha * alpha)
    {
        Initialization();
    }

	void displayTest() const
	{
        double currentTime = tr_.simulator_handle()->current_time().compute_double(0.01);
		std::cout<<std::endl<<"------------------"<<std::endl<<"Vertices Coordinates"<<std::endl<<"------------------"<<std::endl;
		for (Base::Finite_vertices_iterator vit = triangulation_.finite_vertices_begin();
			vit != triangulation_.finite_vertices_end(); ++vit)
		{
			
			std::cout<<"---------------"<<vit->point()<<"------------------"<<std::endl;
			std::cout<<"X :"<<point(vit->point()).x().value_at(currentTime)<<std::endl;
			std::cout<<"Y :"<<point(vit->point()).y().value_at(currentTime)<<std::endl;
			std::cout<<"Z :"<<point(vit->point()).z().value_at(currentTime)<<std::endl;

		}

		std::cout<<std::endl<<"---------------"<<std::endl<<"Edges"<<std::endl<<"------------------"<<std::endl;
        for (All_edges_iterator eit = triangulation_.all_edges_begin();
			eit != triangulation_.all_edges_end(); ++eit) 
		{
			std::cout<<eit->first->vertex(eit->second)->point()<<
                       eit->first->vertex(eit->third )->point()<< std::endl;
        }	

		std::cout<<std::endl<<"---------------"<<std::endl<<"HidenEdges"<<std::endl<<"------------------"<<std::endl;
        for (std::set<Edge>::iterator heit = hiddenEdgeList.begin();
			heit != hiddenEdgeList.end(); ++heit) 
		{ 
			std::cout<<heit->first->vertex(heit->second)->point()<<
                       heit->first->vertex(heit->third )->point()<< std::endl;
        }

		std::cout<<std::endl<<"---------------"<<std::endl<<"HidenFacet"<<std::endl<<"------------------"<<std::endl;
        for (std::set<Facet>::iterator hfit = hiddenFaceList.begin();
			hfit != hiddenFaceList.end(); ++hfit) 
		{ 
			for(int i=0; i<4; i++)
				if(i != hfit->second)
					std::cout<<hfit->first->vertex(i)->point();
			std::cout<<std::endl;
        }

		std::cout<<std::endl<<"---------------"<<std::endl<<"HidenCell"<<std::endl<<"------------------"<<std::endl;
        for (std::set<Cell_handle>::iterator hcit = hiddenCellList.begin();
			hcit != hiddenCellList.end(); ++hcit) 
		{ 
			for(int i=0; i<4; i++)
				std::cout<<(*hcit)->vertex(i)->point();
			std::cout<<std::endl;
        }
					   
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
				hiddenFaceList.erase(f);
			else
				hiddenFaceList.erase(mirror);
		}
		else
			hiddenFaceList.insert(f);

		typename Simulator::Event_key failed = facetsList[f];

		Base::Certificate core = extract_root_stack(failed);

		 if (core.will_fail()) {
			typename Simulator::Time t= core.failure_time();
		    core.pop_failure_time();
			typename Simulator::Event_key k= simulator()->new_event(t, 
					facetShortEvent(core, f, tr_.wrapper_handle()));
      
			facetsList[f] = k;
		} else {
			facetsList[f] = simulator()->null_event();
		}
	}

	void hideShowFace(Edge e){

		//Look for the other edge representations as well
		CellCirculator circ = triangulation_.incident_cells(e);
        CellCirculator done = circ;
        bool removed = false;
        do
        {
            //Find the indexes of the edge in the current cell
            int i = circ->index(e.first->vertex(e.second));
            int j = circ->index(e.first->vertex(e.third));

            Edge currentEdge(circ, i, j);
            
		    //Check if the edge is hidden
		    int edgeCount = hiddenEdgeList.count(currentEdge);

            if(edgeCount > 0)
            {
                removed = true; 
                hiddenEdgeList.erase(currentEdge);
                break;
            }
            circ++;
        }
        while (circ!=done);

        if(!removed)
            hiddenEdgeList.insert(e);

		typename Simulator::Event_key failed = edgesList[e];

		Base::Certificate core = extract_root_stack(failed);

		 if (core.will_fail()) {
			typename Simulator::Time t= core.failure_time();
		    core.pop_failure_time();
			typename Simulator::Event_key k= simulator()->new_event(t, 
										  edgeShortEvent(core, e, tr_.wrapper_handle()));
			edgesList[e] = k;
		} else {
			edgesList[e] = simulator()->null_event();
		}
	}

	void hideShowFace(Cell_handle c){
		//Check if the cell is hidden
		int cell = hiddenCellList.count(c);

		if (cell > 0)
			hiddenCellList.erase(c);
		else
			hiddenCellList.insert(c);
		
		typename Simulator::Event_key failed = cellsList[c];

		Base::Certificate core = extract_root_stack(failed);

		 if (core.will_fail()) {
			typename Simulator::Time t= core.failure_time();
		    core.pop_failure_time();
			typename Simulator::Event_key k= simulator()->new_event(t, 
										  cellShortEvent(core, c, tr_.wrapper_handle()));
			cellsList[c] = k;
		} else {
			cellsList[c] = simulator()->null_event();
		}
	}

#pragma endregion The tree functions for hiding and revealing faces.

    ///Sets the has_certificates_ property and creates/destroys 
    ///the structure's certificates if needed
    void set_has_certificates(bool b)
    {
        if (!has_certificates_ && b) {
            if (triangulation().dimension() == 3) {
	            //Add assertions for our certificates
                /*for (All_edges_iterator eit = triangulation_.all_edges_begin();
	                 eit != triangulation_.all_edges_end(); ++eit) {
	              CGAL_assertion(!has_event(*eit));
	            }
	            for (All_facets_iterator eit = triangulation_.all_facets_begin();
	                 eit != triangulation_.all_facets_end(); ++eit) {
	              CGAL_assertion(!has_event(*eit));
	            }*/
	            create_all_certificates();
        	    has_certificates_=true;
            }
        } else if (has_certificates_ && !b) {
            destroy_all_certificates();
            has_certificates_=false;
        }
    }

protected:
    
    Edge GetEdgeInCell(Cell_handle cc, Edge e)
    {
        int i = cc->index(e.first->vertex(e.second));
        int j = cc->index(e.first->vertex(e.third));

        return Edge(cc, i, j);
    }

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
			    if(!CheckShortCell(cit) && cit->is_valid())
			    {
				    hiddenCellList.insert(cit);
			    }
		    }
    }

    bool CheckShortCell(const Cell_handle cell)
    {
	    std::vector<Point_key> ids;
	    cellPoint(cell, std::back_insert_iterator<std::vector<Point_key> >(ids));
        typename Kinetic_kernel::Function_kernel::Construct_function cf;
        if (ids.size() == 4)
	    {
            CGAL::Sign certSign = s4C3.sign_at(
		        point(ids[0]),
		        point(ids[1]),
		        point(ids[2]),
		        point(ids[3]),
                cf(squared_alpha),
                simulator()->current_time());

	        return certSign != CGAL::NEGATIVE;
        }

        return true;
    }

    bool CheckShortFacet(Facet facet)
    {
	    std::vector<Point_key> ids;
	    facetPoint(facet, std::back_insert_iterator<std::vector<Point_key> >(ids));
        typename Kinetic_kernel::Function_kernel::Construct_function cf;
        if (ids.size() == 3)
        {
	        CGAL::Sign certSign = sTC3.sign_at(
		        point(ids[0]),
		        point(ids[1]),
		        point(ids[2]),
                cf(squared_alpha),
                simulator()->current_time());

            return certSign != CGAL::NEGATIVE;
        }
        //The facet is infinite we dont't put it in too the hiden facets
	    return true;
    }

    bool CheckShortEdge(Edge edge)
    {
	    std::vector<Point_key> ids;
        edgePoint(edge, std::back_insert_iterator<std::vector<Point_key> >(ids));
        typename Kinetic_kernel::Function_kernel::Construct_function cf;
		if(ids.size() == 2)
        {
	        CGAL::Sign certSign = sEC3.sign_at(
		        point(ids[0]),
		        point(ids[1]),
                cf(squared_alpha),
                simulator()->current_time());
        
	        return certSign != CGAL::NEGATIVE;
        }
        //The edge is infinite we dont't put it in too the hiden edges
        return true;
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
		//std::cout<<" cell point:";
	    for (unsigned int i=0; i<4; ++i)
	    {
		    Point_key k= c->vertex(i)->point();
			//std::cout<<k;
		    if(k.is_valid())
		    {
			    *out = k;
			    out++;
		    }	
	    }
		//std::cout<<std::endl;
    }

#pragma endregion Functions for extracting the points out of facets, cells and edges

    void removeShortCertificate(Cell_handle cell)
    {
        Event_key renewed_key = cellsList[cell];
        if(renewed_key != Simulator::Event_key())
		{
			simulator()->delete_event(renewed_key);
		}
		
        cellsList.erase(cell);
    }

    void removeShortCertificate(Edge edge)
    {
        CellCirculator cc = triangulation_.incident_cells(edge);
        CellCirculator done = cc;
        do
        {
            Edge current = GetEdgeInCell(cc, edge);

            if (edgesList.count(current) > 0)
            {
                Event_key for_removal = edgesList[current];
                simulator()->delete_event(for_removal);
                edgesList.erase(current);
                return;
            }

            cc++;
        } while(cc != done);
    }

    void removeShortCertificate(Facet facet)
    {
        if(facetsList.count(facet) > 0)
		{
            Event_key remove_key = facetsList[facet];
			simulator()->delete_event(remove_key);
            facetsList.erase(facet);
		}
        else
        {
            Facet mirror = triangulation_.mirror_facet(facet);
            if (facetsList.count(mirror) > 0)
            {
                Event_key mirror_key = facetsList[mirror];
                simulator()->delete_event(mirror_key);
                facetsList.erase(mirror);
            }
        }

    }

    Certificate cellRootStack(const Cell_handle &c,
			     const typename Simulator::Time &st) const
    {
        typename Kinetic_kernel::Function_kernel::Construct_function cf;
        std::vector<Point_key> ids;
        cellPoint(c, std::back_insert_iterator<std::vector<Point_key> >(ids));
	    
        if (ids.size()==4) 
        {
            return s4C3(point(ids[0]),
		        point(ids[1]),
		        point(ids[2]),
		        point(ids[3]),
			    cf(squared_alpha),
		        st,
		        simulator()->end_time());
        }

        //CGAL_postcondition(0);
        return Certificate();
    }

    Certificate edgeRootStack(const Edge &e,
			     const typename Simulator::Time &st) const
    {
        typename Kinetic_kernel::Function_kernel::Construct_function cf;
        std::vector<Point_key> ids;
        edgePoint(e, std::back_insert_iterator<std::vector<Point_key> >(ids));
	
        if (ids.size()==2) 
        {
            return sEC3(point(ids[0]),
		        point(ids[1]),
				cf(squared_alpha),
		        st,
		        simulator()->end_time());
        }

        //CGAL_postcondition(0);
        return Certificate();
    }

    Certificate facetRootStack(const Facet &f,
			      const typename Simulator::Time &st) const
    {
        typename Kinetic_kernel::Function_kernel::Construct_function cf;
        std::vector<Point_key> ids;
        facetPoint(f, std::back_insert_iterator<std::vector<Point_key> >(ids));
	
	    if (ids.size()==3) 
	    {
     
          return sTC3(point(ids[0]),
		     point(ids[1]),
		     point(ids[2]),
			 cf(squared_alpha),
		     st,
		     simulator()->end_time());
        }

        //CGAL_postcondition(0);
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
            facetsList[f] = k;
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
	      edgesList[e] = k;
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
		  cellsList[c] = k;
        }
    }

    void makeShortCertificate( const Cell_handle &c) 
    {
         makeShortCertificate(c,
		          simulation_traits_object().simulator_handle()->current_time());
    }


    void create_all_certificates() 
    {
        std::cout<<"Creating certificates...";
        CGAL_precondition(!has_certificates_);
 
        for (All_edges_iterator eit = triangulation_.all_edges_begin();
	     eit != triangulation_.all_edges_end(); ++eit) {
			makeShortCertificate(*eit);
			if (is_degree_3(*eit) && !has_degree_4_vertex(*eit)) {
				make_certificate(*eit);
          }
        }

        for (All_facets_iterator fit = triangulation_.all_facets_begin();
	     fit != triangulation_.all_facets_end(); ++fit) {
			makeShortCertificate(*fit);
			if (!has_degree_3_edge(*fit)) {
				make_certificate(*fit);
          }
        }

        for (All_cells_iterator cit= triangulation_.all_cells_begin(); 
	     cit != triangulation_.all_cells_end(); ++cit) {
          v_.create_cell(cit);
	      makeShortCertificate(cit);
        }
      }

    bool hasShortCertificate(const Edge& e)
    {
        CellCirculator cc = triangulation_.incident_cells(e);
        CellCirculator sent = cc;
    
        int count = 0;

        do
        {
            int i = cc->index(e.first->vertex(e.second));
            int j = cc->index(e.first->vertex(e.third));
            count = edgesList.count(Edge(cc, i, j));
            if (count > 0) return true;
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

    NT squared_alpha;
    
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