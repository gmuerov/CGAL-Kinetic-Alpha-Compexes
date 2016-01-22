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
    typedef typename Base::Triangulation Triang;

	typedef typename Base::Facet       Facet;
	typedef typename Base::Cell_handle Cell_handle;
	typedef typename Base::Edge		   Edge;

	typedef typename Base::Moving_object_table::Key Point_key;
    typedef std::pair<Point_key, Point_key> StoredEdge;

	typedef typename TraitsT::ShortTetrahedronCheck S4C3;
	typedef typename TraitsT::ShortTriangleCheck STC3;
	typedef typename TraitsT::ShortEdgeCheck SEC3;

	typedef typename TraitsT::eventShortEdge edgeShortEvent; 
	typedef typename TraitsT::eventShortFacet facetShortEvent; 
	typedef typename TraitsT::eventShortCell cellShortEvent; 
	typedef typename TraitsT::Simulator Simulator;

	int nrOfEdgeFlips,nrOfFacetFlips,nrOfShortEdge,nrOfShortFacet,nrOfShortCell;

	KineticAlphaComplexTriangulationBase(TraitsT tr, NT alpha, Visitor v = Visitor()):
        Base(tr, v), squared_alpha(alpha * alpha)
    {
		nrOfEdgeFlips  = 0;
		nrOfFacetFlips = 0;
		nrOfShortEdge  = 0;
		nrOfShortFacet = 0;
		nrOfShortCell  = 0; 
        Initialization();
    }

    Facet flip(const Edge &e)
	{
		nrOfEdgeFlips++;
        Cell_handle deletedCell = e.first;

		for (int i = 0; i < 4; i++)
        {
            removeShortCertificate(Facet(deletedCell, i));
            for (int j = i + 1; j< 4; j++)
            {
                StoredEdge e(deletedCell->vertex(i)->point(),
                             deletedCell->vertex(j)->point());
                removeShortCertificate(e);
            }
        }
        removeShortCertificate(deletedCell);

        Facet returned = Base::flip(e);
        
        if (returned == Facet())
            return returned;

        std::vector<Cell_handle> cells;
        cells.push_back(returned.first);
        cells.push_back(triangulation_.mirror_facet(returned).first);
        
		for(std::vector<Cell_handle>::iterator cit = cells.begin();
			cit != cells.end(); ++cit)
			RenewCertificates(*cit);

        return returned;
	}

    Edge flip(const Facet &f)
    {
		nrOfFacetFlips++;
        Cell_handle oldCell = f.first;

        Edge returned = Base::flip(f);

        if (returned == Edge())
            return returned;

        CellCirculator edgeCirc = triangulation_.incident_cells(returned);
        CellCirculator done = edgeCirc;
        do
        {
            RenewCertificates(edgeCirc);
            edgeCirc++;
			
        }
        while(edgeCirc != done);
		return returned;
    }

    void audit() const
    {
        printf("Beginning audit.\n");
        Base::audit();

        auditHiddenFaces();

        if (!has_certificates_)
        {
            for (Base::All_edges_iterator eit = triangulation_.all_edges_begin();
	            eit != triangulation_.all_edges_end(); ++eit)
	                CGAL_assertion(edgesList.count(convertEdge(*eit)) <= 0);
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
				if(edgesList.at(convertEdge(*eit)) != Simulator::Event_key() )	
					simulator()->audit_event(edgesList.at(convertEdge(*eit)));
				
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

    void auditHiddenFaces() const
    {
        for (Base::All_cells_iterator cit = triangulation_.all_cells_begin();
	            cit != triangulation_.all_cells_end(); ++cit)
                if (hiddenCellList.count(cit) > 0)
                    CGAL_assertion(CheckShortCell(cit));
                else
                    CGAL_assertion(!CheckShortCell(cit));

        for (Base::Finite_edges_iterator eit = triangulation_.finite_edges_begin();
	                eit != triangulation_.finite_edges_end(); ++eit)
			{
                StoredEdge edge = convertEdge(*eit);
                if(hasShortCertificate(edge))
                    CGAL_assertion( CheckShortEdge(edge));
                else
                    CGAL_assertion(!CheckShortEdge(edge));
				
			}     

        for (Base::Finite_facets_iterator fit = triangulation_.finite_facets_begin();
	            fit != triangulation_.finite_facets_end(); ++fit)
                if(hasShortCertificate(*fit))
                    CGAL_assertion(CheckShortFacet(*fit));
                else
                    CGAL_assertion(!CheckShortFacet(*fit));
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

    typename Simulator::Event_key GetEventKey(const StoredEdge e)
    {
        if (edgesList.count(e) > 0)
            return edgesList[e];
        
        StoredEdge mirror(e.second, e.first);

        if (edgesList.count(mirror) > 0)
            return edgesList[mirror];

        CGAL_ERROR("Couldn't find key for edge");
        return Simulator::Event_key();
    }

    typename Simulator::Event_key GetEventKey(const Cell_handle c)
    {
        Simulator::Event_key key = cellsList[c];
        if (key == Simulator::Event_key())
            CGAL_ERROR("Key not found for cell");

        return cellsList[c];
    }
    

		//----------------------------------displa bigin--------------------
	void displayTest(std::ofstream& outputFile) const
	{	
		outputFile << "Frame"<<std::endl;
        double currentTime = tr_.simulator_handle()->current_time().compute_double(0.01);
		outputFile << "Vertices"<<std::endl;
		for (Base::Finite_vertices_iterator vit = triangulation_.finite_vertices_begin();
			vit != triangulation_.finite_vertices_end(); ++vit)
		{
			outputFile << point(vit->point()).x().value_at(currentTime)<<" ";
			outputFile << point(vit->point()).y().value_at(currentTime)<<" ";
			outputFile << point(vit->point()).z().value_at(currentTime)<<std::endl;

		}

		outputFile << "Edges"<<std::endl;		
        for (All_edges_iterator eit = triangulation_.all_edges_begin();
			eit != triangulation_.all_edges_end(); ++eit) 
		{
			
			if(eit->first->vertex(eit->second)->point().is_valid() && eit->first->vertex(eit->third)->point().is_valid())
			{
				 bool hiddenEdge = false;
				 for (std::set<StoredEdge>::iterator heit = hiddenEdgeList.begin();
					 heit != hiddenEdgeList.end(); ++heit) 
				 { 

					if(eit->first->vertex(eit->second)->point() == heit->first &&
						eit->first->vertex(eit->third)->point() == heit->second ||
						eit->first->vertex(eit->third)->point() == heit->first &&
						eit->first->vertex(eit->second)->point() == heit->second)
					{
						hiddenEdge = true;
					}

				 }

				 if(!hiddenEdge)
				 {
					outputFile << eit->first->vertex(eit->second)->point()<<
							 eit->first->vertex(eit->third )->point()<< std::endl;
				 }

			}
			
        }	
		
		
		outputFile << "Facet"<<std::endl;
		for (All_facets_iterator fit = triangulation_.all_facets_begin();
	                 fit != triangulation_.all_facets_end(); ++fit)
		{ 
			bool pointsValid = true;
			
			for(int i=0; i<4; i++)
				if(i != fit->second)
					if(!fit->first->vertex(i)->point().is_valid())
						pointsValid = false;

			if(pointsValid)
			{
				//check if the Facet is contained in the set
				int face = hiddenFaceList.count(*fit);

				//check for the mirror Facet within the set
				Facet mirror = triangulation_.mirror_facet(*fit);
				int mirroredFace = hiddenFaceList.count(mirror);

				if(face == 0 && mirroredFace == 0)
				{
					for(int i=0; i<4; i++)
						if(i != fit->second)
							outputFile << fit->first->vertex(i)->point();
					outputFile <<std::endl;
				}
			}
			
		}

		outputFile << "Cell"<<std::endl;
		for (Base::All_cells_iterator cit = triangulation_.all_cells_begin();
			cit != triangulation_.all_cells_end(); ++cit)
		{
			
			bool pointsValid = true;
			
			for(int i=0; i<4; i++)
					if(!cit->vertex(i)->point().is_valid())
						pointsValid = false;

			if(pointsValid)
			{

				if(hiddenCellList.find(cit) == hiddenCellList.end())
				{
					for(int i=0; i<4; i++)
						outputFile <<cit->vertex(i)->point();
					outputFile <<std::endl;
				}
				
			}
			
		}
   
	}
	//----------------------------------displa end--------------------
    typename Triang::Vertex_handle insert(Point_key k)
    {
        Triang::Vertex_handle baseVertex = Base::insert(k);

        std::vector<Cell_handle> incidents;
        if(triangulation_.dimension() == 3)
        {
            triangulation_.incident_cells(baseVertex, back_inserter(incidents));

            for(std::vector<Cell_handle>::iterator it = incidents.begin();
                it != incidents.end(); it++)
                RenewCertificates(*it);
        }
        return baseVertex;
    }

    typename Triang::Vertex_handle change_vertex(Point_key k)
    {
        if(!has_certificates_) 
            return NULL;

        Triang::Vertex_handle changedVertex = Base::change_vertex(k);

        std::vector<Cell_handle> cellsAround;

        if(triangulation_.dimension() == 3)
        {
            triangulation_.incident_cells(changedVertex, back_inserter(cellsAround));

            for(std::vector<Cell_handle>::iterator it = cellsAround.begin();
                it != cellsAround.end(); ++it)
            {
                RenewCertificates(*it);
            }

        }

        return changedVertex;
    }

#pragma region Hide/Show functions

	void hideShowFace(Facet f)
    {
        nrOfShortFacet++;
		//check if the Facet is contained in the set
		bool face = hiddenFaceList.find(f) != 
                    hiddenFaceList.end();

		//check for the mirror Facet within the set
		Facet mirror = triangulation_.mirror_facet(f);
		bool mirroredFace = hiddenFaceList.find(mirror) != 
                            hiddenFaceList.end();

		if(face || mirroredFace)
		{
			if (face)
				hiddenFaceList.erase(f);
			else
				hiddenFaceList.erase(mirror);
		}
		else
			hiddenFaceList.insert(f);

		typename Simulator::Event_key failed = facetsList[f];

		Base::Certificate core = extract_root_stack(failed);

		if (core.will_fail()) 
        {
		    typename Simulator::Time t= core.failure_time();
		    core.pop_failure_time();
		    typename Simulator::Event_key k= simulator()->new_event(t, 
				    facetShortEvent(core, f, tr_.wrapper_handle()));
      
		    facetsList[f] = k;
		}
        else
            facetsList.erase(f);
	}

	void hideShowFace(StoredEdge e)
    {
        nrOfShortEdge++;
        if (hiddenEdgeList.count(e) <= 0)
        {
            StoredEdge mirror(e.second, e.first);
            if (hiddenEdgeList.count(mirror) <= 0)
                hiddenEdgeList.insert(e);
            else
                hiddenEdgeList.erase(mirror);
        }
        else
            hiddenEdgeList.erase(e);

		typename Simulator::Event_key failed = edgesList[e];

		Base::Certificate core = extract_root_stack(failed);

		if (core.will_fail()) 
        {
			typename Simulator::Time t= core.failure_time();
		    core.pop_failure_time();
			typename Simulator::Event_key k= simulator()->new_event(t, 
										  edgeShortEvent(core, e, tr_.wrapper_handle()));
			edgesList[e] = k;
		}
        else
            edgesList.erase(e);
	}

	void hideShowFace(Cell_handle c)
    {
        nrOfShortCell++;
		//Check if the cell is hidden
		int cell = hiddenCellList.count(c);

		if (cell > 0)
			hiddenCellList.erase(c);
		else
			hiddenCellList.insert(c);
		
		typename Simulator::Event_key failed = cellsList[c];

		Base::Certificate core = extract_root_stack(failed);

		if (core.will_fail()) 
        {
			typename Simulator::Time t= core.failure_time();
		    core.pop_failure_time();
			typename Simulator::Event_key k= simulator()->new_event(t, 
										  cellShortEvent(core, c, tr_.wrapper_handle()));
			cellsList[c] = k;
		}
         else
             cellsList.erase(c);
	}

#pragma endregion The tree functions for hiding and revealing faces.

    ///Sets the has_certificates_ property and creates/destroys 
    ///the structure's certificates if needed
    void set_has_certificates(bool b)
    {
        if (!has_certificates_ && b) {
            if (triangulation().dimension() == 3) {
	            //Add assertions for our certificates
                for (All_edges_iterator eit = triangulation_.all_edges_begin();
	                 eit != triangulation_.all_edges_end(); ++eit) {
	              CGAL_assertion(!has_event(*eit));
                  CGAL_assertion(!hasShortCertificate(convertEdge(*eit)));
	            }
	            for (All_facets_iterator eit = triangulation_.all_facets_begin();
	                 eit != triangulation_.all_facets_end(); ++eit) {
	              CGAL_assertion(!has_event(*eit));
                  CGAL_assertion(!hasShortCertificate(*eit));
	            }
                for (Base::All_cells_iterator cit = triangulation_.all_cells_begin();
                    cit != triangulation_.all_cells_end(); ++cit){
                        CGAL_assertion(!hasShortCertificate(cit));
                }
	            create_all_certificates();
        	    has_certificates_=true;
            }
        } else if (has_certificates_ && !b) {
            destroy_all_certificates();
            has_certificates_=false;
        }
    }

    typename Base::Moving_object_table* moving_object_table()
    {
        return tr_.active_points_3_table_handle();
    }
	
	
	int getNrOfEdgeFlips()
	{
		return nrOfEdgeFlips;
	}
	
	int getNrOfFacetFlips()
	{
		return nrOfFacetFlips;
	}
	
	int getNrOfShortEdge()
	{
		return nrOfShortEdge;
	}
	
	int getNrOfShortFacet()
	{
		return nrOfShortFacet;
	}
	
	int getNrOfShortCell()
	{
		return nrOfShortCell;
	}

protected:

    void Initialization()
    {
        for (Base::All_edges_iterator eit = triangulation_.all_edges_begin();
		    eit != triangulation_.all_edges_end(); ++eit)
	    {
            StoredEdge convert = convertEdge(*eit);
		    if(!CheckShortEdge(convert) && 
                hiddenEdgeList.count(StoredEdge(convert.second, convert.first))<=0)
		    {
			    hiddenEdgeList.insert(convert);
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

    bool CheckShortCell(const Cell_handle cell) const
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

            return certSign == CGAL::POSITIVE;
        }
        //If the cell is infinite we don't want to put it with the hidden ones
        return true;
    }

    bool CheckShortFacet(Facet facet) const
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

            return certSign == CGAL::POSITIVE;
        }
        //The facet is infinite we don't put it into the hiden facets
	    return true;
    }

    bool CheckShortEdge(StoredEdge edge) const
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
        
	        return certSign == CGAL::POSITIVE;
        }
        //The edge is infinite we don't put it into the hiden edges
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
    void edgePoint(const StoredEdge &e, Oit out) const
    {
        Point_key k = e.first;

		if(k.is_valid())
		{
			*out = k;
			++out;
		}

        Point_key k1 = e.second;
			
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
            Vertex_handle vh = c->vertex(i);
            if (vh != NULL)
            {
		        Point_key k= vh->point();
		        if(k.is_valid())
		        {
			        *out = k;
			        out++;
		        }
            }
	    }
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

        if (hiddenCellList.count(cell) > 0)
            hiddenCellList.erase(cell);
    }

    void removeShortCertificate(StoredEdge edge, bool mirroredEdge = false)
    {
        if(edgesList.find(edge) != edgesList.end())
        {
            Event_key for_removal = edgesList[edge];
            if(for_removal != Simulator::Event_key())
		    {
			    simulator()->delete_event(for_removal);
		    }
            edgesList.erase(edge);
            return;
        }
        
        if(hiddenEdgeList.find(edge) != hiddenEdgeList.end())
            hiddenEdgeList.erase(edge);

        if (!mirroredEdge)
            removeShortCertificate(StoredEdge(edge.second, edge.first), true);
    }

    void removeShortCertificate(Facet facet)
    {

        Facet mirror = triangulation_.mirror_facet(facet);
        bool removed = false;

        if(facetsList.find(facet) != facetsList.end())
		{
            Event_key remove_key = facetsList[facet];
			simulator()->delete_event(remove_key);
            facetsList.erase(facet);
            removed = true;
		}
        if (facetsList.find(mirror) != facetsList.end())
        {
            Event_key mirror_key = facetsList[mirror];
            simulator()->delete_event(mirror_key);
            facetsList.erase(mirror);
            removed = true;
        }

        /*if (!removed)
        {
            std::cout<<"Warning: Facet ";
            printFacet(facet);
            std::cout<<" had no short certificate."<<endl;
        }*/
            

        if(hiddenFaceList.find(facet ) != hiddenFaceList.end())
            hiddenFaceList.erase(facet);
        if(hiddenFaceList.find(mirror) != hiddenFaceList.end())
            hiddenFaceList.erase(mirror);
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

    Certificate edgeRootStack(const StoredEdge &e,
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

    void makeShortCertificate( const StoredEdge &e,
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

    void makeShortCertificate( const StoredEdge &e) {
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
    
    StoredEdge convertEdge(const Edge& e) const
    {
        return StoredEdge(e.first->vertex(e.second)->point(),
                          e.first->vertex(e.third )->point());
    }

    void create_all_certificates() 
    {
        CGAL_precondition(!has_certificates_);
 
        for (All_edges_iterator eit = triangulation_.all_edges_begin();
	     eit != triangulation_.all_edges_end(); ++eit) {
			makeShortCertificate(convertEdge(*eit));
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

    bool hasShortCertificate(const StoredEdge& e) const
    {
        StoredEdge copy(e.second, e.first);
        int count = edgesList.count(e) + edgesList.count(copy);
        return count > 0;
    }

    bool hasShortCertificate(const Facet& f) const
    {
        Facet mirrored = triangulation_.mirror_facet(f);
        bool faceIn = facetsList.find(f) != facetsList.end();
        bool mirrorIn = facetsList.find(mirrored) != facetsList.end();

        return faceIn || mirrorIn;
    }

    bool hasShortCertificate(const Cell_handle& c) const
    { 
        return cellsList.count(c) > 0;
    }

    void RenewCertificates(Cell_handle edgeCirc)
    {
        if (has_certificates_)
        {
            removeShortCertificate(edgeCirc);
            makeShortCertificate(edgeCirc);
        }

        bool cellShort = CheckShortCell(edgeCirc);

        if(!cellShort && edgeCirc->is_valid())
            hiddenCellList.insert(Cell_handle(edgeCirc));

        for(int i = 0; i < 4; i++)
        {
            Facet f(edgeCirc, i);
            
            if (has_certificates_)
            {
                removeShortCertificate(f);
                makeShortCertificate(f);
            }

            if (cellShort && CheckShortFacet(f) &&
                hiddenFaceList.find(triangulation_.mirror_facet(f)) == hiddenFaceList.end())
                hiddenFaceList.insert(f);

            for(int j = i + 1; j < 4; j++)
            { 
                StoredEdge e(edgeCirc->vertex(i)->point(),
                             edgeCirc->vertex(j)->point());
                if(has_certificates_)
                {
                    removeShortCertificate(e);
                    makeShortCertificate(e);
                }

                if (cellShort && CheckShortEdge(e) && 
                     hiddenEdgeList.count(StoredEdge(e.second, e.first))<=0)
                    hiddenEdgeList.insert(e);
            }
        }
    }

    void printFacet(Facet f) const
    {
        for(size_t i = 0; i < 4; i++)
        {
            if (i != f.second)
                std::cout<<f.first->vertex(i)->point();
        }
    }


    NT squared_alpha;
    
    S4C3 s4C3;
    STC3 sTC3;
    SEC3 sEC3;

    std::set<Cell_handle> hiddenCellList;
    std::set<Facet>		  hiddenFaceList;
    std::set<StoredEdge>  hiddenEdgeList;

    std::map<Cell_handle, typename Simulator::Event_key> cellsList;
    std::map<Facet,       typename Simulator::Event_key> facetsList;
    std::map<StoredEdge,  typename Simulator::Event_key> edgesList;

};

#endif