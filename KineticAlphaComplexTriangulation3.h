#ifndef KINETIC_ALPHA_COMPLEX_TRIANGULATION_3_H
#define KINETIC_ALPHA_COMPLEX_TRIANGULATION_3_H

#include <CGAL/Kinetic/basic.h>

#include "KineticAlphaComplexTriangulationBase3.h"

#include <CGAL/Kinetic/Listener.h>
#include <CGAL/Kinetic/Ref_counted.h>

// Triangulations
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

// Local helpers
#include <CGAL/Kinetic/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Kinetic/listeners.h>
#include <CGAL/Kinetic/Delaunay_triangulation_visitor_base_3.h>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4355) // complaint about using 'this' to
#endif                          // initialize a member

template <class Traits>
struct AlphaComplex3Types
{
  typedef typename Traits::Active_points_3_table MPT;
  typedef typename Traits::Kinetic_kernel KK;
  typedef CGAL::Kinetic::Delaunay_triangulation_cell_base_3<Traits> CFBI;

  typedef CGAL::Triangulation_vertex_base_3<typename Traits::Instantaneous_kernel> CVB;
  typedef CGAL::Triangulation_data_structure_3<CVB, CFBI> TDS;

  typedef CGAL::Delaunay_triangulation_3<typename Traits::Instantaneous_kernel, TDS> Default_triangulation;

};

//! A 3D kinetic Delaunay triangulation.
template <class TraitsT,
	  class Visitor= CGAL::Kinetic::Delaunay_triangulation_visitor_base_3,
	  class TriangulationT= typename AlphaComplex3Types<TraitsT>::Default_triangulation>
class KineticAlphaComplexTriangulation3: 
    public CGAL::Kinetic::Ref_counted<KineticAlphaComplexTriangulation3<TraitsT, Visitor, TriangulationT> > 
{

private:
  typedef KineticAlphaComplexTriangulation3<TraitsT, Visitor, TriangulationT> This_AC3;
  typedef KineticAlphaComplexTriangulation3<TraitsT, Visitor, TriangulationT> This;

public:
  typedef typename TraitsT::Kinetic_kernel::Side_of_oriented_sphere_3::result_type Root_stack;

  typedef typename TriangulationT::Cell_handle Cell_handle;
  typedef typename TriangulationT::Facet Facet;
  typedef typename TriangulationT::Edge Edge;

  typedef typename TraitsT Traits;
  typedef typename TraitsT::Simulator::Event_key Event_key;
  typedef typename TraitsT::Active_points_3_table::Key Point_key;
  typedef typename TraitsT::Simulator Simulator;

private:
  struct Base_traits: public TraitsT {
    typedef TriangulationT Triangulation;
    typedef typename TraitsT::Kinetic_kernel::Side_of_oriented_sphere_3 Side_of_oriented_sphere_3;
    typedef typename TraitsT::Kinetic_kernel::Orientation_3 Orientation_3;
    
    typedef typename TraitsT::Kinetic_kernel::ShortEdgeCheck3        ShortEdgeCheck;
    typedef typename TraitsT::Kinetic_kernel::ShortTriangleCheck3    ShortTriangleCheck;
    typedef typename TraitsT::Kinetic_kernel::ShortTetrahedronCheck3 ShortTetrahedronCheck;

    // The next typedef is needed for VC++ as it does not pick the definition from ten lines above
    typedef typename TraitsT::Kinetic_kernel::Side_of_oriented_sphere_3::result_type Root_stack;

    typedef CGAL::Kinetic::internal::Delaunay_3_edge_flip_event<This_AC3, Root_stack> Edge_flip;
    typedef typename CGAL::Kinetic::internal::Delaunay_3_facet_flip_event<This_AC3, Root_stack> Facet_flip;

	typedef typename EventShortEdge<This_AC3, Root_stack> eventShortEdge;
	typedef typename EventShortFacet<This_AC3, Root_stack> eventShortFacet;
	typedef typename EventShortCell<This_AC3, Root_stack> eventShortCell;
	

    Side_of_oriented_sphere_3 side_of_oriented_sphere_3_object() const
    {
      return TraitsT::kinetic_kernel_object().side_of_oriented_sphere_3_object();
    }

    Orientation_3 orientation_3_object() const
    {
      return TraitsT::kinetic_kernel_object().orientation_3_object();
    }

    Base_traits(This_AC3 *t, const TraitsT &tr): TraitsT(tr), wr_(t) {}

    This_AC3* wrapper_handle() {
      return wr_;
    }
    const This_AC3* wrapper_handle() const
    {
      return wr_;
    }

    This_AC3 *wr_;
  };

 
  friend class CGAL::Kinetic::internal::Delaunay_event_base_3<This, Root_stack>;  

  friend class CGAL::Kinetic::internal::Delaunay_3_edge_flip_event<This, Root_stack>;

  friend class CGAL::Kinetic::internal::Delaunay_3_facet_flip_event<This, Root_stack>;

  

  typedef KineticAlphaComplexTriangulationBase<Base_traits, Visitor> ACBase;
  CGAL_KINETIC_DECLARE_LISTENERS(typename TraitsT::Simulator,
				 typename TraitsT::Active_points_3_table)

public:
	
    
    typedef typename ACBase::StoredEdge StoredEdge;

	typename Simulator::Event_key GetEventKey(StoredEdge e)
	{
		return kdel_.GetEventKey(e);
	}

	typename Simulator::Event_key GetEventKey(Facet f)
	{
		return kdel_.GetEventKey(f);
	}

	typename Simulator::Event_key GetEventKey(Cell_handle c)
	{
		return kdel_.GetEventKey(c);
	}
	
	void hideShowFace(StoredEdge e)
	{
		kdel_.hideShowFace(e);
	}

	void hideShowFace(Facet f)
	{
		kdel_.hideShowFace(f);
	}

	void hideShowFace(Cell_handle c)
	{
		kdel_.hideShowFace(c);
	}

	typename ACBase::Point pointEx(Point_key p)
	{
		return kdel_.point(p);
	}

	typedef typename ACBase::Finite_vertices_iterator Finite_vertices_iterator;
	void WriteVerticesAndEdges(std::ofstream& outputFile)
	{
		kdel_.displayTest(outputFile);
	}
	
	
	int getNrOfEdgeFlips()
	{
		return kdel_.getNrOfEdgeFlips();
	}
	
	int getNrOfFacetFlips()
	{
		return kdel_.getNrOfFacetFlips();
	}
	
	int getNrOfShortEdge()
	{
		return kdel_.getNrOfShortEdge();
	}
	
	int getNrOfShortFacet()
	{
		return kdel_.getNrOfShortFacet();
	}
	
	int getNrOfShortCell()
	{
		return kdel_.getNrOfShortCell();
	}
  //! Initialize it.
  KineticAlphaComplexTriangulation3(TraitsT tr,typename TraitsT::Simulator::NT alpha, Visitor v= Visitor()): 
    kdel_(Base_traits(this, tr), alpha, v) {
    CGAL_KINETIC_INITIALIZE_LISTENERS(tr.simulator_handle(),
				      tr.active_points_3_table_handle());
  }

  //! The type of the underlying triangulation
  typedef TriangulationT Triangulation;
  //! access the underlying triangulation
  const Triangulation& triangulation() const
  {
    return kdel_.triangulation();
  }

  Visitor& visitor() {
    return kdel_.visitor();
  }

  const Visitor& visitor() const
  {
    return kdel_.visitor();
  }



  void write(std::ostream &out) const
  {
    kdel_.write(out);
  }


  //! make the structure have or not have certificates
  void set_has_certificates(bool tf) {
    kdel_.set_has_certificates(tf);
  }
  
  void audit() const
  {
    kdel_.audit();
  }

  //! true if the structure has certificates
  bool has_certificates() const
  {
    return kdel_.has_certificates();
  }

  void erase(Point_key k) {
    kdel_.delete_vertex(k);
    on_geometry_changed();
  }

  void set(Point_key k) {
    kdel_.change_vertex(k);
  }

  void insert(Point_key k) {

    kdel_.insert(k);
    on_geometry_changed();
  }

  void flip(const typename ACBase::Edge &edge) {
    kdel_.flip(edge);
    on_geometry_changed();
  }

  void flip(const typename ACBase::Facet &flip_facet) {
    kdel_.flip(flip_facet);
    on_geometry_changed();
  }

  typename ACBase::Simulator* simulator()
  {
      return kdel_.simulator();
  }

  typename ACBase::Moving_object_table* moving_object_table()
  {
      return kdel_.moving_object_table();
  }

  CGAL_KINETIC_LISTENER1(TRIANGULATION)

 void on_geometry_changed() {
    CGAL_KINETIC_NOTIFY(TRIANGULATION);
  }

  ACBase kdel_;
};

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif