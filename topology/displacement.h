//
// Created by kusch on 12/09/19.
//

#ifndef NEST_DISPLACEMENT_H
#define NEST_DISPLACEMENT_H

// C++ includes:
#include <bitset>
#include <fstream>

// Includes from topology:
#include "position.h"

// redefine some  variable in order to avoid scope
#define EQUAL_TMP EQUAL
#define AT_TMP AT
#undef EQUAL
#undef AT

// Include from CGAL library (manage mesh calculation)
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/Surface_mesh_shortest_path/Surface_mesh_shortest_path.h>

namespace nest
{

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                             Point;
typedef CGAL::Surface_mesh<Point>                   TriMesh;
typedef CGAL::AABB_face_graph_triangle_primitive<TriMesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, TriMesh> Traits_surface;
typedef CGAL::Surface_mesh_shortest_path<Traits_surface> Surface_mesh_shortest_path;
typedef CGAL::Surface_mesh_shortest_path< Traits_surface>::Source_point_iterator Source_point_it;
typedef Surface_mesh_shortest_path::Face_location Face_location;
template < int D, class T = double >
class Displacement
{
public:
  /**
   * Constructor.
   */
  Displacement( const Position<D,T>& p_source, const Position<D,T>& p_target,
          const std::bitset< D >& periodic, const Position< D,T >& extent);

  /**
   * Default constructor.
   */
  Displacement();

  /**
   * Copy constructor.
   */
  Displacement( const Displacement& other );

  /**
   * transforme the displacement in list of Position vector
   * @return array of Position
   */
  Position<D> to_vector() const;

  /**
   * The distance between two points.
   * @returns Euclidian norm of the vector or mesh dependant distance.
   */
  T length() const;

  /**
   * Print position to output stream.
   *
   * Format: Only as many coordinates as dimensions,
   *         separated by spaces [default], no trailing space.
   *
   * @param out output stream
   * @param sep separator character
   */
  void print( std::ostream& out, char sep = ' ' ) const;

  /**
   * Unary minus.
   * @returns opposite sens of displacament.
   */
  Displacement operator-() const;

  /**
   * Shift of displacement (use only for shift in the case of anchor).
   * @returns shifting displacement.
   */
  Displacement operator-( const Position< D, T >& anchor ) const;

  /**
   * Define the mesh if the distance is dependant of it
   */
  void define_mesh(const std::string& path_mesh);

protected:
   Position<D,T> p_source_;
   Position<D,T> p_target_;
   std::bitset< D > periodic_;
   Position< D,T > extent_;
   Position<D,T> displ_;
private:
    static bool mesh_use; // Boolean to know if the distance is mesh dependant or not
    static TriMesh primal; // The mesh
};

template < int D, class T >
TriMesh Displacement<D,T>::primal; // the mesh in CGAL

template < int D, class T >
bool Displacement<D,T>::mesh_use = false;

template < int D, class T >
inline Displacement< D, T >::Displacement( const Position<D,T>& p_source, const Position<D,T>& p_target,
        const std::bitset< D >& periodic, const Position< D,T >& extent)
{
  p_source_ = Position< D,T>(p_source);
  p_target_ = Position< D,T>(p_target);
  periodic_ = std::bitset< D >(periodic);
  extent_ = Position< D,T >(extent);
  displ_ = Position< D,T>(p_target - p_source);
  // take in count the boundary condition
  for ( int i = 0; i < D; ++i )
  {
    if ( periodic[ i ] )
    {
      displ_[ i ] = -0.5 * extent[ i ]
        + std::fmod( displ_[ i ] + 0.5 * extent[ i ], extent[ i ] );
      if ( displ_[ i ] < -0.5 * extent[ i ] )
      {
        displ_[ i ] += extent[ i ];
      }
    }
  }
}

template < int D, class T >
inline Displacement< D, T >::Displacement()
{
    displ_ = Position< D > ();
    p_target_ = Position< D > ();
    p_source_ = Position< D > ();
    periodic_ = std::bitset< D >();
    extent_ = Position< D,T >();
}

template < int D, class T >
inline Displacement< D, T >::Displacement( const Displacement< D, T >& other )
{
    displ_ = Position< D > (other.displ_);
    p_target_ = Position< D > (other.p_target_);
    p_source_ = Position< D > (other.p_source_);
    periodic_ = std::bitset< D >(other.periodic_);
    extent_ = Position< D,T >(other.extent_);
}

template < int D, class T >
inline Position<D> Displacement< D, T >::to_vector() const {
    return Position<D>(displ_);
}

template < int D, class T >
inline T Displacement< D, T >::length() const {
    if (mesh_use)
    {
        // create the object for computing path
        Surface_mesh_shortest_path surface_mesh(primal);
        // Find the face of the starting and ending point of the path
        Face_location face_1 = surface_mesh.locate<Traits>(Point(p_source_[0],p_source_[1],p_source_[2]));
        Face_location face_2 = surface_mesh.locate<Traits>(Point(p_target_[0],p_target_[1],p_target_[2]));
        // clear the object for the path and add the starting point
        surface_mesh.remove_all_source_points();
        surface_mesh.add_source_point(face_1);
        surface_mesh.build_sequence_tree();
        // find the path until the end
        Surface_mesh_shortest_path::Shortest_path_result result = surface_mesh.shortest_distance_to_source_points(face_2.first,face_2.second);
      return result.first;
    } else {
        //std::cout<<" no mesh"<<std::endl;
        return displ_.length(); // Euclidian distance
    }
}

template < int D, class T >
inline void Displacement< D, T >::print( std::ostream& out, char sep ) const{
    displ_.print(out,sep);
}

template < int D, class T >
inline Displacement< D, T > Displacement< D, T >::operator-() const
{
  Displacement< D > d;
  d=*this;
  Position<D> p_tmp(d.p_source_);
  d.p_source_ = d.p_target_;
  d.p_target_ = p_tmp;
  d.displ_ = -d.displ_;
  return d;
}

template < int D, class T >
inline Displacement< D, T > Displacement< D, T >::operator-( const Position< D, T >& anchor ) const
{
  Displacement d = Displacement(this->p_source_-anchor,
                                this->p_target_,
                                this->periodic_,
                                this->extent_);
  return d;
}

template < int D, class T >
inline void Displacement< D, T >::define_mesh(const std::string& path_mesh) {
    mesh_use =true;
    TriMesh primal_real;
    primal = primal_real;
    std::ifstream input(path_mesh);
    input >> primal;
    input.close();
}
} // namespace nest
#undef EQUAL
#undef AT
#define EQUAL EQUAL_TMP
#define AT AT_TMP

#endif //NEST_DISPLACEMENT_H
