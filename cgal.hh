#pragma once

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Aff_transformation_3.h>

namespace hnl
{

namespace cgal
{

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Polyhedron = CGAL::Polyhedron_3<Kernel>;

using Primitive = CGAL::AABB_face_graph_triangle_primitive<Polyhedron>;
using Traits = CGAL::AABB_traits<Kernel, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;
using TriangleMeshSide = CGAL::Side_of_triangle_mesh<Polyhedron, Kernel>;
using Transformation = CGAL::Aff_transformation_3<Kernel>;

auto triangle_mesh_side(Polyhedron const& polyhedron)
{
    Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
    tree.accelerate_distance_queries();
    return TriangleMeshSide(tree);
}

bool is_inside(Point const& point, TriangleMeshSide const& triangle_mesh_side)
{
    return triangle_mesh_side(point) == CGAL::ON_BOUNDED_SIDE;
}

bool is_inside(Point const& point, Polyhedron const& polyhedron)
{
    return is_inside(point, triangle_mesh_side(polyhedron));
}

}

}
