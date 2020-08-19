#pragma once

#include "cgal.hh"
#include "container.hh"

namespace hnl
{

namespace mapp
{

using namespace cgal;

auto get_points()
{
    std::vector<Point> points;
    points.emplace_back(3.27, -2, 52.83); //0
    points.emplace_back(4, -2, 61.39); //1
    points.emplace_back(4, 1, 61.39); //2
    points.emplace_back(3.27, 1, 52.83); //3
    points.emplace_back(12.24, -2, 33.63); //4
    points.emplace_back(16.53, -2, 35.45); //5
    points.emplace_back(16.53, 1, 35.45); //6
    points.emplace_back(12.24, 1, 33.63); //7
    return points;
}

void make_face(Polyhedron& polyhedron, Point const& one, Point const& two, Point const& three, Point const& four)
{
    polyhedron.make_triangle(one, two, three);
    polyhedron.make_triangle(one, four, three);
}

auto get_polyhedron(std::vector<Point> const& points)
{
    Polyhedron polyhedron;
    make_face(polyhedron, points[0], points[1], points[2], points[3]);
    make_face(polyhedron, points[4], points[5], points[1], points[0]);
    make_face(polyhedron, points[7], points[6], points[5], points[4]);
    make_face(polyhedron, points[5], points[6], points[2], points[1]);
    make_face(polyhedron, points[6], points[7], points[3], points[2]);
    make_face(polyhedron, points[7], points[4], points[0], points[3]);
    return polyhedron;
}

auto get_polyhedron()
{
    auto poly_points = get_points();
    return get_polyhedron(poly_points);
}

Transformation rotation_z(double angle)
{
    auto cosa = std::cos(angle);
    auto sina = std::sin(angle);
    return {cosa, -sina, 0, sina, cosa, 0, 0, 0, 1};
}

Transformation reflection_z()
{
    auto cosa = std::cos(M_PI);
    auto sina = std::sin(M_PI);
    return {1, 0, 0, 0, cosa, -sina, 0, sina, cosa};
}

auto rotated_poly(double angle, bool reflect)
{
    auto rotation = rotation_z(angle);
    auto reflection = reflect ? reflection_z() : Transformation(CGAL::Identity_transformation());
    std::vector<Point> transformed = transform(get_points(), [&rotation, &reflection](Point const & point) -> Point {
        return reflection(rotation(point));
    });
    return get_polyhedron(transformed);
}

auto get_polyhedrons()
{
    std::vector<Polyhedron> polyhedrons;
    for (int i = 0; i < 8; ++i) {
        polyhedrons.emplace_back(rotated_poly(M_PI_4 * i, false));
        polyhedrons.emplace_back(rotated_poly(M_PI_4 * i, true));
    }
    return polyhedrons;
}

struct Analysis {
    Analysis(std::vector<cgal::Polyhedron> const& polyhedrons_) :
        name("MAPP"),
        polyhedrons(polyhedrons_),
        trees(get_trees(polyhedrons))
    {
        for (auto& tree : trees) tree.accelerate_distance_queries();
        detectors = get_detectors(trees);
    }
    std::deque<cgal::Tree> get_trees(std::vector<cgal::Polyhedron> const& polyhedrons_) const
    {
        std::deque<cgal::Tree> trees_;
        for (auto const& polyhedron : polyhedrons_) trees_.emplace_back(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
        return trees_;
    }
    std::vector<cgal::TriangleMeshSide> get_detectors(std::deque<cgal::Tree> const& trees_) const
    {
        return transform(trees_, [](cgal::Tree const & tree) -> cgal::TriangleMeshSide {
            return cgal::TriangleMeshSide{tree};
        });
    }
    bool is_inside(Point const& point) const {
        return boost::algorithm::any_of(detectors, [&point](auto const& detector) {
            return mapp::is_inside(point, detector);
        });
    }
    std::string name;
    std::vector<cgal::Polyhedron> polyhedrons;
    std::deque<cgal::Tree> trees;
    std::vector<cgal::TriangleMeshSide> detectors;
};

auto analysis() -> Analysis
{
    return {get_polyhedrons()};
}

}

}
