#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/algorithm/min_element.hpp>
#include <iostream>

#define CGAL_USE_BASIC_VIEWER
#include "geometry.hh"

#include "draw_polyhedron.h"
#include "ui_ImageInterface.h"

#include "hepmc.hh"

#include "generic.hh"


namespace neutrino{

namespace gcal{

using Viewer = CGAL::SimplePolyhedronViewerQt<Polyhedron, CGAL::DefaultColorFunctorPolyhedron>;

Viewer& add_pipe(Viewer& viewer)
{
    Point pipe_one(0, 0, 10);
    Point pipe_two(0, 0, -65);
    viewer.add_segment(pipe_one, pipe_two);
    return viewer;
}

Viewer& draw(Viewer& viewer, std::vector<Point> const& points)
{
    Point origin(0, 0, 0);
    for (auto const& point : points) viewer.add_segment(origin, point, CGAL::red());
    return viewer;
}

void execute(QApplication& application, Viewer& viewer)
{
    viewer.show();
    application.exec();
}

hep::FourVector four_vector(Point const& point){
    return {point.x(), point.y(), point.z(), 0};
}

void print_max_eta_phi(){
    auto poly_points = get_points();
        auto min_eta = four_vector(poly_points.front()).eta();
    auto max_eta = four_vector(poly_points.front()).eta();
    auto min_phi = four_vector(poly_points.front()).phi();
    auto max_phi = four_vector(poly_points.front()).phi();

    for (auto const& point : poly_points) {
        auto fv = four_vector(point);
        auto eta = fv.eta();
        min_eta = std::min(eta, min_eta);
        max_eta = std::max(eta, max_eta);
        auto phi = fv.phi();
        min_phi = std::min(phi, min_phi);
        max_phi = std::max(phi, max_phi);
        print("eta", eta, "phi", phi);
    }

    print("eta", min_eta, max_eta);
    print("phi", min_phi, max_phi);

    auto res1 = boost::range::max_element(poly_points, [](auto const& one, auto const& two){
        return four_vector(one).eta() < four_vector(two).eta();
    });

    auto res2 = boost::range::min_element(poly_points, [](auto const& one, auto const& two){
        return four_vector(one).eta() < four_vector(two).eta();
    });

    print(four_vector(*res1), four_vector(*res2));
    print(four_vector(*res1).eta(), four_vector(*res2).eta());

    auto res3 = boost::range::min_element(poly_points, [](auto const& one, auto const& two){
        return four_vector(one).phi() < four_vector(two).phi();
    });

    auto res4 = boost::range::max_element(poly_points, [](auto const& one, auto const& two){
        return four_vector(one).phi() < four_vector(two).phi();
    });

    print(four_vector(*res3), four_vector(*res4));
    print(four_vector(*res3).phi(), four_vector(*res4).phi());
}


}

}

int main(int argc, char** argv)
{

    using namespace neutrino;
    using namespace gcal;

    print_max_eta_phi();
    auto polyhedron = get_polyhedron();
    std::vector<Point> test_points;
    test_points.emplace_back(15, -1, -40);
    test_points.emplace_back(10, -1, -40);

    for (auto const& test_point : test_points) print("Is inside:", is_inside(test_point, polyhedron));

    QApplication application(argc, argv);
    Viewer viewer(application.activeWindow(), polyhedron);
    add_pipe(viewer);
    draw(viewer, test_points);
    execute(application, viewer);
    return 0;
}