
#include "mapp.hh"
#define CGAL_USE_BASIC_VIEWER
#include "ui_ImageInterface.h"
#include "draw_polyhedron.h"

namespace hnl {

namespace cgal {

using Viewer = CGAL::SimplePolyhedronViewerQt<Polyhedron, CGAL::DefaultColorFunctorPolyhedron>;

Viewer& add_pipe(Viewer& viewer) {
    Point pipe_one(0, 0, -65);
    Point pipe_two(0, 0, 65);
    viewer.add_segment(pipe_one, pipe_two);
    return viewer;
}

Viewer& draw(Viewer& viewer, std::vector<Point> const& points) {
    Point origin(0, 0, 0);
    for (auto const& point : points) viewer.add_segment(origin, point, CGAL::red());
    return viewer;
}

void execute(QApplication& application, Viewer& viewer) {
    viewer.show();
    application.exec();
}

}

}

int main(int argc, char** argv) {
    using namespace hnl;
    using namespace mapp;
    using namespace cgal;
    QApplication application(argc, argv);
    Viewer viewer(application.activeWindow());
//     for (auto const& polyhedron : get_polyhedrons()) viewer.add(polyhedron);
    auto polyhedron = get_polyhedron();
    viewer.add(polyhedron);
//     add_pipe(viewer);
  // draw(viewer, test_points);
    execute(application, viewer);
    return 0;
}
