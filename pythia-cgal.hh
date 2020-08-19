#include "pythia.hh"
#include "cgal.hh"

namespace hnl {

inline auto to_cgal(Pythia8::Particle const& particle) -> cgal::Point {
    return {particle.xProd() / 1000, particle.yProd() / 1000, particle.zProd() / 1000}; // convert from mm to m
}

}
