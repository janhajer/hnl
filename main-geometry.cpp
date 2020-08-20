#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/algorithm/min_element.hpp>

#include "HepMC/SimpleVector.h"
#include "HepMC/GenParticle.h"

#include "mapp.hh"
#include "io.hh"
#include "hepmc.hh"

namespace hnl {

namespace cgal {

hep::FourVector four_vector(Point const& point) {
    return {point.x(), point.y(), point.z(), 0};
}

void print_max_eta_phi()
{
    auto poly_points = mapp::get_points();
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

    auto res1 = boost::range::max_element(poly_points, [](auto const & one, auto const & two) {
        return four_vector(one).eta() < four_vector(two).eta();
    });

    auto res2 = boost::range::min_element(poly_points, [](auto const & one, auto const & two) {
        return four_vector(one).eta() < four_vector(two).eta();
    });

//     print(four_vector(*res1), four_vector(*res2));
    print(four_vector(*res1).eta(), four_vector(*res2).eta());

    auto res3 = boost::range::min_element(poly_points, [](auto const & one, auto const & two) {
        return four_vector(one).phi() < four_vector(two).phi();
    });

    auto res4 = boost::range::max_element(poly_points, [](auto const & one, auto const & two) {
        return four_vector(one).phi() < four_vector(two).phi();
    });

//     print(four_vector(*res3), four_vector(*res4));
    print(four_vector(*res3).phi(), four_vector(*res4).phi());
}


}

}

int main() {

   using namespace hnl;
   using namespace cgal;

    auto analysis = mapp::analysis();
// hep::FourVector vector(0,0,0,0);
// hep::Particle particle(vector);
    cgal::Point point(0, 0, 0);
    auto r = analysis.is_inside(point);
    print("inside", r);
    cgal::print_max_eta_phi();
}
