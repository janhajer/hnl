// Form factors against values published with the fits, and the kinematic constraints at q^2 = 0.

#include <cmath>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "FormFactors.hh"
#include "math.hh"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using hnl::sqr;

namespace {

auto f_plus(int from, int to, double q2) {
    return hnl::pseudoscalar_form_factors(from, to, q2, 0., 0.)[0];
}

auto f_0(int from, int to, double q2) {
    return hnl::pseudoscalar_form_factors(from, to, q2, 0., 0.)[1];
}

}

TEST_CASE("D -> pi reproduces FNAL/MILC 22, table XIV", "[form factors]") {
    auto q2_max = sqr(1.86483 - 0.1349768);
    CHECK_THAT(f_plus(421, 211, 0.), WithinRel(0.6300, 1e-4));
    CHECK_THAT(f_plus(421, 211, q2_max), WithinRel(3.119, 5e-3));
    CHECK_THAT(f_0(421, 211, q2_max), WithinRel(1.2783, 5e-3));
}

TEST_CASE("D -> K reproduces the FLAG 2024 f_+(0), eq. (138)", "[form factors]") {
    CHECK_THAT(f_plus(421, 321, 0.), WithinRel(0.7430, 3e-3)); // within its uncertainty 0.0027
}

TEST_CASE("B_s -> D_s reproduces HPQCD 19, appendix A", "[form factors]") {
    auto q2_max = sqr(5.3669 - 1.9683);
    CHECK_THAT(f_plus(531, 431, 0.), WithinRel(0.666, 1e-4));
    CHECK_THAT(f_plus(531, 431, q2_max), WithinRel(1.209, 1e-2));
    CHECK_THAT(f_0(531, 431, q2_max), WithinRel(0.917, 1e-2));
}

TEST_CASE("pseudoscalar form factors satisfy f_+(0) = f_0(0)", "[form factors]") {
    for (auto [from, to] : {std::pair{130, 211}, {321, 111}, {411, 311}, {421, 321}, {411, 111}, {421, 211}, {431, 221},
                            {511, 211}, {521, 111}, {531, 321}, {511, 411}, {521, 421}, {531, 431}}) {
        INFO(from << " -> " << to);
        CHECK_THAT(f_plus(from, to, 0.), WithinRel(f_0(from, to, 0.), 1e-10));
    }
}

TEST_CASE("B -> pi f_0 has no B* pole", "[form factors]") {
    // the BCL fit of FLAG 2024 (appendix B.3.2) has no pole in f_0; with the former pole f_0(q2_max) was 5.7
    auto q2_max = sqr(5.27934 - 0.1349768);
    CHECK(f_0(511, 211, q2_max) < 1.5);
}

TEST_CASE("B -> D* reproduces F(1) of FNAL/MILC 21, eq. (216) of FLAG 2024", "[form factors]") {
    auto m_B = 5.280, m_V = 2.010, r = m_V / m_B; // the masses of the BGL fit
    auto A_1 = hnl::vector_form_factors(511, 413, sqr(m_B - m_V), m_B, m_V)[2];
    auto h_A1 = A_1 * (1. + r) / (2. * std::sqrt(r)); // A_1(w = 1) = R* h_A1(1), R* = 2 sqrt(r) / (1 + r)
    CHECK_THAT(h_A1, WithinRel(0.909, 2e-2)); // within its uncertainty 0.017
}

TEST_CASE("vector form factors satisfy A_0(0) = A_3(0) with a small shift", "[form factors]") {
    struct Transition {
        int from, to;
        double m_from, m_to;
    };
    for (auto t : {Transition{411, 313, 1.86966, 0.89555}, {511, 413, 5.27966, 2.01026}, {511, 213, 5.27966, 0.77526},
                   {531, 433, 5.36692, 2.1122}, {531, 323, 5.36692, 0.89167}, {541, 443, 6.27447, 3.0969}}) {
        INFO(t.from << " -> " << t.to);
        auto ff = hnl::vector_form_factors(t.from, t.to, 0., t.m_from, t.m_to);
        auto A_3 = ((t.m_from + t.m_to) * ff[2] - (t.m_from - t.m_to) * ff[3]) / (2. * t.m_to);
        CHECK_THAT(ff[1], WithinRel(A_3, 1e-10));
        // the fits satisfy the constraint to their precision: the shift is a few % of A_0 at most
        CHECK(std::abs(hnl::A_0_shift(t.from, t.to, t.m_from, t.m_to)) < 0.04 * ff[1]);
    }
}
