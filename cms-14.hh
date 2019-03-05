#pragma once

#include "detector.hh"

namespace neutrino
{

namespace cms14
{

auto tracker() noexcept {
    auto detector = neutrino::tracker();
    detector.transversal = {5_mm, 1100_mm};
    detector.longitudinal = {10_cm, 2800_mm};
    detector.eta = Range<Angle>{ -4._rad, 4._rad};
    return detector;
}

auto chamber() noexcept {
    auto detector = neutrino::chamber();
    detector.transversal = {4_m, 7_m};
    detector.longitudinal = {7_m, 11_m};
    detector.eta = { -4._rad, 4._rad};
    return detector;
}

auto trigger() noexcept -> std::map<Id, Energy> {
    return {{Id::electron, 24_GeV}, {Id::muon, 24_GeV}, {Id::tau, 75_GeV}};
}

auto trigger_2() noexcept -> PairTrigger { // FIXME ATLAS values
    return {
        {{Id::electron, Id::electron}, {{18_GeV, 18_GeV}, boost::none}},
        {{Id::electron, Id::muon}, {{8_GeV, 25_GeV}, std::make_pair(18_GeV, 15_GeV)}},
        {{Id::electron, Id::tau}, {{18_GeV, 30_GeV}, boost::none}},
        {{Id::muon, Id::muon}, {{15_GeV, 15_GeV}, std::make_pair(23_GeV, 9_GeV)}},
        {{Id::muon, Id::tau}, {{15_GeV, 30_GeV}, boost::none}},
        {{Id::tau, Id::tau}, {{40_GeV, 30_GeV}, boost::none}}
    };
}

auto track() noexcept -> std::map<Id, Energy> {
    return {{Id::electron, 5_GeV}, {Id::muon, 5_GeV}, {Id::tau, 5_GeV}};
}

template<typename Particle>
auto electron_efficiency(Particle const& particle) noexcept {
    auto pt = neutrino::pt(particle);
    if (pt <= 0.2_GeV) return 0.;
    auto eta = abs(neutrino::eta(particle));
    if (eta <= 1.2_rad) return pt <= 1_GeV ? value(pt) * 0.96 : 0.97;
    if (eta <= 2.5_rad)
    {
        if (pt <= 1.0_GeV) return value(pt) * 0.85;
        if (pt <= 10_GeV) return 0.82 + value(pt) * 0.01;
        return 0.90;
    }
    if (eta <= 4_rad)
    {
        if (pt <= 1.0_GeV) return value(pt) * 0.8;
        if (pt <= 10_GeV) return 0.8 + value(pt) * 0.01;
        return 0.85;
    }
    return 0.;
}

template<typename Particle>
auto muon_efficiency(Particle const& particle) noexcept {
    auto pt = neutrino::pt(particle);
    if (pt <= 0.2_GeV) return 0.;
    auto eta = abs(neutrino::eta(particle));
//     if (eta <= 1.2_rad) return pt <= 1_GeV ? value(pt) : 1.;
    if (eta <= 2.8_rad) return pt <= 1_GeV ? value(pt) : 1.;
    if (eta <= 4.0_rad) return pt <= 1_GeV ? value(pt) * 0.95 : 0.95;
    return 0.;
}

template<typename Particle>
auto jet_efficiency(Particle const& particle) noexcept {
    auto pt = neutrino::pt(particle);
    if (pt <= 0.2_GeV) return 0.;
    auto eta = abs(neutrino::eta(particle));
    if (eta <= 1.2_rad) return pt <= 1_GeV ? value(pt) * 0.96 : 0.97;
    if (eta <= 2.5_rad) return pt <= 1_GeV ? value(pt) * 0.85 : 0.87;
    if (eta <= 4_rad) return pt <= 1_GeV ? value(pt) * 0.8 : 0.82;
    return 0.;
}

auto efficiency(Particle const& particle) noexcept {
    switch (abs_id(particle))
    {
    case to_underlying(Id::electron) : return electron_efficiency(particle);
    case to_underlying(Id::muon) : return muon_efficiency(particle);
    default : return jet_efficiency(particle);
    }
}

// auto displaced_efficiency(Detector const& detector, hep::Particle const& particle) noexcept {
//     if (before(particle, detector)) return 1.;
//     auto dist = remaining(box(detector), ray(particle));
//     return dist > 0_m ? drop_off(dist, detector) : 0.;
// }

auto analysis() noexcept -> Analysis {
    return {"CMS", tracker(), chamber(), trigger(), trigger_2(), track(), efficiency};
}

}

}
