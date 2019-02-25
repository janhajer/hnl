#pragma once

#include "detector.hh"

namespace neutrino
{

auto cms_14_tracker() noexcept {
    auto detector = tracker();
    detector.transversal = {5_mm, 1100_mm};
    detector.longitudinal = {10_cm, 2800_mm};
    detector.eta = { -4_rad, 4_rad};
    return detector;
}

auto cms_14_chamber() noexcept {
    auto detector = chamber();
    detector.transversal = {4_m, 7_m};
    detector.longitudinal = {7_m, 11_m};
    detector.eta = { -4_rad, 4_rad};
    return detector;
}

auto cms_14_trigger() noexcept -> std::map<Id, Energy> {
    return {{Id::electron, 24_GeV}, {Id::muon, 24_GeV}, {Id::tau, 75_GeV}};
}

auto cms_14_trigger_2() noexcept -> PairTrigger { // FIXME ATLAS values
    return {
        {{Id::electron, Id::electron}, {{18_GeV, 18_GeV}, boost::none}},
        {{Id::electron, Id::muon}, {{8_GeV, 25_GeV}, std::make_pair(18_GeV, 15_GeV)}},
        {{Id::electron, Id::tau}, {{18_GeV, 30_GeV}, boost::none}},
        {{Id::muon, Id::muon}, {{15_GeV, 15_GeV}, std::make_pair(23_GeV, 9_GeV)}},
        {{Id::muon, Id::tau}, {{15_GeV, 30_GeV}, boost::none}},
        {{Id::tau, Id::tau}, {{40_GeV, 30_GeV}, boost::none}}
    };
}

auto cms_14_track() noexcept -> std::map<Id, Energy> {
    return {{Id::electron, 5_GeV}, {Id::muon, 5_GeV}, {Id::tau, 5_GeV}};
}

template<typename Particle>
auto electron_efficiency_cms_14(Particle const& particle) noexcept {
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
auto muon_efficiency_cms_14(Particle const& particle) noexcept {
    auto pt = neutrino::pt(particle);
    if (pt <= 0.2_GeV) return 0.;
    auto eta = abs(neutrino::eta(particle));
//     if (eta <= 1.2_rad) return pt <= 1_GeV ? value(pt) : 1.;
    if (eta <= 2.8_rad) return pt <= 1_GeV ? value(pt) : 1.;
    if (eta <= 4.0_rad) return pt <= 1_GeV ? value(pt) * 0.95 : 0.95;
    return 0.;
}

template<typename Particle>
auto jet_efficiency_cms_14(Particle const& particle) noexcept {
    auto pt = neutrino::pt(particle);
    if (pt <= 0.2_GeV) return 0.;
    auto eta = abs(neutrino::eta(particle));
    if (eta <= 1.2_rad) return pt <= 1_GeV ? value(pt) * 0.96 : 0.97;
    if (eta <= 2.5_rad) return pt <= 1_GeV ? value(pt) * 0.85 : 0.87;
    if (eta <= 4_rad) return pt <= 1_GeV ? value(pt) * 0.8 : 0.82;
    return 0.;
}

auto efficiency_cms_14(Particle const& particle) noexcept {
    switch (abs_id(particle))
    {
    case to_underlying(Id::electron) : return electron_efficiency_cms_14(particle);
    case to_underlying(Id::muon) : return muon_efficiency_cms_14(particle);
    default : return jet_efficiency_cms_14(particle);
    }
}

auto analysis_cms_14() noexcept -> Analysis {
    return {"CMS", cms_14_tracker(), cms_14_chamber(), cms_14_trigger(), cms_14_trigger_2(), cms_14_track(), efficiency_cms_14};
}

}
