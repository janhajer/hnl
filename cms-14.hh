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

auto electron_efficiency_cms_14(Property const& vector) noexcept {
    auto pt = neutrino::pt(vector);
    if (pt <= 0.2_GeV) return 0.;
    auto eta = abs(neutrino::eta(vector));
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

auto muon_efficiency_cms_14(Property const& vector) noexcept {
    auto pt = neutrino::pt(vector);
    if (pt <= 0.2_GeV) return 0.;
    auto eta = abs(neutrino::eta(vector));
//     if (eta <= 1.2_rad) return pt <= 1_GeV ? value(pt) : 1.;
    if (eta <= 2.8_rad) return pt <= 1_GeV ? value(pt) : 1.;
    if (eta <= 4.0_rad) return pt <= 1_GeV ? value(pt) * 0.95 : 0.95;
    return 0.;
}

auto jet_efficiency_cms_14(Property const& vector) noexcept {
    auto pt = neutrino::pt(vector);
    if (pt <= 0.2_GeV) return 0.;
    auto eta = abs(neutrino::eta(vector));
    if (eta <= 1.2_rad) return pt <= 1_GeV ? value(pt) * 0.96 : 0.97;
    if (eta <= 2.5_rad) return pt <= 1_GeV ? value(pt) * 0.85 : 0.87;
    if (eta <= 4_rad) return pt <= 1_GeV ? value(pt) * 0.8 : 0.82;
    return 0.;
}

auto efficiency_cms_14(Property const& vector) noexcept {
    switch (vector.id)
    {
    case Id::electron : return electron_efficiency_cms_14(vector);
    case Id::muon : return muon_efficiency_cms_14(vector);
    case Id::tau : return jet_efficiency_cms_14(vector);
    default : return 0.;
    }
}

auto analysis_cms_14() noexcept -> Analysis {
    auto second_pt = 15_GeV;
    auto vertex_pt = 5_GeV;
    return {"CMS", cms_14_tracker(), cms_14_chamber(), cms_14_trigger(), efficiency_cms_14, second_pt, vertex_pt};
}

}
