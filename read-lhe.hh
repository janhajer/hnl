#pragma once

#include "pythia.hh"
#include "geometry.hh"

namespace hnl {

auto to_cgal(Pythia8::Particle const& particle) -> cgal::Point {
    return {particle.xProd() / 1000, particle.yProd() / 1000, particle.zProd() / 1000}; // convert from mm to m
}

void read_lhe() {
    Pythia8::Pythia pythia;

    pythia.readString("Beams:frameType = 4");
    pythia.readString("Beams:LHEF = ttbar.lhe");
    pythia.init();

    auto analysis = mapp::analysis();
    for (auto event_number = 0; ; ++event_number) {
        if (!pythia.next()) {
            if(pythia.info.atEndOfFile()) break;
            continue;
        }

        for (auto line = 0; line < pythia.event.size(); ++line) {
            auto const& particle = pythia.event[line];
            if (!particle.isFinal() || !particle.isCharged()) continue;
            if (!analysis.is_inside(to_cgal(particle))) continue;
            if (debug) print("Hooray!", particle.vProd().pAbs());
            ++good;
            break;
        }
    }
    pythia.stat();
}

}
