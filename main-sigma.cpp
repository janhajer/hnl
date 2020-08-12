// Simple illustration how to provide (a) your own resonance-width class, and (b) your own cross-section class, with instances handed in to Pythia.
// The hypothetical scenario is that top would have been so long-lived that a toponium resonance Theta could form.
// Then production could proceed via q qbar -> gamma*/Z* -> Theta, with decay either to a fermion pair or (dominantly) to three gluons.
// The implementation is not physically correct in any number of ways, but should exemplify the strategy needed for realistic cases.

#include "generic.hh"
#include "Pythia8/Pythia.h"
#include "Sigma.hh"
#include "ResonanceWidths.hh"

namespace hnl {

const int heavy_neutrino = 9900012;

auto neutrino_coupling = [](double factor) {
    return [factor](int id_heavy, int id_light) -> double {
        if (id_light != 12 && id_light != 14 && id_light != 16) {
            hnl::print(id_light, "is not a light neutrino");
            return 0;
        }
        if (id_heavy != 9900012 && id_heavy != 9900014 && id_heavy != 9900016) {
            hnl::print(id_heavy, "is not a heavy neutrino");
            return 0;
        }
        auto id = id_heavy - 9900000 - id_light;
        return id == 0 ? factor : 0;
        return id == 0 || std::abs(id) == 2 || std::abs(id) == 4 ? factor : 0;
    };
};

double gen_point(double mass) {
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);

    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 14000.");
    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Init:showProcesses = off");
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showMultipartonInteractions = off");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.particleData.m0(heavy_neutrino, mass);
    pythia.particleData.mMin(heavy_neutrino, mass / 2);
    pythia.particleData.mMax(heavy_neutrino, mass * 2);
    pythia.setResonancePtr(new NeutrinoResonance(pythia, neutrino_coupling(1), heavy_neutrino));
    Pythia8::SigmaProcess* sigma = new Sigma(heavy_neutrino, 1.);
    pythia.setSigmaPtr(sigma);

    pythia.init();
//   Pythia8::Hist mTheta("Theta mass", 100, 200., 300.);
    for (int event_number = 0; event_number < 20; ++event_number) {
        if (!pythia.next()) print("Error in event", event_number);
//     mTheta.fill( pythia.process[5].m() );
//     for (int i = 0; i <= 5; ++i) print(pythia.process[i].m());
    }
//     pythia.stat();
//   std::cout << mTheta;
    delete sigma;
    return pythia.info.sigmaGen();
}

}

int main() {

    using namespace hnl;

    std::vector<double> vector;
    for (int i = 1; i <= 10; ++i) {
        double mass = 25 * i;
        auto res = gen_point(mass);
        vector.emplace_back(res);
        print(mass, res);
    }
    for (auto i : vector) print(i);
}
