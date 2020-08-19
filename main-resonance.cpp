
#include "ResonanceWidths.hh"
#include "pythia.hh"
#include "container.hh"
#include "id.hh"

using namespace hnl;

std::vector<std::string> get_table(std::function<double (int id_heavy, int id_light)> const& coupling, double mass, int id) {
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 14000.");
    pythia.readString("Bottomonium:all = on");
    pythia.setResonancePtr(new NeutrinoResonance(pythia, coupling, mass, heavy_neutrino));
    pythia.init();
    auto& particle = *pythia.particleData.particleDataEntryPtr(id);
    std::vector<std::string> vector({std::to_string(id) + ":all = nu void " + std::to_string(2) + " " + std::to_string(0) + " " + std::to_string(0) + " " + std::to_string(mass) + " " + std::to_string(particle.mWidth()) + " " + std::to_string(mass / 2) + " " + std::to_string(mass * 2) + " " + std::to_string(particle.tau0())});
    for (auto i : irange(particle.sizeChannels())) {
        auto c = particle.channel(i);
        vector.emplace_back(std::to_string(id) + ":addChannel = " + std::to_string(1) + " " + std::to_string(c.bRatio()) + " " + std::to_string(101) + " " + std::to_string(c.product(0)) + " " + std::to_string(c.product(1)) + " " + std::to_string(c.product(2)));
    }
    return vector;
}

int main() {
    print(get_table(neutrino_coupling(1), 2., heavy_neutrino));
}
