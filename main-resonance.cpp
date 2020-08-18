
#include "ResonanceWidths.hh"
#include "pythia.hh"
#include "container.hh"
#include "id.hh"

int main() {
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 14000.");
    pythia.readString("Bottomonium:all = on");
    using namespace hnl;
    pythia.init();
    NeutrinoResonance res(pythia, neutrino_coupling(1), 2., heavy_neutrino);
    for (auto channel : irange(pythia.particleData.particleDataEntryPtr(heavy_neutrino)->sizeChannels())) print(res.width(channel, 2.));
}
