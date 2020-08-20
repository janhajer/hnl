#pragma once

#include "Pythia8Plugins/HepMC2.h"

namespace hnl {

namespace hep {

using FourVector = HepMC::FourVector;
using Particle = HepMC::GenParticle;

std::ostream& operator<<(std::ostream& stream, FourVector const& vector) {
    return stream << "(" << vector.t() << ", " << vector.x() << ", " << vector.y() << ", " << vector.z() << ")";
}

}

Pythia8::Vec4 to_pythia(HepMC::FourVector const& vector) {
    return {vector.x(), vector.y(), vector.z(), vector.t()};
}

}

