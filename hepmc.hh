#pragma once

#include <ostream>

#include "HepMC3/Attribute.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/WriterAscii.h"
#include "Pythia8Plugins/HepMC3.h"

namespace hnl {

namespace hep {

using FourVector = HepMC3::FourVector;
using Particle = HepMC3::GenParticle;

inline std::ostream& operator<<(std::ostream& stream, FourVector const& vector) {
    return stream << "(" << vector.t() << ", " << vector.x() << ", " << vector.y() << ", " << vector.z() << ")";
}

}

inline Pythia8::Vec4 to_pythia(HepMC3::FourVector const& vector) {
    return {vector.x(), vector.y(), vector.z(), vector.t()};
}

}
