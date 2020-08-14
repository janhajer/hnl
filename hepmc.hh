#pragma once

#include "Pythia8Plugins/HepMC2.h"

namespace hnl {

Pythia8::Vec4 to_pythia(HepMC::FourVector const& vector) {
    return {vector.x(), vector.y(), vector.z(), vector.t()};
}

}

