#pragma once

#include "Pythia8/Pythia.h"

namespace hnl
{

struct NeutrinoResonanceWidth : public Pythia8::ResonanceWidths {
    NeutrinoResonanceWidth(int idResIn) : thetaWRat(), mW()
    {
        initBasic(idResIn);
    }
private:
    double thetaWRat;
    double mW;
    virtual void initConstants(); // Initialize constants.
    virtual void calcPreFac(bool = false); // Calculate various common prefactors for the current mass.
    virtual void calcWidth(bool = false); // Caclulate width for currently considered channel.
};

void NeutrinoResonanceWidth::initConstants()   // Locally stored properties and couplings: righthanded W mass.
{
    thetaWRat = 1. / (768. * M_PI * Pythia8::pow2(couplingsPtr->sin2thetaW()));
    mW = particleDataPtr->m0(24);

}

void NeutrinoResonanceWidth::calcPreFac(bool) // Calculate various common prefactors for the current mass.
{
    // Common coupling factors.
    alpEM = couplingsPtr->alphaEM(mHat * mHat);
    alpS = couplingsPtr->alphaS(mHat * mHat);
    colQ = 3. * (1. + alpS / M_PI);
    preFac = Pythia8::pow2(alpEM) * thetaWRat * Pythia8::pow5(mHat) / Pythia8::pow4(std::max(mHat, mW));
}

void NeutrinoResonanceWidth::calcWidth(bool) // Calculate width for currently considered channel.
{
    // Check that above threshold.
    if (mHat < mf1 + mf2 + mf3 + MASSMARGIN) return;
// Coupling part of widths to l- q qbar', l- l'+ nu_lR' and c.c.
    widNow = (id2Abs < 9 && id3Abs < 9) ? preFac * colQ * couplingsPtr->V2CKMid(id2, id3) : preFac;
// Phase space corrections in decay. Must have y < 1.
    auto x = (mf1 + mf2 + mf3) / mHat;
    auto x2 = x * x;
    auto fx = 1. - 8. * x2 + 8. * Pythia8::pow3(x2) - Pythia8::pow4(x2) - 24. * Pythia8::pow2(x2) * std::log(x);
    auto y = std::min(0.999, Pythia8::pow2(mHat / mW));
    auto fy = (12. * (1. - y) * log(1. - y) + 12. * y - 6. * y * y - 2.* Pythia8::pow3(y)) / Pythia8::pow4(y);
    widNow *= fx * fy;
}

}
