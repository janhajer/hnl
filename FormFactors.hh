#pragma once

#include <array>

namespace hnl {

// Form factors of the charged-current transition of a pseudoscalar meson into a lighter meson, in the conventions of
// arXiv:1805.08567 (appendix A). Sources: FormFactors.cpp and the physics inputs table in Readme.md.

// {f_+, f_0} for a pseudoscalar daughter; m_from and m_to are the meson masses, used only for B_c -> eta_c
std::array<double, 2> pseudoscalar_form_factors(int id_from, int id_to, double q2, double m_from, double m_to);

// {V, A_0, A_1, A_2} for a vector daughter. A_0 is shifted by the constant A_0_shift so that A_0(0) = A_3(0) =
// ((m_from + m_to) A_1(0) - (m_from - m_to) A_2(0)) / (2 m_to) holds exactly with these masses: the fits satisfy it
// only to their precision, and a_- of 1805.08567 contains (A_0 - A_3) / q^2.
std::array<double, 4> vector_form_factors(int id_from, int id_to, double q2, double m_from, double m_to);

double A_0_shift(int id_from, int id_to, double m_from, double m_to);

}
