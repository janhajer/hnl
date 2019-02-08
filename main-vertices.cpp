#include "root-plot.hh"
#include "scan-hep-plot.hh"
#include "cms-14.hh"

using namespace neutrino;

int main(int argc, char** argv)
{
    auto function = [](auto & reader, auto predicate) noexcept {
        return accumulate(reader, xy_histogram({-3.,3.}, {-2.,2.}), [&](auto const& event) noexcept {
            return predicate(event.particles(), event.vertices());
        });
    };
    scan_vertices(arguments(argc, argv), analysis_cms_14(), function);
}
