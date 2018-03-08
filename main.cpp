#include <iostream>

#include <boost/optional.hpp>

#include <boost/range/irange.hpp>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/transform.hpp>
#include <boost/range/numeric.hpp>

#include "TClonesArray.h"

#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"

template<typename Object>
auto sqr(Object const& object)
{
    return object * object;
}

template<typename Element, template <typename, typename = std::allocator<Element>> class Container>
auto & operator<<(std::ostream& stream, Container<Element> const& container)
{
    for (auto const& element : boost::adaptors::index(container)) stream << '\n' << element.index() << ": " << element.value();
    return stream;
}

void print()
{
    std::cout << std::endl;
}

template<typename Object, typename ... Arguments>
void print(Object const& object, Arguments ... arguments)
{
    std::cout << std::boolalpha << object << ' ';
    print(arguments ...);
}

template<typename Integer>
auto range(Integer integer)
{
    return boost::irange(static_cast<Integer>(0), integer);
}

template<typename Container, typename Function>
auto transform(Container const& container, Function const& function)
{
    return container | boost::adaptors::transformed(function);
}

template<typename Particle>
auto transverse_distance(Particle const& particle)
{
    return std::sqrt(sqr(particle.X) + sqr(particle.Y));
}

auto particle_distance(TClonesArray const& particle_branch, int number)
{
    auto& particle = static_cast<GenParticle&>(*particle_branch.At(number));
    return std::abs(particle.PID) == 13 ? transverse_distance(particle) : 0;
}

auto muon_distance(TClonesArray const& muon_branch, int number)
{
    auto& muon = static_cast<Muon&>(*muon_branch.At(number));
    auto& particle = static_cast<GenParticle&>(*muon.Particle.GetObject());
    return transverse_distance(particle);
}

template<typename Function>
auto displacement(ExRootTreeReader& tree_reader, int entry, TClonesArray const& branch, Function const& function)
{
    tree_reader.ReadEntry(entry);
    std::vector<double> result;
    for (auto number : range(branch.GetEntriesFast())) function(number);
    return !result.empty();
}

auto analyse_events(ExRootTreeReader& tree_reader)
{
    auto& muon_branch = *tree_reader.UseBranch("Muon");
    return boost::accumulate(range(tree_reader.GetEntries()), 0., [&](auto & sum, auto entry) {
        return sum += displacement(tree_reader, entry, muon_branch, [&](auto number) {
            return particle_distance(muon_branch, number);
        }) ? 1 : 0;
    }) / tree_reader.GetEntries();
}

auto AnalyseEvents(ExRootTreeReader& tree_reader)
{
    auto& muon_branch = *tree_reader.UseBranch("Muon");
    auto& particle_branch = *tree_reader.UseBranch("Particle");
    // Loop over all events
    auto number = 0;
    for (auto entry : range(tree_reader.GetEntries())) {
        // Load selected branches with data from specified event
        tree_reader.ReadEntry(entry);
        std::vector<double> result;
        auto entries = muon_branch.GetEntriesFast();
        print("entries", entries);
        if (entries != 2) continue;
        for (auto number : range(entries)) {
            auto& muon = static_cast<Muon&>(*muon_branch.At(number));
//             print("got muon");
            auto& particle = static_cast<GenParticle&>(*muon.Particle.GetObject());
//             print("got particle");
            auto distance = transverse_distance(particle);
//                 if(distance > 0)
            if (distance > 0 && particle.D0 > 0) {
                print(distance, particle.D0);
//             if(distance > 0) print("dist", distance);
//             if (distance > 10 && distance < 200)
                result.emplace_back(distance);
            }
        }
//         if(!result.empty()) print("result:", result.size());
        if (!result.empty()) ++number;
    }
    print("displaced", number);
    return static_cast<double>(number) / tree_reader.GetEntries();
}

auto AnalyseEvents2(ExRootTreeReader& tree_reader)
{
    auto& particle_branch = *tree_reader.UseBranch("Particle");
    auto number_displaced = 0;
    // loop over all events
    for (auto entry : range(tree_reader.GetEntries())) {
        tree_reader.ReadEntry(entry);
        std::vector<double> result;
        auto number_muons = 0;
        // loop over all particles
        for (auto position : range(particle_branch.GetEntriesFast())) {
            auto& particle = static_cast<GenParticle&>(*particle_branch.At(position));
            if (std::abs(particle.PID) != 13) continue;
            print("ID, Status", particle.PID, particle.Status);
            ++number_muons;
            auto distance = transverse_distance(particle);
            if (distance > 10 && distance < 200) {
                result.emplace_back(distance);
            }
        }
//         if (number_muons != 2)
        print("number of muons", number_muons, result.size());
        if (!result.empty()) ++number_displaced;
    }
    print("displaced", number_displaced);
    return static_cast<double>(number_displaced) / tree_reader.GetEntries();
}

auto AnalyseEvents3(ExRootTreeReader& tree_reader)
{
    auto& particle_branch = *tree_reader.UseBranch("Particle");
    auto number_of_muons = 0;
    for (auto entry : range(tree_reader.GetEntries())) {
        tree_reader.ReadEntry(entry);
        std::vector<double> result;
        for (auto position : range(particle_branch.GetEntriesFast())) {
            auto& particle = static_cast<GenParticle&>(*particle_branch.At(position));
            auto distance = transverse_distance(particle);
            if (distance > 100 && distance < 2000)  result.emplace_back(distance);
        }
        if (!result.empty()) ++number_of_muons;
    }
    print("displaced", number_of_muons);
    return static_cast<double>(number_of_muons) / tree_reader.GetEntries();
}

struct File {
    File(std::string const& file_name) : chain("Delphes"), tree_reader(&chain)
    {
        chain.Add(file_name.c_str());
    }
    TChain chain;
    ExRootTreeReader tree_reader;
};

auto analyze(std::string const& file_name)
{
//     TChain chain("Delphes");
//     chain.Add(file_name.c_str());
//     ExRootTreeReader tree_reader(&chain);
    File file(file_name);
    return AnalyseEvents2(file.tree_reader);
}

auto file_name(int number)
{
    auto path = "~/scratch/2.6.2_heavyion/scan/Events/";
    auto name = (number < 10 ? "0" : "") + std::to_string(number);
    auto folder = "run_" + name + "_decayed_1/";
    auto file = "tag_1_delphes_events.root";
    auto result = path + folder + file;
    print(result);
    return result;
}

int main()
{
//     auto range = boost::irange(1, 49);
    auto range = boost::irange(1, 2);
    auto result = transform(range, [](auto number) {
        return analyze(file_name(number));
    });
    for (auto res : result) print(res, '\n');
}
