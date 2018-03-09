#include <iostream>
#include <fstream>

#include <boost/optional.hpp>

#include <boost/algorithm/string.hpp>

#include <boost/range/irange.hpp>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/transform.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/numeric.hpp>

#include "TClonesArray.h"

#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"


    using namespace std::string_literals;

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

std::ostream& operator<<(std::ostream& stream, GenParticle const& particle)
{
    return stream << "(" << particle.PID << ", " << particle.Status << ")";
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
    return boost::accumulate(range(branch.GetEntriesFast()), 0, [&](auto & sum, auto number) {
        return sum += function(number) > 10 ? 1 : 0 ;
    }) > 0;
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
    auto displaced_number = 0;
    for (auto entry : range(tree_reader.GetEntries())) {
        // Load selected branches with data from specified event
        tree_reader.ReadEntry(entry);
        std::vector<double> result;
        auto entries = muon_branch.GetEntriesFast();
        for (auto position : range(entries)) {
            auto& muon = static_cast<Muon&>(*muon_branch.At(position));
            auto& particle = static_cast<GenParticle&>(*muon.Particle.GetObject());
            auto distance = transverse_distance(particle);
//             if (distance > 10 && distance < 200){
            if (distance > 100) {
//                 print(distance);
                result.emplace_back(distance);
            }
        }
        if (!result.empty()) ++displaced_number;
    }
//     print("displaced", displaced_number);
    return static_cast<double>(displaced_number) / tree_reader.GetEntries();
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
            auto distance = transverse_distance(particle);
            if (particle.Status == 1) print("ID, Status", particle, distance);
            ++number_muons;
//             if (distance > 0) {
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

auto AnalyseEvents4(ExRootTreeReader& tree_reader)
{
    auto& particle_branch = *tree_reader.UseBranch("Particle");
    for (auto entry : range(tree_reader.GetEntries())) {
        tree_reader.ReadEntry(entry);
        std::vector<double> result;
        for (auto position : range(particle_branch.GetEntriesFast())) {
            auto& particle = static_cast<GenParticle&>(*particle_branch.At(position));
            (particle.PID == 9900012 || particle.PID == 9900014 || particle.PID == 9900016 || std::abs(particle.PID) == 13) ? print(position, ":      ", particle, transverse_distance(particle)) : print(position, ": ", particle, transverse_distance(particle));
        }
        print("");
    }
    return 0.;
}

auto base_path()
{
    return "~/scratch/2.6.2_heavyion/";
}

auto x_sec_file_name(std::string const& run)
{
    return base_path() + run + "/cross_sections";
}

struct File {
    File(std::string const& name) : file(name) {}
    ~File()
    {
        file.close();
    }
    std::ifstream file;
};

class Line
{
    std::string string;
public:
    friend std::istream& operator>>(std::istream& stream, Line& line)
    {
        std::getline(stream, line.string);
        return stream;
    }
    operator std::string() const
    {
        return string;
    }
};

auto to_string(int number)
{
    return (number < 10 ? "0" : "") + std::to_string(number);
}

auto to_folder(int number)
{
    return "run_" + to_string(number) + "_decayed_1";
}

auto get_xsec(std::string const& run, int number)
{
    print(x_sec_file_name(run));
    File file1(x_sec_file_name(run));
    std::vector<std::string> lines;
    std::copy(std::istream_iterator<Line>(file1.file), std::istream_iterator<Line>(), std::back_inserter(lines));
    for (auto const& line : lines) {
        std::vector<std::string> strings;
        boost::split(strings, line, [](char c) {
            return c == ' ';
        });
        print(strings.at(0), to_folder(number));
        if (strings.at(0) == to_folder(number)) return strings.at(2);
    }
    return "Not found"s;
}

auto file_name(std::string const& run, int number)
{
    return base_path() + run + "/Events/" + to_folder(number) + "/tag_1_delphes_events.root";
}

struct Tree {
    Tree(std::string const& file_name) : chain("Delphes"), tree_reader(&chain)
    {
        chain.Add(file_name.c_str());
    }
    TChain chain;
    ExRootTreeReader tree_reader;
};

int main()
{
    auto run = "plain_scan"s;
//     auto run = "lead_scan"s;


    print(get_xsec(run,1));
    print(get_xsec(run,40));

    print("starting from", file_name(run, 1));
//     auto range = boost::irange(1, 49);
    auto range = boost::irange(1, 3);
    auto result = transform(range, [&run](auto number) {
        Tree tree(file_name(run, number));
        return AnalyseEvents(tree.tree_reader);
    });
    for (auto i : result) print(i);

}
