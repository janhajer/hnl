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

#include "ExRootAnalysis/ExRootTreeReader.h"

#include "classes/DelphesClasses.h"

using namespace std::string_literals;

template<typename Object>
auto sqr(Object const& object)
{
    return object * object;
}

template<typename Particle>
auto transverse_distance(Particle const& particle)
{
    return std::sqrt(sqr(particle.X) + sqr(particle.Y));
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


template<typename Key_, typename Value_>
auto& operator<<(std::ostream& stream, std::pair<Key_, Value_> const& pair)
{
    return stream << '(' << pair.first << ", " << pair.second << ')';
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

auto base_path()
{
    return "/home/ucl/cp3/hajer/scratch/2.6.2_heavyion/";
}

auto x_sec_file_name(std::string const& process)
{
    return base_path() + process + "/cross_sections";
}

auto to_string(int number)
{
    return (number < 10 ? "0" : "") + std::to_string(number);
}

auto join_name(std::string const& first, std::string const& second)
{
    return first + "_" + second;
}

auto run_name(int number)
{
    return join_name("run", to_string(number));
}

auto to_folder(int number)
{
    return join_name(run_name(number), "decayed_1");
}

auto file_name(std::string const& process, int number)
{
    return base_path() + process + "/Events/" + to_folder(number) + "/tag_1_delphes_events.root";
}

auto banner_name(std::string const& process, int number)
{
    auto name = run_name(number);
    return base_path() + process + "/Events/" + name + "/" + join_name(name, "tag_1_banner.txt");
}

struct Tree {
    Tree(std::string const& file_name) : chain("Delphes"), reader(&chain)
    {
        chain.Add(file_name.c_str());
    }
    TChain chain;
    ExRootTreeReader reader;
};

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

template<typename Predicate>
auto read_file(std::string const& file_name, Predicate predicate, int pos, std::string const& name)
{
    File file(file_name);
    std::vector<std::string> lines;
    std::copy(std::istream_iterator<Line>(file.file), std::istream_iterator<Line>(), std::back_inserter(lines));
    for (auto& line : lines) {
        std::vector<std::string> strings;
        boost::trim_if(line, boost::is_any_of("\t "));
        boost::split(strings, line, [](char c) {
            return c == ' ';
        }, boost::token_compress_on);
        if (predicate(strings)) return strings.at(pos);
    }
    return name + " value not found"s;
}

auto get_xsec(std::string const& process, int number)
{
    return read_file(x_sec_file_name(process), [number](std::vector<std::string> const & strings) {
        return strings.size() >= 2 && strings.at(0) == to_folder(number);
    }, 2, "cross section");
}

auto get_mass(std::string const& process, int number)
{
    return read_file(banner_name(process, number), [number](std::vector<std::string> const & strings) {
        return strings.size() > 2 && strings.at(0) == std::to_string(9900012) && strings.at(2) == "#" && strings.at(3) == "mn1";
    }, 1, "mass");
}

auto get_coupling(std::string const& process, int number)
{
    return read_file(banner_name(process, number), [number](std::vector<std::string> const & strings) {
        return strings.size() > 3 && strings.at(0) == std::to_string(4) && strings.at(2) == "#" && strings.at(3) == "vmun1";
    }, 1, "coupling");
}

auto get_width(std::string const& process, int number)
{
    return read_file(banner_name(process, number), [number](std::vector<std::string> const & strings) {
        return strings.size() > 2 && strings.at(0) == "DECAY" && strings.at(1) == std::to_string(9900012);
    }, 2, "width");
}

// auto get_xsec(std::string const& run, int number)
// {
//     File file(x_sec_file_name(run));
//     std::vector<std::string> lines;
//     std::copy(std::istream_iterator<Line>(file.file), std::istream_iterator<Line>(), std::back_inserter(lines));
//     for (auto const& line : lines) {
//         std::vector<std::string> strings;
//         boost::split(strings, line, [](char c) {
//             return c == ' ';
//         });
//         if (strings.size() >= 2 && strings.at(0) == to_folder(number)) return strings.at(2);
//     }
//     return "XSec value not found"s;
// }


// auto get_mass(std::string const& run, int number)
// {
//     File file(banner_name(run, number));
//     std::vector<std::string> lines;
//     std::copy(std::istream_iterator<Line>(file.file), std::istream_iterator<Line>(), std::back_inserter(lines));
//     for (auto& line : lines) {
//         std::vector<std::string> strings;
//         boost::trim_if(line, boost::is_any_of("\t "));
//         boost::split(strings, line, [](char c) {
//             return c == ' ';
//         }, boost::token_compress_on);
//         if (strings.size() > 2 && strings.at(0) == std::to_string(9900012) && strings.at(2) == "#" && strings.at(3) == "mn1") return strings.at(1);
//     }
//     return "Mass value not found"s;
// }

// auto get_coupling(std::string const& run, int number)
// {
//     File file(banner_name(run, number));
//     std::vector<std::string> lines;
//     std::copy(std::istream_iterator<Line>(file.file), std::istream_iterator<Line>(), std::back_inserter(lines));
//     for (auto& line : lines) {
//         std::vector<std::string> strings;
//         boost::trim_if(line, boost::is_any_of("\t "));
//         boost::split(strings, line, [](char c) {
//             return c == ' ';
//         }, boost::token_compress_on);
//         if (strings.size() > 3 && strings.at(0) == std::to_string(4) && strings.at(2) == "#" && strings.at(3) == "vmun1") return strings.at(1);
//     }
//     return "Coupling value not found"s;
// }

// auto get_width(std::string const& run, int number)
// {
//     File file(banner_name(run, number));
//     std::vector<std::string> lines;
//     std::copy(std::istream_iterator<Line>(file.file), std::istream_iterator<Line>(), std::back_inserter(lines));
//     for (auto& line : lines) {
//         std::vector<std::string> strings;
//         boost::trim_if(line, boost::is_any_of("\t "));
//         boost::split(strings, line, [](char c) {
//             return c == ' ';
//         }, boost::token_compress_on);
//         if (strings.size() > 2 && strings.at(0) == "DECAY" && strings.at(1) == std::to_string(9900012)) return strings.at(2);
//     }
//     return "Width value not found"s;
// }

// auto particle_distance(TClonesArray const& particle_branch, int number)
// {
//     auto& particle = static_cast<GenParticle&>(*particle_branch.At(number));
//     return std::abs(particle.PID) == 13 ? transverse_distance(particle) : 0;
// }
//
// auto muon_distance(TClonesArray const& muon_branch, int number)
// {
//     auto& muon = static_cast<Muon&>(*muon_branch.At(number));
//     auto& particle = static_cast<GenParticle&>(*muon.Particle.GetObject());
//     return transverse_distance(particle);
// }
//
// template<typename Function>
// auto displacement(ExRootTreeReader& tree.reader, int entry, TClonesArray const& branch, Function const& function)
// {
//     tree.reader.ReadEntry(entry);
//     return boost::accumulate(range(branch.GetEntriesFast()), 0, [&](auto & sum, auto number) {
//         return sum += function(number) > 10 ? 1 : 0 ;
//     }) > 0;
// }
//
// auto analyse_events(ExRootTreeReader& tree.reader)
// {
//     auto& muon_branch = *tree.reader.UseBranch("Muon");
//     return boost::accumulate(range(tree.reader.GetEntries()), 0., [&](auto & sum, auto entry) {
//         return sum += displacement(tree.reader, entry, muon_branch, [&](auto number) {
//             return particle_distance(muon_branch, number);
//         }) ? 1 : 0;
//     }) / tree.reader.GetEntries();
// }

auto AnalyseEvents(std::string const& process, int number)
{
    Tree tree(file_name(process, number));
    auto& muon_branch = *tree.reader.UseBranch("Muon");
    auto& particle_branch = *tree.reader.UseBranch("Particle");
    // Loop over all events
    auto displaced_number = 0;
    for (auto entry : range(tree.reader.GetEntries())) {
        // Load selected branches with data from specified event
        tree.reader.ReadEntry(entry);
        std::vector<double> result;
        auto entries = muon_branch.GetEntriesFast();
        for (auto position : range(entries)) {
            auto& muon = static_cast<Muon&>(*muon_branch.At(position));
            auto& particle = static_cast<GenParticle&>(*muon.Particle.GetObject());
            auto distance = transverse_distance(particle);
            if (distance > 100) result.emplace_back(distance);
        }
        if (!result.empty()) ++displaced_number;
    }
    return static_cast<double>(displaced_number) / tree.reader.GetEntries();
}

auto AnalyseEvents2(std::string const& process, int number)
{
    Tree tree(file_name(process, number));
    auto& particle_branch = *tree.reader.UseBranch("Particle");
    auto number_displaced = 0;
    // loop over all events
    for (auto entry : range(tree.reader.GetEntries())) {
        tree.reader.ReadEntry(entry);
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
    return static_cast<double>(number_displaced) / tree.reader.GetEntries();
}

// auto AnalyseEvents3(ExRootTreeReader& tree.reader)
// {
//     auto& particle_branch = *tree.reader.UseBranch("Particle");
//     auto number_of_muons = 0;
//     for (auto entry : range(tree.reader.GetEntries())) {
//         tree.reader.ReadEntry(entry);
//         std::vector<double> result;
//         for (auto position : range(particle_branch.GetEntriesFast())) {
//             auto& particle = static_cast<GenParticle&>(*particle_branch.At(position));
//             auto distance = transverse_distance(particle);
//             if (distance > 100 && distance < 2000)  result.emplace_back(distance);
//         }
//         if (!result.empty()) ++number_of_muons;
//     }
//     print("displaced", number_of_muons);
//     return static_cast<double>(number_of_muons) / tree.reader.GetEntries();
// }
//
// auto AnalyseEvents4(ExRootTreeReader& tree.reader)
// {
//     auto& particle_branch = *tree.reader.UseBranch("Particle");
//     for (auto entry : range(tree.reader.GetEntries())) {
//         tree.reader.ReadEntry(entry);
//         std::vector<double> result;
//         for (auto position : range(particle_branch.GetEntriesFast())) {
//             auto& particle = static_cast<GenParticle&>(*particle_branch.At(position));
//             (particle.PID == 9900012 || particle.PID == 9900014 || particle.PID == 9900016 || std::abs(particle.PID) == 13) ? print(position, ":      ", particle, transverse_distance(particle)) : print(position, ": ", particle, transverse_distance(particle));
//         }
//         print("");
//     }
//     return 0.;
// }

template<typename Result>
void save_result(Result const& result, std::string const& process)
{
    for (auto i : result) print(i);
    std::ofstream file("./" + process + ".dat");
    std::ostream_iterator<std::string> iterator(file, "\n");
    boost::copy(result, iterator);
}

int main(int argc, char** argv)
{
    std::vector<std::string> arguments(argv, argv + argc);
//     auto process = "proton_scan"s;
    auto process = "lead_scan"s;

    print("starting from", file_name(process, 1));
    auto range = boost::irange(1, 49);
//     auto range = boost::irange(1, 3);
    auto result = transform(range, [&process](auto number) {
        return get_mass(process, number) + " " + get_coupling(process, number) + " " + std::to_string(AnalyseEvents(process, number)) + " " +  get_xsec(process, number) + " " + get_width(process, number);
    });
    save_result(result, process);
}
