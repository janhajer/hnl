#include <iostream>
#include <fstream>

#include <boost/optional.hpp>

#include <boost/algorithm/string.hpp>

#include <boost/range/irange.hpp>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/transform.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/count_if.hpp>
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

template<typename Object>
struct Branch {
    Branch(TClonesArray* input) : array(input), range(boost::irange(0, array->GetEntriesFast())) {}
    Branch(TClonesArray& input) : array(&input), range(boost::irange(0, array->GetEntriesFast())) {}
    auto begin()
    {
        return range.begin();
    }
    auto end()
    {
        return range.end();
    }
    void update()
    {
        range = boost::irange(0, array->GetEntriesFast());
    }
    auto& at(int position) const
    {
        return static_cast<Object&>(*array->At(position));
    }
    TClonesArray* array;
    boost::integer_range<int> range;
};

struct Tree {
    Tree(std::string const& file_name) : chain("Delphes"), reader(&chain), range(0, 0)
    {
        chain.Add(file_name.c_str());
        range = boost::irange(0ll, reader.GetEntries());
    }
    template<typename Object>
    Branch<Object> use_branch(std::string const& name)
    {
        return reader.UseBranch(name.c_str());
    }
    auto begin()
    {
        return range.begin();
    }
    auto end()
    {
        return range.end();
    }
    TChain chain;
    ExRootTreeReader reader;
    boost::integer_range<long long> range;
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
auto read_file(std::string const& file_name, Predicate predicate, int pos, std::string const& name = "")
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

auto& get_particle(Muon& muon)
{
    return static_cast<GenParticle&>(*muon.Particle.GetObject());
}

template<typename Particle>
void read_entry(Tree& tree, Branch<Particle>& branch, int entry)
{
    tree.reader.ReadEntry(entry);
    branch.update();
}

auto muon_distance(Branch<Muon> const& muons, int position)
{
    auto& muon = muons.at(position);
    auto& particle = get_particle(muon);
    auto distance = transverse_distance(particle);
    if (distance < 10) return 0.f;
    if (std::abs(particle.PID) != 13) print("background");
    return distance;
}

auto number_of_dispalced(Branch<Muon> const& branch)
{
    return boost::count_if(range(branch.array->GetEntriesFast()), [&](auto position) {
        return muon_distance(branch, position) > 10.;
    });
}

auto analyse_events(std::string const& process, int number)
{
    Tree tree(file_name(process, number));
    auto muons = tree.use_branch<Muon>("Muon");
    tree.use_branch<GenParticle>("Particle");
    return boost::count_if(range(tree.reader.GetEntries()), [&](auto entry) {
        read_entry(tree, muons, entry);
        return number_of_dispalced(muons) > 0;
    }) / tree.reader.GetEntries();
}

auto AnalyseEvents(std::string const& process, int number)
{
    Tree tree(file_name(process, number));
    auto muons = tree.use_branch<Muon>("Muon");
    tree.use_branch<GenParticle>("Particle");
    auto displaced_number = 0;
    for (auto entry : tree) {
        read_entry(tree, muons, entry);
        std::vector<double> result;
        for (auto position : muons) {
            auto& muon = muons.at(position);
            auto& particle = get_particle(muon);
            auto distance = transverse_distance(particle);
            if (distance < 10) continue;
            if (std::abs(particle.PID) == 13) result.emplace_back(distance);
            else print("background");
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
        return get_mass(process, number) + " " + get_coupling(process, number) + " " + std::to_string(analyse_events(process, number)) + " " +  get_xsec(process, number) + " " + get_width(process, number);
    });
    save_result(result, process);
}
