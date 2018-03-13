#include <iostream>
#include <fstream>

#include <boost/algorithm/string.hpp>

#include <boost/range/irange.hpp>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/transformed.hpp>
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
    return stream << boost::accumulate(boost::adaptors::index(container), [](std::ostream & stream, auto const & element) {
        return stream << '\n' << element.index() << ": " << element.value();
    });
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

auto to_string(int point)
{
    return (point < 10 ? "0" : "") + std::to_string(point);
}

auto join_name(std::string const& first, std::string const& second)
{
    return first + "_" + second;
}

auto run_name(int point)
{
    return join_name("run", to_string(point));
}

auto to_folder(int point)
{
    return join_name(run_name(point), "decayed_1");
}

auto file_name(std::string const& process, int point)
{
    return base_path() + process + "/Events/" + to_folder(point) + "/tag_1_delphes_events.root";
}

auto banner_name(std::string const& process, int point)
{
    auto name = run_name(point);
    return base_path() + process + "/Events/" + name + "/" + join_name(name, "tag_1_banner.txt");
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

template<typename Predicate>
auto read_file(std::string const& file_name, Predicate predicate, int pos, std::string const& name = "")
{
    File file(file_name);
    std::vector<std::string> lines;
    std::copy(std::istream_iterator<Line>(file.file), std::istream_iterator<Line>(), std::back_inserter(lines));
    for (auto& line : lines) {
        boost::trim_if(line, boost::is_any_of("\t "));
        std::vector<std::string> strings;
        boost::split(strings, line, [](char c) {
            return c == ' ';
        }, boost::token_compress_on);
        if (predicate(strings)) return strings.at(pos);
    }
    return name + " value not found"s;
}

auto get_xsec(std::string const& process, int point)
{
    return read_file(x_sec_file_name(process), [point](std::vector<std::string> const & strings) {
        return strings.size() >= 2 && strings.at(0) == to_folder(point);
    }, 2, "cross section");
}

auto get_mass(std::string const& process, int point)
{
    return read_file(banner_name(process, point), [point](std::vector<std::string> const & strings) {
        return strings.size() > 2 && strings.at(0) == std::to_string(9900012) && strings.at(2) == "#" && strings.at(3) == "mn1";
    }, 1, "mass");
}

auto get_coupling(std::string const& process, int point)
{
    return read_file(banner_name(process, point), [point](std::vector<std::string> const & strings) {
        return strings.size() > 3 && strings.at(0) == std::to_string(4) && strings.at(2) == "#" && strings.at(3) == "vmun1";
    }, 1, "coupling");
}

auto get_width(std::string const& process, int point)
{
    return read_file(banner_name(process, point), [point](std::vector<std::string> const & strings) {
        return strings.size() > 2 && strings.at(0) == "DECAY" && strings.at(1) == std::to_string(9900012);
    }, 2, "width");
}

template<typename Object>
auto& get_particle(Object const& object)
{
    return static_cast<GenParticle&>(*object.Particle.GetObject());
}


template<typename Object>
auto& get(TClonesArray const& array, int position)
{
    return static_cast<Object&>(*array.At(position));
}

template<typename Object>
auto& get_particle(TClonesArray const& muons, int position)
{
    return get_particle(get<Object>(muons, position));
}

template<typename Muon>
auto& get_mother(TClonesArray const& muons, TClonesArray const& particles, int position)
{
    return get<GenParticle>(particles, get_particle<Muon>(muons, position).M1);
}

// template<typename Muon>
// auto& get_grand_mother(TClonesArray const& muons, TClonesArray const& particles, int position)
// {
//     return get<GenParticle>(particles, get_mother<Muon>(muons, particles, position).M1);
// }
//
// template<typename Muon>
// auto& get_grand_grand_mother(TClonesArray const& muons, TClonesArray const& particles, int position)
// {
//     return get<GenParticle>(particles, get_grand_mother<Muon>(muons, particles, position).M1);
// }

auto const neutrino_ID = 9900012;
auto const muon_ID = 13;

auto check_origin(TClonesArray const& particles, int position, int check_id)
{
    auto id = check_id;
    while (std::abs(id) == check_id) {
        auto& mother = get<GenParticle>(particles, position);
        id = mother.PID;
        position = mother.M1;
    };
    return id;
}

template<typename Muon>
auto check_origin(TClonesArray const& muons, TClonesArray const& particles, int position, int check_id)
{
    auto& particle = get_particle<Muon>(muons, position);
    return std::abs(particle.PID) == check_id ? check_origin(particles, particle.M1, check_id) : particle.PID;
}

auto secondary_vertex(TClonesArray const& muons, int position)
{
    auto& muon =  get<Muon>(muons, position);
    auto& particle = get_particle(muon);
    if (std::abs(particle.PID) != 13) print("background");
    return transverse_distance(particle);
}

auto number_of_displaced(TClonesArray const& muons, TClonesArray const& particles)
{
    return boost::count_if(range(muons.GetEntriesFast()), [&muons, &particles](auto position) {
        auto cut = secondary_vertex(muons, position) > 100.;
        if (!cut) return cut;
//         auto& mother = get_mother<Muon>(muons, particles, position);
//         auto& grand_mother = get_grand_mother<Muon>(muons, particles, position);
//         auto& grand_grand_mother = get_grand_grand_mother<Muon>(muons, particles, position);
//         if (mother.PID != neutrino_ID && grand_mother.PID != neutrino_ID && grand_grand_mother.PID != neutrino_ID) print("Muon from", mother.PID, "and", grand_mother.PID, "and", grand_grand_mother.PID);
        print(check_origin<Muon>(muons, particles, position, muon_ID));
        return cut;
    });
}

auto analyse_events(std::string const& process, int point)
{
    TChain chain("Delphes");
    chain.Add(file_name(process, point).c_str());
    ExRootTreeReader reader(&chain);
    auto& muons = *reader.UseBranch("Muon");
    auto& particles = *reader.UseBranch("Particle");
    return boost::count_if(range(reader.GetEntries()), [&reader, &muons, &particles](auto entry) {
        reader.ReadEntry(entry);
        auto number = number_of_displaced(muons, particles);
//         if (number > 1) print(number, "displaced muons");
        return number > 0;
    }) / static_cast<double>(reader.GetEntries());
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
    if (argc < 2) return 0;
    std::vector<std::string> arguments(argv, argv + argc);
    auto process = arguments.at(1);
    print("starting from", file_name(process, 1));
    auto points = boost::irange(1, 49);
    auto result = transform(points, [&process](auto point) {
        return get_mass(process, point) + " " + get_coupling(process, point) + " " + std::to_string(analyse_events(process, point)) + " " +  get_xsec(process, point) + " " + get_width(process, point);
    });
    save_result(result, process);
}
