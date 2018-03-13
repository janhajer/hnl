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

template<typename Object>
auto& get_particle(Object const& object)
{
    return static_cast<GenParticle&>(*object.Particle.GetObject());
}

auto secondary_vertex(TClonesArray const& muons, int position)
{
    auto& muon =  static_cast<Muon&>(*muons.At(position));
    auto& particle = get_particle(muon);
    if (std::abs(particle.PID) != 13) print("background");
    return transverse_distance(particle);
}

auto number_of_dispalced(TClonesArray const& branch)
{
    return boost::count_if(range(branch.GetEntriesFast()), [&](auto position) {
        return secondary_vertex(branch, position) > 100.;
    });
}

auto analyse_events(std::string const& process, int number)
{
    TChain chain("Delphes");
    chain.Add(file_name(process, number).c_str());
    ExRootTreeReader reader(&chain);
    auto& muons = *reader.UseBranch("Muon");
    reader.UseBranch("Particle");
    return boost::count_if(range(reader.GetEntries()), [&](auto entry) {
        reader.ReadEntry(entry);
        return number_of_dispalced(muons) > 0;
    }) / double(reader.GetEntries());
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
    if(argc < 1) return 0;
    auto process = arguments.at(1);
    print("string",process);
    print("starting from", file_name(process, 1));
    auto range = boost::irange(1, 49);
    auto result = transform(range, [&process](auto number) {
        return get_mass(process, number) + " " + get_coupling(process, number) + " " + std::to_string(analyse_events(process, number)) + " " +  get_xsec(process, number) + " " + get_width(process, number);
    });
    save_result(result, process);
}
