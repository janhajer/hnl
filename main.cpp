#include <iostream>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/count_if.hpp>
#include <boost/range/numeric.hpp>
#include <boost/filesystem.hpp>

#include "TClonesArray.h"

#include "ExRootAnalysis/ExRootTreeReader.h"

#include "classes/DelphesClasses.h"

#include "alphanum.hpp"

using namespace std::string_literals;

auto const neutrino_ID = 9900012;
auto const muon_ID = 13;

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

std::ostream& operator<<(std::ostream& stream, boost::filesystem::path const& path)
{
    return stream << path.string();
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

template<typename Container>
void print_line(Container const& container)
{
    for (auto const& element : container) std::cout << element << ", ";
    std::cout << std::endl;
}

std::string join_folder(std::string const& string)
{
    return string;
}

template<typename ... Arguments>
std::string join_folder(std::string const& string, Arguments ... arguments)
{
    return string + "/" + join_folder(arguments...);
}

auto base_path()
{
    return "/home/ucl/cp3/hajer/scratch/2.6.2_heavyion/results";
}

auto is_regular_file = static_cast<bool (*)(boost::filesystem::path const&)>(&boost::filesystem::is_regular_file);

auto has_ending(std::string const& string, std::string const& ending)
{
    return string.length() >= ending.length() ? 0 == string.compare(string.length() - ending.length(), ending.length(), ending) : false;
}

auto event_folder(std::string const& process)
{
    return join_folder(base_path(), process, "Events");
}

template<typename Function>
auto get_paths(boost::filesystem::path const& path, Function function)
{
    auto range = function(boost::make_iterator_range(boost::filesystem::directory_iterator(path), {}));
    //workaround as ranges of paths can not be sorted
    std::vector<boost::filesystem::path> paths;
    boost::range::copy(range, std::back_inserter(paths));
    return paths;
}

template<typename Function>
auto get_file(boost::filesystem::path const& path, Function function)
{
    auto paths = get_paths(path, function);
    if (paths.size() == 1) return paths.front();
    paths.size() > 1 ? print("Too many potential files:", paths) : print("No file");
    return boost::filesystem::path();
}

auto is_delphes = [](boost::filesystem::path const& path)
{
    return has_ending(path.filename().string(), "_delphes_events.root");
};

auto delphes_file(boost::filesystem::path const& path)
{
    return get_file(path, [](auto const & files) {
        return files | boost::adaptors::filtered(is_regular_file) | boost::adaptors::filtered(is_delphes);
    });
}

auto is_banner = [](boost::filesystem::path const& path)
{
    return has_ending(path.filename().string(), "_banner.txt");
};

auto banner_file(boost::filesystem::path const& path)
{
    return get_file(path, [](auto const & files) {
        return files | boost::adaptors::filtered(is_regular_file) | boost::adaptors::filtered(is_banner);
    });
}

auto is_directory = static_cast<bool (*)(boost::filesystem::path const&)>(&boost::filesystem::is_directory);

auto is_decayed_folder = [](boost::filesystem::path const& path)
{
    return has_ending(path.filename().string(), "_decayed_1");
};

auto decayed_folder(std::string const& path_name)
{
    boost::filesystem::path path(path_name);
    if (!boost::filesystem::is_directory(path)) print("Path:", path_name, "does not exist");
    auto paths = get_paths(path, [](auto const & range) {
        return range | boost::adaptors::filtered(is_directory) | boost::adaptors::filtered(is_decayed_folder);
    });
    auto sorted = boost::range::sort(paths, [](auto const & one, auto const & two) {
        return doj::alphanum_comp(one.string(), two.string()) < 0;
    });
    return sorted;
}

struct File {
//     File(std::string const& name) : file(name) {}
    File(boost::filesystem::path const& path) : file(path.string()) {}
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
auto read_file(boost::filesystem::path const& path, Predicate predicate, int pos, std::string const& name = "")
{

    File file(path);
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

auto get_xsec(boost::filesystem::path const& path)
{
    return read_file(banner_file(path), [](auto const & strings) {
        return strings.size() > 4 && strings.at(0) == "#" && strings.at(1) == "Integrated" && strings.at(2) == "weight" && strings.at(3) == "(pb)" && strings.at(4) == ":";
//         return strings.size() >= 2 && strings.at(0) == to_folder(point);
    }, 5, "cross section");
}

auto get_mass(boost::filesystem::path const& path)
{
    return read_file(banner_file(path), [](auto const & strings) {
        return strings.size() > 2 && strings.at(0) == std::to_string(neutrino_ID) && strings.at(2) == "#" && strings.at(3) == "mn1";
    }, 1, "mass");
}

auto get_coupling(boost::filesystem::path const& path)
{
    return read_file(banner_file(path), [](auto const & strings) {
        return strings.size() > 3 && strings.at(0) == std::to_string(4) && strings.at(2) == "#" && strings.at(3) == "vmun1";
    }, 1, "coupling");
}

auto get_width(boost::filesystem::path const& path)
{
    return read_file(banner_file(path), [](auto const & strings) {
        return strings.size() > 2 && strings.at(0) == "DECAY" && strings.at(1) == std::to_string(neutrino_ID);
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

auto origin(TClonesArray const& particles, int position, int check_id)
{
    std::vector<int> ids;
    while (position != -1) {
        auto& mother = get<GenParticle>(particles, position);
        if (std::abs(mother.PID) == check_id) return std::vector<int> {mother.PID};
        ids.emplace_back(mother.PID);
        position = mother.M1;
    };
    return ids;
}

template<typename Muon>
auto origin(TClonesArray const& muons, TClonesArray const& particles, int position, int check_id)
{
    auto& particle = get_particle<Muon>(muons, position);
    return std::abs(particle.PID) == check_id ? std::vector<int> {particle.PID} : origin(particles, particle.M1, check_id);
}

auto secondary_vertex(TClonesArray const& muons, int position)
{
    auto& muon = get<Muon>(muons, position);
    auto& particle = get_particle(muon);
    if (std::abs(particle.PID) != muon_ID) print("Misidentified muon");
    auto dist = transverse_distance(particle);
//     if (dist > 1.) print(dist);
    return dist;
}

auto number_of_displaced(TClonesArray const& muons, TClonesArray const& particles)
{
    return boost::count_if(range(muons.GetEntriesFast()), [&muons, &particles](auto position) {
        auto hit = secondary_vertex(muons, position) > 1.;
        if (!hit) return hit;
        auto ids = origin<Muon>(muons, particles, position, neutrino_ID);
        if (std::abs(ids.front()) != neutrino_ID) print_line(ids);
        return hit;
    });
}

auto analyse_events(boost::filesystem::path const& path)
{
    TChain chain("Delphes");
    chain.Add(delphes_file(path).c_str());
    ExRootTreeReader reader(&chain);
    auto& muons = *reader.UseBranch("Muon");
    auto& particles = *reader.UseBranch("Particle");
    auto entries = reader.GetEntries();
    return entries == 0 ? 0. : boost::count_if(range(reader.GetEntries()), [&reader, &muons, &particles](auto entry) {
        reader.ReadEntry(entry);
        return number_of_displaced(muons, particles) > 0;
    }) / static_cast<double>(entries);
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
    if (argc < 2) {
        print("Please pass the process name as argument");
        return 0;
    }
    std::vector<std::string> arguments(argv, argv + argc);
    auto process = arguments.at(1);
    std::vector<std::string> results;
    for (auto const& folder : decayed_folder(event_folder(process))) {
        auto result = get_mass(folder) + " " + get_coupling(folder) + " " + std::to_string(analyse_events(folder)) + " " +  get_xsec(folder) + " " + get_width(folder);
        print(result);
        results.emplace_back(result);
    }
    save_result(results, process);
}
