#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/size.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/count_if.hpp>
#include <boost/range/algorithm/find_if.hpp>

#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

#include "classes/DelphesClasses.h"

#include "alphanum.hpp"

using namespace std::string_literals;

auto const neutrino_ID = 9900012;
auto const electron_ID = 11;
auto const muon_ID = 13;
auto const tau_ID = 15;

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

std::ostream& operator<<(std::ostream& stream, boost::filesystem::path const& path)
{
    return stream << path.string();
}

std::ostream& operator<<(std::ostream& stream, GenParticle const& particle)
{
    return stream << "(" << particle.PID << ", " << particle.Status << ")";
}

std::ostream& operator<<(std::ostream& stream, Muon const& muon)
{
    return stream << "(" << muon.PT << ", " << muon.Eta << ")";
}

std::ostream& operator<<(std::ostream& stream, Electron const& muon)
{
    return stream << "(" << muon.PT << ", " << muon.Eta << ")";
}

std::ostream& operator<<(std::ostream& stream, Jet const& muon)
{
    return stream << "(" << muon.PT << ", " << muon.Eta << ")";
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

template<typename Lepton>
auto get_id()
{
    print("code gone wrong");
    return 0;
}

template <>
auto get_id<Electron>()
{
    return electron_ID;
}

template <>
auto get_id<Muon>()
{
    return muon_ID;
}

template <>
auto get_id<Jet>()
{
    return tau_ID;
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
    return ""s + getenv("HOME") + "/scratch/results";
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

auto has_ending(std::string const& string, std::string const& ending)
{
    return string.length() >= ending.length() ? 0 == string.compare(string.length() - ending.length(), ending.length(), ending) : false;
}

auto is_directory = static_cast<bool (*)(boost::filesystem::path const&)>(&boost::filesystem::is_directory);

auto is_decayed_folder = [](boost::filesystem::path const& path)
{
    return has_ending(path.filename().string(), "_decayed_1");
};

auto decayed_folders(std::string const& path_name)
{
    boost::filesystem::path path(path_name);
    if (!boost::filesystem::is_directory(path)) print("Path:", path_name, "does not exist");
    auto paths = get_paths(path, [](auto const & range) {
        return range | boost::adaptors::filtered(is_directory) | boost::adaptors::filtered(is_decayed_folder);
    });
    return boost::range::sort(paths, [](auto const & one, auto const & two) {
        return doj::alphanum_comp(one.string(), two.string()) < 0;
    });
}

template<typename Function>
auto get_file(boost::filesystem::path const& path, Function function)
{
    auto paths = get_paths(path, function);
    if (paths.size() == 1) return paths.front();
    paths.size() > 1 ? print("Too many potential files:", paths) : print("No file");
    return boost::filesystem::path();
}

auto is_regular_file = static_cast<bool (*)(boost::filesystem::path const&)>(&boost::filesystem::is_regular_file);

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
struct File {
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

auto split_line(std::string const& line)
{
    std::vector<std::string> strings;
    boost::split(strings, line, [](char c) {
        return c == ' ';
    }, boost::token_compress_on);
    return strings;
}

template<typename Predicate>
auto read_file(boost::filesystem::path const& path, Predicate predicate, int pos)
{
    File file(path);
    std::vector<std::string> lines;
    std::copy(std::istream_iterator<Line>(file.file), std::istream_iterator<Line>(), std::back_inserter(lines));
    auto found = boost::range::find_if(lines, [&predicate](auto & line) {
        boost::trim_if(line, boost::is_any_of("\t "));
        return predicate(split_line(line));
    });
    return found == lines.end() ? "value not found"s : split_line(*found).at(pos);
}

auto get_xsec(boost::filesystem::path const& path)
{
    return read_file(banner_file(path), [](auto const & strings) {
        return strings.size() > 4 && strings.at(0) == "#" && strings.at(1) == "Integrated" && strings.at(2) == "weight" && strings.at(3) == "(pb)" && strings.at(4) == ":";
    }, 5);
}

auto get_mass(boost::filesystem::path const& path)
{
    return read_file(banner_file(path), [](auto const & strings) {
        return strings.size() > 2 && strings.at(0) == std::to_string(neutrino_ID) && strings.at(2) == "#" && strings.at(3) == "mn1";
    }, 1);
}

auto get_coupling(boost::filesystem::path const& path, int pos, std::string const& name)
{
    return read_file(banner_file(path), [&name, pos](auto const & strings) {
        return strings.size() > 3 && strings.at(0) == std::to_string(pos) && strings.at(2) == "#" && strings.at(3) == name;
    }, 1);
}

auto get_e_coupling(boost::filesystem::path const& path)
{
    return get_coupling(path, 1, "ven1");
}

auto get_mu_coupling(boost::filesystem::path const& path)
{
    return get_coupling(path, 4, "vmun1");
}

auto get_tau_coupling(boost::filesystem::path const& path)
{
    return get_coupling(path, 7, "vtan1");
}

auto get_width(boost::filesystem::path const& path)
{
    return read_file(banner_file(path), [](auto const & strings) {
        return strings.size() > 2 && strings.at(0) == "DECAY" && strings.at(1) == std::to_string(neutrino_ID);
    }, 2);
}

template<typename Lepton>
auto& get_particle(Lepton const& lepton)
{
    return static_cast<GenParticle&>(*lepton.Particle.GetObject());
}

auto get_particles(Jet const& jet)
{
    std::vector<GenParticle> particles;
    auto* iterator = static_cast<TRefArrayIter*>(jet.Particles.MakeIterator());
    if (!iterator) return particles;
    while (auto* object = iterator->Next()) particles.emplace_back(*static_cast<GenParticle*>(object));
    return particles;
}

auto origin(TTreeReaderArray<GenParticle> const& particles, int position, int check_id)
{
    std::vector<int> ids;
    while (position != -1) {
        auto& mother = particles.At(position);
        if (std::abs(mother.PID) == check_id) return std::vector<int> {mother.PID};
        ids.emplace_back(mother.PID);
        position = mother.M1;
    };
    return ids;
}

template<typename Lepton>
auto origin(Lepton const& lepton, TTreeReaderArray<GenParticle> const& particles, int check_id)
{
    auto& particle = get_particle(lepton);
    return std::abs(particle.PID) == check_id ? std::vector<int> {particle.PID} : origin(particles, particle.M1, check_id);
}

template<typename One, typename Two>
auto insert(One& one, Two const& two)
{
    one.insert(one.end(), two.begin(), two.end());
}

template<>
auto origin(Jet const& lepton, TTreeReaderArray<GenParticle> const& gen_particles, int check_id)
{
    std::vector<int> result;
    for (auto const& particle : get_particles(lepton)) std::abs(particle.PID) == check_id ? result.emplace_back(particle.PID) : insert(result, origin(gen_particles, particle.M1, check_id));
    return result;
}

template<typename Lepton>
auto secondary_vertex(Lepton const& lepton)
{
    auto& particle = get_particle(lepton);
    if (std::abs(particle.PID) != get_id<Lepton>()) print("Misidentified lepton");
    return transverse_distance(particle);
}

template<>
auto secondary_vertex(Jet const& lepton)
{
    using namespace boost::accumulators;
    accumulator_set<float, stats<tag::mean>> distances;
//     auto hit = false;
    auto particles = get_particles(lepton);
    for (auto const& particle : particles) {
//         if (std::abs(particle.PID) == get_id<Jet>()) hit = true;
        distances(transverse_distance(particle));
    }
//     if (!hit) print("Misidentified tau", boost::adaptors::transform(particles, [](auto const& particle){
//         return std::to_string(particle.PID) + " ";
//     }));
    auto d = mean(distances);
//     if (d > 1) print("displaced tau", d);
    return d;
}

template<typename Lepton>
auto disp(){
    return 10.;
}

template<>
auto disp<Jet>(){
    return 30.;
}

template<typename Lepton>
auto hard(){
    return 25.;
}

template<>
auto hard<Jet>(){
    return 50.;
}

template<typename Lepton>
auto is_hard(Lepton const& lepton)
{
    return secondary_vertex(lepton) < disp<Lepton>() && lepton.PT > hard<Lepton>();
}

template<typename Leptons>
auto number_of_displaced(Leptons const& leptons, TTreeReaderArray<GenParticle> const& particles)
{
    return boost::count_if(leptons, [&particles](auto lepton) {
        auto distance = secondary_vertex(lepton);
        auto hit = distance > disp<typename Leptons::iterator::value_type>();
//         if (!hit) return hit;
//         auto ids = origin(lepton, particles, neutrino_ID);
//         if (std::abs(ids.front()) != neutrino_ID) {
//             print_line(ids);
//             print(distance);
//         };
        return hit;
    });
}

template<typename Leptons>
auto number_of_hard(Leptons const& leptons)
{
    return boost::count_if(leptons, [](auto lepton) {
        return is_hard(lepton);
    });
}

template<typename Predicate>
auto count_if(TTreeReader& reader, Predicate predicate)
{
    auto counter = 0;
    while (reader.Next()) if (predicate()) ++counter;
    return counter;
}

auto get_taus(TTreeReaderArray<Jet> const& jets)
{
    return boost::adaptors::filter(jets, [](auto const & jet) {
        return jet.TauTag;
    });
}

auto get_signal(TTreeReader& reader)
{
    TTreeReaderArray<Electron> electrons(reader, "Electron");
    TTreeReaderArray<Muon> muons(reader, "Muon");
    TTreeReaderArray<Jet> jets(reader, "Jet");
    TTreeReaderArray<GenParticle> particles(reader, "Particle");
    return count_if(reader, [&]() {
        particles.IsEmpty();
        auto taus = get_taus(jets);
//         print("number of taus", boost::size(taus));
//         auto displaced = number_of_displaced(muons, particles);
//         auto hard = number_of_hard(muons);
        auto displaced = number_of_displaced(electrons, particles) + number_of_displaced(muons, particles) + number_of_displaced(taus, particles);
        auto hard = number_of_hard(electrons) + number_of_hard(muons) + number_of_hard(taus);
//         print(displaced, hard,  displaced > 0 && hard > 0);
        return displaced > 0 && hard > 0;
    });
}

auto get_efficiency(boost::filesystem::path const& path)
{
    TFile file(delphes_file(path).c_str(), "read");
    TTreeReader reader("Delphes", &file);
    auto entries = reader.GetEntries(false);
    if (entries == 0) {
        print("No events");
        return 0.;
    }
    return get_signal(reader) / static_cast<double>(entries);
}

template<typename Result>
void save_result(Result const& result, std::string const& process)
{
    print_line(result);
    std::ofstream file("./" + process + ".dat");
    std::ostream_iterator<std::string> iterator(file, "\n");
    boost::copy(result, iterator);
}

auto get_result(boost::filesystem::path const& folder)
{
    auto result = get_mass(folder);
    result += " " + get_e_coupling(folder);
    result += " " + get_mu_coupling(folder);
    result += " " + get_tau_coupling(folder);
    result += " " + std::to_string(get_efficiency(folder));
    result += " " + get_xsec(folder);
    result += " " + get_width(folder);
    print(result);
    return result;
}

auto get_header()
{
    std::string header = "# mass";
    header += " e_coupling";
    header += " mu_coupling";
    header += " tau_coupling";
    header += " efficiency";
    header += " crosssection";
    header += " width";
    print(header);
    return header;
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        print("Please pass the process name as argument");
        return 0;
    }
    std::vector<std::string> arguments(argv, argv + argc);
    auto process = arguments.at(1);
    std::vector<std::string> results {get_header()};
    for (auto const& folder : decayed_folders(event_folder(process))) results.emplace_back(get_result(folder));
    save_result(results, process);
}
