#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
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

template<typename Element>
auto operator+(std::vector<Element> const& one, std::vector<Element> const& two)
{
    auto copy = one;
    copy.insert(copy.end(), two.begin(), two.end());
    return copy;
}

template<typename Element>
auto& operator+=(std::vector<Element>& one, std::vector<Element> const& two)
{
    one.insert(one.end(), two.begin(), two.end());
    return one;
}

template<typename Element, typename Predicate>
auto find_erase(std::vector<Element>& container, Predicate predicate) -> boost::optional<Element> {
    auto found = boost::range::find_if(container, predicate);
    if (found == container.end()) return boost::none;
    auto element = *found;
    container.erase(found);
    return element;
}

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
    return stream << "(" << particle.PID << ", " << transverse_distance(particle) << ")";
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

struct Line {
    friend std::istream& operator>>(std::istream& stream, Line& line)
    {
        std::getline(stream, line.string);
        return stream;
    }
    operator std::string() const
    {
        return string;
    }
private:
    std::string string;
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
auto get_particles(Lepton const& lepton) -> std::vector<GenParticle> {
    return {static_cast<GenParticle&>(*lepton.Particle.GetObject())};
}

template<>
auto get_particles(Jet const& jet) -> std::vector<GenParticle> {
    std::vector<GenParticle> particles;
    auto* iterator = static_cast<TRefArrayIter*>(jet.Particles.MakeIterator());
    if (!iterator) return particles;
    while (auto* object = iterator->Next()) particles.emplace_back(*static_cast<GenParticle*>(object));
    return particles;
}

auto origin(TTreeReaderArray<GenParticle> const& particles, int position, int check_id) -> boost::optional<GenParticle> {
    while (position != -1)
    {
        auto& mother = particles.At(position);
        if (std::abs(mother.PID) == check_id) return mother;
        position = mother.M1;
    };
    return boost::none;
}

template<typename Lepton>
auto origin(Lepton const& lepton, TTreeReaderArray<GenParticle> const& gen_particles, int check_id) -> boost::optional<GenParticle> {
    for (auto const& particle : get_particles(lepton))
    {
        if (std::abs(particle.PID) == check_id) return particle;
        if (auto mother = origin(gen_particles, particle.M1, check_id)) return *mother;
    }
    return boost::none;
}

template<typename Lepton>
auto id()
{
    print("never end up here");
    return 0;
}

template<>
auto id<Electron>()
{
    return electron_ID;
}

template<>
auto id<Muon>()
{
    return muon_ID;
}

template<>
auto id<Jet>()
{
    return tau_ID;
}

template<typename Lepton>
auto name()
{
    print("never end up here");
    return "never end up here";
}

template<>
auto name<Electron>()
{
    return "electron";
}

template<>
auto name<Muon>()
{
    return "muon";
}

template<>
auto name<Jet>()
{
    return "tau";
}

template<typename Lepton>
auto no_particle(Lepton const&)
{
    print("no", name<Lepton>());
    GenParticle particle;
    particle.X = 0;
    particle.Y = 0;
    return particle;
}

auto constituents(Jet const& jet)
{
    TLorentzVector momentum;
    for (auto pos = 0; pos < jet.Constituents.GetEntriesFast(); ++pos) {
        auto* object = jet.Constituents.At(pos);
        if (!object) continue;
        if (object->IsA() == GenParticle::Class()) momentum += static_cast<GenParticle&>(*object).P4();
        else if (object->IsA() == Track::Class()) momentum += static_cast<Track&>(*object).P4();
        else if (object->IsA() == Tower::Class()) momentum += static_cast<Tower&>(*object).P4();
    }
    return momentum;
}

template<typename Lepton>
auto number_of_tracks(Lepton const&)
{
    return 100;
}

template<>
auto number_of_tracks(Jet const& jet)
{
    auto number = 0;
    for (auto pos = 0; pos < jet.Constituents.GetEntriesFast(); ++pos) {
        auto* object = jet.Constituents.At(pos);
        if (!object) continue;
        if (object->IsA() == Track::Class()) ++number;
    }
    return number;
}

struct Lepton {
    template<typename Input>
    Lepton(Input const& lepton, TTreeReaderArray<GenParticle> const& gen_particles) :
        lorentz_vector(lepton.P4())
        , tracks(number_of_tracks(lepton))
    {
        if (auto mother = origin(lepton, gen_particles, id<Input>())) particle = *mother;
        else particle = no_particle(lepton);
        if(tracks != 100) print(tracks);
    }
    TLorentzVector lorentz_vector;
    GenParticle particle;
    int tracks;
};

auto has_secondary_vertex(Lepton const& lepton)
{
    auto d = transverse_distance(lepton.particle);
    return d > 10. && d < 100.;
}

auto is_hard(Lepton const& lepton)
{
    return lepton.lorentz_vector.Pt() > 25.;
}

auto back_to_back(Lepton const& one, Lepton const& two)
{
    return one.lorentz_vector.DeltaR(two.lorentz_vector) > 4.;
}

auto is_signal(std::vector<Lepton>& leptons)
{
    auto displaced = find_erase(leptons, [](auto const & lepton) {
        return has_secondary_vertex(lepton);
    });
    if (!displaced) return false;
    auto hard = find_erase(leptons, [](auto const & lepton) {
        return is_hard(lepton);
    });
    if (!hard) return false;
    return !back_to_back(*displaced, *hard);
}

template<typename Container>
auto get_leptons(Container const& container, TTreeReaderArray<GenParticle> const& particles)
{
    std::vector<Lepton> leptons;
    for (auto const& lepton : container) leptons.emplace_back(lepton, particles);
    return leptons;
}

auto filter_taus(TTreeReaderArray<Jet> const& jets)
{
    return boost::adaptors::filter(jets, [](auto const & jet) {
        return jet.TauTag;
    });
}

template<typename Predicate>
auto count_if(TTreeReader& reader, Predicate predicate)
{
    auto counter = 0;
    while (reader.Next()) if (predicate()) ++counter;
    return counter;
}

auto count_signals(TTreeReader& reader)
{
    TTreeReaderArray<Electron> electrons(reader, "Electron");
    TTreeReaderArray<Muon> muons(reader, "Muon");
    TTreeReaderArray<Jet> jets(reader, "Jet");
    TTreeReaderArray<GenParticle> particles(reader, "Particle");
    return count_if(reader, [&]() {
        particles.IsEmpty();
        auto leptons = get_leptons(electrons, particles) + get_leptons(muons, particles) + get_leptons(filter_taus(jets), particles);
        return is_signal(leptons);
    });
}

auto efficiency(boost::filesystem::path const& path)
{
    TFile file(delphes_file(path).c_str(), "read");
    TTreeReader reader("Delphes", &file);
    auto entries = reader.GetEntries(false);
    if (entries == 0) {
        print("No events");
        return 0.;
    }
    return count_signals(reader) / static_cast<double>(entries);
}

template<typename Result>
void save_result(Result const& result, std::string const& process)
{
    print(result);
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
    result += " " + std::to_string(efficiency(folder));
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
