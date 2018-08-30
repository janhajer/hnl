#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/cxx11/any_of.hpp>
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
#include <boost/optional/optional_io.hpp>

#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

#include "classes/DelphesClasses.h"

#include "alphanum.hpp"

using namespace std::string_literals;

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
    return stream << "ID: " << particle.PID << ", D0: " << transverse_distance(particle) << " mm" << " Status: " << particle.Status;
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

auto const neutrino_ID = 9900012;
auto const electron_ID = 11;
auto const muon_ID = 13;
auto const tau_ID = 15;

template<typename Lepton>
auto ids() -> std::vector<int> {
    print("never end up here");
    return {0};
}

template<>
auto ids<Electron>() -> std::vector<int> {
    return {electron_ID};
}

template<>
auto ids<Muon>() -> std::vector<int> {
    return {muon_ID};
}

template<>
auto ids<Jet>() -> std::vector<int> {
    return {tau_ID};
}

template<>
auto ids<Track>() -> std::vector<int> {
    return {tau_ID, muon_ID, electron_ID};
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

auto has_tag(std::string const& string, std::string const& tag)
{
    return string.find(tag) != std::string::npos;
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

struct Path {
    boost::filesystem::path base;
    std::string tag;
};

auto is_regular_file = static_cast<bool (*)(boost::filesystem::path const&)>(&boost::filesystem::is_regular_file);

auto is_delphes = [](Path const& paths)
{
    return has_ending(paths.base.filename().string(), "_delphes_events.root") && has_tag(paths.base.filename().string(), paths.tag);
};

auto delphes_file(Path const& paths)
{
    return get_file(paths.base, [&paths](auto const & files) {
        return files | boost::adaptors::filtered(is_regular_file) | boost::adaptors::filtered([&paths](auto const & path) {
            return is_delphes({path, paths.tag});
        });
    });
}

auto is_banner = [](Path const& paths)
{
    return has_ending(paths.base.filename().string(), "_banner.txt") && has_tag(paths.base.filename().string(), paths.tag);
};

auto banner_file(Path const& paths)
{
    return get_file(paths.base, [&paths](auto const & files) {
        return files | boost::adaptors::filtered(is_regular_file) | boost::adaptors::filtered([&paths](auto const & path) {
            return is_banner({path, paths.tag});
        });
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

auto get_xsec(Path const& paths)
{
    return read_file(banner_file(paths), [](auto const & strings) {
        return strings.size() > 4 && strings.at(0) == "#" && strings.at(1) == "Integrated" && strings.at(2) == "weight" && strings.at(3) == "(pb)" && strings.at(4) == ":";
    }, 5);
}

auto get_mass(Path const& paths)
{
    return read_file(banner_file(paths), [](auto const & strings) {
        return strings.size() > 2 && strings.at(0) == std::to_string(neutrino_ID) && strings.at(2) == "#" && strings.at(3) == "mn1";
    }, 1);
}

auto get_coupling(Path const& paths, int pos, std::string const& name)
{
    return read_file(banner_file(paths), [&name, pos](auto const & strings) {
        return strings.size() > 3 && strings.at(0) == std::to_string(pos) && strings.at(2) == "#" && strings.at(3) == name;
    }, 1);
}

auto get_e_coupling(Path const& paths)
{
    return get_coupling(paths, 1, "ven1");
}

auto get_mu_coupling(Path const& paths)
{
    return get_coupling(paths, 4, "vmun1");
}

auto get_tau_coupling(Path const& paths)
{
    return get_coupling(paths, 7, "vtan1");
}

auto get_width(Path const& paths)
{
    return read_file(banner_file(paths), [](auto const & strings) {
        return strings.size() > 2 && strings.at(0) == "DECAY" && strings.at(1) == std::to_string(neutrino_ID);
    }, 2);
}

template<typename Lepton>
auto get_particles(Lepton const& lepton) -> std::vector<GenParticle> {
    std::vector<GenParticle> particles;
    if (auto* particle = static_cast<GenParticle*>(lepton.Particle.GetObject())) particles.emplace_back(*particle);
    return particles;
}

template<typename Function>
auto for_each(TRefArray const& ref_array, Function function)
{
    if (auto* iterator = static_cast<TRefArrayIter*>(ref_array.MakeIterator())) while (auto* object = iterator->Next()) function(*object);
}

template<typename Particle>
auto cast(TObject const& object) -> boost::optional<Particle> {
return object.IsA() == Particle::Class() ? boost::optional<Particle>{static_cast<Particle const&>(object)} : boost::optional<Particle>{};
}

template<typename Particle, typename Function>
auto for_each(TRefArray const& ref_array, Function function)
{
    for_each(ref_array, [function](auto const & object) {
        if (auto optional = cast<Particle>(object)) function(*optional);
    });
}

template<>
auto get_particles(Jet const& jet) -> std::vector<GenParticle> {
    std::vector<GenParticle> particles;
    for_each<GenParticle>(jet.Particles, [&particles](auto & particle)
    {
        particles.emplace_back(particle);
    });
    return particles;
}

auto origin(TTreeReaderArray<GenParticle> const& particles, int position, std::vector<int> const& check_ids) -> boost::optional<GenParticle> {
    while (position >= 0 && position < particles.GetSize())
    {
        auto& particle = particles.At(position);
//         print(particle);
        if (boost::algorithm::any_of_equal(check_ids, std::abs(particle.PID))) return particle;
//         if (particle.M2 >= 0) if(auto mother_2 = origin(particles, particle.M2, check_ids)) return mother_2;
        position = particle.M1;
    };
    return boost::none;
}

template<typename Lepton>
auto origin(Lepton const& lepton, TTreeReaderArray<GenParticle> const& particles, std::vector<int> const& check_ids) -> boost::optional<GenParticle> {
    for (auto const& particle : get_particles(lepton))
    {
        if (boost::algorithm::any_of_equal(check_ids, std::abs(particle.PID))) return particle;
        if (auto mother = origin(particles, particle.M1, check_ids)) return mother;
    }
    return boost::none;
}

auto min_displacement()
{
    return 5.;
}

auto is_particle(TTreeReaderArray<GenParticle> const& particles, GenParticle const& particle, std::vector<int> const& check_ids)
{
    return boost::algorithm::any_of_equal(check_ids, std::abs(particle.PID)) && std::abs(particles.At(particle.M1).PID) == neutrino_ID;
}

auto origin2(TTreeReaderArray<GenParticle> const& particles, int position, std::vector<int> const& check_ids) -> boost::optional<GenParticle> {
    while (position >= 0 && position < particles.GetSize())
    {
        auto& particle = particles.At(position);
        if (is_particle(particles, particle, check_ids)) return particle;
//         if (particle.M2 >= 0) if(auto mother_2 = origin2(particles, particle.M2, check_ids)) return mother_2;
        position = particle.M1;
    };
    return boost::none;
}

template<typename Lepton>
auto origin2(Lepton const& lepton, TTreeReaderArray<GenParticle> const& particles, std::vector<int> const& check_ids) -> boost::optional<GenParticle> {
    for (auto const& particle : get_particles(lepton))
    {
        auto& M1 = particles.At(particle.M1);
        if (is_particle(particles, particle, check_ids)) return particle;
        if (auto mother = origin2(particles, particle.M1, check_ids)) return mother;
    }
    return boost::none;
}

auto tree(TTreeReaderArray<GenParticle> const& particles, int position) -> std::vector<GenParticle> {
    std::vector<GenParticle> vector;
    while (position >= 0 && position < particles.GetSize())
    {
        auto& particle = particles.At(position);
        vector.emplace_back(particle);
        position = particle.M1;
    };
    return vector;
}

template<typename Lepton>
auto tree(Lepton const& lepton, TTreeReaderArray<GenParticle> const& particles) -> std::vector<GenParticle> {
//     print("New tree");
    std::vector<GenParticle> vector;
    for (auto const& particle : get_particles(lepton))
    {
        vector.emplace_back(particle);
        vector += tree(particles, particle.M1);
    }
    return vector;
}

template<typename Function>
auto for_each_constituent(TRefArray const& ref_array, Function function)
{
    for_each(ref_array, [function](auto & object) {
        if (auto optional = cast<GenParticle>(object)) return function(*optional);
        if (auto optional = cast<Track>(object)) return function(*optional);
        if (auto optional = cast<Tower>(object)) return function(*optional);
        print("Unexpected Constituent");
    });
}

auto constituents(Jet const& jet)
{
    TLorentzVector momentum;
    for_each_constituent(jet.Constituents, [&momentum](auto & object) {
        momentum += object.P4();
    });
    return momentum;
}

struct Lepton {
    template<typename Input>
    Lepton(Input const& lepton, TTreeReaderArray<GenParticle> const& particles) :
        lorentz_vector(lepton.P4()),
        particle(origin2(lepton, particles, ids<Input>())),
        mother(origin(lepton, particles, {neutrino_ID})),
           charge(lepton.Charge),
           tree(::tree(lepton, particles))
    {
//         print(*this);
    }
    TLorentzVector lorentz_vector;
    boost::optional<GenParticle> particle;
    boost::optional<GenParticle> mother;
    int charge;
    std::vector<GenParticle> tree;
};

std::ostream& operator<<(std::ostream& stream, TLorentzVector const& lepton)
{
    return stream << "Pt: " << lepton.Pt() << " GeV";
}

std::ostream& operator<<(std::ostream& stream, Lepton const& lepton)
{
    return stream << "particle: " << lepton.particle << " mother: " << lepton.mother << ", " << lepton.lorentz_vector << ", Charge: " << lepton.charge;
}

auto has_secondary_vertex(Lepton const& lepton)
{
    auto d = lepton.particle ? transverse_distance(*lepton.particle) : 0;
    return d > min_displacement() && d < 100.;
//     return d > min_displacement() && d < 100. && lepton.mother && std::abs(lepton.mother->PID) == neutrino_ID && lepton.lorentz_vector.Pt() > 5.;
}

auto is_hard(Lepton const& lepton)
{
    return lepton.lorentz_vector.Pt() > 25.;
//     && (lepton.particle ? transverse_distance(*lepton.particle) : 0) < 5. && lepton.mother && lepton.mother->PID != std::abs(neutrino_ID);
}

auto back_to_back(Lepton const& one, Lepton const& two)
{
    return one.lorentz_vector.DeltaR(two.lorentz_vector) > 4.;
}

auto is_displaced_signal(std::vector<Lepton>& leptons)
{
    auto displaced = find_erase(leptons, [](auto const & lepton) {
        return has_secondary_vertex(lepton);
    });
    if (!displaced) return false;
    auto hard = find_erase(leptons, [](auto const & lepton) {
        return is_hard(lepton);
    });
    if (!hard) return false;
    auto good = !back_to_back(*displaced, *hard);
    if (!good) return false;
//     print("Displaced: (", *displaced, ")  Hard: (", *hard, ")");
//     print(displaced->tree);
    return true;
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

template<typename Predicate>
auto count_events_if(TTreeReader& reader, Predicate predicate)
{
//     TTreeReaderArray<Electron> electrons(reader, "Electron");
//     TTreeReaderArray<Muon> muons(reader, "MuonLoose");
//     TTreeReaderArray<Jet> jets(reader, "Jet");
    TTreeReaderArray<GenParticle> particles(reader, "Particle");
    TTreeReaderArray<Track> tracks(reader, "Track");
    return count_if(reader, [&]() {
        particles.IsEmpty();
        auto leptons = get_leptons(tracks, particles);
//         auto leptons = get_leptons(electrons, particles) + get_leptons(muons, particles) + get_leptons(filter_taus(jets), particles);
        return predicate(leptons);
    });
}

auto count_displaced_events(TTreeReader& reader)
{
    return count_events_if(reader, [](auto & leptons) {
        return is_displaced_signal(leptons);
    });
}

template<typename Predicate>
auto events(Path const& paths, Predicate predicate)
{
    TFile file(delphes_file(paths).c_str(), "read");
    TTreeReader reader("Delphes", &file);
    return predicate(reader);
}

auto displaced_events(Path const& paths)
{
    return events(paths, [](auto & reader) {
        return count_displaced_events(reader);
    });
}

auto same_sign(Lepton const& one, Lepton const& two)
{
    return one.charge == two.charge && one.charge > 0;
}

auto is_prompt_signal(std::vector<Lepton>& leptons)
{
    auto hard = find_erase(leptons, [](auto const & lepton) {
        return is_hard(lepton);
    });
    if (!hard) return false;
    auto second = find_erase(leptons, [](auto const & lepton) {
        return is_hard(lepton);
    });
    if (!second) return false;
    return !back_to_back(*second, *hard) && same_sign(*second, *hard);
}

auto count_prompt_events(TTreeReader& reader)
{
    return count_events_if(reader, [](auto & leptons) {
        return is_prompt_signal(leptons);
    });
}

auto prompt_events(Path const& paths)
{
    return events(paths, [](auto & reader) {
        return count_prompt_events(reader);
    });
}

auto all_events(Path const& paths)
{
    return events(paths, [](auto & reader) {
        return reader.GetEntries(false);
    });
}

template<typename Result>
void save_result(Result const& result, std::string const& process)
{
//     print(result);
    std::ofstream file("./" + process + ".dat");
    std::ostream_iterator<std::string> iterator(file, "\n");
    boost::copy(result, iterator);
}

auto get_result(Path const& paths)
{
//     if (get_mass(paths) != "5.000000e+01" && get_mass(paths) != "1.000000e+00") return ""s;
//     if (get_mu_coupling(paths) != "1.000000e-01") return ""s;
    auto result = get_mass(paths);
    result += " " + get_e_coupling(paths);
    result += " " + get_mu_coupling(paths);
    result += " " + get_tau_coupling(paths);
    result += " " + std::to_string(all_events(paths));
    result += " " + std::to_string(displaced_events(paths));
    result += " " + std::to_string(prompt_events(paths));
    result += " " + get_xsec(paths);
    result += " " + get_width(paths);
    print(result);
    return result;
}

auto get_header()
{
    std::string header = "# mass";
    header += " e_coupling";
    header += " mu_coupling";
    header += " tau_coupling";
    header += " total";
    header += " displaced";
    header += " prompt";
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
    auto tag = arguments.size() > 2 ? arguments.at(2) : "tag_1"s;
    std::vector<std::string> results {get_header()};
    for (auto const& folder : decayed_folders(event_folder(process))) results.emplace_back(get_result({folder, tag}));
    save_result(results, process);
}
