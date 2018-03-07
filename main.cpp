#include <iostream>
#include <array>

#include <boost/range/irange.hpp>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/transform.hpp>

// #include "TH1.h"
// #include "TSystem.h"
// #include "THStack.h"
// #include "TLegend.h"
// #include "TPaveText.h"
#include "TClonesArray.h"
// #include "TApplication.h"

#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
// #include "ExRootAnalysis/ExRootResult.h"

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

// struct MyPlots {
//     TH1* fJetPT[2];
//     TH1* fMissingET;
//     TH1* fElectronPT;
// };
//
// auto BookHistograms(ExRootResult& result)
// {
//
//     MyPlots plots;
//     // book 2 histograms for PT of 1st and 2nd leading jets
//     plots.fJetPT[0] = result.AddHist1D(
//                           "jet_pt_0", "leading jet P_{T}",
//                           "jet P_{T}, GeV/c", "number of jets",
//                           50, 0.0, 100.0);
//
//     plots.fJetPT[1] = result.AddHist1D(
//                           "jet_pt_1", "2nd leading jet P_{T}",
//                           "jet P_{T}, GeV/c", "number of jets",
//                           50, 0.0, 100.0);
//
//     plots.fJetPT[0]->SetLineColor(kRed);
//     plots.fJetPT[1]->SetLineColor(kBlue);
//
//     // book 1 stack of 2 histograms
//     auto* stack = result.AddHistStack("jet_pt_all", "1st and 2nd jets P_{T}");
//     stack->Add(plots.fJetPT[0]);
//     stack->Add(plots.fJetPT[1]);
//
//     // book legend for stack of 2 histograms
//     auto* legend = result.AddLegend(0.72, 0.86, 0.98, 0.98);
//     legend->AddEntry(plots.fJetPT[0], "leading jet", "l");
//     legend->AddEntry(plots.fJetPT[1], "second jet", "l");
//
//     // attach legend to stack (legend will be printed over stack in .eps file)
//     result.Attach(stack, legend);
//
//     // book more histograms
//     plots.fElectronPT = result.AddHist1D(
//                             "electron_pt", "electron P_{T}",
//                             "electron P_{T}, GeV/c", "number of electrons",
//                             50, 0.0, 100.0);
//
//     plots.fMissingET = result.AddHist1D(
//                            "missing_et", "Missing E_{T}",
//                            "Missing E_{T}, GeV", "number of events",
//                            60, 0.0, 30.0);
//
//     // book general comment
//
//     auto* comment = result.AddComment(0.64, 0.86, 0.98, 0.98);
//     comment->AddText("demonstration plot");
//     comment->AddText("produced by analyze.C");
//
//     // attach comment to single histograms
//
//     result.Attach(plots.fJetPT[0], comment);
//     result.Attach(plots.fJetPT[1], comment);
//     result.Attach(plots.fElectronPT, comment);
//
//     // show histogram statisics for MissingET
//     plots.fMissingET->SetStats();
//     return plots;
// }

auto distance(GenParticle const& particle)
{
    return std::sqrt(sqr(particle.X) + sqr(particle.Y));
}

void AnalyseEvents(ExRootTreeReader& tree_reader)
{
    auto& muon_branch = *tree_reader.UseBranch("Muon");
    // Loop over all events
    for (auto entry : range(tree_reader.GetEntries())) {
        // Load selected branches with data from specified event
        tree_reader.ReadEntry(entry);
        std::vector<double> result;
        for (auto number : range(muon_branch.GetEntriesFast())) {
            auto& muon = static_cast<Muon&>(*muon_branch.At(number));
            auto& particle = static_cast<GenParticle&>(*muon.Particle.GetObject());
            auto dist = distance(particle);
            if (dist > 10 && dist < 200) result.emplace_back(dist);
        }
        print(result);
    }
}

// void PrintHistograms(ExRootResult& result)
// {
//     result.Print("png");
// }

auto analyze(std::string const& inputFile)
{
    TChain chain("Delphes");
    chain.Add(inputFile.c_str());

//     ExRootResult result;
//     auto plots = BookHistograms(result);

    ExRootTreeReader tree_reader(&chain);
    AnalyseEvents(tree_reader);

//     PrintHistograms(result, plots);

//     result.Write("results.root");
    return 1.5;
}

auto file_name(int number)
{
    auto path = "scan/Events/";
    auto name = (number < 10 ? "0" : "") + std::to_string(number);
    auto folder = "run_" + name + "_decayed_1/";
    auto file = "tag_1_delphes_events.root";
    return path + folder + file;
}

int main()
{
    auto result = transform(range(2), [](auto number) {
        return analyze(file_name(number));
    });
//     print(result);
}
