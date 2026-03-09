#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

#include <string>
#include <thread>
#include <vector>

#include "./Plotter.cxx"
#include "./Simulation_S2008.cxx"
// what
// if simu = runs simulation
// if plot = plots results
// standalone (only applies to simu setting):
// if true, runs only first item in Ex vector and plots in-simulation results
// if false, runs all Ex simulations but doesn't plot
void Runner(TString what = "plot", bool standalone = true)
{
    // Settings
    // Names of particles
    std::string beam {"20Mg"};
    std::string target {"p"};
    std::string light {"p"};
    int neutronPS {0};
    int protonPS {0}; // number of protons in final state
    double T1 {4.24}; // Mev/u

    bool isPS {true};
    std::vector<double> Eexs {0};

    ROOT::EnableThreadSafety();
    std::vector<std::thread> threads;
    auto worker {[](TString str) { return gSystem->Exec(str); }};

    if(what.Contains("simu"))
    {
        Simulation_S2008(beam, target, light, neutronPS, protonPS, T1, Eexs.front(), standalone);
        // if(standalone)
        // {
        // }
        // else
        // {
        // TString haddlist {};
        // TString haddout {};
        // if(isPS)
        // {
        //     haddout = gSelector->GetSimuFile(Eexs.front(), neutronPS, protonPS);
        //     // List of files generated per thread
        //     // Number of threads = 6
        //     int nthreads {6};
        //     // 1e8 events each
        //     for(int i = 1; i <= nthreads; i++)
        //     {
        //         auto str {TString::Format(
        //             "root -l -b -x -q \'Simulation_E796.cpp(\"%s\",\"%s\",\"%s\",%d,%d,%f,%f,%d,%d)\'",
        //             beam.c_str(), target.c_str(), light.c_str(), neutronPS, protonPS, T1, Eexs.front(),
        //             standalone, i)};
        //         gSelector->SetTag(std::to_string(i));
        //         haddlist += gSelector->GetSimuFile(Eexs.front(), neutronPS, protonPS) + " ";
        //         threads.emplace_back(worker, str);
        //     }
        // }
        // else
        // {
        //     for(const auto& Eex : Eexs)
        //     {
        //         auto str {TString::Format(
        //             "root -l -b -x -q \'Simulation_E796.cpp(\"%s\",\"%s\",\"%s\",%d,%d,%f,%f,%d)\'",
        //             beam.c_str(), target.c_str(), light.c_str(), neutronPS, protonPS, T1, Eex, standalone)};
        //         threads.emplace_back(worker, str);
        //     }
        // }
        // for(auto& thread : threads)
        //     thread.join();
        // // Once finished
        // if(isPS)
        // {
        //     std::cout << "Output file: " << haddout << '\n';
        //     std::cout << "Input files: " << haddlist << '\n';
        //     // Merge
        //     gSystem->Exec(TString::Format("hadd -f %s %s", haddout.Data(), haddlist.Data()));
        //     // Remove
        //     gSystem->Exec(TString::Format("rm %s", haddlist.Data()));
        // }
        // }
    }
    if(what.Contains("plot"))
    {
        Plotter(beam, target, light, Eexs.front(), T1, neutronPS, protonPS);
    }
}
