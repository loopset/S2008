#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActTypes.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"

#include <fstream>
void gateBadEvents()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};

    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    auto gated {df.Filter(
                      [](ActRoot::MergerData& mer)
                      {
                          if(mer.fLight.GetNLayers() == 1)
                              return (mer.fLight.GetLayer(0) == "f0");
                          else
                              return false;
                      },
                      {"MergerData"})
                    .Define("ESil0", [](ActRoot::MergerData& mer) { return mer.fLight.fEs.front(); }, {"MergerData"})};

    auto hPID {gated.Histo2D({"hGasSil", ";E_{Sil} [MeV];#Delta E_{gas} [arb. units]", 450, 0, 70, 600, 0, 3000},
                             "ESil0", "fLight.fQave")};

    ActRoot::CutsManager<int> cut;
    cut.ReadCut(0, "./Inputs/bad_f0.root");
    auto cutnode {gated.Filter([&](ActRoot::MergerData& mer)
                               { return cut.IsInside(0, mer.fLight.fEs.front(), mer.fLight.fQave); }, {"MergerData"})};
    std::ofstream streamer {"./Outputs/bad_f0.dat"};
    cutnode.Foreach([&](ActRoot::MergerData& mer) { mer.Stream(streamer); }, {"MergerData"});
    streamer.close();

    auto* c0 {new TCanvas {"c0", "Bad event canvas"}};
    hPID->DrawClone("colz");
    cut.DrawAll();
}
