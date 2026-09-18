#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSilData.h"
#include "ActSilSpecs.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

#include "TCanvas.h"

#include <fstream>
#include <memory>
#include <vector>

void EstimateFromRange()
{
    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EReadSilMod};
    dataman.SetRuns(31, 33);
    auto chain {dataman.GetChain()};
    auto chainTPC {dataman.GetChain(ActRoot::ModeType::EReadTPC)};
    chain->AddFriend(chainTPC.get());
    auto chainMerger {dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chainMerger.get());

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame d {*chain};
    auto df {d.Define("GATCONF", [](ActRoot::ModularData& mod) { return mod.Get("GATCONF"); }, {"ModularData"})};

    auto specs {std::make_shared<ActPhysics::SilSpecs>()};
    specs->ReadFile("../../configs/silspecs.conf");

    // Gate on front events
    auto dfFilter {df.Filter(
        [&](float& gatconf, ActRoot::TPCData& tpc, ActRoot::MergerData& mer, ActRoot::SilData& sil)
        {
            // if(tpc.fClusters.size() == 0)
            //     return false;
            sil.ApplyFinerThresholds(specs);
            if(gatconf == 4 && sil.fSiE["f0"].size() == 1 && sil.fSiE["f1"].size() == 0)
                // if(sil.fSiE["f0"].front() >= 10)
                return true;
            return false;
        },
        {"GATCONF", "TPCData", "MergerData", "SilData"})};

    // Find last voxel along beam direction
    auto def {dfFilter
                  .Define("LastVoxel",
                          [](ActRoot::TPCData& tpc)
                          {
                              double lastX {-1111};
                              for(const auto& cl : tpc.fClusters)
                              {
                                  for(const auto& v : cl.GetVoxels())
                                  {
                                      const auto& p {v.GetPosition()};
                                      // Only for beam row
                                      if(61 <= p.Y() && p.Y() <= 67)
                                      {
                                          if(p.X() > lastX)
                                              lastX = p.X();
                                      }
                                  }
                              }
                              return lastX;
                          },
                          {"TPCData"})
                  .Define("Ef0", [](ActRoot::SilData& sil) { return sil.fSiE["f0"].front(); }, {"SilData"})};

    // Plot
    auto hLastESil {
        def.Histo2D({"hLastESil", "Last beam X vs E_{sil};Last X [pad];E_{Sil} [MeV]", 128, 0, 128, 400, 0, 40},
                    "LastVoxel", "Ef0")};

    // Draw
    auto* c0 {new TCanvas {"c0", "eff estimation canvas"}};
    c0->DivideSquare(4);
    c0->cd(1);
    hLastESil->DrawClone("colz");
}
