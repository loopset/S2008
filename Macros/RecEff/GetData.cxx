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

#include <fstream>
#include <memory>
#include <vector>

void GetData()
{
    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EReadSilMod};
    dataman.SetRuns(31, 33);
    auto chain {dataman.GetChain()};
    auto chainTPC {dataman.GetChain(ActRoot::ModeType::EReadTPC)};
    chain->AddFriend(chainTPC.get());
    auto chainMerger {dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chainMerger.get());

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
                         {"GATCONF", "TPCData", "MergerData", "SilData"})
                       .Define("X",
                               [](ActRoot::TPCData& tpc)
                               {
                                   std::vector<int> ret;
                                   // Clusters
                                   for(const auto& cl : tpc.fClusters)
                                   {
                                       for(const auto& v : cl.GetVoxels())
                                       {
                                           const auto& p {v.GetPosition()};
                                           ret.push_back(p.X());
                                       }
                                   }
                                   // Noise
                                   for(const auto& v : tpc.fRaw)
                                   {
                                       const auto& p {v.GetPosition()};
                                       ret.push_back(p.X());
                                   }
                                   return ret;
                               },
                               {"TPCData"})
                       .Define("Y",
                               [](ActRoot::TPCData& tpc)
                               {
                                   std::vector<int> ret;
                                   // Clusters
                                   for(const auto& cl : tpc.fClusters)
                                   {
                                       for(const auto& v : cl.GetVoxels())
                                       {
                                           const auto& p {v.GetPosition()};
                                           ret.push_back(p.Y());
                                       }
                                   }
                                   // Noise
                                   for(const auto& v : tpc.fRaw)
                                   {
                                       const auto& p {v.GetPosition()};
                                       ret.push_back(p.Y());
                                   }
                                   return ret;
                               },
                               {"TPCData"})
                       .Define("Z",
                               [](ActRoot::TPCData& tpc)
                               {
                                   std::vector<int> ret;
                                   // Clusters
                                   for(const auto& cl : tpc.fClusters)
                                   {
                                       for(const auto& v : cl.GetVoxels())
                                       {
                                           const auto& p {v.GetPosition()};
                                           ret.push_back(p.Z());
                                       }
                                   }
                                   // Noise
                                   for(const auto& v : tpc.fRaw)
                                   {
                                       const auto& p {v.GetPosition()};
                                       ret.push_back(p.Z());
                                   }
                                   return ret;
                               },
                               {"TPCData"})
                       .Define("Ef0",
                               [&](ActRoot::SilData& sil)
                               {
                                   sil.ApplyFinerThresholds(specs);
                                   return sil.fSiE["f0"].front();
                               },
                               {"SilData"})};

    std::cout << "Entries : " << dfFilter.Count().GetValue() << '\n';
    dfFilter.Snapshot("SimpleTree", "./Outputs/simple_tree.root", {"X", "Y", "Z", "Ef0", "fRun", "fEntry"});
    // std::ofstream streamer {"./Inputs/events_f0.dat"};
    // dfFilter.Foreach([&](ActRoot::MergerData& mer) { mer.Stream(streamer); }, {"MergerData"});
    // streamer.close();
}
