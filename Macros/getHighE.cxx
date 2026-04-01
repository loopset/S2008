#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSilData.h"
#include "ActSilSpecs.h"

#include "ROOT/RDataFrame.hxx"

#include <fstream>
#include <memory>
void getHighE()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EReadSilMod};
    dataman.SetRuns(31, 32);
    auto chain {dataman.GetChain()};
    auto chainMerger {dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chainMerger.get());

    ROOT::RDataFrame d {*chain};
    auto df {d.Define("GATCONF", [](ActRoot::ModularData& mod) { return mod.Get("GATCONF"); }, {"ModularData"})};

    auto specs {std::make_shared<ActPhysics::SilSpecs>()};
    specs->ReadFile("../configs/silspecs.conf");

    // Gate on front events poorly reconstructed
    auto dfFilter {df.Filter(
        [&](float gatconf, ActRoot::MergerData& mer, ActRoot::SilData& sil)
        {
            sil.ApplyFinerThresholds(specs);
            if(gatconf == 4 && sil.fSiE["f0"].size() == 1 && mer.fLightIdx == -1)
            {
                return true;
            }
            return false;
        },
        {"GATCONF", "MergerData", "SilData"})};

    std::ofstream streamer {"./Outputs/high_E.dat"};
    dfFilter.Foreach([&](ActRoot::MergerData& mer) { mer.Stream(streamer); }, {"MergerData"});
    streamer.close();
}
