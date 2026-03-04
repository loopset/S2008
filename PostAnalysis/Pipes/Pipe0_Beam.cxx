#include "ActCluster.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"

#include <atomic>
#include <utility>

void Pipe0_Beam(const std::string& beam)
{
    std::string dataconf {"./../configs/data.conf"};

    ROOT::EnableImplicitMT();
    // Read data
    ActRoot::DataManager datman {dataconf, ActRoot::ModeType::EReadSilMod};
    auto chain {datman.GetJoinedData()};
    auto chain2 {datman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain2.get());
    ROOT::RDataFrame df {*chain};

    // Get GATCONF values
    auto defGat {df.Define("GATCONF", [](ActRoot::ModularData& mod)
                           { return static_cast<int>(mod.fLeaves["GATCONF"]); }, {"ModularData"})};

    // Book histograms
    auto hGATCONF {defGat.Histo1D("GATCONF")};

    // And cound CFA triggers
    std::atomic<unsigned long int> cfa {};
    std::atomic<unsigned long int> cfaOther {};
    defGat.Foreach(
        [&](const int& gatconf, ActRoot::MergerData& merger)
        {
            if(gatconf == 64)
            {
                if(merger.fRun == 41 || merger.fRun == 42)
                    cfaOther++;
                else
                    cfa++;
            }
        },
        {"GATCONF", "MergerData"});

    // Draw
    auto* c0 {new TCanvas {"c00", "Pipe 0 canvas 0"}};
    c0->DivideSquare(4);
    c0->cd(1);
    hGATCONF->DrawClone();

    // Print report
    std::cout << "===== GATCONF report =====" << '\n';
    std::cout << "-> CFA/div = " << cfa << '\n';
    std::cout << "-> CFA/div other = " << cfaOther << '\n';
    std::cout << "==========================" << '\n';
}
