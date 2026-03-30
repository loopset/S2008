#include "ActDataManager.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "TFile.h"

void get()
{
    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EReadTPC};
    dataman.SetRuns(31, 31);
    auto chain {dataman.GetChain()};

    ActRoot::TPCData* tpc {};
    chain->SetBranchAddress("TPCData", &tpc);
    chain->GetEntry(26);
    tpc->Print();

    // Write
    auto* fout {new TFile {"./Inputs/event_ransac_r31_entry26.root", "recreate"}};
    fout->WriteObject(tpc, "TPCData");
    fout->Close();
}
