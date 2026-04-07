#include "TString.h"
#include "TSystem.h"

#include <iostream>
#include <vector>

void RunItMinVoxels()
{
    std::vector<int> minvoxels {5, 10, 15, 20, 11};
    int idx {};
    for(const auto& min : minvoxels)
    {
        std::cout << "Run for : " << min << " min voxels" << '\n';
        // 1-> Modify ma.conf
        gSystem->cd("/media/Data/S2008/configs/");
        auto sed0 {TString::Format(
            "sed -i '/\\[RecRANSAC\\]/,/^[[]/ s/MinVoxels: [0-9]*/MinVoxels: %d/' multiaction.conf", min)};
        gSystem->Exec(sed0);
        if(idx == 0) // in first it, minvoxels of ransac is < cleandeltas; change the latter to min
        {
            auto sed1 {TString::Format(
                "sed -i '/\\[CleanDeltas\\]/,/^[[]/ s/MaxVoxels: [0-9]*/MaxVoxels: %d/' multiaction.conf", min)};
            gSystem->Exec(sed1);
        }

        // 2-> Rerun filter and merger
        gSystem->cd("/media/Data/S2008/");
        gSystem->Exec("actroot -f && actroot -m");

        if(idx == 0) // undo change in CleanDeltas for first it
        {
            gSystem->cd("/media/Data/S2008/configs/");
            // ensure the hardcoded value agrees with the wanted one in ma.conf
            auto sed1 {TString::Format(
                "sed -i '/\\[CleanDeltas\\]/,/^[[]/ s/MaxVoxels: [0-9]*/MaxVoxels: %d/' multiaction.conf", 10)};
            gSystem->Exec(sed1);
        }

        // 3-> Run pipelines
        gSystem->cd("/media/Data/S2008/PostAnalysis/");
        gSystem->Exec("root -l -b -x -q 'Runner.cxx(\"12\")'");

        // 4-> Copy outfiles to dir
        gSystem->cd("/media/Data/S2008/");
        auto cp {TString::Format("cp /media/Data/S2008/PostAnalysis/Outputs/tree_ex_20Mg_p_p_sil.root "
                                 "/media/Data/S2008/Macros/RANSAC/Inputs/tree_sil_minvoxels_%02d.root",
                                 min)};
        gSystem->Exec(cp);

        idx++;
    }
}
