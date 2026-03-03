#include "ActDataManager.h"
#include "ActLine.h"
#include "ActMergerData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TFile.h"
#include "TH1.h"
#include "TStopwatch.h"
#include "TString.h"

#include "Math/Point3Dfwd.h"

#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "../../PostAnalysis/Gates.cxx"

void DistRun()
{
    ROOT::EnableImplicitMT();

    ActRoot::DataManager data {"../../configs/data_all.conf", ActRoot::ModeType::EMerge};
    auto chain {data.GetJoinedData()};

    ROOT::RDataFrame df {*chain};

    // Filter side events
    auto gated {df.Filter(S2008::isFront, {"MergerData"})};

    // Define distances in mm
    double base {256.};
    std::vector<double> dists;
    for(double d = 95; d < 110; d += 1.5)
        dists.push_back(base + d);

    int xbins {200};
    std::pair<double, double> xlims {-20, 300};
    int zbins {200};
    std::pair<double, double> zlims {-20, 300};
    // Save
    auto f {std::make_unique<TFile>("./Outputs/histos.root", "recreate")};
    f->WriteObject(&dists, "dists");
    TStopwatch timer {};
    timer.Start();
    for(const auto& dist : dists)
    {
        std::cout << "Distance : " << dist << '\n';
        // Redefine
        auto node {gated.Define("NewSP",
                                [&, dist](ActRoot::MergerData& d)
                                {
                                    // Reconstruct line from BP and SP
                                    auto p {d.fBP};
                                    auto dir {(d.fSP - d.fBP)};
                                    ActRoot::Line line {p, dir, 0};
                                    return line.MoveToX(dist);
                                },
                                {"MergerData"})};
        // Fill histograms!
        ROOT::TThreadedObject<TH2D> hSP {ROOT::TNumSlots {node.GetNSlots()},
                                         "hSP",
                                         TString::Format("Side %.2f mm;Y [mm];Z [mm]", dist),
                                         xbins,
                                         xlims.first,
                                         xlims.second,
                                         zbins,
                                         zlims.first,
                                         zlims.second};
        std::map<int, ROOT::TThreadedObject<TH1D>> pxs, pzs;
        std::vector<int> idxs {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
        for(const auto& idx : idxs)
        {
            pxs.emplace(std::piecewise_construct, std::forward_as_tuple(idx),
                        std::forward_as_tuple(ROOT::TNumSlots {node.GetNSlots()}, TString::Format("px%d", idx),
                                              TString::Format("%.2f mm X proj %d;Y [mm]", dist, idx), xbins,
                                              xlims.first, xlims.second));
            pzs.emplace(std::piecewise_construct, std::forward_as_tuple(idx),
                        std::forward_as_tuple(ROOT::TNumSlots {node.GetNSlots()}, TString::Format("pz%d", idx),
                                              TString::Format("%.2f mm Z proj %d;Z [mm]", dist, idx), zbins,
                                              zlims.first, zlims.second));
        }
        node.ForeachSlot(
            [&](unsigned int slot, ActRoot::MergerData& data, ROOT::Math::XYZPointF& sp)
            {
                hSP.GetAtSlot(slot)->Fill(sp.Y(), sp.Z());
                auto idx {data.fSilNs.front()};
                if(pxs.count(idx))
                {
                    pxs[idx].GetAtSlot(slot)->Fill(sp.Y());
                    pzs[idx].GetAtSlot(slot)->Fill(sp.Z());
                }
            },
            {"MergerData", "NewSP"});

        // Write data
        f->cd();
        auto path {TString::Format("d_%.1f_mm/", dist)};
        auto* dir {f->mkdir(path)};
        dir->cd();
        hSP.Merge()->Write();
        for(auto m : {&pxs, &pzs})
            for(auto& [_, h] : *m)
                h.Merge()->Write();
        f->cd();
    }
    timer.Stop();
    timer.Print();
}
