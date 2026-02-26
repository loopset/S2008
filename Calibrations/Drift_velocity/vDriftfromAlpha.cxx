#include "ActCluster.h"
#include "ActCutsManager.h"
#include "ActTPCData.h"
#include "ActVoxel.h"

#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RResultPtr.hxx"
#include "ROOT/TThreadedObject.hxx"
#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TMathBase.h"
#include "TPaveText.h"
#include "TString.h"

#include "Math/Point3Dfwd.h"

#include <algorithm>
#include <numeric>
#include <utility>
#include <vector>

struct RetValues
{
    TGraph* fGraphLinear {};
    std::pair<double, double> fDrift {};
};

RetValues run_for(int run = 0, bool draw = false)
{
    auto df {ROOT::RDataFrame("GETTree", TString::Format("../../RootFiles/Cluster/Clusters_Run_%04d.root", run))};
    // Gate on events with only one cluster
    auto gated {df.Filter(
        [](ActRoot::TPCData& data)
        {
            auto size {data.fClusters.size() == 1};
            if(!size)
                return false;
            return data.fClusters.front().GetSizeOfVoxels() >= 30;
        },
        {"TPCData"})};

    // Define points
    auto def {gated
                  .Define("LastPoint",
                          [](ActRoot::TPCData& d)
                          {
                              auto cluster {d.fClusters.front()};
                              auto line {cluster.GetRefToLine()};
                              cluster.SortAlongDir();
                              auto lastVoxel {cluster.GetRefToVoxels().back()};
                              auto proj {line.ProjectionPointOnLine(lastVoxel.GetPosition())};
                              return proj;
                          },
                          {"TPCData"})
                  .Define("OtherPoint",
                          [](ActRoot::TPCData& d)
                          {
                              auto cluster {d.fClusters.front()};
                              auto line {cluster.GetRefToLine()};
                              auto otherPoint {line.MoveToX(-50)};
                              return otherPoint;
                          },
                          {"TPCData"})};


    // Source position
    std::pair<double, double> sourceAt {-29.4, 39.2}; // pads

    auto defDrift {def.Define("DeltaT",
                              [&](ActRoot::TPCData& data, ROOT::Math::XYZPointF& last)
                              {
                                  // Get Z at source position
                                  auto& line {data.fClusters.front().GetLine()};
                                  auto source {line.MoveToX(sourceAt.first)};
                                  auto deltaT {last.Z() - source.Z()};
                                  // Convert from btb to micros
                                  deltaT *= 0.32; // 1 btb = 4tb; 1tb = 12.5 MHz
                                  return deltaT;
                              },
                              {"TPCData", "LastPoint"})
                       .Define("Lxy",
                               [&](ROOT::Math::XYZPointF& last)
                               {
                                   auto lastInPlane {last};
                                   last.SetZ(0);
                                   ROOT::Math::XYZPointF source {(float)sourceAt.first, (float)sourceAt.second, 0};
                                   return (last - source).R() * 2 / 10;
                               },
                               {"LastPoint"})
                       .Define("DeltaT2", "DeltaT * DeltaT")
                       .Define("Lxy2", "Lxy * Lxy")

    };

    // Read cuts
    ActRoot::CutsManager<int> cuts;
    cuts.ReadCut(0, TString::Format("./Cuts/low_%d.root", run).Data());
    cuts.ReadCut(1, TString::Format("./Cuts/middle_%d.root", run).Data());
    cuts.ReadCut(2, TString::Format("./Cuts/up_%d.root", run).Data());

    // And apply filters
    std::vector<ROOT::RDF::RNode> nodes;
    std::vector<ROOT::RDF::RResultPtr<TGraph>> gs;
    for(int i = 0; i < 3; i++)
    {
        auto node {defDrift.Filter([i, &cuts](float& deltat2, float& lxy2) { return cuts.IsInside(i, deltat2, lxy2); },
                                   {"DeltaT2", "Lxy2"})};
        auto gnode {node.Graph("DeltaT2", "Lxy2")};
        nodes.push_back(node);
        gs.push_back(gnode);
    }
    // Fits
    std::vector<TF1*> fits;
    std::vector<double> slopes {};
    std::vector<double> uslopes {};
    int idx {};
    for(auto& g : gs)
    {
        g->Fit("pol1", "0QM");
        auto fit {g->GetFunction("pol1")};
        if(fit)
        {
            fits.push_back((TF1*)fit->Clone());
            auto& f {fits.back()};
            f->SetTitle(TString::Format("fit %d", idx));
            f->SetLineColor(6 + idx);
            f->ResetBit(TF1::kNotDraw);
            std::cout << "Idx : " << idx << " Slope : " << fits.back()->GetParameter(1) << '\n';
            slopes.push_back(f->GetParameter(1));
            uslopes.push_back(f->GetParError(1));
        }
        idx++;
    }

    // Plot the DeltaZ and lxy
    auto gDrift {defDrift.Graph("DeltaT", "Lxy")};
    gDrift->SetTitle("Drift;#Delta T [#mu s];L_{xy} [cm]");
    auto gDriftLinear {defDrift.Graph("DeltaT2", "Lxy2")};
    gDriftLinear->SetTitle("Drift;#Delta T^{2} [#mu s^{2}];L_{xy}^{2} [cm^{2}]");

    // Determine crossing point
    ROOT::TThreadedObject<TH2D> h2d {"h2d", "Trajectories;X [pad];Y [pad]", 600, -50, 100, 600, -50, 100};
    def.ForeachSlot(
        [&](unsigned int slot, ROOT::Math::XYZPointF& last, ROOT::Math::XYZPointF& other)
        {
            auto u {(last - other).Unit()};
            auto dist {(last - other).R()};
            double step {0.25};
            for(double d = 0; d < dist; d += step)
            {
                auto p {other + u * d};
                h2d->Fill(p.X(), p.Y());
            }
        },
        {"LastPoint", "OtherPoint"});

    RetValues ret {};
    ret.fGraphLinear = (TGraph*)gDriftLinear->Clone();
    // Add functions to graph
    for(auto& fit : fits)
    {
        ret.fGraphLinear->GetListOfFunctions()->Add(fit);
    }
    // Compute mean
    std::vector<double> weights(uslopes.size());
    std::transform(uslopes.begin(), uslopes.end(), weights.begin(), [](double& u) { return 1. / (u * u); });
    auto mean {TMath::Mean(slopes.begin(), slopes.end(), weights.begin())};
    auto umean {TMath::Sqrt(1. / std::reduce(weights.begin(), weights.end()))};
    // And this is vd^2. Undo linearisation
    auto vd {TMath::Sqrt(TMath::Abs(mean))};
    auto uvd {umean / TMath::Sqrt(4 * TMath::Abs(mean))};
    ret.fDrift = {vd, uvd};


    if(!draw)
        return ret;

    auto* c0 {new TCanvas {"c0", "Drift canvas"}};
    c0->DivideSquare(4);
    c0->cd(1);
    h2d.Merge()->DrawClone("colz");
    c0->cd(2);
    gDrift->DrawClone("ap");
    c0->cd(3);
    gDriftLinear->DrawClone("ap");

    return ret;
}

void vDriftfromAlpha()
{
    ROOT::EnableImplicitMT();

    std::vector<int> runs {19, 21, 22, 23, 24};
    std::vector<double> fields {(6000 - 490), (5500 - 490), (4500 - 490), (4000 - 490), (3000 - 490)};

    std::vector<RetValues> rets {};
    auto* g {new TGraphErrors};
    g->SetTitle("Drift velocity;E [V/cm];v_{drift} [cm/#mus]");
    auto idx {0};
    for(const auto& run : runs)
    {
        rets.push_back(run_for(run, false));
        rets.back().fGraphLinear->SetTitle(TString::Format("Run %d", run));
        // reduced field
        auto redE {fields[idx] / 23.5}; // vertical length of ACTAR
        g->AddPoint(redE, rets.back().fDrift.first);
        g->SetPointError(idx, 0, rets.back().fDrift.second);

        idx++;
    }

    // Fit
    g->Fit("pol2", "0Q");
    auto fit {g->GetFunction("pol2")};
    fit->ResetBit(TF1::kNotDraw);
    // Eval at the actual field settings
    auto set0 {(5130 - 520) / 23.5};
    auto set1 {(5130 - 530) / 23.5};
    auto v0 {fit->Eval(set0)};
    auto v1 {fit->Eval(set1)};
    auto* text {new TPaveText {0.6, 0.2, 0.8, 0.4, "NDC"}};
    text->AddText(TString::Format("v_{d,0} = %.3f mm/#mus", v0 * 10));
    text->AddText(TString::Format("v_{d,1} = %.3f mm/#mus", v1 * 10));
    g->GetListOfFunctions()->Add(text);

    // Draw
    auto* c0 {new TCanvas {"c0", "Graphs canvas"}};
    c0->DivideSquare(runs.size());
    for(int i = 0; i < runs.size(); i++)
    {
        c0->cd(i + 1);
        rets[i].fGraphLinear->GetXaxis()->SetLimits(0, 500);
        rets[i].fGraphLinear->SetMaximum(250);
        rets[i].fGraphLinear->Draw("ap");
    }

    auto* c1 {new TCanvas {"c1", "Drift canvas"}};
    g->SetMarkerStyle(24);
    g->Draw("ap");
}
