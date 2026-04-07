#include "ActMergerData.h"

#include "ROOT/RDF/HistoModels.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"

#include "TCanvas.h"
#include "THStack.h"
#include "TString.h"
#include "TVirtualPad.h"

#include <vector>

#include "../../PostAnalysis/Gates.cxx"

void PlotItMinVoxels()
{
    ROOT::EnableImplicitMT();

    // Histo model
    ROOT::RDF::TH1DModel mECM {"hRecECM", "Rec ECM;E_{CM} [MeV];Counts / 100 keV", 50, 0, 5};
    ROOT::RDF::TH1DModel mEVertex {"hEVertex", "E_{vertex};E_{vertex} [MeV];Counts", 50, 0, 30};

    // Vector with values
    std::vector<int> minvoxels {5, 10, 11, 15, 20};

    std::vector<TH1D*> hsRecECM, hsEVertex;
    auto* stack {new THStack};
    stack->SetTitle(mECM.fTitle);

    auto* stackE {new THStack};
    stackE->SetTitle(mEVertex.fTitle);

    auto* effECM {new THStack};
    effECM->SetTitle("Eff E_{CM};E_{CM} [MeV];Counts[thresh] / Counts[thresh = 5]");

    auto* effEVertex {new THStack};
    effEVertex->SetTitle("Eff E_{vertex};E_{vertex} [MeV];Counts[thresh] / Counts[thresh = 5]");

    for(const auto& min : minvoxels)
    {
        // Get filename
        auto filename {TString::Format("./Inputs/tree_sil_minvoxels_%02d.root", min)};
        ROOT::RDataFrame df {"Final_Tree", filename};

        // Gate on layers
        auto front {df.Filter([](ActRoot::MergerData& mer) { return S2008::isFront(mer); }, {"MergerData"})};
        // And in RANSAC events
        auto ransac {front.Filter("IsRANSAC == true")};

        // Book histograms
        auto hRecECM {ransac.Histo1D(mECM, "Rec_ECM")};
        auto hEVertex {ransac.Histo1D(mEVertex, "EVertex")};

        // Add to vectors
        hsRecECM.push_back((TH1D*)hRecECM->Clone());
        hsRecECM.back()->SetTitle(TString::Format("%02d", min));
        stack->Add(hsRecECM.back(), "hist");
        // ECM clone
        auto* clone {(TH1D*)hsRecECM.back()->Clone()};
        clone->Divide(hsRecECM[0]);
        effECM->Add(clone, "hist");

        ///////////////////////

        hsEVertex.push_back((TH1D*)hEVertex->Clone());
        hsEVertex.back()->SetTitle(TString::Format("%02d", min));
        stackE->Add(hsEVertex.back(), "hist");
        clone = (TH1D*)hsEVertex.back()->Clone();
        clone->Divide(hsEVertex[0]);
        effEVertex->Add(clone, "hist");
    }

    // Draw
    auto* c0 {new TCanvas {"c0", "RANSAC canvas"}};
    c0->DivideSquare(4);
    c0->cd(1);
    stack->Draw("nostack plc pmc");
    gPad->BuildLegend();
    c0->cd(2);
    stackE->Draw("nostack plc pmc");
    gPad->BuildLegend();
    c0->cd(3);
    effECM->Draw("nostack plc pmc");
    gPad->BuildLegend();
    c0->cd(4);
    effEVertex->Draw("nostack plc pmc");
    gPad->BuildLegend();
}
