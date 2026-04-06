#include "ActMergerData.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"

#include "../../PostAnalysis/Gates.cxx"
#include "../../PostAnalysis/HistConfig.h"

void PlotResults()
{
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_20Mg_p_p_sil.root"};

    // Gate on RANSAC events for front layer
    auto gated {df.Filter([](ActRoot::MergerData& mer) { return S2008::isFront(mer); }, {"MergerData"})
                    .Filter("IsRANSAC == true")};

    // Book histograms
    auto hSP {gated.Histo2D(HistConfig::SP, "fSP.fCoordinates.fY", "fSP.fCoordinates.fZ")};
    auto hThetaPhi {gated.Histo2D({"hThetaPhi", "Angles;#theta [#circ];#phi [#circ]", 150, 0, 90, 180, -180, 180},
                                  "fThetaLight", "fPhiLight")};
    auto hPhi {gated.Histo1D({"hPhi", "Phi;#phi [#circ]", 150, -180, 180}, "fPhiLight")};
    auto hThetaPhiN {
        gated.Profile2D({"hThetaPhiN", "NVoxels vs #theta;#theta [#circ];#phi [#circ]", 60, 0, 60, 180, -180, 180},
                        "fThetaLight", "fPhiLight", "NVoxelsLight")};
    auto hEPhiN {gated.Profile2D({"hEPhiN", "NVoxels vs E;E_{vertex} [MeV];#phi [#circ]", 120, 0, 60, 180, -180, 180},
                                 "EVertex", "fPhiLight", "NVoxelsLight")};
    auto hE {gated.Histo1D({"hE", "E for RANSAC;E [MeV]", 120, 0, 60}, "EVertex")};

    // Plot
    auto* c0 {new TCanvas {"c0", "RANSAC results"}};
    c0->DivideSquare(6);
    c0->cd(1);
    hThetaPhi->DrawClone("colz");
    c0->cd(2);
    hThetaPhiN->DrawClone("colz");
    c0->cd(3);
    hPhi->DrawClone();
    c0->cd(4);
    hEPhiN->DrawClone("colz");
    c0->cd(5);
    hE->DrawClone();
}
