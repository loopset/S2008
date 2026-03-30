
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TExec.h"

// #include "../PostAnalysis/Gates.cxx"
#include "../PostAnalysis/HistConfig.h"

void kinRecRANSAC()
{
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", "../PostAnalysis/Outputs/tree_ex_20Mg_p_p_front.root"};

    auto hno {df.Filter("IsRANSAC == false").Histo2D(HistConfig::Kin, "fThetaLight", "EVertex")};
    auto hyes {df.Filter("IsRANSAC == true").Histo2D(HistConfig::Kin, "fThetaLight", "EVertex")};

    auto* c0 {new TCanvas {"c0", "Rec ransac canvas"}};
    c0->DivideSquare(4);
    c0->cd(1);
    hno->DrawClone("colz");
    auto* exec {new TExec{"ex1", "gStyle->SetPalette(kBird);"}};
    exec->Draw();
    hyes->DrawClone("colz same");
    // c0->cd(2);
    // hyes->SetTitle("Kin recovered with RANSAC");
}
