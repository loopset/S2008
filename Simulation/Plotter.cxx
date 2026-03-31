#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"

#include <string>

#include "/media/Data/S2008/PostAnalysis/HistConfig.h"
void Plotter(const std::string& beam, const std::string& target, const std::string& light, double Ex, double T1,
             int neutronPS, int protonPS)
{

    auto filename {TString::Format("/media/Data/S2008/Simulation/Outputs/simu_%s_%s_%s.root", beam.c_str(),
                                   target.c_str(), light.c_str())};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"SimulationTTree", filename};
    auto hKin {df.Histo2D(HistConfig::Kin, "theta3Lab", "EVertex")};
    // Open file
    auto* file {new TFile {filename}};
    auto hEff {(TH2D*)file->Get<TH2D>("hEff2D")};
    auto hSPf0 {(TH2D*)file->Get<TH2D>("hSPf0")};
    auto hSPl0 {(TH2D*)file->Get<TH2D>("hSPl0")};
    auto hSPr0 {(TH2D*)file->Get<TH2D>("hSPr0")};
    auto hRes {(TH2D*)file->Get<TH2D>("hECMRes")};
    // In Direct frame
    auto hEffDir {(TH2D*)file->Get<TH2D>("hEffDir2D")};
    auto hResDir {(TH2D*)file->Get<TH2D>("hT1DirRes")};

    auto* c0 {new TCanvas {"c0", "CM Simulation canvas"}};
    c0->DivideSquare(6);
    c0->cd(1);
    hKin->DrawClone("colz");
    c0->cd(2);
    hEff->DrawClone("colz");
    c0->cd(3);
    hRes->DrawClone("colz");
    c0->cd(4);
    hSPf0->DrawClone("colz");
    c0->cd(5);
    hSPl0->DrawClone("colz");
    c0->cd(6);
    hSPr0->DrawClone("colz");

    auto* c1 {new TCanvas {"c1", "CM Simulation canvas"}};
    c1->DivideSquare(4);
    c1->cd(1);
    hEffDir->DrawClone("colz");
    c1->cd(2);
    hResDir->DrawClone("colz");
}
