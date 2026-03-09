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
    auto hSP {(TH2D*)file->Get<TH2D>("hSP")};
    auto hRes {(TH2D*)file->Get<TH2D>("hECMRes")};

    auto* c0 {new TCanvas {"c0", "Simulation canvas"}};
    c0->DivideSquare(4);
    c0->cd(1);
    hKin->DrawClone("colz");
    c0->cd(2);
    hEff->DrawClone("colz");
    c0->cd(3);
    hSP->DrawClone("colz");
    c0->cd(4);
    hRes->DrawClone("colz");
}
