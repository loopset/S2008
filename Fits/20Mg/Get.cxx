#include "ActKinematics.h"
#include "ActSRIM.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TF1.h"
#include "TFile.h"
#include "TROOT.h"

#include <functional>
#include <memory>

#include "../Classes/DoubleXS.cxx"
#include "../Classes/DoubleXS.h"


void Get()
{
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_20Mg_p_p_front.root"};


    // Read efficiency
    auto file {std::make_unique<TFile>("../../Simulation/Outputs/simu_20Mg_p_p.root")};
    auto heff {file->Get<TH2D>("hEff2D")};
    heff->SetTitle("Simu eff");
    heff->SetDirectory(nullptr);
    file->Close();


    auto h2d {
        df.Histo2D({"h20Mg", "20Mg;#theta_{CM} [#circ];E_{CM} [MeV]", 36, 0, 180, 50, 0, 5}, "Rec_ThetaCM", "Rec_ECM")};
    h2d->SetTitle("Counts");

    // Read srim
    auto srim {new ActPhysics::SRIM};
    srim->ReadTable("beam", "../../Calibrations/SRIM/20Mg_800mbar_95-5.txt");

    // Kinematics
    auto kin {new ActPhysics::Kinematics {"20Mg(p,p)@84"}};

    // Number of beams
    double Nbeams {201311 * 300}; // counter with GATCONF * div factor

    // Density of target
    double rho {4.743e19};

    DoubleXS xs {h2d.GetPtr(), heff, srim, Nbeams, rho, kin};
    xs.Draw();
    xs.Project(40);
    xs.DrawProjectionsThetaCM(
        [](TH1* p)
        {
            p->SetLineColor(9);
            p->GetXaxis()->SetRangeUser(100, 180);
        });
    xs.DrawProjectionsECM([](TH1* p) { p->SetLineColor(46); });

    // // Get projection
    // auto thetaCMmin {140};
    // auto thetaCMmax {145};
    // auto* p {xs.GetProjectionECM(thetaCMmin, thetaCMmax)};
    //
    // auto* c0 {new TCanvas {"c0", "Get canvas"}};
    // c0->DivideSquare(4);
    // c0->cd(1);
    // p->Draw("histe");

    auto outfile {std::make_unique<TFile>("./Outputs/20Mg_preliminary.root", "recreate")};
    for(auto* c : *(gROOT->GetListOfCanvases()))
        c->Write();
    outfile->Close();
}
