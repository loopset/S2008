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

void GetDirect()
{
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_20Mg_p_p_sil.root"};


    // Read efficiency
    auto file {std::make_unique<TFile>("../../Simulation/Outputs/simu_20Mg_p_p.root")};
    auto heff {file->Get<TH2D>("hEffDir2D")};
    heff->SetTitle("Simu eff");
    heff->SetDirectory(nullptr);
    file->Close();


    auto h2d {df.Histo2D({"h20Mg", "20Mg;#theta_{Lab,3}^{Dir} [#circ];T_{1}^{Dir} [MeV]", 36, 0, 180, 50, 0, 5},
                         "Rec_ThetaLabDir", "Rec_EBeamDir")};
    h2d->SetTitle("Counts");

    // Read srim
    auto srim {new ActPhysics::SRIM};
    srim->ReadTable("beam", "../../Calibrations/SRIM/20Mg_800mbar_95-5.txt");

    // Number of beams
    double Nbeams {201311 * 300}; // counter with GATCONF * div factor

    // Density of target
    double rho {4.743e19};

    // Inverse kinematics
    ActPhysics::Kinematics inverse {"20Mg(p,p)@84.85"};
    // Direct kinematics
    auto direct {inverse.GetOtherKinematics()};


    DoubleXS xs {h2d.GetPtr(), heff, srim, Nbeams, rho, direct.get()};
    xs.Draw();

    // Write
    xs.GetHist()->SaveAs("./Outputs/preliminary_xs.root");

    xs.Project(40);
    xs.DrawProjectionsThetaCM(
        [](TH1* p)
        {
            p->SetLineColor(9);
            p->GetXaxis()->SetRangeUser(100, 180);
        });
    xs.DrawProjectionsECM([](TH1* p) { p->SetLineColor(46); });

    // Write one
    xs.WriteInAzureFormat(10, "./Azure/Inputs/lab_1425.dat");
    xs.WriteInAzureFormat(11, "./Azure/Inputs/lab_1475.dat");
}
