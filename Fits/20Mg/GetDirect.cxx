#include "ActKinematics.h"
#include "ActSRIM.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TMathBase.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TROOT.h"
#include "TString.h"

#include <functional>
#include <memory>
#include <vector>

void GetDirect()
{
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_20Mg_p_p_front.root"};


    // Read efficiency
    auto file {std::make_unique<TFile>("../../Simulation/Outputs/simu_20Mg_p_p.root")};
    auto heff {file->Get<TH2D>("hEffDir2D")};
    heff->SetTitle("Simu eff");
    heff->SetDirectory(nullptr);
    file->Close();


    auto h2d {df.Histo2D({"h20Mg", "20Mg;#theta_{Lab,3}^{Dir} [#circ];T_{1}^{Dir} [MeV]", 36, 0, 180, 100, 0, 5}, "Rec_ThetaLabDir",
                         "Rec_EBeamDir")};
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

    h2d->DrawClone("colz");
}
