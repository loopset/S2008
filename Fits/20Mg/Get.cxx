#include "ActKinematics.h"
#include "ActSRIM.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TROOT.h"

#include <functional>
#include <memory>

#include "../Classes/DoubleXS.cxx"
#include "../Classes/DoubleXS.h"


void Get()
{
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_20Mg_p_p_sil.root"};


    // Read efficiency
    auto file {std::make_unique<TFile>("../../Simulation/Outputs/simu_20Mg_p_p.root")};
    auto heff {file->Get<TH2D>("hEff2D")};
    heff->SetTitle("Simu eff");
    heff->SetDirectory(nullptr);
    file->Close();


    auto h2d {
        df.Histo2D({"h20Mg", "20Mg;#theta_{CM} [#circ];E_{CM} [MeV]", 180, 0, 180, 100, 0, 5}, "Rec_ThetaCM", "Rec_ECM")};
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

    DoubleXS xs {h2d.GetPtr(), heff, srim, Nbeams, rho, kin, "CM"};
    xs.Draw();
    xs.Project(40);
    xs.DrawProjectionsThetaCM(
        [](TH1* p)
        {
            p->SetLineColor(8);
            p->GetXaxis()->SetRangeUser(100, 180);
        });
    // Get full projection
    auto* projE {xs.GetHist()->ProjectionY("projE")};

    // Save them
    xs.GetHist()->SaveAs("./Outputs/preliminary_xs_cm.root");

    // Function to return graph with integral in E range for each theta bin
    auto integralInE {[&xs](const std::pair<double, double>& bounds)
                      {
                          auto* g {new TGraphErrors};
                          int idx = 0;
                          for(const auto& p : xs.GetProjsE())
                          {
                              auto low {p->FindBin(bounds.first)};
                              auto up {p->FindBin(bounds.second) - 1};
                              double error {0};
                              auto integral {p->IntegralAndError(low, up, error)};
                              auto angle {(xs.GetIvsE()[idx].first + xs.GetIvsE()[idx].second) / 2};
                              angle = TMath::Cos(angle * TMath::DegToRad());
                              g->AddPoint(angle, integral);
                              g->SetPointError(g->GetN() - 1, 0, error);
                              idx++;
                          }
                          return g;
                      }};
    std::pair<double, double> id {0.8, 1.2};
    std::pair<double, double> is {1.2, 1.8};

    // 5/2+ integral
    auto* gd {integralInE(id)};
    gd->SetTitle("5/2^{+};cos #theta_{CM};#sum d#sigma/d#Omega");
    auto* gs {integralInE(is)};
    gs->SetTitle("1/2^{+};cos #theta_{CM};#sum d#sigma/d#Omega");

    // Build legendre
    auto* ld {new TF1 {"ld", "ROOT::Math::legendre(2, x)", -1, 1}};
    ld->SetTitle("l = 2");
    auto* ls {new TF1 {"ls", "ROOT::Math::legendre(0, x)", -1, 1}};
    ls->SetTitle("l = 0");
    ls->SetLineColor(8);


    auto* c0 {new TCanvas {"c0", "Ang canvas"}};
    c0->DivideSquare(6);
    c0->cd(1);
    projE->Draw("histe");
    gPad->Update();
    for(auto& v : {id.first, id.second, is.second})
    {
        auto* l {new TLine {v, gPad->GetUymin(), v, gPad->GetUymax()}};
        l->SetLineColor(kRed);
        l->Draw();
    }
    c0->cd(2);
    gd->SetMarkerStyle(24);
    gd->Draw("ap");
    c0->cd(3);
    gs->SetMarkerStyle(25);
    gs->Draw("ap");
    c0->cd(4);
    ld->Draw();
    ld->GetXaxis()->SetTitle("cos #theta_{CM}");
    ld->GetYaxis()->SetTitle("P(x)");
    ls->Draw("same");
    gPad->BuildLegend();

    // Draw ECM projs with lines
    xs.DrawProjectionsECM(
        [&id, &is](TH1* p)
        {
            p->SetLineColor(9);

            for(auto& v : {id.first, id.second, is.second})
            {
                auto* l {new TLine {v, gPad->GetUymin(), v, gPad->GetUymax()}};
                l->SetLineColor(kRed);
                l->Draw();
            }
        });
}
