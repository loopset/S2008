#include "ActKinematics.h"
#include "ActSRIM.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"

#include "Math/SpecFuncMathMore.h"

#include <functional>
#include <memory>

#include "../Classes/DoubleXS.cxx"
#include "../Classes/DoubleXS.h"

/*
A very basic class for MultiDecomposition Analysis (MDA)
*/
class MDA
{
private:
    std::map<int, double> fAmps;
    double fRutherford {0};

public:
    MDA(const std::vector<int>& ls, bool withRutherford = true)
    {
        for(const auto& l : ls)
            fAmps[l] = 1;
        if(withRutherford)
            fRutherford = 1;
    }

    int GetNPar() { return fAmps.size() + (fRutherford ? 1 : 0); }

    double operator()(double* x, double* p)
    {
        // Assume x is in CM degrees
        double ret {};
        int idx {0};
        for(auto& [l, amp] : fAmps)
        {
            amp = p[idx];
            ret += amp * ROOT::Math::legendre(l, TMath::Cos(x[0] * TMath::DegToRad()));
            idx++;
        }
        if(fRutherford)
        {
            fRutherford = p[idx];
            ret += fRutherford / TMath::Power(TMath::Sin(x[0] * TMath::DegToRad() / 2), 4);
        }
        return ret;
    }

    TF1* InitTF1()
    {
        // Start at theta = 1 to avoid zero
        auto f {new TF1 {"mda", this, 1, 180, GetNPar()}};
        // Parameter name
        int idx {0};
        for(const auto& [l, amp] : fAmps)
        {
            f->SetParName(idx, TString::Format("l%d", l));
            idx++;
        }
        if(fRutherford)
            f->SetParName(idx, "Ruth");
        // Parameter initialisation
        for(int i = 0; i < GetNPar(); ++i)
        {
            f->SetParameter(i, 1);
            f->SetParLimits(i, 0, 1e4);
        }

        return f;
    }

    void Print() const
    {
        std::cout << "····· MDA ·····" << '\n';
        for(auto& [l, amp] : fAmps)
            std::cout << "  l : " << l << ", amp : " << amp << '\n';
        if(fRutherford)
            std::cout << "  Rutherford : " << fRutherford << '\n';
        std::cout << "······································" << '\n';
    }
};

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
        df.Histo2D({"h20Mg", "20Mg;#theta_{CM} [#circ];E_{CM} [MeV]", 45, 0, 180, 50, 0, 5}, "Rec_ThetaCM", "Rec_ECM")};
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
                              g->AddPointError(angle, integral, 0, error);
                              idx++;
                          }
                          return g;
                      }};
    std::pair<double, double> id {0.8, 1.2};
    std::pair<double, double> is {1.35, 1.8};

    // 5/2+ integral
    auto* gd {integralInE(id)};
    gd->SetTitle("5/2^{+};#theta_{CM};#sum d#sigma/d#Omega");
    auto* gs {integralInE(is)};
    gs->SetTitle("1/2^{+};#theta_{CM};#sum d#sigma/d#Omega");

    // d state
    auto* md {new MDA {{0, 2}, true}};
    auto* fmd {md->InitTF1()};
    gd->Fit(fmd, "R");
    md->Print();

    // s state
    auto* ms {new MDA {{0, 2}, true}};
    auto* fms {ms->InitTF1()};
    gs->Fit(fms, "R");
    ms->Print();



    gStyle->SetOptFit();
    auto* c0 {new TCanvas {"c0", "Ang canvas"}};
    c0->DivideSquare(6);
    c0->cd(1);
    projE->Draw("histe");
    gPad->Update();
    for(auto& v : {id.first, id.second, is.first, is.second})
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
    fmd->Draw();
    c0->cd(5);
    fms->Draw();

    // Draw ECM projs with lines
    xs.DrawProjectionsECM(
        [&id, &is](TH1* p)
        {
            p->SetLineColor(9);

            for(auto& v : {id.first, id.second, is.first, is.second})
            {
                auto* l {new TLine {v, gPad->GetUymin(), v, gPad->GetUymax()}};
                l->SetLineColor(kRed);
                l->Draw();
            }
        });
}
