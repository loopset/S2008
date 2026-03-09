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

class DoubleXS
{
private:
    std::vector<TH1D*> fProjsThetaCM {};
    std::vector<TH1D*> fProjsECM {};
    TH2* fOriginal {};
    TH2* fHist {};
    TH2* fEff {};
    ActPhysics::SRIM* fsrim {};
    TProfile2D* fEffProf {};
    TF1* fFuncOmega {};
    TH2* fThick {};
    ActPhysics::Kinematics* fKin {};
    double fNbeams {};
    double fDensity {};

public:
    DoubleXS(TH2* hData, TH2* hEff, ActPhysics::SRIM* srim, ActPhysics::Kinematics* kin, double nb, double rho);

    void Draw();

    void Project(int cThresch = 50);
    void DrawProjectionsThetaCM(const std::function<void(TH1* h)>& apply = nullptr);
    void DrawProjectionsECM(const std::function<void(TH1* h)>& apply = nullptr);
    TH1D* GetProjectionECM(double thetamin, double thetamax);

private:
    void ApplyEff();
    void ApplySolidAngle();
    void ApplyThickness();
    void ApplyNormalisation();
    double TransformCMtoLab(double ecm);
};

DoubleXS::DoubleXS(TH2* hData, TH2* hEff, ActPhysics::SRIM* srim, ActPhysics::Kinematics* kin, double nb, double rho)
    : fHist((TH2*)hData->Clone("fHist")),
      fEff(hEff),
      fsrim(srim),
      fKin(kin),
      fNbeams(nb),
      fDensity(rho)
{
    fOriginal = (TH2*)hData->Clone("hOriginal");
    fHist->SetTitle("Cross section");
    fHist->GetZaxis()->SetTitle("d#sigma/d#Omega [mb/sr]");
    fHist->GetZaxis()->SetMaxDigits(2);
    fHist->Sumw2();
    ApplyEff();
    ApplySolidAngle();
    ApplyThickness();
    ApplyNormalisation();
}

void DoubleXS::ApplyEff()
{
    // Build TProfile
    fEffProf = new TProfile2D {"hEffProf",
                               "Mean efficiency;#theta_{CM} [#circ];E_{CM} [MeV];#epsilon",
                               fHist->GetNbinsX(),
                               fHist->GetXaxis()->GetXmin(),
                               fHist->GetXaxis()->GetXmax(),
                               fHist->GetNbinsY(),
                               fHist->GetYaxis()->GetXmin(),
                               fHist->GetYaxis()->GetXmax()};
    // Fill it
    for(int x = 1; x <= fEff->GetNbinsX(); x++)
    {
        for(int y = 1; y <= fEff->GetNbinsY(); y++)
        {
            auto content {fEff->GetBinContent(x, y)};
            if(content <= 0)
                continue;
            auto thetacm {fEff->GetXaxis()->GetBinCenter(x)};
            auto ecm {fEff->GetYaxis()->GetBinCenter(y)};
            fEffProf->Fill(thetacm, ecm, content);
        }
    }
    // Divide
    fHist->Divide(fEffProf);
}

void DoubleXS::ApplySolidAngle()
{

    auto binOmega {fHist->GetXaxis()->GetBinUpEdge(1)}; // in deg
    fFuncOmega = new TF1 {
        "fFuncOmega",
        TString::Format("TMath::TwoPi() * TMath::Sin(x * TMath::DegToRad()) * (%f) * TMath::DegToRad()", binOmega),
        fHist->GetXaxis()->GetXmin(), fHist->GetXaxis()->GetXmax()};
    fFuncOmega->SetTitle("Solid angle func;#theta_{CM} [#circ];2#pi sin(#theta_{CM}) #Delta #theta_{CM} [sr]");
    fHist->Divide(fFuncOmega);
}

double DoubleXS::TransformCMtoLab(double ecm)
{
    auto m1 {fKin->GetParticle(1).GetMass()};
    auto m2 {fKin->GetParticle(2).GetMass()};
    return (m1 + m2) / m2 * ecm;
}

void DoubleXS::ApplyThickness()
{
    fThick = new TH2D {"Thick",
                       "SRIM thickness per E_{CM} bin;#theta_{CM} [#circ];E_{CM} [MeV];Thick [cm]",
                       fHist->GetNbinsX(),
                       fHist->GetXaxis()->GetXmin(),
                       fHist->GetXaxis()->GetXmax(),
                       fHist->GetNbinsY(),
                       fHist->GetYaxis()->GetXmin(),
                       fHist->GetYaxis()->GetXmax()};
    if(!fsrim)
        return;
    for(int y = 1; y <= fThick->GetNbinsY(); y++)
    {
        auto ecmlow {fThick->GetYaxis()->GetBinLowEdge(y)};
        auto ecmup {fThick->GetYaxis()->GetBinUpEdge(y)};
        // And eval srim ranges in LAB
        auto elablow {TransformCMtoLab(ecmlow)};
        auto rlow {fsrim->EvalRange("beam", elablow)};
        auto elabup {TransformCMtoLab(ecmup)};
        auto rup {fsrim->EvalRange("beam", elabup)};
        // And DeltaX
        auto deltax {TMath::Abs(rup - rlow)};

        // Fill that value for all X bins
        for(int x = 1; x <= fThick->GetNbinsX(); x++)
        {
            fThick->SetBinContent(x, y, deltax / 10); // mm to cm (later on, density is in atoms/cm3)
            fThick->SetBinError(x, y, 0);             // this has no error (in principle...)
        }
    }
    fHist->Divide(fThick);
}

void DoubleXS::ApplyNormalisation()
{
    // Divide by density
    fHist->Scale(1. / fDensity);
    // Divide by Nbeams
    fHist->Scale(1. / fNbeams);
    // Apply convert to [mb/sr]
    fHist->Scale(1e27);
}

void DoubleXS::Draw()
{
    auto* c0 {new TCanvas {"cxsCalc", "Double XS calculations"}};
    c0->DivideSquare(6);
    c0->cd(1);
    fOriginal->Draw("colz");
    c0->cd(2);
    fEff->Draw("colz");
    c0->cd(3);
    fEffProf->Draw("colz1");
    c0->cd(4);
    fFuncOmega->Draw();
    c0->cd(5);
    fThick->Draw("colz");
    c0->cd(6);
    fHist->Draw("colz");
}

void DoubleXS::Project(int cThresh)
{
    // Project onto thetaCM axis
    int idx {};
    for(int y = 1; y <= fHist->GetNbinsY(); y++)
    {
        auto proj {fHist->ProjectionX(TString::Format("pTheta%d", idx), y, y)};
        if(proj->GetEntries() < cThresh)
        {
            delete proj;
            continue;
        }
        else
        {
            auto low {fHist->GetYaxis()->GetBinLowEdge(y)};
            auto up {fHist->GetYaxis()->GetBinUpEdge(y)};
            proj->SetTitle(
                TString::Format("E_{CM} #in [%.2f,%.2f) MeV;#theta_{CM} [#circ];d#sigma/d#Omega [mb/sr]", low, up));
            fProjsThetaCM.push_back(proj);
            idx++;
        }
    }
    // Project onto ECM axis
    idx = 0;
    for(int x = 1; x <= fHist->GetNbinsX(); x++)
    {
        auto proj {fHist->ProjectionY(TString::Format("pE%d", idx), x, x)};
        if(proj->GetEntries() < cThresh)
        {
            delete proj;
            continue;
        }
        else
        {
            auto low {fHist->GetXaxis()->GetBinLowEdge(x)};
            auto up {fHist->GetXaxis()->GetBinUpEdge(x)};
            proj->SetTitle(
                TString::Format("#theta_{CM} #in [%.2f,%.2f)#circ;E_{CM} [MeV];d#sigma/d#Omega [mb/sr]", low, up));
            fProjsECM.push_back(proj);
            idx++;
        }
    }
}

void DoubleXS::DrawProjectionsThetaCM(const std::function<void(TH1*)>& apply)
{
    auto* c {new TCanvas {"cxsTheta", "ThetaCM projection canvas"}};
    c->DivideSquare(fProjsThetaCM.size());
    for(int i = 0; i < fProjsThetaCM.size(); i++)
    {
        c->cd(i + 1);
        auto& h {fProjsThetaCM[i]};
        if(apply)
            apply(h);
        h->Draw("histe");
    }
}

void DoubleXS::DrawProjectionsECM(const std::function<void(TH1*)>& apply)
{
    auto* c {new TCanvas {"cxsE", "ECM projection canvas"}};
    c->DivideSquare(fProjsECM.size());
    for(int i = 0; i < fProjsECM.size(); i++)
    {
        c->cd(i + 1);
        auto& h {fProjsECM[i]};
        if(apply)
            apply(h);
        h->Draw("histe");
    }
}

TH1D* DoubleXS::GetProjectionECM(double thetamin, double thetamax)
{
    auto bmin {fHist->GetXaxis()->FindBin(thetamin)};
    auto min {fHist->GetXaxis()->GetBinCenter(bmin)};
    auto bmax {fHist->GetXaxis()->FindBin(thetamax) - 1};
    auto max {fHist->GetXaxis()->GetBinCenter(bmax)};
    auto* p {fHist->ProfileY("pInRange", bmin, bmax)};
    p->SetTitle(TString::Format("#theta_{CM} #in [%.2f,%.2f)", min, max));
    return p;
}


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


    auto h2d {df.Histo2D({"h20Mg", "20Mg;#theta_{CM} [#circ];E_{CM} [MeV]", 36, 0, 180, 50, 0, 5}, "Rec_ThetaCM",
                         "Rec_ECM")};
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

    DoubleXS xs {h2d.GetPtr(), heff, srim, kin, Nbeams, rho};
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
