#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TString.h"

#include <memory>
#include <utility>
#include <vector>

#include "../../PostAnalysis/HistConfig.h"

class DoubleXS
{
private:
    std::vector<std::pair<double, double>> fIvs {};
    std::vector<TH1*> fProjs {};
    std::vector<TProfile*> fEffs {};
    TH2* fHistOrig {};
    TH2* fHist {};
    TH2* fEff {};
    TProfile2D* fEffProf {};

public:
    DoubleXS(TH2* hData, TH2* hEff);

    void Draw();

private:
    void ApplyEff();
    void InitIntervals(double emin, double emax, double estep);
    TH1* GetProjection(int idx, double low, double up);
    TProfile* GetMeanEfficiency(int idx, double low, double up);
};

DoubleXS::DoubleXS(TH2* hData, TH2* hEff)
    : fHistOrig(hData),
      fHist((TH2*)hData->Clone()),
      fEff(hEff)
{
    ApplyEff();
}

void DoubleXS::ApplyEff()
{
    // Build TProfile
    fEffProf = new TProfile2D {"hEffProf",
                               "#theta_{CM} [#circ];E_{CM} [MeV]",
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
            auto thetacm {fEff->GetXaxis()->GetBinCenter(x)};
            auto ecm {fEff->GetYaxis()->GetBinCenter(y)};
            fEffProf->Fill(thetacm, ecm, content);
        }
    }
    // Divide
    fHist->Divide(fEffProf);
}

void DoubleXS::Draw()
{
    auto* c0 {new TCanvas {"cxs0", "Double XS canvas 0"}};
    c0->DivideSquare(6);
    c0->cd(1);
    fHist->Draw("colz");
    c0->cd(2);
    fEff->Draw("colz");
    c0->cd(3);
    fEffProf->Draw("colz1");
}

void DoubleXS::InitIntervals(double emin, double emax, double estep)
{
    int idx {};
    for(double e = emin; e < emax; e += estep)
    {
        // Calculate interval
        std::pair<double, double> ivs {e, e + estep};
        // Push projection
        fProjs.push_back(GetProjection(idx, ivs.first, ivs.second));
        // Push efficiency
        fEffs.push_back(GetMeanEfficiency(idx, ivs.first, ivs.second));
        idx++;
    }
}

TH1* DoubleXS::GetProjection(int idx, double low, double up)
{
    auto blow {fHistOrig->FindBin(low)};
    auto bup {fHistOrig->FindBin(up)};
    auto elow {fHistOrig->GetXaxis()->GetBinCenter(blow)};
    auto eup {fHistOrig->GetXaxis()->GetBinCenter(bup)};
    auto* proj {fHistOrig->ProjectionX(TString::Format("proj%d", idx), blow, bup)};
    proj->SetTitle(TString::Format("E #in [%.2f,%.2f)", elow, eup));
    return proj;
}

TProfile* DoubleXS::GetMeanEfficiency(int idx, double low, double up)
{
    auto blow {fHistOrig->FindBin(low)};
    auto bup {fHistOrig->FindBin(up)};
    auto elow {fHistOrig->GetXaxis()->GetBinCenter(blow)};
    auto eup {fHistOrig->GetXaxis()->GetBinCenter(bup)};
    auto* prof {fHistOrig->ProfileX(TString::Format("prof%d", idx), blow, bup)};
    prof->SetTitle(TString::Format("E #in [%.2f,%.2f)", elow, eup));
    return prof;
}

void Get()
{
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", "../../PostAnalysis/Outputs/tree_ex_20Mg_p_p_front.root"};


    // Read efficiency
    auto file {std::make_unique<TFile>("../../Simulation/Outputs/simu_20Mg_p_p.root")};
    auto heff {file->Get<TH2D>("hEff2D")};
    heff->SetDirectory(nullptr);
    file->Close();

    auto h2d {df.Histo2D({"h20Mg", "20Mg;#theta_{CM} [#circ];E_{CM} [MeV]", 180, 0, 180, 50, 0, 5}, "Rec_ThetaCM", "Rec_ECM")};
    h2d->SetTitle("Efficiency");

    DoubleXS xs {h2d.GetPtr(), heff};
    xs.Draw();


    // Draw
    auto* c0 {new TCanvas {"c0", "Get dxs"}};
    c0->DivideSquare(4);
    c0->cd(1);
    h2d->DrawClone("colz");
    c0->cd(2);
    heff->DrawClone("colz");
}
