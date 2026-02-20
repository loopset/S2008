#include "ActMergerData.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/TThreadedObject.hxx>

#include "TCanvas.h"
#include "TColor.h"
#include "TExec.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TVirtualPad.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

#include <vector>

void plotECM_intervalsThetaCM()
{
    // Read data
    auto filename {TString::Format("../PostAnalysis/Outputs/tree_ex_20Mg_p_p.root")};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", filename};
    auto nodeSil {df.Filter("IsL1 == false")};
    auto nodeL1 {df.Filter("IsL1 == true")};

    // Model of histogram
    ROOT::RDF::TH1DModel mECM {"hECM", "E_{CM};E_{CM} [MeV];Counts / 10 keV", 500, 0, 5};
    // Get histograms of Ecm on intervals of RP.x()
    std::vector<ROOT::TThreadedObject<TH1D>*> hsSil, hsL1;
    double step {10}; // deg
    double Thetamin {40};
    double Thetamax {180};
    int idx {};
    for(double theta = Thetamin; theta < Thetamax; theta += step)
    {
        hsSil.push_back(new ROOT::TThreadedObject<TH1D>(
            TString::Format("hECM%d", idx),
            TString::Format("#theta_{CM} [%.2f, %.2f);E_{CM} [MeV];Counts / 10 keV", theta, theta + step),
            mECM.fNbinsX, mECM.fXLow, mECM.fXUp));
        hsL1.push_back(new ROOT::TThreadedObject<TH1D>(
            TString::Format("hECM%d", idx),
            TString::Format("#theta_{CM} [%.2f, %.2f);E_{CM} [MeV];Counts / 10 keV", theta, theta + step),
            mECM.fNbinsX, mECM.fXLow, mECM.fXUp));
        idx++;
    }
    // Initialize slot 0 to not crash
    for(auto& h : hsSil)
        h->GetAtSlot(0);
    for(auto& h : hsL1)
        h->GetAtSlot(0);

    // Fill histograms
    // auto hECM2d {df.Histo2D(HistConfig::ThetaCMECM, "Rec_ThetaCM", "Rec_ECM")};
    // Silicons
    nodeSil.ForeachSlot(
        [&](unsigned int slot, double thetaCM, double ecm)
        {
            // get the hstogram we have to fill
            for(size_t i = 0; i < hsSil.size(); i++)
            {
                double theta = Thetamin + i * step;
                if(thetaCM >= theta && thetaCM < theta + step)
                {
                    hsSil[i]->GetAtSlot(slot)->Fill(ecm);
                }
            }
        },
        {"Rec_ThetaCM", "Rec_ECM"});
    // L1
    nodeL1.ForeachSlot(
        [&](unsigned int slot, double thetaCM, double ecm)
        {
            // get the hstogram we have to fill
            for(size_t i = 0; i < hsSil.size(); i++)
            {
                double theta = Thetamin + i * step;
                if(thetaCM >= theta && thetaCM < theta + step)
                {
                    hsL1[i]->GetAtSlot(slot)->Fill(ecm);
                }
            }
        },
        {"Rec_ThetaCM", "Rec_ECM"});


    // Styling options
    gStyle->SetPalette(kRainBow);

    // Silicon stack
    auto* stack {new THStack};
    stack->SetTitle("Silicon E_{CM};E_{CM} [MeV];Counts / 20 keV");
    // L1 stack
    auto* lstack {new THStack};
    lstack->SetTitle("L1 E_{CM};E_{CM} [MeV];Counts / 10 keV");


    // Plot them in canvas
    auto* c0 {new TCanvas("c0", "Silicon ECM")};
    c0->DivideSquare(hsSil.size());
    int p {1};
    for(auto& h : hsSil)
    {
        c0->cd(p);
        auto merged {h->Merge()};
        // merged->GetXaxis()->SetRangeUser(0.6, 4.1);
        merged->DrawClone();
        auto* clone {(TH1D*)merged->Clone()};
        clone->SetLineWidth(2);
        clone->Rebin(2);
        stack->Add(clone);
        p++;
    }

    // L1 canvas
    auto* c1 {new TCanvas("c1", "L1 ECM")};
    c1->DivideSquare(hsL1.size());
    p = 1;
    for(auto& h : hsL1)
    {
        c1->cd(p);
        auto merged {h->Merge()};
        merged->DrawClone();
        auto* clone {(TH1D*)merged->Clone()};
        clone->SetLineWidth(2);
        // clone->Rebin(2);
        lstack->Add(clone);
        p++;
    }

    auto* c2 {new TCanvas {"c2", "Ecm 2d canvas"}};
    c2->Divide(1, 2);
    c2->cd(1);
    stack->Draw("plc pmc");
    auto* leg = gPad->BuildLegend(0.8, 0.1, 0.95, 0.85);
    c2->cd(2);
    lstack->Draw("plc pmc");
    leg = gPad->BuildLegend(0.8, 0.1, 0.95, 0.85);

    // Create lines and texts
    // These are for 21Mg!
    // std::vector<double> resonances {0.717, 0.85, 0.98, 1.28, 1.885};
    // No data for 21Al
    std::vector<double> resonances {};
    std::vector<std::string> comments {"?", "New?", "", "", "New!"};
    auto pad {c2->cd(1)};
    pad->Update();
    for(const auto& res : resonances)
    {
        auto* line {new TLine {res, pad->GetUymin(), res, pad->GetUymax()}};
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw();
    }
    pad = c2->cd(2);
    pad->Update();
    for(const auto& res : resonances)
    {
        auto* line {new TLine {res, pad->GetUymin(), res, pad->GetUymax()}};
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw();
    }

    // And now text
    pad = c2->cd(1);
    idx = 0;
    for(const auto& res : resonances)
    {
        auto* text {new TPaveText {res + 0.05, 2500, res + 0.25, 2900, "NB"}};
        text->SetFillStyle(0);
        text->SetBorderSize(0);
        text->AddText(TString::Format("%.2f MeV", res));
        // if(comments[idx].size())
        //     text->AddText(comments[idx].c_str());
        // Style
        text->SetTextFont(22);
        text->SetTextSize(0.04);
        text->Draw();
        idx++;
    }
}
