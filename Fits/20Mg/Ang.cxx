#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH2.h"
#include "TLine.h"
#include "TMath.h"
#include "TVirtualPad.h"

#include <memory>
#include <utility>
#include <vector>
void Ang()
{
    auto file {std::make_unique<TFile>("./Outputs/preliminary_xs.root")};
    auto* h {file->Get<TH2D>("fHist")};
    h->SetDirectory(nullptr);
    file->Close();

    // Get full projection
    auto* proj {h->ProjectionY()};

    // And integrate in regions
    int cThresh {40};
    int idx = 0;
    std::vector<TH1D*> projsE {};
    std::vector<std::pair<double, double>> ivs {};
    for(int x = 1; x <= h->GetNbinsX(); x++)
    {
        auto proj {h->ProjectionY(TString::Format("pE%d", idx), x, x)};
        if(proj->GetEntries() < cThresh)
        {
            delete proj;
            continue;
        }
        else
        {
            auto low {h->GetXaxis()->GetBinLowEdge(x)};
            auto up {h->GetXaxis()->GetBinUpEdge(x)};
            proj->SetTitle(
                TString::Format("#theta_{Lab} #in [%.2f,%.2f)#circ;E_{Lab} [MeV];d#sigma/d#Omega [mb/sr]", low, up));
            projsE.push_back(proj);
            ivs.push_back({low, up});
            idx++;
        }
    }

    // Function to return graph with integral in E range for each theta bin
    auto integralInE {[&ivs, &projsE](const std::pair<double, double>& bounds)
                      {
                          auto* g {new TGraphErrors};
                          g->SetTitle(";#theta_{Lab} [#circ];#sum xs");
                          int idx = 0;
                          for(const auto& p : projsE)
                          {
                              auto low {p->FindBin(bounds.first)};
                              auto up {p->FindBin(bounds.second)};
                              double error {0};
                              auto integral {p->IntegralAndError(low, up, error)};
                              g->AddPointError((ivs[idx].first + ivs[idx].second) / 2, integral);
                              g->SetPointError(g->GetN() - 1, 0, error);
                              idx++;
                          }
                          return g;
                      }};
    std::pair<double, double> id {0.8, 1.3};
    std::pair<double, double> is {1.3, 1.9};

    // 5/2+ integral
    auto* gd {integralInE(id)};
    gd->SetTitle("5/2^{+}");
    auto* gs {integralInE(is)};
    gs->SetTitle("1/2^{+}");

    // Build legendre
    auto* ld {new TF1 {"ld", "ROOT::Math::legendre(2, x)", -1, 1}};
    auto* ls {new TF1 {"ls", "ROOT::Math::legendre(0, x)", -1, 1}};
    ls->SetLineColor(8);

    auto* c0 {new TCanvas {"c0", "Ang canvas"}};
    c0->DivideSquare(6);
    c0->cd(1);
    h->Draw("colz");
    c0->cd(2);
    proj->Draw("histe");
    gPad->Update();
    for(auto& v : {id.first, id.second, is.second})
    {
        auto* l {new TLine {v, gPad->GetUymin(), v, gPad->GetUymax()}};
        l->SetLineColor(kRed);
        l->Draw();
    }
    c0->cd(3);
    gd->SetMarkerStyle(24);
    gd->Draw("ap");
    c0->cd(4);
    gs->SetMarkerStyle(25);
    gs->Draw("ap");
    c0->cd(5);
    ld->Draw();
    ls->Draw("same");
    gPad->BuildLegend();

    auto* c1 {new TCanvas {"c1", "Proj canvas"}};
    c1->DivideSquare(projsE.size());
    for(int i = 0; i < projsE.size(); i++)
    {
        c1->cd(i + 1);
        projsE[i]->Draw("histe");
        gPad->Update();
        for(auto& v : {id.first, id.second, is.second})
        {
            auto* l {new TLine {v, gPad->GetUymin(), v, gPad->GetUymax()}};
            l->SetLineColor(kRed);
            l->Draw();
        }
    }
}
