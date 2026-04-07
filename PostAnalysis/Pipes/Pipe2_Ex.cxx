#ifndef Pipe2_Ex_cxx
#define Pipe2_Ex_cxx

#include "ActCutsManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
#include "Rtypes.h"

#include "TAttLine.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"
#include "TVirtualPad.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../HistConfig.h"

void Pipe2_Ex(const std::string& beam, const std::string& target, const std::string& light)
{
    // Read data
    auto filename {TString::Format("./Outputs/tree_pid_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"PID_Tree", filename};

    // Init SRIM
    auto* srim {new ActPhysics::SRIM};
    // Correct SRIM names
    std::string srimName {};
    if(light == "d")
        srimName = "2H";
    else if(light == "p")
        srimName = "1H";
    else if(light == "t")
        srimName = "3H";
    else if(light == "3He")
        srimName = "3He";
    else if(light == "4He")
        srimName = "4He";
    int pressure {800}; // 20Me beam
    if(beam == "20Ne" || beam == "20Na")
        pressure = 950;
    srim->ReadTable(light,
                    TString::Format("../Calibrations/SRIM/%s_%dmbar_95-5.txt", srimName.c_str(), pressure).Data());
    srim->ReadTable(beam, TString::Format("../Calibrations/SRIM/%s_%dmbar_95-5.txt", beam.c_str(), pressure).Data());


    // Read cuts
    ActRoot::CutsManager<std::string> cuts;
    cuts.ReadCut("ep_range", "./Cuts/elastic_ep_range.root");
    cuts.ReadCut("debug", "./Cuts/debug_l1.root");
    cuts.ReadCut("debug_ep_range", "./Cuts/debug_ep_range.root");

    // Build energy at vertex
    auto dfVertex = df.Define("EVertex",
                              [&](const ActRoot::MergerData& d)
                              {
                                  double ret {};
                                  if(d.fLight.IsFilled() && std::isfinite(d.fLight.fTL))
                                      ret = srim->EvalInitialEnergy(light, d.fLight.fEs.front(), d.fLight.fTL);
                                  else if(d.fLight.IsL1()) // L1 trigger
                                      ret = srim->EvalEnergy(light, d.fLight.fTL);
                                  return ret;
                              },
                              {"MergerData"});

    // Init particles
    ActPhysics::Particle pb {beam};
    auto mbeam {pb.GetMass()};
    ActPhysics::Particle pt {target};
    auto mtarget {pt.GetMass()};
    ActPhysics::Particle pl {light};
    // Beam energies
    // First set of data, runs [31-40]
    double EBeamFirst {4.24}; // MeV/u from SRIM interpolation of exp. 20Mg range
    // 2nd set of data, runs [47, onwards] -> is the SAME
    double EBeamSecond {4.24}; // MeV/u from SRIM
    std::map<int, double> EBeams;
    for(int run = 31; run <= 40; run++)
        EBeams[run] = EBeamFirst;
    for(int run = 47; run <= 140; run++)
        EBeams[run] = EBeamSecond;
    // Allegedly wrong LISE calculation below
    // double EBeamIni {4.05235}; // AMeV at X = 0 of pad plane; energy meassure 5.44 before cfa
    // Energy of unidentified beam on Sunday 12
    // double EBeamIni {4.17}; // AMeV at X = 0. WARNING: SUSPICTION THIS IS NOT 20Mg

    // All the above beam energies include energy losses in CFA, window, etc...
    // They're given at X = 0 of the pad plane
    ActPhysics::Kinematics kin {pb, pt, pl, EBeamFirst * pb.GetAMU()};
    auto T1thresh {kin.GetT1Thresh()};
    // Vector of kinematics as one object is needed per
    // processing slot (since we are changing EBeam in each entry)
    std::vector<ActPhysics::Kinematics> vkins {df.GetNSlots()};
    for(auto& k : vkins)
        k = kin;

    // Beam energy calculation and ECM
    auto def {
        dfVertex
            .Define("EBeam",
                    [&](const ActRoot::MergerData& d)
                    {
                        double EBeam {};
                        if(EBeams.count(d.fRun))
                            EBeam = EBeams[d.fRun];
                        else
                            throw std::runtime_error("Defining EBeam: no initial beam energy for run " +
                                                     std::to_string(d.fRun));
                        return srim->Slow(beam, EBeam * pb.GetAMU(), d.fRP.X());
                    },
                    {"MergerData"})
            .DefineSlot("Rec_EBeam", // assuming Ex = 0 using outgoing light particle kinematics
                        [&](unsigned int slot, double EVertex, const ActRoot::MergerData& d)
                        {
                            // no need for slots here but for the sake of consistency with next calculations...
                            return vkins[slot].ReconstructBeamEnergyFromLabKinematics(EVertex, d.fThetaLight *
                                                                                                   TMath::DegToRad());
                        },
                        {"EVertex", "MergerData"})
            .Define("ECM", [&](double EBeam) { return (mtarget / (mbeam + mtarget)) * EBeam; }, {"EBeam"})
            .Define("Rec_ECM", [&](double rec_EBeam) { return (mtarget / (mbeam + mtarget)) * rec_EBeam; },
                    {"Rec_EBeam"})
            .Filter("fRP.fCoordinates.fX <= 200") // Mask decays by position... for 20Na; for 20Mg ~ 205 mm
            .Filter(                    // filter events poorly reconstructed with NaN when evaluating momentum
                [&pb](double rec_ebeam) // these are few events but since we throw an exception, halts all calculations
                {
                    // For some events p is imaginary... poor reconstruction of who knows
                    double E {rec_ebeam + pb.GetMass()};
                    double p {TMath::Sqrt(E * E - pb.GetMass() * pb.GetMass())};
                    return std::isfinite(p);
                },
                {"Rec_EBeam"})};

    def =
        def.DefineSlot("Ex",
                       [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
                       {
                           vkins[slot].SetBeamEnergy(EBeam);
                           return vkins[slot].ReconstructExcitationEnergy(EVertex, (d.fThetaLight) * TMath::DegToRad());
                       },
                       {"MergerData", "EVertex", "EBeam"})
            .DefineSlot("Rec_Ex",
                        [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
                        {
                            vkins[slot].SetBeamEnergy(EBeam);
                            return vkins[slot].ReconstructExcitationEnergy(EVertex,
                                                                           (d.fThetaLight) * TMath::DegToRad());
                        },
                        {"MergerData", "EVertex", "Rec_EBeam"})
            .DefineSlot("ThetaCM",
                        [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
                        {
                            vkins[slot].SetBeamEnergy(EBeam);
                            return vkins[slot].ReconstructTheta3CMFromLab(EVertex,
                                                                          (d.fThetaLight) * TMath::DegToRad()) *
                                   TMath::RadToDeg();
                        },
                        {"MergerData", "EVertex", "EBeam"})
            .DefineSlot("Rec_ThetaCM",
                        [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
                        {
                            vkins[slot].SetBeamEnergy(EBeam);
                            return vkins[slot].ReconstructTheta3CMFromLab(EVertex,
                                                                          (d.fThetaLight) * TMath::DegToRad()) *
                                   TMath::RadToDeg();
                        },
                        {"MergerData", "EVertex", "Rec_EBeam"})
            .DefineSlot("Rec_ThetaLabDir",
                        [&](unsigned int slot, double rec_ebeam, double rec_thetacm)
                        {
                            vkins[slot].SetBeamEnergy(rec_ebeam);
                            return vkins[slot].ComputeOtherInLab(rec_thetacm * TMath::DegToRad()).first *
                                   TMath::RadToDeg();
                        },
                        {"Rec_EBeam", "Rec_ThetaCM"})
            .DefineSlot("Rec_EBeamDir",
                        [&](unsigned int slot, double rec_ebeam)
                        {
                            vkins[slot].SetBeamEnergy(rec_ebeam);
                            return vkins[slot].ComputeEquivalentOtherT1(rec_ebeam);
                        },
                        {"Rec_EBeam"});

    // Define range of heavy particle
    def = def.Define("RangeHeavy", [&](ActRoot::MergerData& d) { return d.fRP.X() + d.fHeavy.fTL; }, {"MergerData"});


    // Create node to gate on different conditions: silicon layer, l1, etc
    // L0 trigger
    auto nodel0 {def.Filter([](ActRoot::MergerData& d) { return d.fLight.IsL1() == false; }, {"MergerData"})};
    // L0 -> side silicons
    auto nodeLat {nodel0.Filter([](ActRoot::MergerData& d)
                                { return d.fLight.fLayers.front() == "l0" || d.fLight.fLayers.front() == "r0"; },
                                {"MergerData"})};
    // L0 -> front silicons
    auto nodeFront {
        nodel0.Filter([](ActRoot::MergerData& d) { return d.fLight.fLayers.front() == "f0"; }, {"MergerData"})};

    // L1 trigger
    auto nodel1 {def.Filter([](ActRoot::MergerData& d) { return d.fLight.IsL1() == true; }, {"MergerData"})};

    // Selection in Ep vs R20Mg plot
    auto nodeEpRSil {nodel0.Filter([&](float range, double elab) { return cuts.IsInside("ep_range", range, elab); },
                                   {"RangeHeavy", "EVertex"})};
    auto nodeEpFront {nodeFront.Filter([&](float range, double elab) { return cuts.IsInside("ep_range", range, elab); },
                                       {"RangeHeavy", "EVertex"})};
    // Side events in Ep vs R20Mg plot
    auto nodeEpSide {nodeLat.Filter([&](float range, double elab) { return cuts.IsInside("ep_range", range, elab); },
                                    {"RangeHeavy", "EVertex"})};

    // Combine nodes
    auto nodeL1GatedSil {def.Filter(
        [&](ActRoot::MergerData& mer, float range, double elab)
        {
            if(mer.fLight.IsL1())
                return true;
            else
            {
                return cuts.IsInside("ep_range", range, elab);
            }
        },
        {"MergerData", "RangeHeavy", "EVertex"})};


    // Kinematics and Ex
    auto hKin {def.Histo2D(HistConfig::KinEl, "fThetaLight", "EVertex")};
    auto hKinCM {def.Histo2D(HistConfig::KinCM, "ThetaCM", "EVertex")};

    // Create vector of nodes and labels
    std::vector<std::string> labels {"All", "Lat", "Front", "L1"};
    std::vector<ROOT::RDF::RNode> nodes {def, nodeLat, nodeFront, nodel1};
    std::vector<ROOT::RDF::RNode> gatedNodes {nodeEpRSil, nodeEpSide, nodeEpFront, nodel1};

    // EBeam
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hsEBeam;
    for(int i = 0; i < labels.size(); i++)
    {
        auto h {nodes[i].Histo1D(HistConfig::EBeam, "EBeam")};
        hsEBeam.push_back(h);
    }
    // Rec EBeam
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hsRecEBeam;
    for(int i = 0; i < labels.size(); i++)
    {
        auto h {nodes[i].Histo1D(HistConfig::EBeam, "Rec_EBeam")};
        hsRecEBeam.push_back(h);
    }

    // Ex
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hsEx;
    for(int i = 0; i < labels.size(); i++)
    {
        auto h {nodes[i].Histo1D(HistConfig::Ex, "Ex")};
        hsEx.push_back(h);
    }
    // ECM
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hsECM;
    for(int i = 0; i < labels.size(); i++)
    {
        auto h {nodes[i].Histo1D(HistConfig::ECM, "ECM")};
        hsECM.push_back(h);
    }
    std::vector<ROOT::RDF::RResultPtr<TH2D>> hsECM2d;
    for(int i = 0; i < labels.size(); i++)
    {
        auto h {nodes[i].Histo2D(HistConfig::ThetaCMECM, "ThetaCM", "ECM")};
        hsECM2d.push_back(h);
    }
    // RPx
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hsRPx;
    for(int i = 0; i < labels.size(); i++)
    {
        auto h {nodes[i].Histo1D(HistConfig::RPx, "fRP.fCoordinates.fX")};
        hsRPx.push_back(h);
    }
    // Rec_ECM
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hsRecECM;
    for(int i = 0; i < labels.size(); i++)
    {
        auto h {nodes[i].Histo1D(HistConfig::ECM, "Rec_ECM")};
        hsRecECM.push_back(h);
    }
    // And now for gated on Ep vs R
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hsGatedRecECM;
    for(int i = 0; i < labels.size(); i++)
    {
        auto h {gatedNodes[i].Histo1D(HistConfig::ECM, "Rec_ECM")};
        hsGatedRecECM.push_back(h);
    }
    auto hThetaBeam {def.Histo2D(HistConfig::ThetaBeam, "fRP.fCoordinates.fX", "fThetaBeam")};
    auto hRP {def.Histo2D(HistConfig::RP, "fRP.fCoordinates.fX", "fRP.fCoordinates.fY")};
    auto hThetaCMLab {def.Histo2D(HistConfig::ThetaCMLab, "fThetaLight", "ThetaCM")};
    // Ex dependences
    auto hExThetaCM {def.Histo2D(HistConfig::ExThetaCM, "ThetaCM", "Ex")};
    auto hExThetaLab {def.Histo2D(HistConfig::ExThetaLab, "fThetaLight", "Ex")};
    auto hExRP {def.Histo2D(HistConfig::ExRPx, "fRP.fCoordinates.fX", "Ex")};
    auto hExZ {nodel0.Histo2D(HistConfig::ExZ, "fSP.fCoordinates.fZ", "Ex")};
    // Heavy histograms
    auto hThetaHLLab {def.Histo2D(HistConfig::ChangeTitle(HistConfig::ThetaHeavyLight, "Lab correlations"),
                                  "fThetaLight", "fThetaHeavy")};
    auto hRecECM {def.Histo1D(HistConfig::ECM, "Rec_ECM")};
    hRecECM->SetTitle("Rec E_{CM} with E_{x} = 0");
    auto hRangeHeavyEVertex {def.Histo2D(
        {"hRangeHeavyEVertex", "R vs E_{vertex};Range heavy [mm];E_{vertex} [MeV]", 150, 0, 300, 200, 0, 60},
        "RangeHeavy", "EVertex")};

    // Histograms for online analysis
    auto hRPxELab {
        nodel1.Filter("70 < fThetaLight && fThetaLight < 80")
            .Histo2D({"hRPxELab", "#theta_{Lab} in [70, 80];RP.X [mm];E_{Vertex} [#circ]", 400, 0, 260, 300, 0, 30},
                     "fRP.fCoordinates.fX", "EVertex")};

    auto hRPxThetaLab {
        nodel1.Histo2D({"hRPxThetaLab", "L1 exlusion zone;RP.X [mm];#theta_{Lab} [#circ]", 400, 0, 260, 300, 0, 90},
                       "fRP.fCoordinates.fX", "fThetaLight")};

    auto hECMRPx {nodeFront.Histo2D(HistConfig::RPxECM, "fRP.fCoordinates.fX", "ECM")};
    auto hRecECMRPx {def.Histo2D(HistConfig::RPxECM, "fRP.fCoordinates.fX", "Rec_ECM")};
    auto hEpRMg {def.Histo2D(HistConfig::EpRMg, "RangeHeavy", "EVertex")};
    // ECM from cuts in Ep vs R20Mg histo
    auto hECMCutSil {nodeEpRSil.Histo1D(HistConfig::ECM, "ECM")};
    auto hECMCutFront {nodeEpFront.Histo1D(HistConfig::ECM, "ECM")};
    auto hECMCutSide {nodeEpSide.Histo1D(HistConfig::ECM, "ECM")};


    // Save only the Ep_Range selection with silicons
    auto countSil {nodeEpRSil.Count()};
    auto countFront {nodeEpFront.Count()};
    auto outfile {TString::Format("./Outputs/tree_ex_%s_%s_%s_sil.root", beam.c_str(), target.c_str(), light.c_str())};
    // nodeL1GatedSil.Define("IsL1", [](ActRoot::MergerData& mer) { return mer.fLight.IsL1(); }, {"MergerData"})
    //     .Snapshot("Final_Tree", outfile);
    nodeEpRSil.Snapshot("Final_Tree", outfile);
    std::cout << "Saving Final_Tree in " << outfile << '\n';
    std::cout << "Counts Sil : " << countSil.GetValue() << '\n';
    std::cout << "Counts Front : " << countFront.GetValue() << '\n';

    // std::ofstream streamer {"./debug_ep_range.dat"};
    // auto nodeStreamer {def.Filter([&](double e, float range) { return cuts.IsInside("debug_ep_range", range, e); },
    //                               {"EVertex", "RangeHeavy"})};
    // auto hKinDebug {nodeStreamer.Histo2D(HistConfig::KinEl, "fThetaLight", "EVertex")};
    // auto hCorrEDebug {
    //     nodel0.Histo2D({"hDebug", ";TL;RP.X", 300, 0, 200, 300, 0, 200}, "RangeHeavy", "fRP.fCoordinates.fX")
    // nodel0.Filter([](ActRoot::MergerData& m) { return m.fLight.IsL1() == false; }, {"MergerData"})
    //     .Filter("120 <= RangeHeavy && RangeHeavy <= 140")
    //     .Define("ESil", [](ActRoot::MergerData& m) { return m.fLight.fEs.front(); }, {"MergerData"})
    //     .Histo2D({"hDebug", "Debug ep range;#theta_{Lab} [#circ];E_{Sil} [MeV]", 300, 0, 100, 300, 0, 15},
    //              "fThetaLight", "ESil")
    // .Histo2D({"hECorr", "Check E sil rec;E_{Sil} [MeV];E_{Vertex} [MeV]", 300, 0, 15, 300, 0, 15}, "ESil",
    //          "EVertex")
    // };
    // nodeStreamer.Foreach([&](ActRoot::MergerData& d) { d.Stream(streamer); }, {"MergerData"});
    // streamer.close();


    // Set styles for histograms
    std::vector<int> colors {-1, 8, 46, 9};
    for(int i = 0; i < hsEBeam.size(); i++)
    {
        auto& col {colors[i]};
        auto title {labels[i].c_str()};
        hsEBeam[i]->SetTitle(title);
        hsEBeam[i]->SetLineColor(col);

        hsRecEBeam[i]->SetTitle(title);
        hsRecEBeam[i]->SetLineColor(col);

        hsEx[i]->SetTitle(title);
        hsEx[i]->SetLineColor(col);

        hsECM[i]->SetTitle(title);
        hsECM[i]->SetLineColor(col);

        hsECM2d[i]->SetTitle(title);

        hsRPx[i]->SetTitle(title);
        hsRPx[i]->SetLineColor(col);

        hsRecECM[i]->SetTitle(title);
        hsRecECM[i]->SetLineColor(col);

        hsGatedRecECM[i]->SetTitle(title);
        hsGatedRecECM[i]->SetLineColor(col);

        if(i == 0)
        {
            hsEBeam[i]->SetTitle("E_{beam}");
            hsRecEBeam[i]->SetTitle("Rec E_{beam}");
            hsEx[i]->SetTitle("E_{x}");
            hsECM[i]->SetTitle("E_{CM}");
            hsECM2d[i]->SetTitle("E_{CM} vs #theta_{CM}");
            hsRPx[i]->SetTitle("RP_{x}");
            hsRecECM[i]->SetTitle("E_{CM} from proton");
            hsGatedRecECM[i]->SetTitle("E_{CM} from proton + EpR cut");
        }
    }


    auto* c22 {new TCanvas("c22", "Pipe2 canvas 2")};
    c22->DivideSquare(6);
    c22->cd(1);
    hRP->DrawClone();
    c22->cd(2);
    for(int i = 0; i < hsRPx.size(); i++)
        hsRPx[i]->DrawClone(i == 0 ? "" : "same");
    gPad->BuildLegend();
    c22->cd(3);
    for(int i = 0; i < hsRPx.size(); i++)
        hsEBeam[i]->DrawClone(i == 0 ? "" : "same");
    gPad->BuildLegend();
    c22->cd(4);
    for(int i = 0; i < hsRPx.size(); i++)
        hsRecEBeam[i]->DrawClone(i == 0 ? "" : "same");
    gPad->BuildLegend();
    // c22->cd(5);
    // hKinDebug->SetTitle("DebugKin in weird Ep_R plot");
    // hKinDebug->DrawClone("colz");
    // c22->cd(6);
    // hCorrEDebug->DrawClone("colz");

    auto* c21 {new TCanvas("c21", "Pipe2 canvas 1")};
    c21->DivideSquare(6);
    c21->cd(1);
    hKin->DrawClone("colz");
    auto* theo {kin.GetKinematicLine3()};
    kin.SetEx(0.6);
    auto* theoIne {kin.GetKinematicLine3()};
    theoIne->SetLineColor(46);
    theoIne->SetLineStyle(kDashed);
    theo->Draw("same");
    theoIne->Draw("same");
    c21->cd(2);
    for(int i = 0; i < hsRPx.size() - 1; i++)
    {
        if(i == 0)
            hsEx[i]->Add(hsEx.back().GetPtr(), -1); // substract L1 ex... plot it separately
        hsEx[i]->DrawClone(i == 0 ? "" : "same");
    }
    gPad->BuildLegend();
    c21->cd(3);
    hsEx.back()->SetTitle("E_{x} with L1");
    hsEx.back()->DrawClone();
    c21->cd(4);
    hExThetaLab->DrawClone("colz");
    c21->cd(5);
    hExZ->DrawClone("colz");
    c21->cd(6);
    hExRP->DrawClone("colz");

    // auto* c23 {new TCanvas {"c23", "Pipe2 canvas 3"}};
    // c23->DivideSquare(4);
    // c23->cd(1);
    // hThetaHLLab->DrawClone("colz");
    // c23->cd(2);
    // hThetaCMLab->DrawClone("colz");
    // c23->cd(3);

    auto* c24 {new TCanvas {"c24", "Pipe2 canvas 4"}};
    c24->DivideSquare(6);
    c24->cd(1);
    for(int i = 0; i < hsRPx.size(); i++)
        hsECM[i]->DrawClone(i == 0 ? "" : "same");
    gPad->BuildLegend();
    c24->cd(2);
    hECMRPx->DrawClone("colz");
    c24->cd(3);
    hRPxThetaLab->DrawClone("colz");
    c24->cd(4);
    hsECM2d.front()->DrawClone("colz");
    // hECM2dL1->DrawClone("colz");
    c24->cd(5);
    hEpRMg->DrawClone("colz");
    cuts.DrawCut("ep_range");
    // cuts.DrawCut("debug_ep_range");
    c24->cd(6);
    hECMCutSil->SetTitle("E_{CM} from RP.X + EpR cut");
    hECMCutSil->SetLineColor(1);
    hECMCutSil->DrawClone();
    hECMCutFront->SetLineColor(colors[2]);
    hECMCutFront->DrawClone("same");
    hECMCutSide->SetLineColor(colors[1]);
    hECMCutSide->DrawClone("same");
    //// hECMRPx->DrawClone("colz");

    auto* c25 {new TCanvas {"c25", "Pipe2 canvas 5"}};
    c25->DivideSquare(4);
    c25->cd(1);
    for(int i = 0; i < hsRPx.size(); i++)
        hsRecECM[i]->DrawClone(i == 0 ? "" : "same");
    gPad->BuildLegend();
    c25->cd(2);
    for(int i = 0; i < hsRPx.size() - 1; i++)
        hsGatedRecECM[i]->DrawClone(i == 0 ? "" : "same");
    gPad->BuildLegend();
    c25->cd(3);
    hsRecECM.back()->SetTitle("L1 E_{CM} from protons");
    hsRecECM.back()->DrawClone();


    // // Save to file
    // auto fout {std::make_unique<TFile>("./Outputs/gated_ecm.root", "update")};
    // hEpRMg->Write();
    // hECMCut->Write();
    // fout->Close();
}
#endif
