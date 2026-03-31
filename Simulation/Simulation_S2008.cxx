#include "ActColors.h"
#include "ActCutsManager.h"
#include "ActKinematicGenerator.h"
#include "ActKinematics.h"
#include "ActParticle.h"
#include "ActRunner.h"
#include "ActSRIM.h"
#include "ActSilMatrix.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"

#include "Rtypes.h"

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"
#include "TProfile2D.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TTree.h"

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "/media/Data/S2008/PostAnalysis/Gates.cxx"
#include "/media/Data/S2008/PostAnalysis/HistConfig.h"
#include "/media/Data/S2008/PostAnalysis/Utils.cxx"

using XYZPoint = ROOT::Math::XYZPoint;
using XYZPointF = ROOT::Math::XYZPointF;
using XYZVector = ROOT::Math::XYZVector;
using XYZVectorF = ROOT::Math::XYZVectorF;

void ApplyNaN(double& val, double thresh = 0, const std::string& comment = "stopped")
{
    if(val <= thresh)
        val = std::nan(comment.c_str());
}

std::pair<XYZPoint, XYZPoint> SampleVertex(double meanZ, double sigmaZ, TH3D* h, double lengthX)
{

    // X is always common for both manners
    double Xstart {0};
    double Xrp {gRandom->Uniform() * lengthX};
    // Y depends completely on the method of calculation
    double Ystart {-1};
    double Yrp {-1};
    // Z of beam at entrance
    double Zstart {gRandom->Gaus(meanZ, sigmaZ)};
    double Zrp {-1};
    // Ystart in this case is sampled from the histogram itself!
    double thetaXY {};
    double thetaXZ {};
    h->GetRandom3(Ystart, thetaXY, thetaXZ);
    // Mind that Y is not centred in the histogram value!
    // Rp values are computed as follows:
    Yrp = Ystart - Xrp * TMath::Tan(thetaXY * TMath::DegToRad());
    Zrp = Zstart - Xrp * TMath::Tan(thetaXZ * TMath::DegToRad());
    XYZPoint start {Xstart, Ystart, Zstart};
    XYZPoint vertex {Xrp, Yrp, Zrp};
    return {std::move(start), std::move(vertex)};
}

void ApplySilRes(double& e, double sigma)
{
    e = gRandom->Gaus(e, sigma * TMath::Sqrt(e / 5.5));
}

double AngleWithNormal(const ROOT::Math::XYZVector& dir, const ROOT::Math::XYZVectorF& normal)
{
    auto dot {dir.Unit().Dot(normal.Unit())};
    return TMath::ACos(dot);
}


void Simulation_S2008(const std::string& beam, const std::string& target, const std::string& light, int neutronPS,
                      int protonPS, double T1, double Ex, bool standalone, int thread = -1)
{
    // set batch mode if not an independent function
    if(!standalone)
        gROOT->SetBatch(true);

    // Resolutions
    const double sigmaSil {0.060 / 2.355};
    const double sigmaPercentBeam {0};
    const double sigmaAngleLight {0.95 / 2.355};
    // Parameters of beam in mm
    // Center in Z is defined from silicon matrices
    const double zVertexSigma {2.84}; // mm

    // number of iterations
    const int iterations {static_cast<int>(standalone ? 1e6 : 3e7)};

    // Which parameters will be activated
    bool stragglingInGas {true};
    bool stragglingInSil {true};
    bool silResolution {true};
    bool thetaResolution {true};

    // TPC basic parameters
    ActRoot::TPCParameters tpc {"Actar"};

    // Vertex sampling
    auto beamfile {std::make_unique<TFile>("/media/Data/E796v2/Macros/Emittance/Outputs/histos.root")};
    auto* hBeam {beamfile->Get<TH3D>("h3d")};
    if(!hBeam)
        throw std::runtime_error("Simulation_S2008(): Could not load beam emittance histogram");
    hBeam->SetDirectory(nullptr);
    beamfile.reset();

    // Kinematics
    ActPhysics::Particle p1 {beam};
    ActPhysics::Particle p2 {target};
    ActPhysics::Particle p3 {light};
    // Automatically compute 4th particle
    ActPhysics::Kinematics kaux {p1, p2, p3};
    ActPhysics::Particle p4 {kaux.GetParticle(4)};
    // Binary kinematics generator
    ActSim::KinematicGenerator kingen {p1, p2, p3, p4, (protonPS > 0 ? protonPS : 0), (neutronPS > 0 ? neutronPS : 0)};
    kingen.Print();
    auto* kin {kingen.GetBinaryKinematics()};

    // Silicon specs
    auto* specs {new ActPhysics::SilSpecs};
    specs->ReadFile("/media/Data/S2008/configs/silspecs.conf");
    // Front silicons
    auto* f0sm {S2008::GetFrontMatrix()};
    auto silCentre = f0sm->GetMeanZ({4, 7});
    specs->GetLayer("f0").ReplaceWithMatrix(f0sm);
    auto beamOffset = -5.40; // mm
    double zVertexMean {silCentre + beamOffset};
    TString secondLayer {"f1"};
    // Left silicons
    auto* l0sm {S2008::GetLeftMatrix()};
    l0sm->MoveZTo(zVertexMean - 2.67, {4});
    specs->GetLayer("l0").ReplaceWithMatrix(l0sm);
    // Right silicons
    auto* r0sm {S2008::GetRightMatrix()};
    r0sm->MoveZTo(zVertexMean - 1.51, {4});
    specs->GetLayer("r0").ReplaceWithMatrix(r0sm);

    // CUTS ON SILICON ENERGY, depending on particle and layer
    std::map<std::string, std::pair<double, double>> eLoss0Cuts {};
    for(const auto& layer : {"f0", "l0", "r0"})
    {
        ActRoot::CutsManager<std::string> cut;
        cut.ReadCut(light, TString::Format("/media/Data/S2008/PostAnalysis/Cuts/pid_%s_%s_%s.root", light.c_str(),
                                           layer, beam.c_str())
                               .Data());
        if(cut.GetCut(light))
        {
            eLoss0Cuts[layer] = cut.GetXRange(light);
            std::cout << BOLDGREEN << "-> ESil in " << layer << " : " << light << ": [" << eLoss0Cuts[layer].first
                      << ", " << eLoss0Cuts[layer].second << "] MeV" << RESET << '\n';
        }
        else
        {
            std::cout << BOLDRED << "Simulation_S2008(): could not read PID cut for " << light
                      << " -> using default eLoss0Cut" << RESET << '\n';
            eLoss0Cuts[layer] = {0, 1000};
        }
    }

    // Histograms
    TH1::AddDirectory(false);
    // To compute a fine-grain efficiency, we require at least a binning width of 0.25 degrees!
    auto hThetaCM {HistConfig::ThetaCM.GetHistogram()};
    auto hThetaCMAll {HistConfig::ChangeTitle(HistConfig::ThetaCM, "ThetaCM all", "All").GetHistogram()};
    // And also efficiency in LAB!
    auto hThetaLabAll {HistConfig::ChangeTitle(HistConfig::ThetaCM, "ThetaLab all", "LabAll").GetHistogram()};
    auto hThetaLab {HistConfig::ChangeTitle(HistConfig::ThetaCM, "ThetaLab", "Lab").GetHistogram()};
    auto hDistF0 {HistConfig::ChangeTitle(HistConfig::TL, "Distance to F0").GetHistogram()};
    auto hKinVertex {HistConfig::ChangeTitle(HistConfig::KinSimu, "Kinematics at vertex").GetHistogram()};
    auto hKinSampled {HistConfig::ChangeTitle(HistConfig::KinSimu, "Sampled kinematics").GetHistogram()};
    std::map<std::string, std::shared_ptr<TH2D>> hsSP, hsSPTheta;
    for(const auto& layer : {"f0", "l0", "r0"})
    {
        hsSP[layer] = HistConfig::SP.GetHistogram();
        hsSP[layer]->SetNameTitle(TString::Format("hSP%s", layer), layer);

        hsSPTheta[layer] = std::make_unique<TProfile2D>(
            "hSPTheta", "SP vs #theta_{CM};Y [mm];Z [mm];#theta_{CM} [#circ]", 75, 0, 300, 75, 0, 300);
        hsSPTheta[layer]->SetTitle(layer);
    }
    auto hEexAfter {HistConfig::ChangeTitle(HistConfig::Ex, "Ex after resolutions").GetHistogram()};
    auto hRP {HistConfig::RP.GetHistogram()};
    auto hRPz {std::make_unique<TH2D>("hRPz", "RP;Y [mm];Z [mm]", 550, 0, 256, 550, 0, 256)};
    // Debug histograms
    auto hDeltaE {
        std::make_unique<TH2D>("hDeltaEE", "#Delta E - E;E_{in} [MeV];#Delta E_{0} [MeV]", 300, 0, 60, 300, 0, 60)};
    auto hELoss0 {std::make_unique<TH2D>("hELoss0", "ELoss0;E_{in} [MeV];#Delta E_{0} [MeV]", 200, 0, 40, 200, 0, 40)};
    auto hThetaLabNormal {std::make_unique<TH2D>("hThetaLabNormal",
                                                 "Theta in sil corr;#theta_{Lab} [#circ];#theta_{Normal Si} [#circ]",
                                                 250, 0, 90, 250, 0, 90)};
    auto hRPxEBeam {HistConfig::EBeamRPx.GetHistogram()};

    // Phi efficiency
    auto hPhiAll {std::make_unique<TH1D>("hPhiAll", "#phi eff;#phi [#circ];", 400, 0, 360)};
    auto hPhiLab {std::make_unique<TH1D>("hPhiLab", "#phi eff;#phi [#circ];", 400, 0, 360)};

    // 3D efficiency
    auto hEffAll {HistConfig::Eff2D.GetHistogram()};
    auto hEffAfter {HistConfig::Eff2D.GetHistogram()};
    hEffAfter->SetTitle("Efficiency");

    auto hEffDirAll {HistConfig::Eff2D.GetHistogram()};
    auto hEffDirAfter {HistConfig::Eff2D.GetHistogram()};
    hEffDirAfter->SetTitle("Efficiency;#theta_{3, Lab}^{dir} [#circ];T_{1}^{Dir} [MeV]");

    // ECM resoltion
    auto hECMRes {HistConfig::ECMECM.GetHistogram()};

    // T1Dir resolution
    auto hT1DirRes {HistConfig::ECMECM.GetHistogram()};
    hT1DirRes->SetTitle("T_{1} direct res;T_{1} sampled [MeV];T_{1} rec [MeV]");

    // Load SRIM tables
    // The name of the file sets particle + medium
    auto* srim {new ActPhysics::SRIM()};
    srim->ReadTable("light", "/media/Data/S2008/Calibrations/SRIM/1H_800mbar_95-5.txt"); // always 1H at 800 mbar?
                                                                                         // pressure will change
    srim->ReadTable("beam",
                    TString::Format("/media/Data/S2008/Calibrations/SRIM/%s_800mbar_95-5.txt", beam.c_str()).Data());
    srim->ReadTable("lightInSil", "/media/Data/S2008/Calibrations/SRIM/1H_silicon.txt");

    // Random generator
    gRandom->SetSeed();
    // Runner: contains utility functions to execute multiple actions
    ActSim::Runner runner(nullptr, nullptr, gRandom, 0);

    // Output from simulation!
    auto aux {thread != -1 ? TString::Format("_%d", thread).Data() : ""};
    auto filename {TString::Format("/media/Data/S2008/Simulation/Outputs/simu_%s_%s_%s%s.root", beam.c_str(),
                                   target.c_str(), light.c_str(), aux)};
    // We only store a few things in the TTree
    // 1-> Excitation energy
    // 2-> Theta in CM frame
    // 3-> Weight of the generator: for three-body reactions (phase spaces) the other two
    // variables need to be weighted by this value. For binary reactions, weight = 1
    // 4-> Energy at vertex
    // 5-> Theta in Lab frame
    auto* outFile {new TFile(filename, standalone ? "read" : "recreate")};
    auto* outTree {new TTree("SimulationTTree", "A TTree containing only our Eex obtained by simulation")};
    if(standalone)
        outTree->SetDirectory(nullptr);
    ROOT::Math::XYZPoint sp_tree {};
    outTree->Branch("SP", &sp_tree);
    double theta3CM_tree {};
    outTree->Branch("theta3CM", &theta3CM_tree);
    double Ex_tree {};
    outTree->Branch("Eex", &Ex_tree);
    double weight_tree {};
    outTree->Branch("weight", &weight_tree);
    double EVertex_tree {};
    outTree->Branch("EVertex", &EVertex_tree);
    double theta3Lab_tree {};
    outTree->Branch("theta3Lab", &theta3Lab_tree);
    double rpx_tree {};
    outTree->Branch("RPx", &rpx_tree);
    int silIdx_tree {};
    outTree->Branch("SilIdx", &silIdx_tree);
    double esil0_tree {};
    outTree->Branch("ESil0", &esil0_tree);
    // Lorentz vectors
    std::vector<ROOT::Math::XYZTVector> lor_tree {};
    outTree->Branch("Lor", &lor_tree);
    XYZVector beta_tree {};
    outTree->Branch("Beta", &beta_tree);

    //---- SIMULATION STARTS HERE
    ROOT::EnableImplicitMT();

    // timer
    TStopwatch timer {};
    timer.Start();
    // print fancy info
    std::cout << BOLDMAGENTA << "Running for Ex = " << Ex << " MeV" << RESET << '\n';
    std::cout << BOLDGREEN;
    const int percentPrint {5};
    int step {iterations / (100 / percentPrint)};
    int nextPrint {step};
    int percent {};
    int lightIn {};
    int heavyIn {};
    for(long int reaction = 0; reaction < iterations; reaction++)
    {
        // Print progress
        if(reaction >= nextPrint)
        {
            percent = 100 * (reaction + 1) / iterations;
            int nchar {percent / percentPrint};
            std::cout << "\r" << std::string((int)(percent / percentPrint), '|') << percent << "%";
            std::cout.flush();
            nextPrint += step;
        }
        // 1-> Sample vertex
        auto [start, vertex] {SampleVertex(zVertexMean, zVertexSigma, hBeam, tpc.X())};

        // 2-> Beam energy according to its sigma
        auto TBeam {runner.RandomizeBeamEnergy(
            T1 * p1.GetAMU(),
            sigmaPercentBeam * T1 * p1.GetAMU())}; // T1 in Mev / u * mass of beam in u = total kinetic energy
        // And slow according to distance travelled
        auto distToVertex {(vertex - start).R()};
        TBeam = srim->SlowWithStraggling("beam", TBeam, distToVertex);
        if(TBeam <= 0)
            continue;
        hRPxEBeam->Fill(vertex.X(), TBeam);

        // 3-> Run kinematics!
        kingen.SetBeamAndExEnergies(TBeam, Ex);
        double theta3Lab {};
        double phi3Lab {};
        double T3Lab {};
        double theta4Lab {};
        double phi4Lab {};
        double weight {1};
        // Uniform phi always and it is the same for CM and Lab
        auto phiCM {gRandom->Uniform(0, TMath::TwoPi())};
        // thetaCM following xs or not
        double thetaCM = TMath::ACos(gRandom->Uniform(-1, 1));
        kin->ComputeRecoilKinematics(thetaCM, phiCM);
        // Set info
        theta3Lab = kin->GetTheta3Lab();
        phi3Lab = phiCM;
        T3Lab = kin->GetT3Lab();
        theta4Lab = kin->GetTheta4Lab();
        phi4Lab = phiCM;
        // Compute ECM
        double ECM {TBeam * (p2.GetMass()) / (p1.GetMass() + p2.GetMass())};
        // Beam energy in direct reaction
        double T1Dir {kin->ComputeEquivalentOtherT1(TBeam)};

        // Efficiencies
        double thetaCMEff {kin->ReconstructTheta3CMFromLab(T3Lab, theta3Lab)};
        hThetaCMAll->Fill(thetaCMEff * TMath::RadToDeg());
        double theta3LabEff {theta3Lab}; // before implementing resolution in angle
        hThetaLabAll->Fill(theta3LabEff * TMath::RadToDeg());
        hPhiAll->Fill(phi3Lab * TMath::RadToDeg());
        hEffAll->Fill(thetaCMEff * TMath::RadToDeg(), ECM);
        // Theta3 in direct reaction
        double theta3LabDir {kin->ComputeOtherInLab(thetaCMEff).first * TMath::RadToDeg()};
        hEffDirAll->Fill(theta3LabDir, T1Dir);
        hRP->Fill(vertex.X(), vertex.Y());

        // 4-> Include thetaLab resolution to compute thetaCM and Ex afterwards
        if(thetaResolution)
            theta3Lab = gRandom->Gaus(theta3Lab, sigmaAngleLight * TMath::DegToRad());

        // 5-> Propagate track from vertex to silicon wall using SilSpecs class
        // And using the angle with the uncertainty already in
        ROOT::Math::XYZVector dirBeamFrame {TMath::Cos(theta3Lab), TMath::Sin(theta3Lab) * TMath::Sin(phi3Lab),
                                            TMath::Sin(theta3Lab) * TMath::Cos(phi3Lab)};
        ROOT::Math::XYZVector heavyBeamFrame {TMath::Cos(theta4Lab), TMath::Sin(theta4Lab) * TMath::Sin(phi4Lab),
                                              TMath::Sin(theta4Lab) * TMath::Cos(phi4Lab)};
        // Declare beam direction
        auto beamDir {(vertex - start).Unit()};
        // Rotate to world = geometry frame
        auto dirWorldFrame {runner.RotateToWorldFrame(dirBeamFrame, beamDir)};
        auto heavyWorldFrame {runner.RotateToWorldFrame(heavyBeamFrame, beamDir)};

        // Light particle
        int silIndex0 {-1};
        std::string firstLayer {};
        XYZPoint silPoint0InMM {};
        std::shared_ptr<ActPhysics::SilMatrix> sm {};
        for(const auto& layer : {"f0", "l0", "r0"})
        {
            auto [index, sp] = specs->FindSPInLayer(layer, vertex, dirWorldFrame);
            if(index != -1)
            {
                silIndex0 = index;
                firstLayer = layer;
                silPoint0InMM = sp;
                sm = specs->GetLayer(layer).GetSilMatrix();
                break;
            }
        }

        // skip tracks that doesn't reach silicons
        if(silIndex0 == -1)
            continue;
        lightIn++;

        // Heavy particle: NOT NEEDED BC HEAVY PARTICLES STOP BF SILICONS
        // // Only for front reactions, since we have verified that for l0 doesnt apply
        // auto [heavyIndex0, heavyPoint0] {specs->FindSPInLayer("f0", vertex, heavyWorldFrame)};
        // if(heavyIndex0 != -1)
        // {
        //     // This means we would measure multiplicity two in the layer -> skip event
        //     // We do not consider ELosses in gas bc heavy particle has always a large energy
        //     // that wont stop it on the gas
        //     // std::cout << "=========================" << '\n';
        //     // std::cout << "Vertex : " << vertex << '\n';
        //     // std::cout << "SP : " << heavyPoint0 << '\n';
        //     // std::cout << "Heavy reached layer theta: " << theta4Lab * TMath::RadToDeg() << '\n';
        //     heavyIn++;
        //     continue;
        // }

        // Apply SilMatrix cut
        if(TString(firstLayer).Contains("f"))
        {
            if(!sm->IsInside(silIndex0, silPoint0InMM.Y(), silPoint0InMM.Z()))
                continue;
        }
        else
        {
            if(!sm->IsInside(silIndex0, silPoint0InMM.X(), silPoint0InMM.Z()))
                continue;
        }

        // Define SP distance
        auto distance0 {(vertex - silPoint0InMM).R()};
        auto T3EnteringSil {srim->SlowWithStraggling("light", T3Lab, distance0)};
        ApplyNaN(T3EnteringSil);
        // nan if stopped in gas
        if(!std::isfinite(T3EnteringSil))
            continue;

        // First layer of silicons
        auto& layer {specs->GetLayer(firstLayer)};
        // Angle with normal
        auto normal {layer.GetNormal()};
        auto angleNormal0 {AngleWithNormal(dirWorldFrame, normal)};
        auto T3AfterSil0 {
            srim->SlowWithStraggling("lightInSil", T3EnteringSil, layer.GetUnit().GetThickness(), angleNormal0)};
        auto eLoss0 {T3EnteringSil - T3AfterSil0};
        // Apply resolution
        if(T3AfterSil0 != 0)
        {
            ApplySilRes(eLoss0, sigmaSil);
            T3AfterSil0 = T3EnteringSil - eLoss0;
        }
        ApplyNaN(eLoss0, layer.GetThresholds().at(1), "thresh"); // assuming common threshold for all
        // nan if bellow threshold
        if(!std::isfinite(eLoss0))
            continue;
        // Debug histograms
        hDeltaE->Fill(T3EnteringSil, eLoss0);
        hThetaLabNormal->Fill(theta3Lab * TMath::RadToDeg(), angleNormal0 * TMath::RadToDeg());

        // 6-> Same but to silicon layer 1 if exists
        double T3AfterInterGas {};
        double distance1 {};
        int silIndex1 {};
        ROOT::Math::XYZPoint silPoint1 {};
        double eLoss1 {};
        double T3AfterSil1 {};
        bool isPunch {};
        if(T3AfterSil0 > 0 && (firstLayer == "f0"))
        {
            // first, propagate in gas
            auto [silIndex1, silPoint1InMM] {specs->FindSPInLayer(secondLayer.Data(), vertex, dirWorldFrame)};
            if(silIndex1 == -1)
                continue;

            distance1 = (silPoint1InMM - silPoint0InMM).R();
            T3AfterInterGas = srim->SlowWithStraggling("light", T3AfterSil0, distance1);
            ApplyNaN(T3AfterInterGas);
            if(!std::isfinite(T3AfterInterGas))
                continue;

            // now, silicon if we have energy left
            if(T3AfterInterGas > 0)
            {
                // For S2008 angleNormal0 = angleNormal1 but this is not general
                auto angleNormal1 {angleNormal0};
                T3AfterSil1 = srim->SlowWithStraggling("lightInSil", T3AfterInterGas,
                                                       specs->GetLayer(secondLayer.Data()).GetUnit().GetThickness(),
                                                       angleNormal1);
                auto eLoss1 {T3AfterInterGas - T3AfterSil1};
                ApplySilRes(eLoss1, sigmaSil);
                T3AfterSil1 = T3AfterInterGas - eLoss1;
                ApplyNaN(eLoss1, specs->GetLayer(secondLayer.Data()).GetThresholds().at(1), "thresh");
                isPunch = true;
            }
        }

        // 7 -> Reconstruct energy at vertex
        // INFO: using reconstruction SRIM
        double EBefSil0 {};
        if(isPunch && T3AfterSil1 == 0 && std::isfinite(eLoss1))
        {
            double EAfterSil0 {srim->EvalInitialEnergy("light", eLoss1, distance1)};
            EBefSil0 = eLoss0 + EAfterSil0;
        }
        else if(!isPunch && T3AfterSil0 == 0)
            EBefSil0 = eLoss0;
        else
            EBefSil0 = -1;

        // 7->
        // we are ready to reconstruct Eex with all resolutions implemented
        //(d,light) is investigated gating on Esil1 = 0!
        // bool cutEAfterSil0 {EBefSil0 != -1 && !isPunch};
        bool cutEAfterSil0 {T3AfterSil0 == 0};
        bool cutELoss0 {eLoss0Cuts[firstLayer].first <= eLoss0 && eLoss0 <= eLoss0Cuts[firstLayer].second};
        if(cutEAfterSil0 && cutELoss0) // fill histograms
        {
            auto T3Recon {srim->EvalInitialEnergy("light", EBefSil0, distance0)};
            auto ExAfter {kin->ReconstructExcitationEnergy(T3Recon, theta3Lab)};
            auto thetaCM {kin->ReconstructTheta3CMFromLab(T3Recon, theta3Lab)};
            // Eval ECM
            auto T1Rec {kin->ReconstructBeamEnergyFromLabKinematics(T3Recon, theta3Lab)};
            auto ECMRec {T1Rec * (p2.GetMass()) / (p1.GetMass() + p2.GetMass())};
            auto T1DirRec {kin->ComputeEquivalentOtherT1(T1Rec)};
            hECMRes->Fill(ECM, ECMRec);
            hT1DirRes->Fill(T1Dir, T1DirRec);

            // fill histograms
            hDistF0->Fill(distance0);
            hKinSampled->Fill(theta3LabEff * TMath::RadToDeg(), T3Lab);
            hKinVertex->Fill(theta3Lab * TMath::RadToDeg(), T3Recon);
            hEexAfter->Fill(ExAfter, weight);
            if(TString(firstLayer).Contains("f"))
            {
                hsSP[firstLayer]->Fill(silPoint0InMM.Y(), silPoint0InMM.Z());
                hsSPTheta[firstLayer]->Fill(silPoint0InMM.Y(), silPoint0InMM.Z(), thetaCM * TMath::RadToDeg());
                hRPz->Fill(vertex.Y(), vertex.Z());
            }
            else
            {
                hsSP[firstLayer]->Fill(silPoint0InMM.X(), silPoint0InMM.Z());
                hsSPTheta[firstLayer]->Fill(silPoint0InMM.X(), silPoint0InMM.Z(), thetaCM * TMath::RadToDeg());
                // hRPz->Fill(vertex.X(), vertex.Z());
            }
            // RP histogram
            hELoss0->Fill(T3EnteringSil, eLoss0);

            // Efficiency computation: passed histogram
            // WARNING: we must use the original thetaCM generated by the kinematic generator
            // otherwise we could be biasing the efficiency: a original thetaCM could be reconstructed shifted
            // (which makes sense after implementing resolutions) and hence contributing to another bin!!!
            // Besides, this could cause errors when making the division: passed counts > 0 / all counts == 0!!
            hThetaCM->Fill(thetaCMEff * TMath::RadToDeg());
            hThetaLab->Fill(theta3LabEff * TMath::RadToDeg());
            hPhiLab->Fill(phi3Lab * TMath::RadToDeg());
            hEffAfter->Fill(thetaCMEff * TMath::RadToDeg(), ECM);
            hEffDirAfter->Fill(theta3LabDir, T1Dir);

            // Save lorentz vectors
            lor_tree.clear();
            for(int lor = 0, size = kingen.GetNt(); lor < size; lor++)
            {
                auto* vec {kingen.GetLorentzVector(lor)};
                lor_tree.emplace_back(*vec);
            }

            // write to TTree
            sp_tree = silPoint0InMM;
            Ex_tree = ExAfter;
            weight_tree = weight;
            theta3CM_tree = thetaCM * TMath::RadToDeg();
            EVertex_tree = T3Recon;
            theta3Lab_tree = theta3Lab * TMath::RadToDeg();
            rpx_tree = vertex.X();
            silIdx_tree = silIndex0;
            esil0_tree = eLoss0;
            beta_tree = *kingen.GetBeta();
            outTree->Fill();
        }
    }
    std::cout << "\r" << std::string(100 / percentPrint, '|') << 100 << "%";
    std::cout.flush();
    std::cout << RESET << '\n';
    std::cout << "HeavyIn/LightIn : " << (double)heavyIn / lightIn * 100 << '\n';

    // Efficiencies as quotient of histograms in TEfficiency class
    auto* eff {new TEfficiency(*hThetaCM, *hThetaCMAll)};
    eff->SetNameTitle("eff", TString::Format("#theta_{CM} eff E_{x} = %.2f MeV", Ex));
    auto* effLab {new TEfficiency(*hThetaLab, *hThetaLabAll)};
    effLab->SetNameTitle("effLab", TString::Format("#theta_{Lab} eff E_{x} = %.2f MeV", Ex));
    auto* effPhi {new TEfficiency(*hPhiLab, *hPhiAll)};
    effPhi->SetNameTitle("effPhi", TString::Format("#phi_{Lab} eff E_{x} = %.2f MeV", Ex));
    // Manual computation of efficiencies
    auto* hEff2D {(TH2D*)hEffAfter->Clone("hEff2D")};
    hEff2D->Divide(hEffAll.get());
    auto* hEffDir2D {(TH2D*)hEffDirAfter->Clone("hEffDir2D")};
    hEffDir2D->Divide(hEffDirAll.get());

    // SAVING
    if(!standalone)
    {
        outFile->cd();
        outTree->Write();
        eff->Write();
        effLab->Write();
        effPhi->Write();
        for(auto& [_, h] : hsSP)
            h->Write();
        hRP->Write("hRP");
        hEff2D->Write("hEff2D");
        hEffDir2D->Write("hEffDir2D");
        hECMRes->Write("hECMRes");
        hT1DirRes->Write("hT1DirRes");
        outFile->Close();
        delete outFile;
        outFile = nullptr;
    }

    // plotting
    if(standalone)
    {
        // draw theoretical kinematics
        ActPhysics::Kinematics theokin {p1, p2, p3, p4, T1 * p1.GetAMU(), Ex};
        auto* gtheo {theokin.GetKinematicLine3()};

        auto* c0 {new TCanvas("c0", "Canvas for inspection 0")};
        c0->DivideSquare(6);
        c0->cd(1);
        hThetaCM->DrawClone();
        c0->cd(2);
        hThetaCMAll->DrawClone();
        c0->cd(3);
        // hDistF0->DrawClone();
        hKinSampled->DrawClone("colz");
        gtheo->Draw("l");
        c0->cd(4);
        // hThetaCMAll->DrawClone();
        hRP->DrawClone("colz");
        c0->cd(5);
        hRPz->DrawClone("colz");
        f0sm->DrawClone();
        c0->cd(6);
        hELoss0->DrawClone("colz");
        // hThetaLabNormal->DrawClone("colz");

        auto* c1 {new TCanvas("cAfter", "Canvas for inspection 1")};
        c1->DivideSquare(6);
        c1->cd(1);
        hKinVertex->DrawClone("colz");
        gtheo->Draw("same");
        c1->cd(2);
        // hSP->DrawClone("col");
        hsSP["f0"]->DrawClone("col");
        f0sm->Draw();
        c1->cd(3);
        hEexAfter->DrawClone("hist");
        c1->cd(4);
        eff->Draw("apl");
        c1->cd(5);
        // hSPTheta->DrawClone("colz");
        // effLab->Draw("apl");
        effPhi->Draw("apl");
        c1->cd(6);
        hDeltaE->DrawClone("colz");

        auto* c2 {new TCanvas {"c2", "2D Eff canvas"}};
        c2->DivideSquare(6);
        TString opt {"colz"};
        c2->cd(1);
        hEff2D->DrawClone(opt);
        c2->cd(2);
        hEffDir2D->DrawClone(opt);
        c2->cd(3);
        hECMRes->DrawClone("colz");
        c2->cd(4);
        hT1DirRes->DrawClone("colz");
        c2->cd(5);
        hsSP["l0"]->DrawClone("col");
        l0sm->Draw();
        c2->cd(6);
        hsSP["r0"]->DrawClone("col");
        r0sm->Draw();
    }

    // deleting news
    delete srim;

    timer.Stop();
    timer.Print();
}
