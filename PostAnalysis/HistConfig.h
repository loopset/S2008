#ifndef HistConfig_h
#define HistConfig_h

#include "ROOT/RDF/HistoModels.hxx"

#include "TString.h"

namespace HistConfig
{
using namespace ROOT::RDF;

const TH2DModel PID {"hPID", "PID;E_{Sil} [MeV];Q_{ave} [mm^{-1}]", 400, 0, 40, 800, 0, 2000};

const TH2DModel PIDTwo {"hPIDTwo", "PID with two silicons;E_{1} [MeV];E_{0} [MeV]", 800, 0, 40, 800, 0, 40};

const TH2DModel SP {"hSP", "SP;X or Y [mm];Z [mm]", 200, -10, 300, 200, -10, 300};

const TH2DModel RP {"hRP", "RP;X [mm];Y [mm]", 200, -10, 300, 200, -10, 300};

const TH1DModel RPx {"hRPx", "RPx;X [mm];Counts", 200, -10, 300};

const TH1DModel TL {"hTL", "Track length; TL [mm]", 300, 0, 600};

const TH2DModel Kin {"hKin", "Kinematics;#theta_{Lab} [#circ];E_{Vertex} [MeV]", 250, 0, 60, 250, 0, 20};

const TH2DModel KinEl {"hKinEl", "Kinematics;#theta_{Lab} [#circ];E_{Vertex} [MeV]", 600, 0, 180, 400, 0, 20};

const TH2DModel KinSimu {"hKin", "Simulation kinematics;#theta_{Lab} [#circ];E_{Vertex} [MeV]", 600, 0, 90, 600, 0, 40};

const TH2DModel KinCM {"hKinCM", "CM kinematics;#theta_{CM} [#circ];E_{Vertex} [MeV]", 400, 0, 180, 400, 0, 20};

const TH1DModel Ex {"hEx", TString::Format("Excitation energy;E_{x} [MeV];Counts / %.f keV", (10. - (-5.)) / 200 * 1e3),
                    200, -5, 10};

const TH1DModel ThetaCM {"hThetaCM", "ThetaCM;#theta_{CM} [#circ]", 600, 0, 180};

const TH2DModel ZThetaZ {"hZThetaZ", "Emittance along Z;Z [mm];#theta_{Z} [#circ]", 600, 0, 270, 600, -10, 10};

const TH2DModel YPhiY {"hYPhiY", "Emittance along Y;Y [mm];#phi_{Y} [#circ]", 600, 0, 270, 600, -10, 10};

const TH2DModel ThetaBeam {
    "hThetaBeam", "#theta_{Beam} against RP.X;RP.X() [mm];#theta_{Beam} [#circ]", 200, -5, 270, 200, -1, 10};

const TH2DModel ExZ {"hExZ", "E_{x} vs SP.Z();SP.Z() [mm];E_{x} [MeV]", 300, -10, 450, 200, -5, 10};

const TH2DModel ExThetaCM {
    "hExThetaCM", "E_{x} vs #theta_{CM};#theta_{CM} [#circ];E_{x} [MeV]", 400, 0, 180, 200, -10, 10};

const TH2DModel ExThetaLab {
    "hExThetaLab", "E_{x} vs #theta_{Lab};#theta_{Lab} [#circ];E_{x} [MeV]", 400, 0, 100, 200, -10, 10};

const TH2DModel ExRPx {"hExRPX", "E_{x} vs RP.X;RP.X() [mm];E_{x} [MeV]", 200, -10, 300, 200, -10, 10};

const TH2DModel ThetaHeavyLight {
    "hThetaHL", "#theta heavy vs light;#theta_{Light} [#circ];#theta_{Heavy} [#circ]", 400, 0, 90, 400, 0, 15};

const TH2DModel ThetaCMLab {
    "hThetaCMLab", "CM vs Lab correlations;#theta_{Lab} [#circ];#theta_{CM} [#circ]", 400, 0, 90, 400, 0, 90};

const TH2DModel RPxThetaCM {
    "hRPxThetaCM", "RP.X vs #theta_{CM} correlations;RP.X [mm];#theta_{CM} [#circ]", 200, 0, 300, 100, 0, 60};

const TH2DModel RPxECM {"hRPxECM", "ECM vs RP.X;RP.X [mm];E_{CM} [MeV]", 200, 0, 300, 150, 0, 15};

const TH1DModel ECM {"hECM", "E_{CM};E_{CM} [MeV];Counts / 10 keV", 500, 0, 5};

const TH2DModel ThetaCMECM {
    "hThetaCMECM", "#theta_{CM} vs E_{CM};#theta_{CM} [#circ];E_{CM} [MeV]", 400, 0, 180, 160, 0, 8};

const TH2DModel EpRMg {"hEpRMg", "Ep vs R Mg; Total path ^{20}Mg [mm];E_{light} [MeV]", 150, 0, 300, 300, 0, 20};

const TH1DModel EBeam {"hEBeam", "Beam energy;E_{Beam} [MeV]", 180, 0, 90};

const TH2DModel EBeamRPx {"hEBeamRPx", "Beam energy;RP.X [mm];E_{Beam} [MeV]", 400, 0, 270, 180, 0, 90};

const TH2DModel Eff2D {"hEff2D", "2D efficiency;#theta_{CM} [#circ];E_{CM} [MeV]", 720, 0, 180, 320, 0, 8};


const TH2DModel ECMECM {"hECMRes", "E_{CM} resolution;#E_{CM}^{nom} [MeV];E_{CM} [MeV]", 500, 0, 5, 500, 0, 5};
template <typename T>
T ChangeTitle(T model, const TString& title, const TString& label = "");
} // namespace HistConfig

template <typename T>
T HistConfig::ChangeTitle(T model, const TString& title, const TString& label)
{
    auto ret {model};
    if(label.Length() > 0)
        ret.fName = model.fName + label;
    TString old {model.fTitle};
    auto end {old.First(';')};
    TString nt {title + old(end, old.Length() - end)};
    ret.fTitle = nt;
    return ret;
}

#endif
