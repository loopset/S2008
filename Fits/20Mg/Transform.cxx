#include "ActKinematics.h"
#include "ActParticle.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"

#include <iostream>


TGraph* kinlabtocm(TGraph* gkin, TGraph* gcmlab)
{
    auto* gret {new TGraph};
    gret->SetTitle(";#theta_{CM} [#circ];E_{lab} [MeV]");
    // Invert function
    TF1 fcmlab {"fcmlab", [&](double* x, double* p) { return gcmlab->Eval(x[0]); }, 0, 180, 0};
    for(int p = 0; p < gkin->GetN(); p++)
    {
        auto lab {gkin->GetPointX(p)};
        auto e {gkin->GetPointY(p)};
        auto cm {fcmlab.GetX(lab)};
        gret->AddPoint(cm, e);
    }
    return gret;
}

void Transform()
{
    // Inverse kinematics
    ActPhysics::Kinematics inv {"20Mg(p,p)@84"};
    auto ginv {inv.GetThetaLabvsThetaCMLine()};
    auto gkininv {inv.GetKinematicLine3()};
    auto gkininvcm {kinlabtocm(gkininv, ginv)};
    inv.ComputeOtherInLab(0);

    // Calculate T2 needed for ECM to be equal in both cases
    ActPhysics::Particle p {"p"};
    ActPhysics::Particle mg {"20Mg"};
    double Tmg {84};
    double Tp {p.GetMass() / mg.GetMass() * Tmg};
    std::cout << "20Mg beam energy : " << Tmg << '\n';
    std::cout << "p beam energy : " << Tp << '\n';

    // Direct kinematics
    ActPhysics::Kinematics dir {"p(20Mg,p)@4.2"};
    auto gdir {dir.GetThetaLabvsThetaCMLine()};
    auto gkindir {dir.GetKinematicLine3()};
    auto gkindircm {kinlabtocm(gkindir, gdir)};

    // Equivalences
    auto* gequiv {new TGraph};
    gequiv->SetTitle("#theta_{lab} trans;Inverse;Direct");
    for(int p = 0; p < ginv->GetN(); p++)
    {
        // Inverse
        auto cm {ginv->GetPointX(p)};
        auto invlab {ginv->GetPointY(p)};
        // Direct
        auto dirlab {gdir->Eval(cm)};
        // Fill
        gequiv->AddPoint(invlab, dirlab);
    }
    auto* gequive {new TGraph};
    gequive->SetTitle("E_{lab} trans;Inverse;Direct");
    for(int p = 0; p < gkininvcm->GetN(); p++)
    {
        // Inverse
        auto cm {gkininvcm->GetPointX(p)};
        auto einv {gkininvcm->GetPointY(p)};
        // Direct
        auto edir {gkindircm->Eval(cm)};
        // Fill
        gequive->AddPoint(einv, edir);
    }

    auto* gmanual {new TGraph};
    auto* gmanuale {new TGraph};
    for(double thetacm = 0; thetacm <= 180; thetacm += 1)
    {
        inv.ComputeRecoilKinematics(thetacm * TMath::DegToRad(), 0);
        auto [theta, e] {inv.ComputeOtherInLab(thetacm * TMath::DegToRad())};
        gmanual->AddPoint(inv.GetTheta3Lab() * TMath::RadToDeg(), theta * TMath::RadToDeg());
        gmanuale->AddPoint(inv.GetT3Lab(), e);

    }

    auto* c0 {new TCanvas {"c0", "transform canvas"}};
    c0->DivideSquare(8);
    c0->cd(1);
    ginv->Draw("al");
    c0->cd(2);
    gkininvcm->Draw("al");
    c0->cd(3);
    gdir->Draw("al");
    c0->cd(4);
    gkindircm->Draw("al");
    c0->cd(5);
    gequiv->Draw("apl");
    c0->cd(6);
    gmanual->Draw("al");
    c0->cd(7);
    gequive->Draw("al");
    c0->cd(8);
    gmanuale->Draw("al");

    auto* c1 {new TCanvas {"c1", "e trans canvas"}};
    c1->DivideSquare(4);
    c1->cd(1);
    gkindir->Draw("al");
    c1->cd(2);
    // inv.GetOtherKinematics()->SetBeamEnergy(Tp);
    inv.GetOtherKinematics()->GetKinematicLine3()->Draw("al");

}
