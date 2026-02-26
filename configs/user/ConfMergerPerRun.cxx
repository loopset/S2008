#include "ConfMergerPerRun.h"

#include "ActMergerDetector.h"

#include <stdexcept>

bool ActAlgorithm::ConfMergerPerRun::Run()
{
    // 20Mg runs
    double vd {}; // drift velocity in mm/us
    if(fMergerData->fRun <= 38)
        vd = 6.756;
    else if(fMergerData->fRun > 38 && fMergerData->fRun <= 64)
        vd = 6.747;
    else
        throw std::runtime_error("ConfMergerPerRun::Run(): 20Na runs not yet implemented");
    // Compute factor
    auto factor {vd * 4 * 0.08}; // conversion factor
    fMergerDet->SetDriftFactor(factor);
    if(fIsVerbose)
    {
        std::cout << "Setting drift factor : " << factor << '\n';
        std::cout << "Checking : " << fMergerDet->GetDriftFactor() << '\n';
    }

    return true;
}

extern "C" ActAlgorithm::ConfMergerPerRun* Create()
{
    return new ActAlgorithm::ConfMergerPerRun;
}
