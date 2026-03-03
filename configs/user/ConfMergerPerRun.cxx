#include "ConfMergerPerRun.h"

#include "ActMergerDetector.h"


bool ActAlgorithm::ConfMergerPerRun::Run()
{
    // 20Mg runs
    double vd {}; // drift velocity in mm/us
    if(fMergerData->fRun <= 38)
        vd = 6.756;
    else // applies for both 20Mg and 20Na
        vd = 6.747;
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
