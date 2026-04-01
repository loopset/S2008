#include "RecRANSAC.h"

#include "ActAFindRP.h"
#include "ActCluster.h"
#include "ActColors.h"
#include "ActInputParser.h"
#include "ActMultiAction.h"
#include "ActRANSAC.h"
#include "ActTPCData.h"
#include "ActVoxel.h"

#include "TRandom.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <memory>
#include <tuple>
#include <vector>

void ActAlgorithm::RecRANSAC::ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block)
{
    fIsEnabled = block->GetBool("IsEnabled");
    if(!fIsEnabled)
        return;
    if(block->CheckTokenExists("NIter"))
        fNIter = block->GetInt("NIter");
    if(block->CheckTokenExists("MinVoxels"))
        fMinVoxels = block->GetInt("MinVoxels");
    if(block->CheckTokenExists("DistThresh"))
        fDistThresh = block->GetDouble("DistThresh");
}

void ActAlgorithm::RecRANSAC::Run()
{
    if(!fIsEnabled)
        return;
    gRandom->SetSeed(1);

    // Get noise
    const auto& noise {fTPCData->fRaw};

    // Determine beam likes in event, calling DetermineBeamLikes from FindRP before
    if(fMultiAction->HasAction("FindRP"))
    {
        auto findrp {std::dynamic_pointer_cast<ActAlgorithm::Actions::FindRP>(fMultiAction->GetAction("FindRP"))};
        if(findrp)
            findrp->ExecInnerAction("DetermineBeamLikes");
    }

    // Mask voxels above beam region, because they have noise in Z
    auto itBeam {std::find_if(fTPCData->fClusters.begin(), fTPCData->fClusters.end(),
                              [&](const ActRoot::Cluster& cl) { return cl.GetIsBeamLike(); })};

    std::vector<ActRoot::Voxel> maskedNoise {};
    if(itBeam != fTPCData->fClusters.end())
    {
        float ymin, ymax;
        std::tie(ymin, ymax) = itBeam->GetYRange();
        if(fIsVerbose)
        {
            std::cout << "Yrange : " << ymin << ", " << ymax << '\n';
        }
        std::copy_if(noise.begin(), noise.end(), std::back_inserter(maskedNoise),
                     [&](const ActRoot::Voxel& v)
                     {
                         const auto& p {v.GetPosition()};
                         auto condY {p.Y() < ymin || p.Y() > ymax};
                         return condY;
                     });
    }
    else
        return; // cannot do anything

    // // Trigger only when all are beam-likes
    // auto trigger {std::all_of(fTPCData->fClusters.begin(), fTPCData->fClusters.end(),
    //                           [](const ActRoot::Cluster& cl) { return cl.GetIsBeamLike(); })};
    // if(!trigger)
    //     return;
    // // if(fTPCData->fClusters.size() > 1)
    // //     return;

    // Check if maskedNoise has enough voxels
    if(maskedNoise.size() < fMinVoxels)
        return;

    ActAlgorithm::RANSAC ransac {fNIter, fMinVoxels, fDistThresh};
    auto [clusters, back] {ransac.Run(maskedNoise)};
    if(fIsVerbose)
    {
        std::cout << BOLDGREEN << "-- RecRANSAC --" << '\n';
        std::cout << "Noise size   : " << noise.size() << '\n';
        std::cout << "Masked noise : " << maskedNoise.size() << '\n';
        std::cout << "New clusters : " << clusters.size() << RESET << '\n';
    }
    if(clusters.size())
    {
        // Get clusters with best chi2
        std::sort(clusters.begin(), clusters.end(), [](const ActRoot::Cluster& a, const ActRoot::Cluster& b)
                  { return a.GetLine().GetChi2() < b.GetLine().GetChi2(); });
        // Push the first one
        auto& cluster {clusters.front()};
        // Set flag
        cluster.SetFlag("IsRANSAC", true);
        // Push to vector
        fTPCData->fClusters.push_back(cluster);
    }
}

void ActAlgorithm::RecRANSAC::Print() const
{
    std::cout << BOLDCYAN << "····· " << GetActionID() << " ·····" << '\n';
    if(!fIsEnabled)
    {
        std::cout << "······························" << RESET << '\n';
        return;
    }
    std::cout << "  NIter      : " << fNIter << '\n';
    std::cout << "  MinVoxels  : " << fMinVoxels << '\n';
    std::cout << "  DistThresh : " << fDistThresh << RESET << '\n';
}

// Create symbol to load class from .so
extern "C" ActAlgorithm::RecRANSAC* CreateUserAction()
{
    return new ActAlgorithm::RecRANSAC;
}
