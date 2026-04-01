#include "TopoSplit.h"

#include "ActCluster.h"
#include "ActColors.h"
#include "ActLine.h"
#include "ActMultiAction.h"
#include "ActTPCData.h"
#include "ActVAction.h"
#include "ActVoxel.h"

#include "TMath.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>

void ActAlgorithm::TopoSplit::ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block)
{
    fIsEnabled = block->GetBool("IsEnabled");
    if(!fIsEnabled)
        return;
    if(block->CheckTokenExists("Chi2Thresh"))
        fChi2Thresh = block->GetDouble("Chi2Thresh");
    if(block->CheckTokenExists("SliceStep", true))
        fSliceStep = block->GetInt("SliceStep");
    if(block->CheckTokenExists("MedianFactor", true))
        fMedianFactor = block->GetDouble("MedianFactor");
    if(block->CheckTokenExists("StdDevThresh", true))
        fStdDevThresh = block->GetDouble("StdDevThresh");
    if(block->CheckTokenExists("BeamWindowY"))
        fWindowY = block->GetDouble("BeamWindowY");
    if(block->CheckTokenExists("BeamWindowZ"))
        fWindowZ = block->GetDouble("BeamWindowZ");
    if(block->CheckTokenExists("Debug", true))
        fDebug = block->GetBool("Debug");
}

void ActAlgorithm::TopoSplit::Run()
{
    if(!fIsEnabled)
        return;
    std::vector<ActRoot::Cluster> toAppend {};

    for(auto it = fTPCData->fClusters.begin(); it != fTPCData->fClusters.end();)
    {
        // 1-> Check whether we meet conditions to execute this
        auto& line {it->GetLine()};
        bool isBadFit {line.GetChi2() > fChi2Thresh};
        if(!isBadFit)
        {
            it++;
            continue;
        }

        // 2-> Get ref to voxels
        auto& voxels {it->GetRefToVoxels()};

        // 3-> Calculate slices
        auto [xmin, xmax] {it->GetXRange()};
        int nSlices {static_cast<int>((xmax - xmin) / fSliceStep)};

        // 4-> Fit in slice
        std::vector<double> maxSigma;
        std::vector<XYZPointF> gps;
        for(int i = 0; i < nSlices; i++)
        {
            double x0 {xmin + i * fSliceStep};
            double x1 {x0 + fSliceStep};
            std::vector<ActRoot::Voxel> sliceVoxels {};
            std::copy_if(voxels.begin(), voxels.end(), std::back_inserter(sliceVoxels), [&](const ActRoot::Voxel& v)
                         { return x0 <= v.GetPosition().X() && v.GetPosition().X() < x1; });
            // Fit NOT WEIGHTED by charge
            ActRoot::Line sliceLine {};
            sliceLine.FitVoxels(sliceVoxels, false, true, true);
            // Store
            maxSigma.push_back(std::max(sliceLine.GetSigmas().Y(), sliceLine.GetSigmas().Z()));
            gps.push_back(sliceLine.GetPoint());
            if(fIsVerbose && fDebug)
            {
                std::cout << "[x0, x1] : [" << x0 << ", " << x1 << "] maxSigma : " << maxSigma.back() << '\n';
            }
        }

        // 5-> Find median
        auto median {TMath::Median(maxSigma.size(), maxSigma.data())};
        // Define threshold
        double threshold {(fStdDevThresh == -1) ? median * fMedianFactor : fStdDevThresh};
        auto maxIt {std::find_if(maxSigma.begin(), maxSigma.end(), [&](double& s) { return s > threshold; })};
        int splitIdx {(maxIt != maxSigma.end()) ? (int)std::distance(maxSigma.begin(), maxIt) : -1};

        // And assign gravity point to index before divergence
        XYZPointF gp {};
        if(splitIdx != -1) // found split
        {
            if(splitIdx > 0)
                gp = gps[splitIdx - 1];
            if(splitIdx == 0)
                gp = gps[0];
        }
        else // found no gp -> call BreakChi2 if possible
        {
            if(!fMultiAction->HasAction("BreakChi2")) // BreakChi2 not specified in ma.conf
                return;
            auto breakchi2 {fMultiAction->GetAction("BreakChi2")};
            auto status {breakchi2->GetIsEnabled()};
            breakchi2->SetIsEnabled(true);
            breakchi2->Run();
            breakchi2->SetIsEnabled(status);
            return;
        }

        // 6-> Find window point
        int xWindowSize {4};
        std::vector<ActRoot::Voxel> windowVoxels {};
        auto xaux {xmin};
        std::copy_if(voxels.begin(), voxels.end(), std::back_inserter(windowVoxels), [&](const ActRoot::Voxel& v)
                     { return xaux <= v.GetPosition().X() && v.GetPosition().X() < xaux + xWindowSize; });
        auto [clusters, _] {fAlgo->Run(windowVoxels)};
        // 7-> WP correspond to gravity point closer to middle of chamber in Y dir
        int yRef {64};
        auto minIt {std::min_element(clusters.begin(), clusters.end(),
                                     [&](const ActRoot::Cluster& a, const ActRoot::Cluster& b)
                                     {
                                         auto ya {a.GetLine().GetPoint().Y()};
                                         auto yb {b.GetLine().GetPoint().Y()};
                                         return std::abs(ya - yRef) < std::abs(yb - yRef);
                                     })};
        XYZPointF wp {-1, -1, -1};
        if(minIt != clusters.end())
            wp = minIt->GetLine().GetPoint();

        if(fIsVerbose)
        {
            std::cout << BOLDGREEN << "····· TopoSplit verbose ·····" << '\n';
            std::cout << "Median: " << median << '\n';
            std::cout << "Thresh: " << threshold << '\n';
            std::cout << "splitIdx: " << splitIdx << '\n';
            std::cout << "x0 : " << xmin + splitIdx * fSliceStep << '\n';
            std::cout << "WP: " << wp << '\n';
            std::cout << "GP: " << gp << RESET << '\n';

            if(fDebug) // just for visualization purposes
            {
                fTPCData->fRPs.push_back(wp);
                fTPCData->fRPs.push_back(gp);
            }
        }

        // 8-> And separate and run continuity (same as BreackChi2)
        // And define cylinder
        ActRoot::Line cylinder {gp, wp};
        auto toMove {std::partition(voxels.begin(), voxels.end(),
                                    [&](const ActRoot::Voxel& voxel)
                                    {
                                        auto pos {voxel.GetPosition()};
                                        pos += XYZVectorF {0.5, 0.5, 0.5};
                                        auto proj {cylinder.ProjectionPointOnLine(pos)};
                                        auto condY {std::abs(proj.Y() - pos.Y()) <= fWindowY};
                                        auto condZ {std::abs(proj.Z() - pos.Z()) <= fWindowZ};
                                        return condY && condZ;
                                    })};
        // Create vector to move to
        std::vector<ActRoot::Voxel> notBeam {};
        notBeam.insert(notBeam.end(), std::make_move_iterator(toMove), std::make_move_iterator(voxels.end()));
        voxels.erase(toMove, voxels.end());

        // Check if it satisfies minimum voxel requirement to be clusterized again
        if(voxels.size() <= fAlgo->GetMinPoints())
        {
            // Delete remaining beam-like cluster
            it = fTPCData->fClusters.erase(it);
        }
        else
        {
            // Reprocess
            it->ReFit();
            // Reset ranges
            it->ReFillSets();
            // Set it is beam-like
            // INFO: update March 2025. Disable this since it can
            // create false positives of beam-likes
            // it->SetBeamLike(true);

            // And of course, add to iterator
            it++;
        }
        // 4-> Run cluster algorithm again in the not beam voxels
        std::vector<ActRoot::Cluster> newClusters;
        std::vector<ActRoot::Voxel> noise;
        std::tie(newClusters, noise) = fAlgo->Run(notBeam);
        // Set flag accordingly
        for(auto& cl : newClusters)
            cl.SetIsBreakBeam(true);
        // Move to vector
        toAppend.insert(toAppend.end(), std::make_move_iterator(newClusters.begin()),
                        std::make_move_iterator(newClusters.end()));
        if(fIsVerbose)
        {
            std::cout << BOLDGREEN << "New clusters   = " << newClusters.size() << '\n';
            std::cout << "-------------------" << RESET << '\n';
        }
    }
    // Write to vector of clusters
    fTPCData->fClusters.insert(fTPCData->fClusters.end(), std::make_move_iterator(toAppend.begin()),
                               std::make_move_iterator(toAppend.end()));
}

void ActAlgorithm::TopoSplit::Print() const
{
    std::cout << BOLDCYAN << "····· " << GetActionID() << " ·····" << '\n';
    if(!fIsEnabled)
    {
        std::cout << "······························" << RESET << '\n';
        return;
    }
    std::cout << "  Chi2Thresh : " << fChi2Thresh << '\n';
    std::cout << "  SliceStep  : " << fSliceStep << '\n';
    std::cout << "  MedianFact : " << fMedianFactor << '\n';
    std::cout << "  StdDevThre : " << fStdDevThresh << '\n';
    std::cout << "  WindowY    : " << fWindowY << '\n';
    std::cout << "  WindowZ    : " << fWindowZ << '\n';
}

extern "C" ActAlgorithm::TopoSplit* CreateUserAction()
{
    return new ActAlgorithm::TopoSplit;
}
