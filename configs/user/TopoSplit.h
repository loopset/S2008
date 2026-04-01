#ifndef TopoSplit_h
#define TopoSplit_h

#include "ActContinuity.h"
#include "ActVAction.h"
#include "ActVoxel.h"

#include <memory>
#include <vector>


namespace ActAlgorithm
{
class TopoSplit : public VAction
{
public:
    // Copy to Continuity but with reduced NVoxels (= 1)
    std::shared_ptr<ActAlgorithm::Continuity> fRedCont {};
    // Parameters of the action
    double fChi2Thresh {};      //!< Chi2 thresh. Same as in BreakChi2
    int fSliceStep {2};         //!< Step in X to compute slices
    double fMedianFactor {2.5}; //!< If using median of max stddevs, median * fMedianFactor is the threshold
    double fStdDevThresh {-1};  //!< Disabled by default
    int fIdxOffset {-1};        //!< Number of slices left/right of threshold position
    double fWindowY {2};        //!< Beam window along Y
    double fWindowZ {2};        //!< Beam window along Z
    double fYRef {64};          //!< Default Y ref to get points closer to beam entrance along Y
    double fDebug {false};      //!< Enhanced fVerbose

public:
    TopoSplit() : VAction("TopoSplit") {}

    void ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block) override;
    void Run() override;
    void Print() const override;

private:
    XYZPointF FindPointCloserToY(const std::vector<ActRoot::Voxel>& covels, double xmin, double width, double yRef);
    void FallbackToBreakChi2();
};
} // namespace ActAlgorithm

#endif
