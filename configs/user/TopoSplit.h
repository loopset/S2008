#ifndef TopoSplit_h
#define TopoSplit_h

#include "ActVAction.h"


namespace ActAlgorithm
{
class TopoSplit : public VAction
{
public:
    // Parameters of the action
    double fChi2Thresh {};      //!< Chi2 thresh. Same as in BreakChi2
    int fSliceStep {2};         //!< Step in X to compute slices
    double fMedianFactor {2.5}; //!< If using median of max stddevs, median * fMedianFactor is the threshold
    double fStdDevThresh {-1};  //!< Disabled by default
    double fWindowY {2};        //!< Beam window along Y
    double fWindowZ {2};        //!< Beam window along Z
    double fDebug {false};      //!< Enhanced fVerbose

public:
    TopoSplit() : VAction("TopoSplit") {}

    void ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block) override;
    void Run() override;
    void Print() const override;
};
} // namespace ActAlgorithm

#endif
