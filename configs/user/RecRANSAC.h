#include "ActVAction.h"

namespace ActAlgorithm
{
class RecRANSAC : public VAction
{
public:
    // Parameters of the action
    int fNIter {125};       //!< Number of iterations
    int fMinVoxels {10};   //!< Min of voxels
    double fDistThresh {3}; //!< Distance threshold

public:
    RecRANSAC() : VAction("RecRANSAC") {}

    void ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block) override;
    void Run() override;
    void Print() const override;
};
} // namespace ActAlgorithm
