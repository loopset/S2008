#ifndef PCASplit_h
#define PCASplit_h

#include "ActVAction.h"


namespace ActAlgorithm
{
class PCASplit : public VAction
{
public:
    // Parameters of the action
public:
    PCASplit() : VAction("PCASplit") {}

    void ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block) override;
    void Run() override;
    void Print() const override;
};
} // namespace ActAlgorithm

#endif
