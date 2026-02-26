#ifndef ConfMergerPerRun_h
#define ConfMergerPerRun_h

#include "ActVTask.h"

namespace ActAlgorithm
{
class ConfMergerPerRun : public VTask
{
public:
    ConfMergerPerRun() : VTask("ConfMergerPerRun") {}

    void ReadConfiguration() override {}
    bool Run() override;
    void Print() override {};
};
} // namespace ActAlgorithm
#endif // !ConfMergerPerRun_h
