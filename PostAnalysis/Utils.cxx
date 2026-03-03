#ifndef Utils_cxx
#define Utils_cxx

#include "ActSilMatrix.h"

namespace S2008
{
ActPhysics::SilMatrix* GetFrontMatrix()
{
    auto* sm {new ActPhysics::SilMatrix {"f0"}};
    sm->Read("/media/Data/S2008/Macros/SMs/Outputs/f0_matrix.root");
    return sm;
}
} // namespace S2008

#endif // !Utils_cxx
