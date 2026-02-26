#ifndef Gates_h
#define Gates_h

#include "ActMergerData.h"
namespace S2008
{
auto isFront {[](ActRoot::MergerData& m)
              {
                  if(m.fLight.GetNLayers() == 1)
                      if(m.fLight.GetLayer(0) == "f0")
                          return true;
                  return false;
              }};
auto isLeft {[](ActRoot::MergerData& m)
             {
                 if(m.fLight.GetNLayers() == 1)
                     if(m.fLight.GetLayer(0) == "l0")
                         return true;
                 return false;
             }};
auto isRight {[](ActRoot::MergerData& m)
              {
                  if(m.fLight.GetNLayers() == 1)
                      if(m.fLight.GetLayer(0) == "r0")
                          return true;
                  return false;
              }};
} // namespace S2008

#endif // !Gates_h
