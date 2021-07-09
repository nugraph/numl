#pragma once

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include <tuple>

namespace numl {

  std::array<float, 3> GlobalToLocal(geo::WireID id, float time);

  void GetDUNE10ktGlobalWireTDC(detinfo::DetectorPropertiesData const& detProp,
                                unsigned int localWire, double localTDC, unsigned int plane, unsigned int tpc,
                                unsigned int& globalWire, unsigned int& globalPlane, double& globalTDC);

} // namespace numl
