#include "numl/NeutrinoML/NeutrinoMLUtils.h"

#include "larcore/Geometry/Geometry.h"

using std::array;
// using std::make_tuple;

using art::ServiceHandle;
using geo::Geometry;
using geo::WireID;
using detinfo::DetectorPropertiesService;

namespace numl {

  //---------------------------------------------------------------------------
  array<float, 3> GlobalToLocal(WireID id, float time) {

    auto geo = ServiceHandle<Geometry>();

    TVector3 point = geo->TPC(id.TPC).Plane(id.Plane).GetCenter() - geo->TPC(0).Plane(id.Plane).GetCenter();
    std::cout << "TPC offset is " << point.X() << ", " << point.Y() << ", " << point.Z() << std::endl;

    std::cout << "Plane is " << id.Plane << ", sin φ_z is " << geo->TPC(id.TPC).Plane(id.Plane).SinPhiZ()
      << ", cos φ_z is " << geo->TPC(id.TPC).Plane(id.Plane).CosPhiZ() << std::endl;

    // unsigned int nWiresTPC = 400;
    // unsigned int wireGap = 4;
    // double driftLen = geo->TPC(tpc).DriftDistance();
    // double apaLen = geo->TPC(tpc).Width() - geo->TPC(tpc).ActiveWidth();
    // double driftVel = detprop->DriftVelocity();
    // unsigned int drift_size = (driftLen / driftVel) * 2; // Time in ticks to cross a TPC 
    // unsigned int apa_size   = 4*(apaLen / driftVel) * 2; // Width of the whole APA in TDC
    // // std::cout << "Drift size is " << drift_size << ", APA size is " << apa_size << std::endl;
    // globalWire = 0;
    // globalPlane = 0;
    // // Collection plane has more wires
    // if(plane == 2){
    //   nWiresTPC = 480;
    //   wireGap = 5;
    //   globalPlane = 2;
    // }
    // bool includeZGap = true;
    // if(includeZGap) nWiresTPC += wireGap;
    // // 10kt has four real TPCs and two dummies in each slice
    // //
    // //                 |--|-----|-----|-----|-----|--| /  /
    // //      y ^        |11| 10  |  9  |  8  |  7  | 6|/  /
    // //        | -| z   |--|-----|-----|-----|-----|--|  /
    // //        | /      | 5|  4  |  3  |  2  |  1  | 0| /
    // //  x <---|/       |--|-----|-----|-----|-----|--|/
    // //                     ^  wires  ^ ^  wires  ^
    // //
    // // We already filtered out the dummies, so we can assume 0->3 as follows:
    // //
    // //                 |-----|-----|-----|-----| /  /
    // //      y ^        |  7  |  6  |  5  |  4  |/  /
    // //        | -| z   |-----|-----|-----|-----|  /
    // //        | /      |  3  |  2  |  1  |  0  | /
    // //  x <---|/       |-----|-----|-----|-----|/
    // //                  ^  wires  ^ ^  wires  ^
    // //
    // size_t tpc_x = (tpc%6) - 1;   // x coordinate in 0->4 range
    // size_t tpc_xy = (tpc%12) - 1; // xy coordinate as 0->3 & 6->9 (converted from 1->4, 7->10)
    // if (tpc_xy > 3) tpc_xy -= 2;  // now subtract 2 so it's in the 0->7 range
    // // Induction views depend on the drift direction
    // if (plane < 2 and tpc%2 == 1) globalPlane = !plane;
    // else globalPlane = plane;
    // int offset = 752; // Offset between upper and lower modules in induction views, from Robert & Dorota's code
    // // Second induction plane gets offset from the back of the TPC
    // if (globalPlane != 1) globalWire += (tpc/12)*nWiresTPC;
    // else globalWire += ((300-tpc)/12)*nWiresTPC;
    // // Reverse wires and add offset for upper modules in induction views
    // if (tpc_xy > 3 and globalPlane < 2) globalWire += geo->Nwires(globalPlane, tpc, 0) + offset - localWire;
    // else globalWire += localWire;
    // if (tpc_x % 2 == 0) globalTDC = localTDC;
    // else globalTDC = (2*drift_size) - localTDC;
    // if (tpc_x > 1) globalTDC += 2 * (drift_size + apa_size);

    return std::array<float, 3> { (float)id.Plane, (float)id.Wire, time };

    // so

    // we want to get tpc position from the geometry service

    // how do we get the geometry service?

    // how do we get the wire pitch?

    // we should copy the old function here too

    // so it builds

  } // function GlobalToLocal

  //---------------------------------------------------------------------------
  void GetDUNE10ktGlobalWireTDC(detinfo::DetectorPropertiesData const& detProp,
                                unsigned int localWire, double localTDC, unsigned int plane, unsigned int tpc,
                                unsigned int& globalWire, unsigned int& globalPlane, double& globalTDC)
  {
    auto geo = ServiceHandle<Geometry>();

    unsigned int nWiresTPC = 400;
    unsigned int wireGap = 4;
    double driftLen = geo->TPC(tpc).DriftDistance();
    double apaLen = geo->TPC(tpc).Width() - geo->TPC(tpc).ActiveWidth();
    double driftVel = detProp.DriftVelocity();
    unsigned int drift_size = (driftLen / driftVel) * 2; // Time in ticks to cross a TPC 
    unsigned int apa_size   = 4*(apaLen / driftVel) * 2; // Width of the whole APA in TDC
    // std::cout << "Drift size is " << drift_size << ", APA size is " << apa_size << std::endl;
    globalWire = 0;
    globalPlane = 0;
    // Collection plane has more wires
    if(plane == 2){
      nWiresTPC = 480;
      wireGap = 5;
      globalPlane = 2;
    }
    bool includeZGap = true;
    if(includeZGap) nWiresTPC += wireGap;
    // 10kt has four real TPCs and two dummies in each slice
    //
    //                 |--|-----|-----|-----|-----|--| /  /
    //      y ^        |11| 10  |  9  |  8  |  7  | 6|/  /
    //        | -| z   |--|-----|-----|-----|-----|--|  /
    //        | /      | 5|  4  |  3  |  2  |  1  | 0| /
    //  x <---|/       |--|-----|-----|-----|-----|--|/
    //                     ^  wires  ^ ^  wires  ^
    //
    // We already filtered out the dummies, so we can assume 0->3 as follows:
    //
    //                 |-----|-----|-----|-----| /  /
    //      y ^        |  7  |  6  |  5  |  4  |/  /
    //        | -| z   |-----|-----|-----|-----|  /
    //        | /      |  3  |  2  |  1  |  0  | /
    //  x <---|/       |-----|-----|-----|-----|/
    //                  ^  wires  ^ ^  wires  ^
    //
    size_t tpc_x = (tpc%6) - 1;   // x coordinate in 0->4 range
    size_t tpc_xy = (tpc%12) - 1; // xy coordinate as 0->3 & 6->9 (converted from 1->4, 7->10)
    if (tpc_xy > 3) tpc_xy -= 2;  // now subtract 2 so it's in the 0->7 range
    // Induction views depend on the drift direction
    if (plane < 2 and tpc%2 == 1) globalPlane = !plane;
    else globalPlane = plane;
    int offset = 752; // Offset between upper and lower modules in induction views, from Robert & Dorota's code
    // Second induction plane gets offset from the back of the TPC
    if (globalPlane != 1) globalWire += (tpc/12)*nWiresTPC;
    else globalWire += ((300-tpc)/12)*nWiresTPC;
    // Reverse wires and add offset for upper modules in induction views
    if (tpc_xy > 3 and globalPlane < 2) globalWire += geo->Nwires(globalPlane, tpc, 0) + offset - localWire;
    else globalWire += localWire;
    if (tpc_x % 2 == 0) globalTDC = localTDC;
    else globalTDC = (2*drift_size) - localTDC;
    if (tpc_x > 1) globalTDC += 2 * (drift_size + apa_size);
  } // function PixelMapProducer::GetDUNE10ktGlobalWireTDC

} // namespace numl
