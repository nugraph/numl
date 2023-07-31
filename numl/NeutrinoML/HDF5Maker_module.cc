////////////////////////////////////////////////////////////////////////
// Class:       HDF5Maker
// Plugin Type: analyzer (art v3_06_03)
// File:        HDF5Maker_module.cc
//
// Generated at Wed May  5 08:23:31 2021 by Jeremy Hewes using cetskelgen
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////


#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"

#include "hep_hpc/hdf5/make_ntuple.hpp"
using std::array;
using std::endl;
using std::setfill;
using std::set;
using std::setw;
using std::string;
using std::vector;

using simb::MCParticle;
using simb::MCTruth;
using sim::TrackIDE;
using recob::Hit;
using recob::OpHit;
using recob::OpFlash;
using recob::SpacePoint;
using recob::Wire;

using mf::LogInfo;

using art::ServiceHandle;
using cheat::BackTrackerService;
using cheat::ParticleInventoryService;
using detinfo::DetectorClocksService;
using detinfo::DetectorPropertiesService;

using hep_hpc::hdf5::Column;
using hep_hpc::hdf5::make_ntuple;
using hep_hpc::hdf5::make_scalar_column;
using hep_hpc::hdf5::make_column;

class HDF5Maker : public art::EDAnalyzer {
public:
  explicit HDF5Maker(fhicl::ParameterSet const& p);
  ~HDF5Maker() noexcept {}; // bare pointers are cleaned up by endSubRun

  HDF5Maker(HDF5Maker const&) = delete;
  HDF5Maker(HDF5Maker&&) = delete;
  HDF5Maker& operator=(HDF5Maker const&) = delete;
  HDF5Maker& operator=(HDF5Maker&&) = delete;

  void beginSubRun(art::SubRun const& sr) override;
  void endSubRun(art::SubRun const& sr) override;
  void analyze(art::Event const& e) override;

private:

  std::string fTruthLabel;
  std::string fHitLabel;
  std::string fHitTruthLabel;
  std::string fSPLabel;

  bool fUseMap;
  std::string fEventInfo;
  std::string fOutputName;

  hep_hpc::hdf5::File fFile;  ///< output HDF5 file

  hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>     // event id (run, subrun, event)
  >* fEventNtuple; ///< event ntuple

  hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
                        hep_hpc::hdf5::Column<int, 1>,    // is cc
			hep_hpc::hdf5::Column<int, 1>,     // pdg code
                        hep_hpc::hdf5::Column<float, 1>,  // nu energy
                        hep_hpc::hdf5::Column<float, 1>,  // lep energy
                        hep_hpc::hdf5::Column<float, 1>,   // nu dir (x, y, z)
			hep_hpc::hdf5::Column<float, 1>, // nu vertex position (x, y, z)
			hep_hpc::hdf5::Column<float, 1>, // nu vertex position, corrected (x, y, z) 
			hep_hpc::hdf5::Column<int, 1>, // nu vertex corrected wire pos
			hep_hpc::hdf5::Column<float, 1> // nu vertex corrected wire time
  >* fEventNtupleNu; ///< event ntuple with neutrino information

  hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
                        hep_hpc::hdf5::Column<int, 1>,    // spacepoint id
                        hep_hpc::hdf5::Column<float, 1>,  // 3d position (x, y, z)
                        hep_hpc::hdf5::Column<int, 1>     // 2d hit (u, v, y)
  >* fSpacePointNtuple; ///< spacepoint ntuple

  hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
                        hep_hpc::hdf5::Column<int, 1>,    // hit id
                        hep_hpc::hdf5::Column<float, 1>,  // hit integral
                        hep_hpc::hdf5::Column<float, 1>,  // hit rms
                        hep_hpc::hdf5::Column<int, 1>,    // tpc id
                        hep_hpc::hdf5::Column<int, 1>,    // global plane
                        hep_hpc::hdf5::Column<float, 1>,  // global wire
                        hep_hpc::hdf5::Column<float, 1>,  // global time
                        hep_hpc::hdf5::Column<int, 1>,    // raw plane
                        hep_hpc::hdf5::Column<float, 1>,  // raw wire
                        hep_hpc::hdf5::Column<float, 1>   // raw time
  >* fHitNtuple; ///< hit ntuple

  hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
                        hep_hpc::hdf5::Column<int, 1>,    // g4 id
                        hep_hpc::hdf5::Column<int, 1>,    // particle type
                        hep_hpc::hdf5::Column<int, 1>,    // parent g4 id
                        hep_hpc::hdf5::Column<float, 1>,  // momentum
                        hep_hpc::hdf5::Column<float, 1>,  // start position (x, y, z)
                        hep_hpc::hdf5::Column<float, 1>,  // end position (x, y, z)
                        hep_hpc::hdf5::Column<std::string, 1>, // start process
                        hep_hpc::hdf5::Column<std::string, 1>  // end process
  >* fParticleNtuple; ///< particle ntuple

  hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
                        hep_hpc::hdf5::Column<int, 1>,    // hit id
                        hep_hpc::hdf5::Column<int, 1>,    // g4 id
                        hep_hpc::hdf5::Column<float, 1>   // deposited energy [ MeV ]
  >* fEnergyDepNtuple; ///< energy deposition ntuple
};

  using ProxyPfpColl_t = decltype(proxy::getCollection<std::vector<recob::PFParticle> >(
                std::declval<art::Event>(),std::declval<art::InputTag>(),
                proxy::withAssociated<larpandoraobj::PFParticleMetadata>(std::declval<art::InputTag>()),
                proxy::withAssociated<recob::Slice>(std::declval<art::InputTag>()),
                proxy::withAssociated<recob::Cluster>(std::declval<art::InputTag>()),
                proxy::withAssociated<recob::Vertex>(std::declval<art::InputTag>()),
                proxy::withAssociated<anab::T0>(std::declval<art::InputTag>()) ));
  using ProxyPfpElem_t = ProxyPfpColl_t::element_proxy_t;

  // proxy to connect cluster to hit
  using ProxyClusColl_t = decltype(proxy::getCollection<std::vector<recob::Cluster>>(
      std::declval<art::Event>(), std::declval<art::InputTag>(),
      proxy::withAssociated<recob::Hit>(std::declval<art::InputTag>())));
  using ProxyClusElem_t = ProxyClusColl_t::element_proxy_t;

  void True2RecoMappingXYZ(float& t, float& x, float& y, float& z);
  void ApplySCEMappingXYZ(float& x, float& y, float& z);
  void AddDaughters(const std::map<unsigned int, unsigned int>& pfpmap, const ProxyPfpElem_t &pfp_pxy, const ProxyPfpColl_t &pfp_pxy_col, std::vector<ProxyPfpElem_t> &slice_v);
  float GetMetaData(const ProxyPfpElem_t &pfp_pxy, string metaDataName);
  int NearWire(const geo::Geometry& geo, const geo::PlaneID &ip, const float x, const float y, const float z);

  struct BtPart {
  public:

    BtPart(const int pdg_, const int category_, const float px_, const float py_, const float pz_, const float e_,
	   const std::vector<unsigned int> &tids_, const float start_x_, const float start_y_, const float start_z_, const float start_t_) :
      pdg(pdg_),category(category_),px(px_),py(py_),pz(pz_),e(e_),
      tids(tids_), start_x(start_x_), start_y(start_y_), start_z(start_z_), start_t(start_t_)
      {}

    BtPart(const int pdg_, const int category_, const float px_, const float py_, const float pz_, const float e_,
	   const unsigned int tid_, const float start_x_, const float start_y_, const float start_z_, const float start_t_) :
      pdg(pdg_), category(category_), px(px_), py(py_), pz(pz_), e(e_),
      start_x(start_x_), start_y(start_y_), start_z(start_z_), start_t(start_t_)
      { tids.push_back(tid_); }

    int pdg;
    int category;
    float px, py, pz, e;
    std::vector<unsigned int> tids;
    int nhits = 0;
    float start_x, start_y, start_z, start_t;
  };

  enum ParticleCategory {
    Pion     = 0,
    Muon     = 1,
    Kaon     = 2,
    Proton   = 3,
    Electron = 4,
    Michel   = 5,
    Delta    = 6,
    OtherNu  = 7,
    Photon   = 8
  };

  std::vector<BtPart> initBacktrackingParticleVec(const std::vector<sim::MCShower> &inputMCShower,
						  const std::vector<sim::MCTrack> &inputMCTrack,
						  const std::vector<recob::Hit> &inputHits,
						  const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart,
						  int nhitcut = 5);

};


HDF5Maker::HDF5Maker(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fTruthLabel(p.get<std::string>("TruthLabel")),
    fHitLabel(  p.get<std::string>("HitLabel")),
    fHitTruthLabel(  p.get<std::string>("HitTruthLabel","")),
    fSPLabel(   p.get<std::string>("SPLabel")),
    fUseMap(    p.get<bool>("UseMap", false)),
    fEventInfo( p.get<std::string>("EventInfo")),
    fOutputName(p.get<std::string>("OutputName"))
{
  if (fEventInfo != "none" && fEventInfo != "nu")
    throw art::Exception(art::errors::Configuration)
      << "EventInfo must be \"none\" or \"nu\", not " << fEventInfo;
}

void HDF5Maker::analyze(art::Event const& e)
{
  const cheat::BackTrackerService* bt = 0;
  if (!fUseMap) {
    art::ServiceHandle<cheat::BackTrackerService> bt_h;
    bt = bt_h.get();
  }

  int run = e.id().run();
  int subrun = e.id().subRun();
  int event = e.id().event();

  std::array<int, 3> evtID { run, subrun, event };

  // Fill event table
  if (fEventInfo == "none") {
    fEventNtuple->insert( evtID.data() );
    mf::LogInfo("HDF5Maker") << "Filling event table"
                             << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                             << ", event " << evtID[2];
  }
  if (fEventInfo == "nu") {
    // Get MC truth
    art::Handle<std::vector<simb::MCTruth>> truthHandle;
    e.getByLabel(fTruthLabel, truthHandle);
    if (!truthHandle.isValid() || truthHandle->size() == 0) {
      throw art::Exception(art::errors::LogicError)
        << "Expected to find exactly one MC truth object!";
    }
    simb::MCNeutrino nutruth = truthHandle->at(0).GetNeutrino();
    array<float, 3> nuDirection {
      (float)nutruth.Nu().Momentum().Vect().Unit().X(),
      (float)nutruth.Nu().Momentum().Vect().Unit().Y(),
      (float)nutruth.Nu().Momentum().Vect().Unit().Z()
    };

    array<float, 3> nuVtx { (float)nutruth.Nu().Vx(), (float)nutruth.Nu().Vy(), (float)nutruth.Nu().Vz() };
    
    float nuT = nutruth.Nu().T();
    array<float, 3> nuVtxCorr { (float)nutruth.Nu().Vx(), (float)nutruth.Nu().Vy(), (float)nutruth.Nu().Vz() };
    True2RecoMappingXYZ(nuT, nuVtxCorr[0], nuVtxCorr[1], nuVtxCorr[2]);

    vector<int> nearwires;
    for (auto p : geo->IteratePlaneIDs()) nearwires.push_back(NearWire(*geo,p,nuVtxCorr[0],nuVtxCorr[1],nuVtxCorr[2]));

    fEventNtupleNu->insert( evtID.data(),
      nutruth.CCNC() == simb::kCC,
      nutruth.Nu().PdgCode(),
      nutruth.Nu().E(),
      nutruth.Lepton().E(),
      nuDirection.data(),
      nuVtx.data(),
      nuVtxCorr.data(),
      nearwires.data(),
      nuT
    );
    LogInfo("HDF5Maker") << "Filling event table"
                         << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                         << ", event " << evtID[2]
                         << "\nis cc? " << (nutruth.CCNC() == simb::kCC)
                         << ", nu energy " << nutruth.Nu().E()
                         << ", lepton energy " << nutruth.Lepton().E()
                         << "\nnu direction x " << nuDirection[0] << ", y "
                         << nuDirection[1] << ", z " << nuDirection[2];
  } // if nu event info

  // Get spacepoints from the event record
  art::Handle<std::vector<recob::SpacePoint>> spListHandle;
  std::vector<art::Ptr<recob::SpacePoint>> splist;
  if (e.getByLabel(fSPLabel, spListHandle))
    art::fill_ptr_vector(splist, spListHandle);

  // Get hits from the event record
  art::Handle<std::vector<recob::Hit>> hitListHandle;
  std::vector<art::Ptr<recob::Hit>> hitlist;
  if (e.getByLabel(fHitLabel, hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  // Get assocations from spacepoints to hits
  art::FindManyP<recob::Hit> fmp(spListHandle, e, fSPLabel);
  std::vector<std::vector<art::Ptr<recob::Hit>>> sp2Hit(splist.size());
  for (size_t spIdx = 0; spIdx < sp2Hit.size(); ++spIdx) {
    sp2Hit[spIdx] = fmp.at(spIdx);
  } // for spacepoint

  // Fill spacepoint table
  for (size_t i = 0; i < splist.size(); ++i) {

    std::array<float, 3> pos {
      (float)splist[i]->XYZ()[0],
      (float)splist[i]->XYZ()[1],
      (float)splist[i]->XYZ()[2]
    };

    std::array<int, 3> hitID { -1, -1, -1 };
    for (size_t j = 0; j < sp2Hit[i].size(); ++j)
      hitID[sp2Hit[i][j]->View()] = sp2Hit[i][j].key();

    fSpacePointNtuple->insert(evtID.data(),
      splist[i]->ID(), pos.data(), hitID.data()
    );

    mf::LogInfo("HDF5Maker") << "Filling spacepoint table"
                             << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                             << ", event " << evtID[2]
                             << "\nspacepoint id " << splist[i]->ID()
                             << "\nposition x " << pos[0] << ", y " << pos[1]
                             << ", z " << pos[2]
                             << "\nhit ids " << hitID[0] << ", " << hitID[1]
                             << ", " << hitID[2];

  } // for spacepoint

  std::set<int> g4id;
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e, clockData);

  std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> hittruth;
  if (fUseMap) {
    hittruth = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> >(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitListHandle, e, fHitTruthLabel));
  }

  // Loop over hits
  for (art::Ptr<recob::Hit> hit : hitlist) {

    // Fill hit table
    geo::WireID wireid = hit->WireID();
    size_t plane = wireid.Plane;
    size_t wire = wireid.Wire;
    double time = hit->PeakTime();
    fHitNtuple->insert(evtID.data(),
      hit.key(), hit->Integral(), hit->RMS(), wireid.TPC,
      plane, wire, time,
      wireid.Plane, wireid.Wire, hit->PeakTime()
    );

    mf::LogInfo("HDF5Maker") << "Filling hit table"
                             << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                             << ", event " << evtID[2]
                             << "\nhit id " << hit.key() << ", integral "
                             << hit->Integral() << ", RMS " << hit->RMS()
                             << ", TPC " << wireid.TPC
                             << "\nglobal plane " << plane << ", global wire "
                             << wire << ", global time " << time
                             << "\nlocal plane " << wireid.Plane
                             << ", local wire " << wireid.Wire
                             << ", local time " << hit->PeakTime();

    // Fill energy deposit table
    if (fUseMap) {
      std::vector<art::Ptr<simb::MCParticle>> particle_vec = hittruth->at(hit.key());
      std::vector<anab::BackTrackerHitMatchingData const *> match_vec = hittruth->data(hit.key());
      //loop over particles
      for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
	g4id.insert(particle_vec[i_p]->TrackId());
	fEnergyDepNtuple->insert(evtID.data(),
		hit.key(), particle_vec[i_p]->TrackId(), match_vec[i_p]->ideFraction
	);
	mf::LogInfo("HDF5Maker") << "Filling energy deposit table"
	        		 << "\nrun " << evtID[0] << ", subrun " << evtID[1]
			         << ", event " << evtID[2]
			         << "\nhit id " << hit.key() << ", g4 id "
			         << particle_vec[i_p]->TrackId() << ", energy fraction "
			         << match_vec[i_p]->ideFraction;
      }
    } else {
      for (const sim::TrackIDE& ide : bt->HitToTrackIDEs(clockData, hit)) {
        if (ide.trackID < 0)
          throw art::Exception(art::errors::LogicError)
            << "Negative track ID (" << ide.trackID << ") found in simulated "
            << "energy deposits! This is usually an indication that you're "
            << "running over simulation from before the larsoft Geant4 "
            << "refactor, which is not supported due to its incomplete MC "
            << "truth record.";
	g4id.insert(ide.trackID);
	fEnergyDepNtuple->insert(evtID.data(),
		hit.key(), ide.trackID, ide.energy);
	mf::LogInfo("HDF5Maker") << "Filling energy deposit table"
                                 << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                                 << ", event " << evtID[2]
			         << "\nhit id " << hit.key() << ", g4 id "
			         << ide.trackID << ", energy "
			         << ide.energy << " MeV";
      } // for energy deposit
    } // if using microboone map method or not
  } // for hit

  art::ServiceHandle<cheat::ParticleInventoryService> pi;
  std::set<int> allIDs = g4id; // Copy original so we can safely modify it

  // Add invisible particles to hierarchy
  for (int id : g4id) {
    const simb::MCParticle* p = pi->TrackIdToParticle_P(abs(id));
    while (p->Mother() != 0) {
      allIDs.insert(abs(p->Mother()));
      p = pi->TrackIdToParticle_P(abs(p->Mother()));
    }
  }

  // Loop over true particles and fill table
  for (int id : allIDs) {
    const simb::MCParticle* p = pi->TrackIdToParticle_P(abs(id));
    if (p==NULL) continue;
    std::array<float, 3> particleStart { (float)p->Vx(), (float)p->Vy(), (float)p->Vz() };
    std::array<float, 3> particleEnd { (float)p->EndX(), (float)p->EndY(), (float)p->EndZ() };
    fParticleNtuple->insert(evtID.data(),
      abs(id), p->PdgCode(), p->Mother(), (float)p->P(),
      particleStart.data(), particleEnd.data(),
      p->Process(), p->EndProcess()
    );
    mf::LogInfo("HDF5Maker") << "Filling particle table"
                             << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                             << ", event " << evtID[2]
                             << "\ng4 id " << abs(id) << ", pdg code "
                             << p->PdgCode() << ", parent " << p->Mother()
                             << ", momentum " << p->P()
                             << "\nparticle start x " << particleStart[0]
                             << ", y " << particleStart[1]
                             << ", z " << particleStart[2]
                             << "\nparticle end x " << particleEnd[0] << ", y "
                             << particleEnd[1] << ", z " << particleEnd[2]
                             << "\nstart process " << p->Process()
                             << ", end process " << p->EndProcess();
  }
} // function HDF5Maker::analyze

void HDF5Maker::beginSubRun(art::SubRun const& sr) {

  struct timeval now;
  gettimeofday(&now, NULL);

  // Open HDF5 output
  std::ostringstream fileName;
  fileName << fOutputName << "_r" << std::setfill('0') << std::setw(5) << sr.run()
    << "_s" << std::setfill('0') << std::setw(5) << sr.subRun() << "_ts" << std::setw(6) << now.tv_usec << ".h5";

  fFile = hep_hpc::hdf5::File(fileName.str(), H5F_ACC_TRUNC);

  if (fEventInfo == "none")
    fEventNtuple = new hep_hpc::hdf5::Ntuple(
      hep_hpc::hdf5::make_ntuple({fFile, "event_table", 1000},
        hep_hpc::hdf5::make_column<int>("event_id", 3)
    ));
  if (fEventInfo == "nu")
    fEventNtupleNu = new hep_hpc::hdf5::Ntuple(
      hep_hpc::hdf5::make_ntuple({fFile, "event_table", 1000},
      hep_hpc::hdf5::make_column<int>("event_id", 3),
      hep_hpc::hdf5::make_scalar_column<int>("is_cc"),
      hep_hpc::hdf5::make_scalar_column<int>("nu_pdg"),
      hep_hpc::hdf5::make_scalar_column<float>("nu_energy"),
      hep_hpc::hdf5::make_scalar_column<float>("lep_energy"),
      hep_hpc::hdf5::make_column<float>("nu_dir", 3),
      hep_hpc::hdf5::make_column<float>("nu_vtx", 3),
      hep_hpc::hdf5::make_column<float>("nu_vtx_corr", 3),
      hep_hpc::hdf5::make_column<int>("nu_vtx_wire_pos", ServiceHandle<geo::Geometry>()->Nviews()),
      hep_hpc::hdf5::make_scalar_column<float>("nu_vtx_wire_time")
    ));

  fSpacePointNtuple = new hep_hpc::hdf5::Ntuple(
    hep_hpc::hdf5::make_ntuple({fFile, "spacepoint_table", 1000},
      hep_hpc::hdf5::make_column<int>("event_id", 3),
      hep_hpc::hdf5::make_scalar_column<int>("spacepoint_id"),
      hep_hpc::hdf5::make_column<float>("position", 3),
      hep_hpc::hdf5::make_column<int>("hit_id", 3)
  ));

  fHitNtuple = new hep_hpc::hdf5::Ntuple(
    hep_hpc::hdf5::make_ntuple({fFile, "hit_table", 1000},
      hep_hpc::hdf5::make_column<int>("event_id", 3),
      hep_hpc::hdf5::make_scalar_column<int>("hit_id"),
      hep_hpc::hdf5::make_scalar_column<float>("integral"),
      hep_hpc::hdf5::make_scalar_column<float>("rms"),
      hep_hpc::hdf5::make_scalar_column<int>("tpc"),
      hep_hpc::hdf5::make_scalar_column<int>("global_plane"),
      hep_hpc::hdf5::make_scalar_column<float>("global_wire"),
      hep_hpc::hdf5::make_scalar_column<float>("global_time"),
      hep_hpc::hdf5::make_scalar_column<int>("local_plane"),
      hep_hpc::hdf5::make_scalar_column<float>("local_wire"),
      hep_hpc::hdf5::make_scalar_column<float>("local_time")
  ));

  fParticleNtuple = new hep_hpc::hdf5::Ntuple(
    hep_hpc::hdf5::make_ntuple({fFile, "particle_table", 1000},
      hep_hpc::hdf5::make_column<int>("event_id", 3),
      hep_hpc::hdf5::make_scalar_column<int>("g4_id"),
      hep_hpc::hdf5::make_scalar_column<int>("type"),
      hep_hpc::hdf5::make_scalar_column<int>("parent_id"),
      hep_hpc::hdf5::make_scalar_column<float>("momentum"),
      hep_hpc::hdf5::make_column<float>("start_position", 3),
      hep_hpc::hdf5::make_column<float>("end_position", 3),
      hep_hpc::hdf5::make_scalar_column<std::string>("start_process"),
      hep_hpc::hdf5::make_scalar_column<std::string>("end_process")
  ));

  fEnergyDepNtuple = new hep_hpc::hdf5::Ntuple(
    hep_hpc::hdf5::make_ntuple({fFile, "edep_table", 1000},
      hep_hpc::hdf5::make_column<int>("event_id", 3),
      hep_hpc::hdf5::make_scalar_column<int>("hit_id"),
      hep_hpc::hdf5::make_scalar_column<int>("g4_id"),
      hep_hpc::hdf5::make_scalar_column<float>("energy")
  ));
}

void HDF5Maker::endSubRun(art::SubRun const& sr) {
  if (fEventInfo == "none") delete fEventNtuple;
  if (fEventInfo == "nu") delete fEventNtupleNu;
  delete fSpacePointNtuple;
  delete fHitNtuple;
  delete fParticleNtuple;
  delete fEnergyDepNtuple;
  fFile.close();
}


// apply the mapping of XYZ true -> XYZ position as it would be recosntructed.
// takes into account SCE, trigger time offset, and wirecell-pandora offset.
// to be applied to truth xyz in order to compare to reconstructed variables
// e.g. used for resolution plots
void HDF5Maker::True2RecoMappingXYZ(float& t, float& x, float& y, float& z)
{
  ApplySCEMappingXYZ(x, y, z);
  auto const &detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  auto const &detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  double g4Ticks = detClocks->TPCG4Time2Tick(t) + detProperties->GetXTicksOffset(0, 0, 0) - detProperties->TriggerOffset();
  float _xtimeoffset = detProperties->ConvertTicksToX(g4Ticks, 0, 0, 0);
  x += _xtimeoffset;
  x += fXOffset;
}

// apply the mapping of XYZ true -> XYZ position after SCE-induced shift.
// to be applied to truth xyz in order to compare to reconstructed variables
// e.g. used for resolution plots
void HDF5Maker::ApplySCEMappingXYZ(float& x, float& y, float& z)
{
  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  if (SCE->EnableSimSpatialSCE() == true)
    {
      auto offset = SCE->GetPosOffsets(geo::Point_t(x, y, z));
      x -= offset.X();
      y += offset.Y();
      z += offset.Z();
    }
}

void HDF5Maker::AddDaughters(const std::map<unsigned int, unsigned int>& _pfpmap,
			     const ProxyPfpElem_t &pfp_pxy,
			     const ProxyPfpColl_t &pfp_pxy_col,
			     std::vector<ProxyPfpElem_t> &slice_v)
{

  auto daughters = pfp_pxy->Daughters();

  slice_v.push_back(pfp_pxy);

  for (auto const &daughterid : daughters)
  {

    if (_pfpmap.find(daughterid) == _pfpmap.end())
    {
      std::cout << "Did not find DAUGHTERID in map! error" << std::endl;
      continue;
    }

    auto pfp_pxy2 = pfp_pxy_col.begin();
    for (size_t j = 0; j < _pfpmap.at(daughterid); ++j) ++pfp_pxy2;

    AddDaughters(_pfpmap,*pfp_pxy2, pfp_pxy_col, slice_v);

  } // for all daughters

  return;
} // AddDaughters

float HDF5Maker::GetMetaData(const ProxyPfpElem_t &pfp_pxy, string metaDataName)
{
  const auto &pfParticleMetadataList = pfp_pxy.get<larpandoraobj::PFParticleMetadata>();
  if (pfParticleMetadataList.size() == 0) return -999.;
  for (unsigned int j = 0; j < pfParticleMetadataList.size(); ++j)
    {
      const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
      auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
      if (!pfParticlePropertiesMap.empty())
	{
	  for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it)
	    {
	      if (it->first == metaDataName)
		return it->second;
	    } // for map elements
	} // if pfp metadata map not empty
    } // for list
  return -999.;
}

int HDF5Maker::NearWire(const geo::Geometry& geo, const geo::PlaneID &ip, const float x, const float y, const float z)
{
  geo::PlaneGeo const& plane = geo.Plane(ip);
  geo::WireID wireID;
  try {
    wireID = plane.NearestWireID(geo::Point_t(x,y,z));
  }
  catch (geo::InvalidWireError const& e) {
    if (!e.hasSuggestedWire()) throw;
    wireID = plane.ClosestWireID(e.suggestedWireID());
  }
  return wireID.Wire;
}

std::vector<HDF5Maker::BtPart> HDF5Maker::initBacktrackingParticleVec(const std::vector<sim::MCShower> &inputMCShower,
								      const std::vector<sim::MCTrack> &inputMCTrack,
								      const std::vector<recob::Hit> &inputHits,
								      const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart,
								      int nhitcut) {

  std::vector<BtPart> btparts_v;
  for (auto mcs : inputMCShower)
  {
    int category = -1;
    if (mcs.Process() == "primary" && mcs.PdgCode()==22) category=ParticleCategory::Photon;//photon
    else if (mcs.MotherPdgCode() == 111 && mcs.Process() == "Decay" && mcs.MotherProcess() == "primary") category=ParticleCategory::Photon;//photons from pi0
    else if (mcs.Process() == "primary" && std::abs(mcs.PdgCode())==11) category=ParticleCategory::Electron;//electron
    else if (std::abs(mcs.MotherPdgCode()) == 13 && mcs.Process() == "Decay") category=ParticleCategory::Michel;//michel electron
    else if (std::abs(mcs.MotherPdgCode()) == 13 && mcs.Process() == "muIoni" && mcs.Start().Momentum().P()>10) category=ParticleCategory::Delta;//delta
    else if (std::abs(mcs.MotherPdgCode()) == 211 && mcs.Process() == "hIoni" && mcs.Start().Momentum().P()>10) category=ParticleCategory::Delta;//delta
    else continue;

    sim::MCStep mc_step_shower_start = mcs.DetProfile();
    btparts_v.push_back(BtPart(mcs.PdgCode(), category, mcs.Start().Momentum().Px() * 0.001, mcs.Start().Momentum().Py() * 0.001,
			       mcs.Start().Momentum().Pz() * 0.001, mcs.Start().Momentum().E() * 0.001, mcs.DaughterTrackID(),
			       mc_step_shower_start.X(), mc_step_shower_start.Y(), mc_step_shower_start.Z(), mc_step_shower_start.T()));

  }
  for (auto mct : inputMCTrack)
  {

    int category = -1;
    if (std::abs(mct.PdgCode())==13) category=ParticleCategory::Muon;  //muon
    else if (std::abs(mct.PdgCode())==2212) category=ParticleCategory::Proton;//proton
    else if (std::abs(mct.PdgCode())==211) category=ParticleCategory::Pion; //charged pion
    else if (std::abs(mct.PdgCode())==321) category=ParticleCategory::Kaon; //charged kaon
    else continue;

    sim::MCStep mc_step_track_start = mct.Start();
    btparts_v.push_back(BtPart(mct.PdgCode(), category, mct.Start().Momentum().Px() * 0.001, mct.Start().Momentum().Py() * 0.001,
			       mct.Start().Momentum().Pz() * 0.001, mct.Start().Momentum().E() * 0.001, mct.TrackID(),
			       mc_step_track_start.X(), mc_step_track_start.Y(), mc_step_track_start.Z(), mc_step_track_start.T()));

  }
  // Now let's fill the nhits member using all input hits
  for (unsigned int ih = 0; ih < inputHits.size(); ih++)
  {
    auto assmcp = assocMCPart->at(ih);
    auto assmdt = assocMCPart->data(ih);
    for (unsigned int ia = 0; ia < assmcp.size(); ++ia)
    {
      auto mcp = assmcp[ia];
      auto amd = assmdt[ia];
      if (amd->isMaxIDE != 1) continue;

      for (auto &btp : btparts_v)
      {
        if (std::find(btp.tids.begin(), btp.tids.end(), mcp->TrackId()) != btp.tids.end())
        {
          btp.nhits++;
        }
      }
    }
  }
  std::vector<BtPart> btparts_v_final;
  for (size_t i=0;i<btparts_v.size();++i) {
    if (btparts_v[i].nhits>nhitcut) {
      btparts_v_final.push_back(btparts_v[i]);
    }
  }
  return btparts_v_final;
}

DEFINE_ART_MODULE(HDF5Maker)
