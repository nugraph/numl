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

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "hep_hpc/hdf5/make_ntuple.hpp"

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
			hep_hpc::hdf5::Column<int, 1>, // nu pdg
                        hep_hpc::hdf5::Column<float, 1>,  // nu energy
                        hep_hpc::hdf5::Column<float, 1>,  // lep energy
                        hep_hpc::hdf5::Column<float, 1>   // nu dir (x, y, z)
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
                        hep_hpc::hdf5::Column<float, 1>,  // deposited energy [ MeV ]
                        hep_hpc::hdf5::Column<float, 1>,  // x position
                        hep_hpc::hdf5::Column<float, 1>,  // y position
                        hep_hpc::hdf5::Column<float, 1>   // z position
  >* fEnergyDepNtuple; ///< energy deposition ntuple
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

    std::array<float, 3> nuMomentum {
      (float)nutruth.Nu().Momentum().Vect().Unit().X(),
      (float)nutruth.Nu().Momentum().Vect().Unit().Y(),
      (float)nutruth.Nu().Momentum().Vect().Unit().Z()
    };

    fEventNtupleNu->insert( evtID.data(),
      nutruth.CCNC() == simb::kCC,
      nutruth.Nu().PdgCode(),
      nutruth.Nu().E(),
      nutruth.Lepton().E(),
      nuMomentum.data()
    );
    mf::LogInfo("HDF5Maker") << "Filling event table"
                             << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                             << ", event " << evtID[2]
                             << "\nis cc? " << (nutruth.CCNC() == simb::kCC)
                             << ", nu energy " << nutruth.Nu().E()
                             << ", lepton energy " << nutruth.Lepton().E()
                             << "\nnu momentum x " << nuMomentum[0] << ", y "
                             << nuMomentum[1] << ", z " << nuMomentum[2];
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
  //   if (fUseMap) {
  //     std::vector<art::Ptr<simb::MCParticle>> particle_vec = hittruth->at(hit.key());
  //     std::vector<anab::BackTrackerHitMatchingData const *> match_vec = hittruth->data(hit.key());
  //     //loop over particles
  //     for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
	// g4id.insert(particle_vec[i_p]->TrackId());
	// fEnergyDepNtuple->insert(evtID.data(),
	// 	hit.key(), particle_vec[i_p]->TrackId(), match_vec[i_p]->ideFraction
	// );
	// mf::LogInfo("HDF5Maker") << "Filling energy deposit table"
	//         		 << "\nrun " << evtID[0] << ", subrun " << evtID[1]
	// 		         << ", event " << evtID[2]
	// 		         << "\nhit id " << hit.key() << ", g4 id "
	// 		         << particle_vec[i_p]->TrackId() << ", energy fraction "
	// 		         << match_vec[i_p]->ideFraction;
  //     }
  //   } else {
      for (const sim::IDE* ide : bt->HitToSimIDEs_Ps(clockData, hit)) {
        if (ide->trackID < 0)
          throw art::Exception(art::errors::LogicError)
            << "Negative track ID (" << ide->trackID << ") found in simulated "
            << "energy deposits! This is usually an indication that you're "
            << "running over simulation from before the larsoft Geant4 "
            << "refactor, which is not supported due to its incomplete MC "
            << "truth record.";
	g4id.insert(ide->trackID);
	fEnergyDepNtuple->insert(evtID.data(),
		hit.key(), ide->trackID, ide->energy, ide->x, ide->y, ide->z);
	mf::LogInfo("HDF5Maker") << "Filling energy deposit table"
                                 << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                                 << ", event " << evtID[2]
			         << "\nhit id " << hit.key() << ", g4 id "
			         << ide->trackID << ", energy "
			         << ide->energy << " MeV, position ("
               << ide->x << ", " << ide->y << ", " << ide->z << ")";
      } // for energy deposit
    // } // if using microboone map method or not
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
        hep_hpc::hdf5::make_column<float>("nu_dir", 3)
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
      hep_hpc::hdf5::make_scalar_column<float>("energy"),
      hep_hpc::hdf5::make_scalar_column<float>("x_position"),
      hep_hpc::hdf5::make_scalar_column<float>("y_position"),
      hep_hpc::hdf5::make_scalar_column<float>("z_position")
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

DEFINE_ART_MODULE(HDF5Maker)
