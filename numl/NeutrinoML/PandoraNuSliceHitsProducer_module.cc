////////////////////////////////////////////////////////////////////////
// Class:       PandoraNuSliceHitsProducer
// Plugin Type: producer (art v3_06_03)
// File:        PandoraNuSliceHitsProducer_module.cc
//
// Generated at Tue May 25 10:39:19 2021 by Giuseppe Cerati using cetskelgen
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include <memory>

#include "canvas/Persistency/Common/Assns.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

class PandoraNuSliceHitsProducer;

using HitParticleAssociations = art::Assns<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData>;

class PandoraNuSliceHitsProducer : public art::EDProducer {
public:
  explicit PandoraNuSliceHitsProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PandoraNuSliceHitsProducer(PandoraNuSliceHitsProducer const&) = delete;
  PandoraNuSliceHitsProducer(PandoraNuSliceHitsProducer&&) = delete;
  PandoraNuSliceHitsProducer& operator=(PandoraNuSliceHitsProducer const&) = delete;
  PandoraNuSliceHitsProducer& operator=(PandoraNuSliceHitsProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  std::string fPfpLabel;
  std::string fSliceLabel;
  std::string fHitLabel;
  std::string fHitTruthLabel;
};


PandoraNuSliceHitsProducer::PandoraNuSliceHitsProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fPfpLabel(p.get<std::string>("PfpLabel","pandora")),
  fSliceLabel(p.get<std::string>("SliceLabel","pandora")),
  fHitLabel(p.get<std::string>("HitLabel","gaushit")),
  fHitTruthLabel(p.get<std::string>("HitTruthLabel","gaushitTruthMatch"))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces<std::vector<recob::Hit>>();
  produces<HitParticleAssociations>();

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void PandoraNuSliceHitsProducer::produce(art::Event& e)
{
  // Implementation of required member function here.

  auto outputHits = std::make_unique<std::vector<recob::Hit> >();
  auto outputHitPartAssns = std::make_unique<HitParticleAssociations>();
  art::PtrMaker<recob::Hit> hitPtrMaker(e);

  art::ValidHandle<std::vector<recob::PFParticle> > inputPfp = e.getValidHandle<std::vector<recob::PFParticle> >(fPfpLabel);
  auto assocPfpSlice = std::unique_ptr<art::FindManyP<recob::Slice> >(new art::FindManyP<recob::Slice>(inputPfp, e, fPfpLabel));

  art::ValidHandle<std::vector<recob::Slice> > inputSlice = e.getValidHandle<std::vector<recob::Slice> >(fSliceLabel);
  auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit> >(new art::FindManyP<recob::Hit>(inputSlice, e, fSliceLabel));

  art::Handle< std::vector< recob::Hit > > hitListHandle;
  e.getByLabel(fHitLabel, hitListHandle);
  std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> hittruth = 
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> >(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitListHandle, e, fHitTruthLabel));

  for (size_t ipfp = 0; ipfp<inputPfp->size(); ipfp++) {

    art::Ptr<recob::PFParticle> pfp(inputPfp,ipfp);
    if (pfp->IsPrimary()==false) continue;
    auto PDG = fabs(pfp->PdgCode());
    if (PDG != 12 && PDG != 14) continue;

    auto assocSlice = assocPfpSlice->at(pfp.key());
    auto sliceHits = assocSliceHit->at(assocSlice[0].key());

    for (size_t ihit = 0; ihit<sliceHits.size(); ++ihit) {
      auto hit = sliceHits.at(ihit);
      outputHits->emplace_back(*hit);

      std::vector<art::Ptr<simb::MCParticle>> particle_vec = hittruth->at(hit.key());
      std::vector<anab::BackTrackerHitMatchingData const *> match_vec = hittruth->data(hit.key());
      const art::Ptr<recob::Hit> ahp = hitPtrMaker(outputHits->size() - 1);
      for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
	outputHitPartAssns->addSingle(particle_vec[i_p],ahp,*match_vec[i_p]);
      }
    }
  }

  e.put(std::move(outputHits));
  e.put(std::move(outputHitPartAssns));

}

DEFINE_ART_MODULE(PandoraNuSliceHitsProducer)
