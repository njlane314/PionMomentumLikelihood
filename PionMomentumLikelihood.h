#ifndef PIONMOMENTUMLIKELIHOOD_H
#define PIONMOMENTUMLIKELIHOOD_H

#include <vector>
#include <string>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"		

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

namespace ubpiontraj 
{
   class PionMomentumLikelihood;
}

class ubpiontraj::PionMomentumLikelihood : public art::EDAnalyzer 
{
   public:
      explicit PionMomentumLikelihood(fhicl::ParameterSet const& p);

      PionMomentumLikelihood(PionMomentumLikelihood const&) = delete;
      PionMomentumLikelihood(PionMomentumLikelihood&&) = delete;
      PionMomentumLikelihood& operator=(PionMomentumLikelihood const&) = delete;
      PionMomentumLikelihood& operator=(PionMomentumLikelihood&&) = delete;

      void analyze(art::Event const& event) override;

      void beginJob() override;
      void endJob() override;

      void beginSubRun(const art::SubRun& sr);
      void endSubRun(const art::SubRun& sr);

   private:
      std::string m_SimParticleLabel;
      
      std::string m_RecoParticleLabel;
      std::string m_RecoTrackLabel;
      std::string m_RecoHitLabel;
      
      std::string m_RecoTrackHitAssocLabel;
      std::string m_RecoHitSimParticleLabel;

      std::string m_RecoCaloLabel;
      bool m_Debug;

      art::Handle<std::vector<simb::MCParticle>> m_SimParticleHandle;

      std::vector<art::Ptr<simb::MCParticle>> m_SimParticles;
      std::map<int, art::Ptr<simb::MCParticle>> m_SimParticleMap;

      art::Handle<std::vector<recob::PFParticle>> m_RecoParticleHandle;
      std::vector<art::Ptr<recob::PFParticle>> m_RecoParticles;

      art::Handle<std::vector<recob::Track>> m_RecoTrackHandle;
      std::vector<art::Ptr<recob::Track>> m_RecoTracks;

      art::Handle<std::vector<recob::Hit>> m_RecoHitHandle;
      std::vector<art::Ptr<recob::Hit>> m_RecoHits;

      art::FindManyP<recob::Track>* m_RecoParticleTrackAssoc;
      art::FindManyP<recob::Hit>* m_RecoTrackHitAssoc;

      art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData>* m_RecoHitSimParticleAssoc;

      art::FindManyP<anab::Calorimetry>* m_RecoTrackCaloAssoc;
      
      TTree* m_SimTree;
      TTree* m_SpTree;

      double m_sim_p, m_sim_w;
      std::vector<double> m_sp_x, m_sp_y, m_sp_z, m_sp_e;
};

#endif // PIONMOMENTUMLIKELIHOOD_H