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

      int count;
      
      TTree* m_Tree;

      double m_px, m_py, m_pz, m_purity, m_completeness;
      int m_pdg;
      std::vector<double> m_x, m_y, m_z, m_e;
};

#endif // PIONMOMENTUMLIKELIHOOD_H