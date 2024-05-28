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

#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"

#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "Objects/Include/Trajectory.h"
#include "Module/Include/PionSimulationAnalyser.h"
#include "Module/Include/PionReconstructionAnalyser.h"

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

      void analyze(art::Event const& e) override;

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
      
      TTree* m_SpTree;
   
      int m_sp_n;
      double m_sp_p;
      std::vector<double> m_sp_x, m_sp_y, m_sp_z, m_sp_e;
};

#endif // PIONMOMENTUMLIKELIHOOD_H