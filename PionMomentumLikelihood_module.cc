#include "PionMomentumLikelihood.h"
//_________________________________________________________________________________________
ubpiontraj::PionMomentumLikelihood::PionMomentumLikelihood(fhicl::ParameterSet const& p) : 
   EDAnalyzer{p},
   m_SimParticleLabel(p.get<std::string>("sim_particle_label", "largeant")),
   m_RecoParticleLabel(p.get<std::string>("reco_particle_label", "pandora")),
   m_RecoTrackLabel(p.get<std::string>("reco_track_label", "pandora")),
   m_RecoHitLabel(p.get<std::string>("reco_hit_label", "pandora")),
   m_RecoTrackHitAssocLabel(p.get<std::string>("reco_trackhitassoc_label", "pandora")),
   m_RecoHitSimParticleLabel(p.get<std::string>("reco_hitsimparticle_label", "gaushitTruthMatch")),
   m_RecoCaloLabel(p.get<std::string>("reco_calo_label", "pandoracali")),
   m_Debug(p.get<bool>("debug", false))
{
   if(m_Debug){
      std::cout << ">>> [PionMomentumLikelihood] Constructing class..." << std::endl;
   }
}
//_________________________________________________________________________________________
void ubpiontraj::PionMomentumLikelihood::beginJob()
{
   if(m_Debug){
      std::cout << ">>> [PionMomentumLikelihood] Initialising trees and branches..." << std::endl;
   }
   art::ServiceHandle<art::TFileService> tfs;

   m_SpTree = tfs->make<TTree>("SpTree", "Space Point Tree");
   m_SpTree->Branch("sp_p", &m_sp_p);
   m_SpTree->Branch("sp_n", &m_sp_n);
   m_SpTree->Branch("sp_x", &m_sp_x);
   m_SpTree->Branch("sp_y", &m_sp_y);
   m_SpTree->Branch("sp_z", &m_sp_z);
   m_SpTree->Branch("sp_e", &m_sp_e);
}
//_________________________________________________________________________________________
void ubpiontraj::PionMomentumLikelihood::endJob()
{
    if(m_Debug){
        std::cout << ">>> [PionMomentumLikelihood] Calling end job..." << std::endl;
    }
}
//_________________________________________________________________________________________
void ubpiontraj::PionMomentumLikelihood::analyze(art::Event const& event)
{
   if(!event.getByLabel(m_SimParticleLabel, m_SimParticleHandle)){
      if(m_Debug){
         std::cerr << ">>> [PionMomentumLikelihood] Error! Could not get simulation handle for label: " << m_SimParticleLabel << std::endl;
      }
      return;
   }

   art::fill_ptr_vector(m_SimParticles, m_SimParticleHandle);

   for(const art::Ptr<simb::MCParticle>& initialparticle : m_SimParticles){
      m_SimParticleMap.insert(std::make_pair(initialparticle->TrackId(), initialparticle));
   }
   if(!event.getByLabel(m_RecoParticleLabel, m_RecoParticleHandle)) 
      throw cet::exception("PionReconstructionAnalyser") << "No PFParticle Data Products Found!" << std::endl;

   if(!event.getByLabel(m_RecoTrackLabel, m_RecoTrackHandle))
      throw cet::exception("PionReconstrucionAnalyser") << "No Track Data Products Found!" << std::endl;

   if(!event.getByLabel(m_RecoHitLabel, m_RecoHitHandle)) 
      throw cet::exception("PionReconstructionAnalyser") << "No Track Data Products Found!" << std::endl;

   art::fill_ptr_vector(m_RecoParticles, m_RecoParticleHandle);
   art::fill_ptr_vector(m_RecoTracks, m_RecoTrackHandle);
   art::fill_ptr_vector(m_RecoHits, m_RecoHitHandle);

   m_RecoParticleTrackAssoc = new art::FindManyP<recob::Track>(m_RecoParticles, event, m_RecoParticleLabel);     
   m_RecoTrackHitAssoc = new art::FindManyP<recob::Hit>(m_RecoTracks, event, m_RecoTrackHitAssocLabel);

   m_RecoHitSimParticleAssoc = new art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData>(m_RecoHitHandle, event, m_RecoHitSimParticleLabel);
   m_RecoTrackCaloAssoc = new art::FindManyP<anab::Calorimetry>(m_RecoTracks, event, m_RecoCaloLabel);

   if(m_Debug){
      std::cout << "///----------------------------------------" << std::endl;
      std::cout << ">>> [PionMomentumLikelihood] Analysing event... " << std::endl;
   }

   for(art::Ptr<simb::MCParticle>& simParticle : m_SimParticles){
      if(abs(simParticle->PdgCode()) == 211 && (simParticle->Process() == "decay" || simParticle->Process() == "primary")){
         m_sp_p = simParticle->P();

         m_sp_n = 0;
         m_sp_x.clear();
         m_sp_y.clear(); 
         m_sp_z.clear();
         m_sp_e.clear();

         for(const art::Ptr<recob::PFParticle>& recoParticle : m_RecoParticles){
            std::vector<art::Ptr<recob::Track>> recoTracks = m_RecoParticleTrackAssoc->at(recoParticle.key());
            if(recoTracks.size() > 1) continue;

            art::Ptr<recob::Track> trk = recoTracks.at(0);
            std::vector<art::Ptr<recob::Hit>> hitVec = m_RecoTrackHitAssoc->at(trk.key());
            std::vector<art::Ptr<anab::Calorimetry>> caloVec = m_RecoTrackCaloAssoc->at(trk.key());

            std::unordered_map<int, double> trackMap;
            int maxHits = -1;
            
            simb::MCParticle const* matchedSimParticle = nullptr;
            std::vector<simb::MCParticle const*> depositingSimParticles;
            std::vector<anab::BackTrackerHitMatchingData const*> matchingHits;

            for(size_t h = 0; h < hitVec.size(); ++h){
               depositingSimParticles.clear();
               matchingHits.clear();
               m_RecoHitSimParticleAssoc->get(hitVec[h].key(), depositingSimParticles, matchingHits);

               for(size_t p = 0; p < depositingSimParticles.size(); ++p){
                  trackMap[depositingSimParticles[p]->TrackId()]++; 

                  if(trackMap[depositingSimParticles[p]->TrackId()] > maxHits){
                     maxHits = trackMap[depositingSimParticles[p]->TrackId()];
                     matchedSimParticle = depositingSimParticles[p];
                  }
               }
            }

            if(matchedSimParticle->TrackId() != simParticle->TrackId()) continue;

            for(size_t p = 0; p < caloVec.size(); p++){
               if(caloVec.at(p)->PlaneID().Plane != 2) continue;

               art::Ptr<anab::Calorimetry> calo = caloVec.at(p);
               if(calo->XYZ().size() < 2) continue; 

               for(size_t sp = 0; sp < calo->XYZ().size(); sp++){
                  const auto& coord = calo->XYZ().at(sp);
                  double dEdx = calo->dEdx().at(sp);

                  m_sp_n = calo->XYZ().size();
                  m_sp_x.push_back(coord.X());
                  m_sp_y.push_back(coord.Y());
                  m_sp_z.push_back(coord.Z());
                  m_sp_e.push_back(dEdx);
               }
            }         
         }

         if(m_Debug){
            std::cout << ">>> [PionMomentumLikelihood] Filling tree for a pion..." << std::endl;
         }

         m_SpTree->Fill();
      }
   }
}
//_________________________________________________________________________________________
void ubpiontraj::PionMomentumLikelihood::beginSubRun(const art::SubRun& sr)
{
   if(m_Debug){
      std::cout << ">>> [PionMomentumLikelihood] Beginning subrun..." << std::endl;
   }
}
//_________________________________________________________________________________________
void ubpiontraj::PionMomentumLikelihood::endSubRun(const art::SubRun& sr)
{
   if(m_Debug){
      std::cout << ">>> [PionMomentumLikelihood] Ending subrun..." << std::endl;
   }
}
//_________________________________________________________________________________________
DEFINE_ART_MODULE(ubpiontraj::PionMomentumLikelihood)