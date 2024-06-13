// lar -c run_PionMomentumLikelihood.fcl -s Data/PhysicsRun-2016_7_29_22_26_32-0007003-01286_20160805T195522_ext_unbiased_20160806T215518_merged_20231124T135742_simmxd_detsim_mix_r1a_r1b_20231124T175_a4687571-abb8-4bf6-8fd9-38e59ea36114.root 
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

   m_Tree = tfs->make<TTree>("Tree", "Tree");
   m_Tree->Branch("px", &m_px);
   m_Tree->Branch("py", &m_py);
   m_Tree->Branch("pz", &m_pz);

   m_Tree->Branch("x", &m_x);
   m_Tree->Branch("y", &m_y);
   m_Tree->Branch("z", &m_z);
   m_Tree->Branch("e", &m_e);
   m_Tree->Branch("pdg", &m_pdg);

   m_Tree->Branch("purity", &m_purity);
   m_Tree->Branch("completeness", &m_completeness);

   count = 0;
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
   art::Handle<std::vector<simb::MCParticle>> simParticleHandle;
   art::Handle<std::vector<recob::PFParticle>> recoParticleHandle;
   art::Handle<std::vector<recob::Track>> recoTrackHandle;
   art::Handle<std::vector<recob::Hit>> recoHitHandle;

   std::vector<art::Ptr<simb::MCParticle>> simParticles;
   std::vector<art::Ptr<recob::PFParticle>> recoParticles;
   std::vector<art::Ptr<recob::Track>> recoTracks;
   std::vector<art::Ptr<recob::Hit>> recoHits;

   art::FindManyP<recob::Track>* m_RecoParticleTrackAssoc;
   art::FindManyP<recob::Hit>* m_RecoTrackHitAssoc;

   art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData>* m_RecoHitSimParticleAssoc;
   art::FindManyP<anab::Calorimetry>* m_RecoTrackCaloAssoc;

   if (!event.getByLabel(m_SimParticleLabel, simParticleHandle)) {
      if (m_Debug) {
         std::cerr << ">>> [PionMomentumLikelihood] Error! Could not get simulation handle for label: " << m_SimParticleLabel << std::endl;
      }
      return;
   }

   if (!event.getByLabel(m_RecoParticleLabel, recoParticleHandle)) 
      throw cet::exception("PionReconstructionAnalyser") << "No PFParticle Data Products Found!" << std::endl;

   if (!event.getByLabel(m_RecoTrackLabel, recoTrackHandle))
      throw cet::exception("PionReconstrucionAnalyser") << "No Track Data Products Found!" << std::endl;

   if (!event.getByLabel(m_RecoHitLabel, recoHitHandle)) 
      throw cet::exception("PionReconstructionAnalyser") << "No Hit Data Products Found!" << std::endl;

   art::fill_ptr_vector(simParticles, simParticleHandle);   
   art::fill_ptr_vector(recoParticles, recoParticleHandle);
   art::fill_ptr_vector(recoTracks, recoTrackHandle);
   art::fill_ptr_vector(recoHits, recoHitHandle);

   m_RecoParticleTrackAssoc = new art::FindManyP<recob::Track>(recoParticles, event, m_RecoParticleLabel);     
   m_RecoTrackHitAssoc = new art::FindManyP<recob::Hit>(recoTracks, event, m_RecoTrackHitAssocLabel);

   m_RecoHitSimParticleAssoc = new art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData>(recoHitHandle, event, m_RecoHitSimParticleLabel);
   m_RecoTrackCaloAssoc = new art::FindManyP<anab::Calorimetry>(recoTracks, event, m_RecoCaloLabel);

   if (m_Debug) {
      std::cout << "///----------------------------------------" << std::endl;
      std::cout << ">>> [PionMomentumLikelihood] Analysing event..." << std::endl;
   }

   for (art::Ptr<simb::MCParticle>& sim : simParticles) {
      if (sim->Mother() != 0)  continue;
      if (abs(sim->PdgCode()) == 211 && (sim->Process() == "decay" || sim->Process() == "primary")) {
         if (m_Debug) {
            std::cout << ">>> [PionMomentumLikelihood] Found pion!" << std::endl;
         }
         
         m_px = 0.; 
         m_py = 0.; 
         m_pz = 0.;

         m_x.clear();
         m_y.clear();
         m_z.clear();
         m_e.clear();

         m_purity = 0.;
         m_completeness = 0.;

         bool foundpion = false;
         for (const art::Ptr<recob::PFParticle>& reco : recoParticles) {
            std::vector<art::Ptr<recob::Track>> recoTracks = m_RecoParticleTrackAssoc->at(reco.key());
            if (recoTracks.size() != 1) continue;
   
            for (const art::Ptr<recob::Track>& trk : recoTracks) {
               std::vector<art::Ptr<recob::Hit>> hitVec = m_RecoTrackHitAssoc->at(trk.key());
               std::vector<art::Ptr<anab::Calorimetry>> caloVec = m_RecoTrackCaloAssoc->at(trk.key());

               std::unordered_map<int, int> trackMap;
               int maxHits = -1;
               
               std::vector<simb::MCParticle const*> depositingSimParticles;
               std::vector<anab::BackTrackerHitMatchingData const*> matchingHits;
               simb::MCParticle const* matchedSimParticle = nullptr;

               int totalSpacePoints = 0;
               int matchedSpacePoints = 0;

               for (size_t h = 0; h < hitVec.size(); ++h) {
                  m_RecoHitSimParticleAssoc->get(hitVec[h].key(), depositingSimParticles, matchingHits);

                  for (size_t p = 0; p < depositingSimParticles.size(); ++p) {
                     trackMap[depositingSimParticles[p]->TrackId()]++; 

                     if (trackMap[depositingSimParticles[p]->TrackId()] > maxHits) {
                        maxHits = trackMap[depositingSimParticles[p]->TrackId()];
                        matchedSimParticle = depositingSimParticles[p];
                     }

                     if (depositingSimParticles[p]->TrackId() == sim->TrackId()) {
                        matchedSpacePoints++;
                     }
                  }
                  totalSpacePoints++;
               }

               if (matchedSimParticle == nullptr) continue; 
               if (matchedSimParticle->TrackId() != sim->TrackId()) continue;
               foundpion = true;
               count += 1; 
               if (m_Debug) {
                  std::cout << ">>> [PionMomentumLikelihood] Found a matching simulation pion!" << std::endl;
               }

               for (size_t p = 0; p < caloVec.size(); p++) {
                  if (caloVec.at(p)->PlaneID().Plane != 2) continue;

                  art::Ptr<anab::Calorimetry> calo = caloVec.at(p);
                  if (calo->XYZ().size() < 2) continue; 

                  std::vector<TVector3> coords;
                  for (const auto& pos : calo->XYZ()) {
                      coords.emplace_back(pos.X(), pos.Y(), pos.Z());
                  }

                  TVector3 vec = coords.back() - coords.front();
                  TVector3 momentum(sim->Px(), sim->Py(), sim->Pz());
                  if (vec.Dot(momentum) < 0) {
                     std::reverse(coords.begin(), coords.end());
                  }

                  for (const auto& coord : coords) {
                     m_x.push_back(coord.X());
                     m_y.push_back(coord.Y());
                     m_z.push_back(coord.Z());
                  }
               }

               if (totalSpacePoints > 0) {
                  m_purity = static_cast<double>(matchedSpacePoints) / totalSpacePoints;
                  std::cout << m_purity << std::endl;
               }

               double totalTrueSpacePoints = sim->NumberTrajectoryPoints();
               if (totalTrueSpacePoints > 0) {
                  m_completeness = static_cast<double>(matchedSpacePoints) / totalTrueSpacePoints;
                  std::cout << m_completeness << std::endl;
               }

               if (foundpion) break;
            }      

            if (foundpion) break;
         }

         if (foundpion) {
            if (m_Debug) {
               std::cout << ">>> [PionMomentumLikelihood] Filling tree for a pion..." << std::endl;
            }

            m_px = sim->Px();
            m_py = sim->Py(); 
            m_pz = sim->Pz();

            m_pdg = sim->PdgCode();
            
            m_Tree->Fill();
         }
      }
   }

   if (m_Debug) {
      std::cout << ">>> [PionMomentumLikelihood] Finished analysing event!" << std::endl;
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
      std::cout << ">>> [PionMomentumLikelihood] Pions found: " << count << std::endl;
   }
}
//_________________________________________________________________________________________
DEFINE_ART_MODULE(ubpiontraj::PionMomentumLikelihood)