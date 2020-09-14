//
// Created by mikhail on 6/16/20.
//

#include "reco_acceptance.h"

namespace AnalysisTree {
void RecoAcceptance::Init(std::map<std::string, void *> &branch_map) {
  reco_header_ = static_cast<EventHeader *>(branch_map.at("event_header"));
  sim_header_ = static_cast<EventHeader *>(branch_map.at("sim_header"));
  sim_tracks_ = static_cast<Particles *>(branch_map.at("sim_tracks"));
  reco_tracks_ = static_cast<Particles *>(branch_map.at("mdc_vtx_tracks"));
  meta_hits_ = static_cast<HitDetector *>(branch_map.at("meta_hits"));
  reco_sim_matching_ =
      static_cast<Matching *>(branch_map.at("mdc_vtx_tracks2sim_tracks"));
  mdc_meta_matching_ = static_cast<Matching *>(branch_map.at("mdc_vtx_tracks2meta_hits"));
  auto sim_event_config = config_->GetBranchConfig("sim_header");
  auto reco_event_config = config_->GetBranchConfig("event_header");
  auto sim_tracks_config = config_->GetBranchConfig("sim_tracks");
  auto reco_tracks_config = config_->GetBranchConfig("mdc_vtx_tracks");
  auto meta_hits_config = config_->GetBranchConfig("meta_hits");

  fields_id_.insert(std::make_pair(
      HITS_TOF_RPC, reco_event_config.GetFieldId("selected_tof_rpc_hits")));
  fields_id_.insert(
      std::make_pair(SIM_GEANT_PID, sim_tracks_config.GetFieldId("geant_pid")));
  fields_id_.insert(
      std::make_pair(IS_PRIMARY, sim_tracks_config.GetFieldId("is_primary")));
  fields_id_.insert(std::make_pair(RECO_GEANT_PID,
                                   reco_tracks_config.GetFieldId("geant_pid")));
  fields_id_.insert(std::make_pair(PSI_RP,
                                   sim_event_config.GetFieldId("reaction_plane")));

  for (int i = 0; i < 8; ++i) {
    int percentile = 2 + i * 5;
    std::string name = "pdg_tracks_prim_" + std::to_string(percentile);
    pdg_tracks_prim_.push_back(new TH2F(name.data(),
                                        "; y-y_{beam};p_{T}, [GeV/c]; conuts",
                                        100, -1.0, 1.0, 100, 0.0, 2.0 ));
    name = "pdg_tracks_sec_" + std::to_string(percentile);
    pdg_tracks_sec_.push_back(new TH2F(name.data(),
                                       ";y-y_{beam};p_{T}, [GeV/c];  conuts",
                                       100, -1.0, 1.0, 100, 0.0, 2.0 ));

    name = "pid_tracks_prim_" + std::to_string(percentile);
    pid_tracks_prim_.push_back(new TH2F(name.data(),
                                        ";y-y_{beam};p_{T}, [GeV/c];  conuts",
                                        100, -1.0, 1.0, 100, 0.0, 2.0 ));
    name = "pid_tracks_sec_" + std::to_string(percentile);
    pid_tracks_sec_.push_back(new TH2F(name.data(),
                                       ";y-y_{beam};p_{T}, [GeV/c];  conuts",
                                       100, -1.0, 1.0, 100, 0.0, 2.0 ));
    name = "pid_tracks_mismatch_" + std::to_string(percentile);
    pid_tracks_mismatch_.push_back(
        new TH2F(name.data(), ";y-y_{beam};p_{T}, [GeV/c];  conuts", 100, -1.0, 1.0, 100, 0.0,2.0));
    name = "pid_reco_" + std::to_string(percentile);
    pid_reco_.push_back(
        new TH2F(name.data(), ";y-y_{beam};p_{T}, [GeV/c];  conuts", 100, -1.0, 1.0, 100, 0.0,2.0));

    name = "pid_prim_phi_pt_midrapidity_" + std::to_string(percentile);
    pid_prim_phi_pt_midrapidity_.push_back(
        new TH2F(name.data(), ";p_{T}, [GeV/c]; #phi, [rad]; conuts", 100, 0.0,2.0,
                 100, -3.15, 3.15));
    name = "pdg_prim_phi_pt_midrapidity_" + std::to_string(percentile);
    pdg_prim_phi_pt_midrapidity_.push_back(
        new TH2F(name.data(), ";p_{T}, [GeV/c]; #phi, [rad]; conuts", 100, 0.0,2.0,
                 100, -3.15, 3.15));

    name = "pid_prim_delta_phi_pt_midrapidity_" + std::to_string(percentile);
    pid_prim_delta_phi_pt_midrapidity_.push_back(
        new TH2F(name.data(), ";p_{T}, [GeV/c]; #phi-#Psi_{RP}, [rad]; conuts", 100, 0.0,2.0,
                 100, -3.5, 3.5));
    name = "pdg_prim_delta_phi_pt_midrapidity_" + std::to_string(percentile);
    pdg_prim_delta_phi_pt_midrapidity_.push_back(
        new TH2F(name.data(), ";p_{T}, [GeV/c]; #phi-#Psi_{RP}, [rad]; conuts", 100, 0.0,2.0,
                 100, -3.5, 3.5));
  }
  momentum_err_ = new TProfile("momentum_err", ";p, [GeV/c]; relative error",
                               100, 0.0, 3.5);
}

void RecoAcceptance::Exec() {
  auto hits_tof_rpc = reco_header_->GetField<int>(fields_id_.at(HITS_TOF_RPC));
  auto psi_rp = sim_header_->GetField<float>( fields_id_.at(PSI_RP) );
  int centrality_class = (int)HadesUtils::Centrality::GetClass(
      hits_tof_rpc, HadesUtils::DATA_TYPE::AuAu_1_23AGeV);
  if (centrality_class > 7) {
    return;
  }
  std::vector<int> sim_matches;
  int n_reco_tracks = reco_tracks_->GetNumberOfChannels();
  for (int i = 0; i < n_reco_tracks; ++i) {
    int sim_id = reco_sim_matching_->GetMatchDirect(i);
    int meta_id = mdc_meta_matching_->GetMatchDirect(i);
    auto r_track = reco_tracks_->GetChannel(i);
    auto r_hit = meta_hits_->GetChannel(meta_id);
    auto s_track = sim_tracks_->GetChannel((sim_id));

    float m_reco = r_track.GetMass();
    float m_sim = s_track.GetMass();
    auto p_reco = r_track.Get4MomentumByMass(m_reco);
    auto p_sim = s_track.Get4MomentumByMass(m_sim);

    if( r_track.GetField<int>(fields_id_.at(RECO_GEANT_PID)) == 14 ){
      pid_reco_.at(centrality_class)->Fill(p_reco.Rapidity() - 0.74, p_reco.Pt());
      if (s_track.GetField<bool>(fields_id_.at(IS_PRIMARY))){
        pid_tracks_prim_.at(centrality_class)
            ->Fill(p_reco.Rapidity() - 0.74, p_reco.Pt());
        if (-0.05 < p_reco.Rapidity() - 0.74 &&
            p_reco.Rapidity() - 0.74 < 0.05) {
          pid_prim_phi_pt_midrapidity_.at(centrality_class)
              ->Fill(p_reco.Pt(), p_reco.Phi());
          auto delta_phi = p_reco.Phi() - psi_rp;
          if (delta_phi < -M_PI)
            delta_phi += 2 * M_PI;
          if (delta_phi > M_PI)
            delta_phi -= 2 * M_PI;
          pid_prim_delta_phi_pt_midrapidity_.at(centrality_class)
              ->Fill(p_reco.Pt(), delta_phi);
        }
      }
      if (!s_track.GetField<bool>(fields_id_.at(IS_PRIMARY))){
        pid_tracks_sec_.at(centrality_class)
            ->Fill(p_reco.Rapidity() - 0.74, p_reco.Pt());
        if (s_track.GetField<int>(fields_id_.at(SIM_GEANT_PID)) == 14)
          pdg_tracks_sec_.at(centrality_class)
              ->Fill(p_sim.Rapidity() - 0.74, p_sim.Pt());
      }
    }
    if( s_track.GetField<int>(fields_id_.at(RECO_GEANT_PID)) == 14 ){
      if (s_track.GetField<bool>(fields_id_.at(IS_PRIMARY))){
        pdg_tracks_prim_.at(centrality_class)
            ->Fill(p_sim.Rapidity() - 0.74, p_sim.Pt());
        if( -0.05 < p_sim.Rapidity() - 0.74 && p_sim.Rapidity() - 0.74 < 0.05 ) {
          pdg_prim_phi_pt_midrapidity_.at(centrality_class)
              ->Fill(p_sim.Pt(), p_sim.Phi());
          auto delta_phi = p_sim.Phi()-psi_rp;
          if ( delta_phi < -M_PI )
            delta_phi+=2*M_PI;
          if (delta_phi > M_PI)
            delta_phi-=2*M_PI;
          pdg_prim_delta_phi_pt_midrapidity_.at(centrality_class)
              ->Fill(p_sim.Pt(), delta_phi);
        }
      }
    }
    if (s_track.GetField<int>(fields_id_.at(SIM_GEANT_PID)) != 14)
      pid_tracks_mismatch_.at(centrality_class)
          ->Fill(p_sim.Rapidity() - 0.74, p_sim.Pt());
    sim_matches.push_back(sim_id);
  }
}

void RecoAcceptance::Finish() {
  momentum_err_->Write();
  for (size_t i = 0; i < pdg_tracks_prim_.size(); ++i) {
    pdg_tracks_prim_.at(i)->Write();
    pdg_tracks_sec_.at(i)->Write();


    pid_tracks_prim_.at(i)->Write();
    pid_tracks_sec_.at(i)->Write();
    pid_reco_.at(i)->Write();

    pid_tracks_mismatch_.at(i)->Write();

    pid_prim_phi_pt_midrapidity_.at(i)->Write();
    pdg_prim_phi_pt_midrapidity_.at(i)->Write();

    pid_prim_delta_phi_pt_midrapidity_.at(i)->Write();
    pdg_prim_delta_phi_pt_midrapidity_.at(i)->Write();
  }
}
} // namespace AnalysisTree