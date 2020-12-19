//
// Created by mikhail on 6/29/20.
//

#include "sim_acceptance.h"

namespace AnalysisTree {
void SimAcceptance::Init(std::map<std::string, void *> &branch_map) {
  sim_header_ = static_cast<EventHeader *>(branch_map.at("sim_header"));
  reco_header_ = static_cast<EventHeader *>(branch_map.at("event_header"));
  sim_tracks_ = static_cast<Particles *>(branch_map.at("sim_tracks"));
  reco_tracks_ = static_cast<Particles *>(branch_map.at("mdc_vtx_tracks"));
  auto sim_event_config = config_->GetBranchConfig("sim_header");
  auto reco_event_config = config_->GetBranchConfig("event_header");
  auto sim_tracks_config = config_->GetBranchConfig("sim_tracks");
  auto reco_tracks_config = config_->GetBranchConfig("mdc_vtx_tracks");

  fields_id_.insert(std::make_pair(
      HITS_TOF_RPC, reco_event_config.GetFieldId("selected_tof_rpc_hits")));
  fields_id_.insert(
      std::make_pair(SIM_GEANT_PID, sim_tracks_config.GetFieldId("geant_pid")));
  fields_id_.insert(
      std::make_pair(IS_PRIMARY, sim_tracks_config.GetFieldId("is_primary")));
  fields_id_.insert(std::make_pair(PSI_RP,
                                   sim_event_config.GetFieldId("reaction_plane")));

  for (int i = 0; i < 8; ++i) {
    int percentile = 2 + i * 5;
    std::string name = "gen_tracks_prim_" + std::to_string(percentile);
    gen_tracks_prim_.push_back(
        new TH2F(name.data(), ";y-y_{beam};p_{T}, [GeV/c]; conuts",
                 100, -1.0, 1.0, 100, 0.0,2.0));
    name = "gen_tracks_sec_" + std::to_string(percentile);
    gen_tracks_sec_.push_back(
        new TH2F(name.data(), ";y-y_{beam};p_{T}, [GeV/c]; conuts",
                 100, -1.0, 1.0, 100, 0.0,2.0));

    name = "gen_prim_phi_pt_midrapidity_" + std::to_string(percentile);
    gen_prim_phi_pt_rapidity_.push_back(
        new TH3F(name.data(), ";y_{cm};p_{T}, [GeV/c]; #phi, [rad]; conuts",
                 15, -0.75, 0.75,
                 20, 0.0,2.0,
                 32, -3.2, 3.2));
    float y_axis[16];
    for(int j=0; j<16; ++j){ y_axis[j]=-0.75f+0.1f* (float) j; }
    float pt_axis[]={0, 0.29375, 0.35625, 0.41875, 0.48125, 0.54375, 0.61875, 0.70625, 0.81875, 1.01875, 2.0};
    float phi_axis[17];
    for(int j=0; j<17; ++j){ phi_axis[j]=-3.2f+0.4f* (float) j; }
    name = "gen_prim_delta_phi_pt_midrapidity_" + std::to_string(percentile);
    gen_prim_delta_phi_pt_rapidity_.push_back(
        new TH3F(name.data(), ";y_{cm};p_{T}, [GeV/c]; #phi-#Psi_{RP}, [rad]; conuts",
                 15, y_axis,
                 10, pt_axis,
                 16, phi_axis));
  }
  entries_vs_pT_y_n_tracks_sector_ = new TH3F( "gen_prim_pT_y_n_tracks_sector",
                                               ";y;p_{T} [GeV/c];N tracks in sector",
                                               100, -1, 1.0,
                                               100, 0.0, 2.0,
                                               30, 0.0, 30.0);
}
void SimAcceptance::Exec() {
  auto hits_tof_rpc = reco_header_->GetField<int>(fields_id_.at(HITS_TOF_RPC));
  int centrality_class = (int) HadesUtils::Centrality::GetClass(hits_tof_rpc, HadesUtils::DATA_TYPE::AuAu_1_23AGeV);
  auto psi_rp = sim_header_->GetField<float>( fields_id_.at(PSI_RP) );
  if (centrality_class > 7) {
    return;
  }
  int n_sim_tracks = sim_tracks_->GetNumberOfChannels();
  std::vector n_tracks_sectors{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int n_reco_tracks = reco_tracks_->GetNumberOfChannels();
  for( int i = 0; i < n_reco_tracks; ++i ){
    auto r_track = reco_tracks_->GetChannel(i);
    auto sector = WhatSector( r_track.GetPhi() );
    n_tracks_sectors.at(sector)++;
  }
  for (int i = 0; i < n_sim_tracks; ++i) {
    auto s_track = (sim_tracks_->GetChannel(i));
    if (s_track.GetPid() != pid_code_)
      continue;
    float m_sim = s_track.GetMass();
    auto p_sim = s_track.Get4MomentumByMass(m_sim);
    if (s_track.GetField<bool>(fields_id_.at(IS_PRIMARY))) {
      gen_tracks_prim_.at(centrality_class)
          ->Fill(p_sim.Rapidity() - 0.74, p_sim.Pt());
      gen_prim_phi_pt_rapidity_.at(centrality_class)
          ->Fill(p_sim.Rapidity() - 0.74, p_sim.Pt(), p_sim.Phi());
      auto delta_phi = p_sim.Phi()-psi_rp;
      if ( delta_phi < -M_PI )
        delta_phi+=2*M_PI;
      if (delta_phi > M_PI)
        delta_phi-=2*M_PI;
      gen_prim_delta_phi_pt_rapidity_.at(centrality_class)
          ->Fill(p_sim.Rapidity() - 0.74, p_sim.Pt(), delta_phi);
      auto sector = WhatSector(p_sim.Phi());
      entries_vs_pT_y_n_tracks_sector_->Fill( p_sim.Rapidity() - 0.74, p_sim.Pt(), n_tracks_sectors.at(sector) );
    }
    if (!s_track.GetField<bool>(fields_id_.at(IS_PRIMARY)))
      gen_tracks_sec_.at(centrality_class)
          ->Fill(p_sim.Rapidity() - 0.74, p_sim.Pt());
  }
}
void SimAcceptance::Finish() {
  entries_vs_pT_y_n_tracks_sector_->Write();
  for (size_t i = 0; i < gen_tracks_prim_.size(); ++i) {
    gen_tracks_prim_.at(i)->Write();
    gen_tracks_sec_.at(i)->Write();
    gen_prim_phi_pt_rapidity_.at(i)->Write();
    gen_prim_delta_phi_pt_rapidity_.at(i)->Write();
  }
}
void SimAcceptance::SetPidCode(int pid_code) { pid_code_ = pid_code; }
} // namespace AnalysisTree