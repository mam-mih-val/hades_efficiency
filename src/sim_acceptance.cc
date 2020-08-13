//
// Created by mikhail on 6/29/20.
//

#include "sim_acceptance.h"

namespace AnalysisTree {
void SimAcceptance::Init(std::map<std::string, void *> &branch_map) {
  sim_header_ = static_cast<EventHeader *>(branch_map.at("sim_header"));
  reco_header_ = static_cast<EventHeader *>(branch_map.at("event_header"));
  sim_tracks_ = static_cast<Particles *>(branch_map.at("sim_tracks"));
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

  for (int i = 0; i < 8; ++i) {
    int percentile = 2 + i * 5;
    std::string name = "gen_acceptance_prim_" + std::to_string(percentile);
    gen_tracks_prim_.push_back(
        new TH2F(name.data(), ";p_{T}, [GeV/c]; y-y_{beam}; conuts", 100, 0.0,
                 2.0, 100, -1.0, 1.0));
    name = "gen_acceptance_sec_" + std::to_string(percentile);
    gen_tracks_sec_.push_back(
        new TH2F(name.data(), ";p_{T}, [GeV/c]; y-y_{beam}; conuts", 100, 0.0,
                 2.0, 100, -1.0, 1.0));
    name = "sim_density_" + std::to_string(percentile);
  }
}
void SimAcceptance::Exec() {
  auto hits_tof_rpc = reco_header_->GetField<int>(fields_id_.at(HITS_TOF_RPC));
  int centrality_class = (int) HadesUtils::Centrality::GetClass(hits_tof_rpc, HadesUtils::DATA_TYPE::AuAu_1_23AGeV);
  if (centrality_class > 7) {
    return;
  }
  int n_sim_tracks = sim_tracks_->GetNumberOfChannels();
  for (int i = 0; i < n_sim_tracks; ++i) {
    auto s_track = (sim_tracks_->GetChannel(i));
    if (s_track.GetField<int>(fields_id_.at(SIM_GEANT_PID)) != 14)
      continue;
    float m_sim = s_track.GetMass();
    auto p_sim = s_track.Get4MomentumByMass(m_sim);
    if (s_track.GetField<bool>(fields_id_.at(IS_PRIMARY)))
      gen_tracks_prim_.at(centrality_class)
          ->Fill(p_sim.Pt(), p_sim.Rapidity() - 0.74);
    if (!s_track.GetField<bool>(fields_id_.at(IS_PRIMARY)))
      gen_tracks_sec_.at(centrality_class)
          ->Fill(p_sim.Pt(), p_sim.Rapidity() - 0.74);
    if (centrality_class > 0)
      continue;
    if (!s_track.GetField<bool>(fields_id_.at(IS_PRIMARY)))
      continue;
  }
}
void SimAcceptance::Finish() {
  for (size_t i = 0; i < gen_tracks_prim_.size(); ++i) {
    gen_tracks_prim_.at(i)->Write();
    gen_tracks_sec_.at(i)->Write();
  }
}
} // namespace AnalysisTree