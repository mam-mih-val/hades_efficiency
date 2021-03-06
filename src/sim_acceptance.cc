//
// Created by mikhail on 6/29/20.
//

#include "sim_acceptance.h"
#include <TDatabasePDG.h>

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

  for (int i = 0; i < 12; ++i) {
    int percentile = 2 + i * 5;
    float y_axis[16];
    for(int j=0; j<16; ++j){ y_axis[j]=-0.75f+0.1f* (float) j; }
    float pt_axis[]={0, 0.29375, 0.35625, 0.41875, 0.48125, 0.54375, 0.61875, 0.70625, 0.81875, 1.01875, 2.0};
    std::string name = "gen_tracks_prim_" + std::to_string(percentile);
    gen_tracks_prim_.push_back(
        new TH2F(name.data(), ";y-y_{beam};p_{T}, [GeV/c]; conuts",
                 180, -0.85, 0.95, 200, 0.0, 2.0));
    name = "gen_tracks_sec_" + std::to_string(percentile);
    gen_tracks_sec_.push_back(
        new TH2F(name.data(), ";y-y_{beam};p_{T}, [GeV/c]; conuts",
                 180, -0.85, 0.95, 200, 0.0, 2.0));
    name = "gen_occupancy_" + std::to_string(percentile);
    gen_occupancy_.push_back(
        new TH2F(name.data(), ";#eta;p, [GeV/c]; conuts",
                 210, 0.0, 2.1, 250, 0.0, 2.5));
  }
}
void SimAcceptance::Exec() {
  auto centrality = reco_header_->GetField<float>(
      config_->GetBranchConfig("event_header").GetFieldId("selected_tof_rpc_hits_centrality") );
  auto centrality_class = (size_t) ( (centrality-2.5)/5.0 );
  if (centrality_class > 11) {
    return;
  }
  int n_sim_tracks = sim_tracks_->GetNumberOfChannels();
  for (int i = 0; i < n_sim_tracks; ++i) {
    auto s_track = (sim_tracks_->GetChannel(i));
    float m_sim = s_track.GetMass();
    auto p_sim = s_track.Get4MomentumByMass(m_sim);
    int charge=0;
    if( TDatabasePDG::Instance()->GetParticle(s_track.GetPid()) )
      charge = TDatabasePDG::Instance()->GetParticle(s_track.GetPid())->Charge() / 3;
    else{
      charge = (int)(s_track.GetPid() / 1E+4) % (int)1e+3;
    }
    if(charge!=0)
      gen_occupancy_.at(centrality_class)->Fill(p_sim.PseudoRapidity(), p_sim.P());
    if (s_track.GetPid() != pid_code_)
      continue;

    if (s_track.GetField<bool>(fields_id_.at(IS_PRIMARY))) {
      gen_tracks_prim_.at(centrality_class)
          ->Fill(p_sim.Rapidity() - 0.74, p_sim.Pt());
    }
    if (!s_track.GetField<bool>(fields_id_.at(IS_PRIMARY)))
      gen_tracks_sec_.at(centrality_class)
          ->Fill(p_sim.Rapidity() - 0.74, p_sim.Pt());
  }
}
void SimAcceptance::Finish() {
  for (size_t i = 0; i < gen_tracks_prim_.size(); ++i) {
    gen_tracks_prim_.at(i)->Write();
    gen_tracks_sec_.at(i)->Write();
    gen_occupancy_.at(i)->Write();
  }
}
void SimAcceptance::SetPidCode(int pid_code) { pid_code_ = pid_code; }
} // namespace AnalysisTree