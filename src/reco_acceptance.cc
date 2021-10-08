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
  y_beam_ = data_header_->GetBeamRapidity();
  fields_id_.insert(std::make_pair(
      HITS_TOF_RPC, reco_event_config.GetFieldId("selected_tof_rpc_hits")));
  fields_id_.insert(
      std::make_pair(SIM_GEANT_PID, sim_tracks_config.GetFieldId("geant_pid")));
  fields_id_.insert(
      std::make_pair(IS_PRIMARY, sim_tracks_config.GetFieldId("is_primary")));
  fields_id_.insert(std::make_pair(RECO_GEANT_PID,
                                   reco_tracks_config.GetFieldId("geant_pid")));
  fields_id_.insert(std::make_pair(LAYERS_0,
                                   reco_tracks_config.GetFieldId("layers_0")));
  fields_id_.insert(std::make_pair(LAYERS_1,
                                   reco_tracks_config.GetFieldId("layers_1")));
  fields_id_.insert(std::make_pair(LAYERS_2,
                                   reco_tracks_config.GetFieldId("layers_2")));
  fields_id_.insert(std::make_pair(LAYERS_3,
                                   reco_tracks_config.GetFieldId("layers_3")));
  fields_id_.insert(std::make_pair(CHI2,
                                     reco_tracks_config.GetFieldId("chi2")));
  fields_id_.insert(std::make_pair(DCA_XY,
                                     reco_tracks_config.GetFieldId("dca_xy")));
  fields_id_.insert(std::make_pair(DCA_Z,
                                     reco_tracks_config.GetFieldId("dca_z")));
  fields_id_.insert(std::make_pair(PSI_RP,
                                   sim_event_config.GetFieldId("reaction_plane")));

  pdg_tracks_cent_ = new TH3F("pdg_tracks_cent", ";y-y_{beam};p_{T}, [GeV/c]; centrality (%)",
                                   100, -0.85, 1.15,
                                   125, 0.0, 2.5,
                                   20, 0.0, 100.5);

  for (int i = 0; i < 12; ++i) {
    int percentile = 2 + i * 5;
    float y_axis[16];
    for(int j=0; j<16; ++j){ y_axis[j]=-0.75f+0.1f* (float) j; }
    float pt_axis[]={0, 0.29375, 0.35625, 0.41875, 0.48125, 0.54375, 0.61875, 0.70625, 0.81875, 1.01875, 2.0};
    std::string name = "pdg_tracks_prim_" + std::to_string(percentile);
    pdg_tracks_prim_.push_back(new TH2F(name.data(),
                                        "; y-y_{beam};p_{T}, [GeV/c]; conuts",
                                        100, -0.85, 1.15, 125, 0.0, 2.5 ));
    name = "pdg_tracks_sec_" + std::to_string(percentile);
    pdg_tracks_sec_.push_back(new TH2F(name.data(),
                                       ";y-y_{beam};p_{T}, [GeV/c];  conuts",
                                       100, -0.85, 1.15, 125, 0.0, 2.5 ));

    name = "pid_tracks_prim_" + std::to_string(percentile);
    pid_tracks_prim_.push_back(new TH2F(name.data(),
                                        ";y-y_{beam};p_{T}, [GeV/c];  conuts",
                                        100, -0.85, 1.15, 125, 0.0, 2.5 ));
    name = "pid_tracks_sec_" + std::to_string(percentile);
    pid_tracks_sec_.push_back(new TH2F(name.data(),
                                       ";y-y_{beam};p_{T}, [GeV/c];  conuts",
                                       100, -0.85, 1.15, 125, 0.0, 2.5 ));
    name = "pid_tracks_mismatch_" + std::to_string(percentile);
    pid_tracks_mismatch_.push_back(
        new TH2F(name.data(), ";y-y_{beam};p_{T}, [GeV/c];  conuts", 100, -0.85, 1.15, 125, 0.0, 2.5));
    name = "pid_reco_" + std::to_string(percentile);
    pid_reco_.push_back(
        new TH2F(name.data(), ";y-y_{beam};p_{T}, [GeV/c];  conuts", 100, -0.85, 1.15, 125, 0.0, 2.5));
    name = "rec_occupancy_" + std::to_string(percentile);
    rec_occupancy_.push_back(
        new TH2F(name.data(), ";#eta;p, [GeV/c];  conuts", 210, 0.0, 2.1, 250, 0.0, 5.0));

    name = "n_tracks_in_sector_" + std::to_string(percentile);
    n_tracks_in_sector_.push_back( new TH1F( name.c_str(), ";N_{tr}", 30, 0, 30 ) );
  }
  momentum_err_ = new TProfile("momentum_err", ";p, [GeV/c]; relative error",
                               100, 0.0, 3.5);
}

std::array<int, 6> RecoAcceptance::CalcRecoSectorsOccupancy(int pid){
  std::array<int, 6> sectors_occupancy{};
  for( auto& n : sectors_occupancy )
    n=0;
  for (int i_track = 0; i_track < reco_tracks_->GetNumberOfChannels(); ++i_track) {
    auto track = reco_tracks_->GetChannel(i_track);
    auto phi = track.GetPhi();
    auto sector = WhatSector(phi);
    auto pdg = track.GetPid();
    if( pid != -1 && pid != pdg )
      continue;
    try {
      sectors_occupancy.at(sector) += 1;
    } catch (std::out_of_range&) {}
  }
  return sectors_occupancy;
}

std::array<int, 6> RecoAcceptance::CalcSimSectorsOccupancy(int pid){
  std::array<int, 6> sectors_occupancy{};
  for( auto& n : sectors_occupancy )
    n=0;
  for (int i_track = 0; i_track < sim_tracks_->GetNumberOfChannels(); ++i_track) {
    auto track = sim_tracks_->GetChannel(i_track);
    auto phi = track.GetPhi();
    auto sector = WhatSector(phi);
    auto pdg = track.GetPid();
    if( pid != -1 && pid != pdg )
      continue;
    try {
      sectors_occupancy.at(sector) += 1;
    } catch (std::out_of_range&) {}
  }
  return sectors_occupancy;
}



void RecoAcceptance::Exec() {
  auto reco_occupancy = CalcRecoSectorsOccupancy();
  auto centrality = reco_header_->GetField<float>(
      config_->GetBranchConfig("event_header").GetFieldId("selected_tof_rpc_hits_centrality") );
  auto centrality_class = (size_t) ( (centrality-2.5)/5.0 );
  if (centrality_class > 11)
    return;
  for( size_t i=0; i<6; ++i ){
    n_tracks_in_sector_.at(centrality_class)->Fill( reco_occupancy.at(i) );
  }
  std::vector<int> sim_matches;
  int n_reco_tracks = reco_tracks_->GetNumberOfChannels();
  for (int i = 0; i < n_reco_tracks; ++i) {
    int sim_id = reco_sim_matching_->GetMatchDirect(i);
    if( sim_id == AnalysisTree::UndefValueInt )
      continue;
    int meta_id = mdc_meta_matching_->GetMatchDirect(i);
    auto r_track = reco_tracks_->GetChannel(i);
    auto r_hit = meta_hits_->GetChannel(meta_id);
    auto s_track = sim_tracks_->GetChannel((sim_id));
    rec_occupancy_.at(centrality_class)->Fill(s_track.GetEta(), s_track.GetP());
    auto chi2 = r_track.GetField<float>(fields_id_.at(CHI2));
    auto dca_xy = r_track.GetField<float>(fields_id_.at(DCA_XY));
    auto dca_z = r_track.GetField<float>(fields_id_.at(DCA_Z));
    if( chi2 > 100.0 )
      continue;
    if( fabsf(dca_xy) > 10.0 )
      continue;
    if( fabsf(dca_z) > 10.0 )
      continue;
    float m_reco = r_track.GetMass();
    float m_sim = s_track.GetMass();
    auto p_reco = r_track.Get4MomentumByMass(m_reco);
    auto p_sim = s_track.Get4MomentumByMass(m_sim);

    if( r_track.GetPid()==pid_code_ ){
      pid_reco_.at(centrality_class)->Fill(p_reco.Rapidity() - y_beam_, p_reco.Pt());
      if (s_track.GetField<bool>(fields_id_.at(IS_PRIMARY))){
        pid_tracks_prim_.at(centrality_class)
            ->Fill(p_reco.Rapidity() - y_beam_, p_reco.Pt());
        if( s_track.GetPid() == pid_code_ ) {
          pdg_tracks_prim_.at(centrality_class)
              ->Fill(p_sim.Rapidity() - y_beam_, p_sim.Pt());
          pdg_tracks_cent_->Fill(p_sim.Rapidity() - y_beam_, p_sim.Pt(), centrality);
        }
      }
      if (!s_track.GetField<bool>(fields_id_.at(IS_PRIMARY))){
        pid_tracks_sec_.at(centrality_class)
            ->Fill(p_reco.Rapidity() - y_beam_, p_reco.Pt());
        if (s_track.GetPid() == pid_code_)
          pdg_tracks_sec_.at(centrality_class)
              ->Fill(p_sim.Rapidity() - y_beam_, p_sim.Pt());
      }
    }
    if (s_track.GetPid()!=pid_code_ && r_track.GetPid()==pid_code_ )
      pid_tracks_mismatch_.at(centrality_class)
          ->Fill(p_sim.Rapidity() - y_beam_, p_sim.Pt());
    sim_matches.push_back(sim_id);
  }
}

void RecoAcceptance::Finish() {
  momentum_err_->Write();
  pdg_tracks_cent_->Write();
  for (size_t i = 0; i < pdg_tracks_prim_.size(); ++i) {
    pdg_tracks_prim_.at(i)->Write();
    pdg_tracks_sec_.at(i)->Write();

    pid_tracks_prim_.at(i)->Write();
    pid_tracks_sec_.at(i)->Write();
    pid_reco_.at(i)->Write();

    pid_tracks_mismatch_.at(i)->Write();

    n_tracks_in_sector_.at(i)->Write();

    rec_occupancy_.at(i)->Write();
  }
}
void RecoAcceptance::SetPidCode(int pid_code) { pid_code_ = pid_code; }
} // namespace AnalysisTree