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
  y_beam_ = data_header_->GetBeamRapidity();

  fields_id_.insert(std::make_pair(
      HITS_TOF_RPC, reco_event_config.GetFieldId("selected_tof_rpc_hits")));
  fields_id_.insert(
      std::make_pair(SIM_GEANT_PID, sim_tracks_config.GetFieldId("geant_pid")));
  fields_id_.insert(
      std::make_pair(IS_PRIMARY, sim_tracks_config.GetFieldId("is_primary")));
  fields_id_.insert(std::make_pair(PSI_RP,
                                   sim_event_config.GetFieldId("reaction_plane")));
  std::vector<double> M0_axis;
  for(int j=0; j<21; ++j){ M0_axis.push_back(5.0f* (float) j); }
  std::vector<double> theta_axis;
  for(int j=0; j<49; ++j){ theta_axis.push_back(0.3f+ 0.025*(double) j); }
  std::vector<double> delta_phi_axis;
  for(int j=0; j<17; ++j){ delta_phi_axis.push_back(-3.2+ 0.4*j); }
  std::vector<double> y_axis;
  std::vector<double> pt_axis;
  if( pid_code_ == 2212 ) {
    pt_axis = {0, 0.29375, 0.35625, 0.41875, 0.48125, 0.54375,0.61875, 0.70625, 0.81875, 1.01875, 2.0};
    for(int j=0; j<16; ++j){ y_axis.push_back(-0.75+0.1* (double) j); }
  }
  if( abs(pid_code_) == 211 ){
    pt_axis = {0, 0.08, 0.105, 0.13, 0.155, 0.18, 0.21, 0.25, 0.315, 0.535, 1.0};
    for(int j=0; j<18; ++j){ y_axis.push_back(-0.65+0.1* (double) j); }
  }

  gen_tracks_prim_cent_ = new TH3F("gen_tracks_prim_cent", ";y-y_{beam};p_{T}, [GeV/c]; centrality (%)",
                                   y_axis.size()-1, y_axis.data(),
                                   pt_axis.size()-1, pt_axis.data(),
                                   M0_axis.size()-1, M0_axis.data());
  theta_centrality_ = new TH2F( "gen_tracks_theta_centrality", ";#theta [rad];centrality (%)",
                               48, 0.3, 1.5,
                               20, 0, 100);
  theta_pT_centrality_ = new TH3F("gen_theta_pT_centrality", ";#theta [rad];p_{T}, [GeV/c];centrality (%)",
                                  theta_axis.size()-1, theta_axis.data(),
                                  pt_axis.size()-1, pt_axis.data(),
                                  M0_axis.size()-1, M0_axis.data());

  pT_delta_phi_centrality_ = new TH3F("gen_pT_delta_phi_centrality", ";p_{T} [GeV/c];#phi-#Psi [rad];centrality (%)",
                                      pt_axis.size()-1, pt_axis.data(),
                                      delta_phi_axis.size()-1, delta_phi_axis.data(),
                                      M0_axis.size()-1, M0_axis.data());
  y_pT_theta_ = new TH3F("y_pT_theta_", ";y-y_{cm};p_{T} [GeV/c];#theta [rad]",
                         y_axis.size()-1, y_axis.data(),
                         pt_axis.size()-1, pt_axis.data(),
                        theta_axis.size()-1, theta_axis.data());
  theta_centrality_all_ = new TH2F( "gen_theta_centrality_all",
                                   ";#theta [rad];centrality (%)",
                                   theta_axis.size()-1, theta_axis.data(),
                                   M0_axis.size()-1, M0_axis.data());
  for (int i = 0; i < 12; ++i) {
    int percentile = 2 + i * 5;
    float y_axis[16];
    for(int j=0; j<16; ++j){ y_axis[j]=-0.75f+0.1f* (float) j; }
    float pt_axis[]={0, 0.29375, 0.35625, 0.41875, 0.48125, 0.54375, 0.61875, 0.70625, 0.81875, 1.01875, 2.0};
    std::string name = "gen_tracks_prim_" + std::to_string(percentile);
    gen_tracks_prim_.push_back(
        new TH2F(name.data(), ";y-y_{beam};p_{T}, [GeV/c]; conuts",
                 100, -0.85, 1.15, 125, 0.0, 2.5));
    name = "gen_tracks_sec_" + std::to_string(percentile);
    gen_tracks_sec_.push_back(
        new TH2F(name.data(), ";y-y_{beam};p_{T}, [GeV/c]; conuts",
                 100, -0.85, 1.15, 125, 0.0, 2.5));
    name = "gen_occupancy_" + std::to_string(percentile);
  }
}
void SimAcceptance::Exec() {
  auto centrality = reco_header_->GetField<float>(
      config_->GetBranchConfig("event_header").GetFieldId("selected_tof_rpc_hits_centrality") );
  auto centrality_class = (size_t) ( (centrality-2.5)/5.0 );
  if( h1_centrality_parameters_ ){
    auto nhits_meta = reco_header_->GetField<int>(
        config_->GetBranchConfig("event_header").GetFieldId("selected_tof_rpc_hits") );
    auto mult_bin = h1_centrality_parameters_->FindBin( nhits_meta );
    centrality_class = (size_t)  h1_centrality_parameters_->GetBinContent( mult_bin )-1;
    centrality = 2.5+5.0*centrality_class;
  }
  if (centrality_class > 11) {
    return;
  }
  auto psi_rp = sim_header_->GetField<float>( config_->GetBranchConfig("sim_header").GetFieldId("reaction_plane") );
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
    if( charge != 0 )
      theta_centrality_all_->Fill( p_sim.Theta(), centrality );
    if (s_track.GetPid() != pid_code_)
      continue;
    if (s_track.GetField<bool>(fields_id_.at(IS_PRIMARY))) {
      gen_tracks_prim_.at(centrality_class)
          ->Fill(p_sim.Rapidity() - y_beam_, p_sim.Pt());
      gen_tracks_prim_cent_->Fill( p_sim.Rapidity() - y_beam_, p_sim.Pt(), centrality );
      theta_centrality_->Fill( p_sim.Theta(), centrality );
      theta_pT_centrality_->Fill( p_sim.Theta(), p_sim.Pt(), centrality );
      y_pT_theta_->Fill( p_sim.Rapidity() - y_beam_, p_sim.Pt(), p_sim.Theta() );
      if( -0.05 < p_sim.Rapidity() - y_beam_ && p_sim.Rapidity() - y_beam_ < 0.05 ){
        auto delta_phi = p_sim.Phi() - psi_rp;
        if( delta_phi < -M_PI )
          delta_phi+=2*M_PI;
        if( delta_phi > M_PI )
          delta_phi-=2*M_PI;
        pT_delta_phi_centrality_->Fill( p_sim.Pt(), delta_phi, centrality );
      }
    }
    if (!s_track.GetField<bool>(fields_id_.at(IS_PRIMARY)))
      gen_tracks_sec_.at(centrality_class)
          ->Fill(p_sim.Rapidity() - y_beam_, p_sim.Pt());
  }
}
void SimAcceptance::Finish() {
  for (size_t i = 0; i < gen_tracks_prim_.size(); ++i) {
    gen_tracks_prim_.at(i)->Write();
    gen_tracks_sec_.at(i)->Write();
  }
  theta_centrality_all_->Write();
  gen_tracks_prim_cent_->Write();
  theta_pT_centrality_->Write();
  theta_centrality_->Write();
  pT_delta_phi_centrality_->Write();
  y_pT_theta_->Write();
}
void SimAcceptance::SetPidCode(int pid_code) { pid_code_ = pid_code; }
void SimAcceptance::SetCentralityParameters(TH1 *h_1_centrality_parameters) {
  h1_centrality_parameters_ = h_1_centrality_parameters;
}
} // namespace AnalysisTree