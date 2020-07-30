//
// Created by mikhail on 6/29/20.
//

#ifndef EFFICIENCY_SRC_SIM_ACCEPTANCE_H_
#define EFFICIENCY_SRC_SIM_ACCEPTANCE_H_

#include <AnalysisTree/EventHeader.hpp>
#include <AnalysisTree/FillTask.hpp>

#include <AnalysisTree/Cuts.hpp>
#include <AnalysisTree/Detector.hpp>
#include <TChain.h>
#include <TProfile2D.h>

#include "centrality.h"


namespace AnalysisTree {
class SimAcceptance : public FillTask {
  void Init(std::map<std::string, void *> &branch_map) override {
    config_->Print();
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
        std::make_pair(PSI_RP, sim_event_config.GetFieldId("reaction_plane")));
    fields_id_.insert(std::make_pair(
        SIM_GEANT_PID, sim_tracks_config.GetFieldId("geant_pid")));
    fields_id_.insert(
        std::make_pair(IS_PRIMARY, sim_tracks_config.GetFieldId("is_primary")));
    gen_flow_y_005_015_ = new TProfile( "gen_flow_c_0_5_y_005_015", ";pt; v_{1}", 16, 0.0, 1.6 );
    gen_flow_y_035_045_ = new TProfile( "gen_flow_c_0_5_y_035_045", ";pt; v_{1}", 16, 0.0, 1.6 );
    gen_flow_y_065_075_ = new TProfile( "gen_flow_c_0_5_y_065_075", ";pt; v_{1}", 16, 0.0, 1.6 );

    for (int i = 0; i < 8; ++i) {
      int percentile = 2 + i * 5;
      std::string name = "gen_acceptance_prim_" + std::to_string(percentile);
      sim_acceptance_prim_.push_back(
          new TH2F(name.data(), ";p_{T}, [GeV/c]; y-y_{beam}; conuts", 100, 0.0,
                   2.0, 100, -1.0, 1.0));
      name = "gen_acceptance_sec_" + std::to_string(percentile);
      sim_acceptance_sec_.push_back(
          new TH2F(name.data(), ";p_{T}, [GeV/c]; y-y_{beam}; conuts", 100, 0.0,
                   2.0, 100, -1.0, 1.0));
      name = "sim_density_" + std::to_string(percentile);
    }
    proton_yield_ = new TH1F("gen_proton_yield", ";N protons; counts", 100, 0, 100);
  }
  void Exec() override {
    auto hits_tof_rpc =
        reco_header_->GetField<int>(fields_id_.at(HITS_TOF_RPC));
    int centrality_class =
        Centrality::GetInstance()->GetCentralityClass5pc(hits_tof_rpc);
    float psi_rp = sim_header_->GetField<float>(fields_id_.at(PSI_RP));
    if (centrality_class > 7) {
      return;
    }
    int n_protons{0};
    proton_yield_->Fill(n_protons);
    int n_sim_tracks = sim_tracks_->GetNumberOfChannels();
    for (int i = 0; i < n_sim_tracks; ++i) {
      auto s_track = (sim_tracks_->GetChannel(i));
      if (s_track.GetField<int>(fields_id_.at(SIM_GEANT_PID)) != 14)
        continue;
      float m_sim = s_track.GetMass();
      auto p_sim = s_track.Get4MomentumByMass(m_sim);
      if (s_track.GetField<bool>(fields_id_.at(IS_PRIMARY)))
        sim_acceptance_prim_.at(centrality_class)->Fill(p_sim.Pt(), p_sim.Rapidity() - 0.74);
      if (!s_track.GetField<bool>(fields_id_.at(IS_PRIMARY)))
        sim_acceptance_sec_.at(centrality_class)->Fill(p_sim.Pt(), p_sim.Rapidity() - 0.74);
      if( centrality_class > 0 )
        continue;
      if( !s_track.GetField<bool>(fields_id_.at(IS_PRIMARY)) )
        continue;

      double y_cm = p_sim.Rapidity() - 0.74;
      double v_n = cos( p_sim.Phi() - psi_rp );
      if( -0.15 <= y_cm && y_cm <= -0.05  ) {
        gen_flow_y_005_015_->Fill(p_sim.Pt(), v_n);
        continue;
      }
      if( -0.45 <= y_cm && y_cm <= -0.35  ) {
        gen_flow_y_035_045_->Fill(p_sim.Pt(), v_n );
        continue;
      }
      if( -0.75 <= y_cm && y_cm <= -0.65  ) {
        gen_flow_y_065_075_->Fill(p_sim.Pt(), v_n );
        continue;
      }
    }
  }
  void Finish() override {
    proton_yield_->Write();
    for (size_t i = 0; i < sim_acceptance_prim_.size(); ++i) {
      sim_acceptance_prim_.at(i)->Write();
      sim_acceptance_sec_.at(i)->Write();
    }
    gen_flow_y_005_015_->Write();
    gen_flow_y_035_045_->Write();
    gen_flow_y_065_075_->Write();
  }

private:
  enum FIELDS {
    HITS_TOF_RPC,
    PSI_RP,
    SIM_GEANT_PID,
    IS_PRIMARY,
  };
  EventHeader *sim_header_{nullptr};
  EventHeader *reco_header_{nullptr};
  Particles *sim_tracks_{nullptr};

  std::map<int, int> fields_id_;
  TH1F *proton_yield_{nullptr};
  std::vector<TH2F *> sim_acceptance_prim_;
  std::vector<TH2F *> sim_acceptance_sec_;
  TProfile* gen_flow_y_005_015_;
  TProfile* gen_flow_y_035_045_;
  TProfile* gen_flow_y_065_075_;
};
} // namespace AnalysisTree
#endif // EFFICIENCY_SRC_SIM_ACCEPTANCE_H_
