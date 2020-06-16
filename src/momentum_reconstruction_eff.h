//
// Created by mikhail on 6/16/20.
//

#ifndef QUALITY_ASSURANCE_SRC_TREE_READER_H_
#define QUALITY_ASSURANCE_SRC_TREE_READER_H_

#include <AnalysisTree/EventHeader.h>
#include <AnalysisTree/FillTask.h>

#include <TChain.h>
#include <TProfile2D.h>

namespace AnalysisTree {
class MomentumReconctructionEff : FillTask{
public:
  MomentumReconctructionEff() = default;
  ~MomentumReconctructionEff() = default;
  void Init( std::map<std::string, void*>& branch_map ) override {
    config_->Print();
    sim_header_ = static_cast<EventHeader*>( branch_map.at( "sim_header" ) );
    reco_header_ = static_cast<EventHeader*>(branch_map.at( "event_header" ));
    sim_tracks_ = static_cast<TrackDetector *>(branch_map.at( "sim_tracks" ));
    reco_tracks_ = static_cast<TrackDetector *>(branch_map.at( "mdc_vtx_tracks" ));
    sim_reco_matching_ = static_cast<Matching *>(branch_map.at( "sim_reco_tracks" ));
    auto sim_event_config = config_->GetBranchConfig("sim_header");
    auto reco_event_config = config_->GetBranchConfig("event_header");
    auto sim_tracks_config = config_->GetBranchConfig("sim_tracks");
    auto reco_tracks_config = config_->GetBranchConfig("mdc_vtx_tracks");

    fields_id_.insert(  std::make_pair( HITS_TOF_RPC, reco_event_config.GetFieldId("selected_hits_tof_rpc") )  );
    fields_id_.insert(  std::make_pair( PSI_RP, sim_event_config.GetFieldId("reaction_plane") )  );

    momentum_efficiency_ = new TProfile2D( "momentum_eff", ";pt, [GeV/c]; y; efficiency", 100, 0.0, 2.5, 100, 0.0, 1.5 );
  }
  void Exec() override {
    short n_tracks = reco_tracks_->GetNumberOfChannels();
    for (short i = 0; i < n_tracks; ++i) {
      short sim_id = sim_reco_matching_->GetMatchInverted(i);
      auto r_track = static_cast<Particle>(reco_tracks_->GetChannel(i));
      auto s_track = static_cast<Particle>(sim_tracks_->GetChannel( (sim_id) ));
      float m_reco = r_track.GetMass();
      float m_sim = s_track.GetMass();

      auto p_reco = r_track.Get4MomentumByMass(m_reco);
      auto p_sim = s_track.Get4MomentumByMass(m_sim);

      double eff = 1.0 - fabs(p_reco.P() - p_sim.P())/p_sim.P();
      momentum_efficiency_->Fill( p_sim.Pt(), p_sim.Rapidity(), eff );
    }
  }
  void Finish() override {
    momentum_efficiency_->Write();
  }
private:
  enum FIELDS{
    HITS_TOF_RPC,
    PSI_RP
  };
  EventHeader* sim_header_;
  EventHeader* reco_header_;
  TrackDetector* sim_tracks_;
  TrackDetector* reco_tracks_;
  Matching* sim_reco_matching_;
  std::map<int, int> fields_id_;
  TProfile2D* momentum_efficiency_{nullptr};
};
} // namespace AnalysisTree
#endif // QUALITY_ASSURANCE_SRC_TREE_READER_H_
