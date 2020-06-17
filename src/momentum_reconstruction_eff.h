//
// Created by mikhail on 6/16/20.
//

#ifndef QUALITY_ASSURANCE_SRC_TREE_READER_H_
#define QUALITY_ASSURANCE_SRC_TREE_READER_H_

#include <AnalysisTree/EventHeader.h>
#include <AnalysisTree/FillTask.h>

#include <AnalysisTree/Cuts.h>
#include <TChain.h>
#include <TH3F.h>
#include <TProfile2D.h>

#include "centrality.h"

namespace AnalysisTree {
class MomentumReconctructionEff : public FillTask{
public:
  MomentumReconctructionEff() = default;
  ~MomentumReconctructionEff() = default;
  void Init( std::map<std::string, void*>& branch_map ) override {
    config_->Print();
    sim_header_ = static_cast<EventHeader*>( branch_map.at( "sim_header" ) );
    reco_header_ = static_cast<EventHeader*>(branch_map.at( "event_header" ));
    sim_tracks_ = static_cast<Particles *>(branch_map.at( "sim_tracks" ));
    reco_tracks_ = static_cast<Particles *>(branch_map.at( "mdc_vtx_tracks" ));
    sim_reco_matching_ = static_cast<Matching *>(branch_map.at( "sim_reco_tracks" ));
    auto sim_event_config = config_->GetBranchConfig("sim_header");
    auto reco_event_config = config_->GetBranchConfig("event_header");
    auto sim_tracks_config = config_->GetBranchConfig("sim_tracks");
    auto reco_tracks_config = config_->GetBranchConfig("mdc_vtx_tracks");

    fields_id_.insert(  std::make_pair( HITS_TOF_RPC, reco_event_config.GetFieldId("selected_tof_rpc_hits") )  );
    fields_id_.insert(  std::make_pair( PSI_RP, sim_event_config.GetFieldId("reaction_plane") )  );
    fields_id_.insert(  std::make_pair(
        SIM_GEANT_PID, sim_tracks_config.GetFieldId("geant_pid") )  );
    fields_id_.insert(  std::make_pair(
        RECO_GEANT_PID, reco_tracks_config.GetFieldId("geant_pid") )  );

    for( int i=0; i<4; ++i ) {
      int percentile = 5+i*10;
      std::string name = "reco_counts_"+std::to_string(percentile);
      reco_momentum_counts_.push_back(new TH3F(name.data(),";p_{T}, [GeV/c]; y; #phi, [rad];conuts",
                                               100, 0.0, 2.5,
                                               100, 0.0, 1.6,
                                               100, -3.15, 3.15));
      name = "sim_counts_"+std::to_string(percentile);
      sim_momentum_counts_.push_back(new TH3F(name.data(),";p_{T}, [GeV/c]; y; #phi, [rad];conuts",
                                              100, 0.0, 2.5,
                                              100, 0.0, 1.6,
                                              100, -3.15, 3.15));
      name = "sim_density_"+std::to_string(percentile);
      sim_track_density_.push_back(new TH2F(name.data(),";#Delta#phi, [rad]; #theta, [rad]; conuts",
                                               100, -3.15, 3.15, 100, 0.0, 1.7));
      name = "reco_density_"+std::to_string(percentile);
      reco_track_density_.push_back(new TH2F(name.data(),";#Delta#phi, [rad]; #theta, [rad]; conuts",
                                               100, -3.15, 3.15, 100, 0.0, 1.7));

    }
    momentum_err_ = new TProfile( "momentum_err", ";p, [GeV/c]; relative error", 100, 0.0, 3.5 );

  }
  void Exec() override {
    auto hits_tof_rpc = reco_header_->GetField<int>(fields_id_.at(HITS_TOF_RPC));
    int centrality_class = Centrality::GetInstance()->GetCentralityClass10pc(hits_tof_rpc);
    if( centrality_class > 3 ) {
      return;
    }
    auto psi_rp = sim_header_->GetField<float>(fields_id_.at(PSI_RP));
    int n_reco_tracks = reco_tracks_->GetNumberOfChannels();
    for (int i = 0; i < n_reco_tracks; ++i) {
      int sim_id = sim_reco_matching_->GetMatchDirect(i);
      auto r_track = (reco_tracks_->GetChannel(i));
      auto s_track = (sim_tracks_->GetChannel( (sim_id) ));
      if( r_track.GetField<int>( fields_id_.at(RECO_GEANT_PID) ) != 14 )
        continue;
      float m_reco = r_track.GetMass();
      float m_sim = s_track.GetMass();

      auto p_reco = r_track.Get4MomentumByMass(m_reco);
      auto p_sim = s_track.Get4MomentumByMass(m_sim);

      double eff = fabs(p_reco.P() - p_sim.P())/p_sim.P();
      momentum_err_->Fill(p_sim.P(), eff);
      reco_momentum_counts_.at(centrality_class)->Fill( p_reco.Pt(), p_reco.Rapidity(), p_reco.Phi() );
      float delta_phi = p_reco.Phi() - psi_rp;
      if( delta_phi < - TMath::Pi() )
        delta_phi += 2*TMath::Pi();
      if( delta_phi >  TMath::Pi() )
        delta_phi -= 2*TMath::Pi();
      reco_track_density_.at(centrality_class)->Fill(delta_phi, p_reco.Theta());
    }
    int n_sim_tracks = sim_tracks_->GetNumberOfChannels();
    for( int i=0; i<n_sim_tracks; ++i ){
      auto s_track = (sim_tracks_->GetChannel( i ));
      if( s_track.GetField<int>( fields_id_.at(SIM_GEANT_PID) ) != 14 )
        continue;
      float m_sim = s_track.GetMass();
      auto p_sim = s_track.Get4MomentumByMass(m_sim);
      sim_momentum_counts_.at(centrality_class)->Fill( p_sim.Pt(), p_sim.Rapidity(), p_sim.Phi() );
      float delta_phi = p_sim.Phi() - psi_rp;
      if( delta_phi < - TMath::Pi() )
        delta_phi += 2*TMath::Pi();
      if( delta_phi >  TMath::Pi() )
        delta_phi -= 2*TMath::Pi();
      sim_track_density_.at(centrality_class)->Fill(delta_phi, p_sim.Theta());
    }
  }
  void Finish() override {
    momentum_err_->Write();
    for( size_t i=0; i<reco_momentum_counts_.size(); ++i){
      reco_momentum_counts_.at(i)->Write();
      sim_momentum_counts_.at(i)->Write();
//      reco_momentum_counts_.at(i)->Divide(sim_momentum_counts_.at(i));
//      std::string name = "efficiency_"+std::to_string(5+i*10);
//      reco_momentum_counts_.at(i)->Write(name.data());

      sim_track_density_.at(i)->Write();
      reco_track_density_.at(i)->Write();
    }
  }
private:
  enum FIELDS{
    HITS_TOF_RPC,
    PSI_RP,
    SIM_GEANT_PID,
    RECO_GEANT_PID
  };
  EventHeader* sim_header_{nullptr};
  EventHeader* reco_header_{nullptr};
  Particles* sim_tracks_{nullptr};
  Particles* reco_tracks_{nullptr};
  Matching* sim_reco_matching_{nullptr};

  std::map<int, int> fields_id_;
  TProfile* momentum_err_{nullptr};
  std::vector<TH3F*> sim_momentum_counts_;
  std::vector<TH3F*> reco_momentum_counts_;
  std::vector<TH2F*> sim_track_density_;
  std::vector<TH2F*> reco_track_density_;
};
} // namespace AnalysisTree
#endif // QUALITY_ASSURANCE_SRC_TREE_READER_H_
