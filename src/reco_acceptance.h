//
// Created by mikhail on 6/16/20.
//

#ifndef QUALITY_ASSURANCE_SRC_TREE_READER_H_
#define QUALITY_ASSURANCE_SRC_TREE_READER_H_

#include <AnalysisTree/EventHeader.hpp>
#include <AnalysisTree/FillTask.hpp>

#include <AnalysisTree/Cuts.hpp>
#include <AnalysisTree/Detector.hpp>
#include <TChain.h>
#include <TH3F.h>
#include <TProfile2D.h>

#include "centrality.h"

namespace AnalysisTree {
class RecoAcceptance : public FillTask{
public:
  RecoAcceptance() = default;
  ~RecoAcceptance() = default;
  void Init( std::map<std::string, void*>& branch_map ) override {
    config_->Print();
    sim_header_ = static_cast<EventHeader*>( branch_map.at( "sim_header" ) );
    reco_header_ = static_cast<EventHeader*>(branch_map.at( "event_header" ));
    sim_tracks_ = static_cast<Particles *>(branch_map.at( "sim_tracks" ));
    reco_tracks_ = static_cast<Particles *>(branch_map.at( "mdc_vtx_tracks" ));
    meta_hits_ = static_cast<HitDetector *>(branch_map.at( "meta_hits" ));
    sim_reco_matching_ = static_cast<Matching *>(branch_map.at( "sim_reco_tracks" ));
    mdc_meta_matching_ = static_cast<Matching *>(branch_map.at( "mdc_meta_match" ));
    auto sim_event_config = config_->GetBranchConfig("sim_header");
    auto reco_event_config = config_->GetBranchConfig("event_header");
    auto sim_tracks_config = config_->GetBranchConfig("sim_tracks");
    auto reco_tracks_config = config_->GetBranchConfig("mdc_vtx_tracks");
    auto meta_hits_config = config_->GetBranchConfig("meta_hits");

    fields_id_.insert(  std::make_pair( HITS_TOF_RPC, reco_event_config.GetFieldId("selected_tof_rpc_hits") )  );
    fields_id_.insert(  std::make_pair( PSI_RP, sim_event_config.GetFieldId("reaction_plane") )  );
    fields_id_.insert(  std::make_pair(
        SIM_GEANT_PID, sim_tracks_config.GetFieldId("geant_pid") )  );
    fields_id_.insert(  std::make_pair(
        IS_PRIMARY, sim_tracks_config.GetFieldId("is_primary") )  );
    fields_id_.insert(  std::make_pair(
        RECO_GEANT_PID, reco_tracks_config.GetFieldId("geant_pid") )  );
    fields_id_.insert(  std::make_pair(
        DCA_XY, reco_tracks_config.GetFieldId("dca_xy") )  );
    fields_id_.insert(  std::make_pair(
        DCA_Z, reco_tracks_config.GetFieldId("dca_z") )  );
    fields_id_.insert(  std::make_pair(META_MASS2, meta_hits_config.GetFieldId("mass2") )  );

    for( int i=0; i<8; ++i ) {
      int percentile = 2+i*5;
      std::string name = "pdg_acceptance_prim_"+std::to_string(percentile);
      pdg_acceptance_prim_.push_back(new TH2F(name.data(),";p_{T}, [GeV/c]; y-y_{beam}; conuts",
                                               100, 0.0, 2.0,
                                               100, -1.0, 1.0));
      name = "pdg_acceptance_sec_"+std::to_string(percentile);
      pdg_acceptance_sec_.push_back(new TH2F(name.data(),";p_{T}, [GeV/c]; y-y_{beam}; conuts",
                                              100, 0.0, 2.0,
                                              100, -1.0, 1.0));

      name = "pid_acceptance_prim_"+std::to_string(percentile);
      pid_acceptance_prim_.push_back(new TH2F(name.data(),";p_{T}, [GeV/c]; y-y_{beam}; conuts",
                                               100, 0.0, 2.0,
                                               100, -1.0, 1.0));
      name = "pid_acceptance_sec_"+std::to_string(percentile);
      pid_acceptance_sec_.push_back(new TH2F(name.data(),";p_{T}, [GeV/c]; y-y_{beam}; conuts",
                                              100, 0.0, 2.0,
                                              100, -1.0, 1.0));

    }
    momentum_err_ = new TProfile( "momentum_err", ";p, [GeV/c]; relative error", 100, 0.0, 3.5 );
    proton_yield_ = new TH1F( "pid_proton_yield", ";N protons; counts", 100, 0, 100 );
    mass_mismatched_ = new TH1F( "mass_mismatched", "mass of mismatched;m, [GeV/c^{2}]", 500, 0.0, 2.5 );
    mass_matched_ = new TH1F( "mass_matched", "mass of matched;m, [GeV/c^{2}]", 500, 0.0, 2.5 );
    mass_all_ = new TH1F( "mass_all", "mass of all;m, [GeV/c^{2}]", 500, 0.0, 2.5 );
  }
  void Exec() override {
    auto hits_tof_rpc = reco_header_->GetField<int>(fields_id_.at(HITS_TOF_RPC));
    int centrality_class = Centrality::GetInstance()->GetCentralityClass5pc(hits_tof_rpc);
    if( centrality_class > 7 ) {
      return;
    }
    int n_reco_tracks = reco_tracks_->GetNumberOfChannels();
    int n_protons{0};
    for (int i = 0; i < n_reco_tracks; ++i) {
      int sim_id = sim_reco_matching_->GetMatchDirect(i);
      int meta_id = mdc_meta_matching_->GetMatchDirect(i);
      auto r_track = reco_tracks_->GetChannel(i);
      auto r_hit = meta_hits_->GetChannel(meta_id);
      auto s_track = sim_tracks_->GetChannel( (sim_id) );
      float m_reco = r_track.GetMass();
      float m_sim = s_track.GetMass();

      auto p_reco = r_track.Get4MomentumByMass(m_reco);
      auto p_sim = s_track.Get4MomentumByMass(m_sim);
//      if( fabsf(r_track.GetField<float>( fields_id_.at(DCA_XY) )) > 10.0f )
//        continue;
//      if( fabsf(r_track.GetField<float>( fields_id_.at(DCA_Z) ) ) > 10.0f )
//        continue;
      if( s_track.GetField<int>( fields_id_.at(SIM_GEANT_PID) ) == 14 &&
          r_track.GetField<int>( fields_id_.at(RECO_GEANT_PID) ) == 14 ){
        if( s_track.GetField<bool>( fields_id_.at(IS_PRIMARY) ) )
          pdg_acceptance_prim_.at(centrality_class)->Fill( p_reco.Pt(), p_reco.Rapidity() - 0.74 );
        if( !s_track.GetField<bool>( fields_id_.at(IS_PRIMARY) ) )
          pdg_acceptance_sec_.at(centrality_class)->Fill( p_reco.Pt(), p_reco.Rapidity() - 0.74 );
        double eff = fabs(p_reco.P() - p_sim.P())/p_sim.P();
        momentum_err_->Fill(p_sim.P(), eff);
      }

      if( r_track.GetField<int>( fields_id_.at(RECO_GEANT_PID) ) == 14 ){
        if( s_track.GetField<bool>( fields_id_.at(IS_PRIMARY) ) )
          pid_acceptance_prim_.at(centrality_class)->Fill( p_reco.Pt(), p_reco.Rapidity() - 0.74 );
        if( !s_track.GetField<bool>( fields_id_.at(IS_PRIMARY) ) )
          pid_acceptance_sec_.at(centrality_class)->Fill( p_reco.Pt(), p_reco.Rapidity() - 0.74 );
        n_protons++;
      }
      if( s_track.GetField<int>( fields_id_.at(SIM_GEANT_PID) ) != 14 )
        continue;
      auto mass = sqrtf(r_hit.GetField<float>(fields_id_.at(META_MASS2)));
      if(r_track.GetField<int>( fields_id_.at(RECO_GEANT_PID) ) != 14 ){
        mass_mismatched_->Fill( mass );
      }if(r_track.GetField<int>( fields_id_.at(RECO_GEANT_PID) ) == 14 ){
        mass_matched_->Fill( mass );
      }
      mass_all_->Fill(mass);
    }
    proton_yield_->Fill(n_protons);
  }
  void Finish() override {
    momentum_err_->Write();
    proton_yield_->Write();
    mass_mismatched_->Write();
    mass_matched_->Write();
    mass_all_->Write();
    for( size_t i=0; i< pdg_acceptance_prim_.size(); ++i){
      pdg_acceptance_prim_.at(i)->Write();
      pdg_acceptance_sec_.at(i)->Write();

      pid_acceptance_prim_.at(i)->Write();
      pid_acceptance_sec_.at(i)->Write();
    }
  }
private:
  enum FIELDS{
    HITS_TOF_RPC,
    PSI_RP,
    SIM_GEANT_PID,
    RECO_GEANT_PID,
    IS_PRIMARY,
    DCA_XY,
    DCA_Z,
    META_MASS2,
  };
  EventHeader* sim_header_{nullptr};
  EventHeader* reco_header_{nullptr};
  Particles* sim_tracks_{nullptr};
  Particles* reco_tracks_{nullptr};
  HitDetector* meta_hits_{nullptr};
  Matching* sim_reco_matching_{nullptr};
  Matching*mdc_meta_matching_{nullptr};

  std::map<int, int> fields_id_;
  TProfile* momentum_err_{nullptr};
  TH1F* proton_yield_{nullptr};
  TH1F*mass_mismatched_{nullptr};
  TH1F*mass_matched_{nullptr};
  TH1F*mass_all_{nullptr};
  std::vector<TH2F*> pdg_acceptance_prim_;
  std::vector<TH2F*> pid_acceptance_prim_;
  std::vector<TH2F*> pdg_acceptance_sec_;
  std::vector<TH2F*> pid_acceptance_sec_;
};
} // namespace AnalysisTree
#endif // QUALITY_ASSURANCE_SRC_TREE_READER_H_
