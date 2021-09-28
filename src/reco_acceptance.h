//
// Created by mikhail on 6/16/20.
//

#ifndef QUALITY_ASSURANCE_SRC_TREE_READER_H_
#define QUALITY_ASSURANCE_SRC_TREE_READER_H_

#include <TChain.h>
#include <TH3F.h>
#include <TProfile2D.h>

#include <AnalysisTree/EventHeader.hpp>
#include <AnalysisTree/FillTask.hpp>
#include <AnalysisTree/Cuts.hpp>
#include <AnalysisTree/Detector.hpp>
#include <AnalysisTree/Matching.hpp>

namespace AnalysisTree {
class RecoAcceptance : public FillTask{
public:
  RecoAcceptance() = default;
  ~RecoAcceptance() override = default;
  void Init( std::map<std::string, void*>& branch_map ) override;
  void Exec() override;
  void Finish() override;
  void SetPidCode(int pid_code);

private:
  std::array<int, 6> CalcRecoSectorsOccupancy(int pid=-1);
  std::array<int, 6> CalcSimSectorsOccupancy(int pid=-1);
  size_t WhatSector( double phi ){
    if( 0.0 < phi && phi < M_PI/3.0 )
      return 0;
    if( M_PI/3.0 < phi && phi < 2.0*M_PI/3.0 )
      return 1;
    if( 2*M_PI/3.0 < phi && phi < M_PI )
      return 2;
    if( -M_PI < phi && phi < -2*M_PI/3.0 )
      return 3;
    if( -2*M_PI/3.0 < phi && phi < -M_PI/3.0 )
      return 4;
    if( -M_PI/3.0 < phi && phi < 0.0 )
      return 5;
    return -1;
  }
  enum FIELDS{
    HITS_TOF_RPC,
    SIM_GEANT_PID,
    RECO_GEANT_PID,
    IS_PRIMARY,
    PSI_RP,
    LAYERS_0,
    LAYERS_1,
    LAYERS_2,
    LAYERS_3,
    CHI2,
    DCA_XY,
    DCA_Z
  };
  EventHeader* reco_header_{nullptr};
  EventHeader* sim_header_{nullptr};
  Particles* sim_tracks_{nullptr};
  Particles* reco_tracks_{nullptr};
  HitDetector* meta_hits_{nullptr};
  Matching*reco_sim_matching_{nullptr};
  Matching* mdc_meta_matching_{nullptr};
  int pid_code_=2212;

  std::map<int, int> fields_id_;
  TProfile* momentum_err_{nullptr};

  std::vector<TH1F*> n_tracks_in_sector_;

  std::vector<TH2F *> rec_occupancy_;
  std::vector<TH2F*> pdg_tracks_prim_; // Gen-PID == Reco-PID, is_primary
  std::vector<TH2F*> pid_tracks_prim_; // Gen-PID == Reco-PID, !is_primary
  std::vector<TH2F*> pdg_tracks_sec_; // is_primary
  std::vector<TH2F*> pid_tracks_sec_; // !is_primary
  std::vector<TH2F*> pid_tracks_mismatch_; // Gen-PID != Reco-PID
  std::vector<TH2F*> pid_reco_; // pid_prim+pid_sec+pid_mismatch
};
} // namespace AnalysisTree
#endif // QUALITY_ASSURANCE_SRC_TREE_READER_H_
