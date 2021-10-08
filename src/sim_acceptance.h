//
// Created by mikhail on 6/29/20.
//

#ifndef EFFICIENCY_SRC_SIM_ACCEPTANCE_H_
#define EFFICIENCY_SRC_SIM_ACCEPTANCE_H_

#include <TChain.h>
#include <TProfile2D.h>

#include <AnalysisTree/EventHeader.hpp>
#include <AnalysisTree/FillTask.hpp>
#include <AnalysisTree/Cuts.hpp>
#include <AnalysisTree/Detector.hpp>
#include <AnalysisTree/DataHeader.hpp>

#include <TH3F.h>

namespace AnalysisTree {
class SimAcceptance : public FillTask {
public:
  void Init(std::map<std::string, void *> &branch_map) override;
  void Exec() override;
  void Finish() override;
  void SetPidCode(int pid_code);

private:
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
  enum FIELDS {
    HITS_TOF_RPC,
    SIM_GEANT_PID,
    IS_PRIMARY,
    PSI_RP
  };
  EventHeader *sim_header_{nullptr};
  EventHeader *reco_header_{nullptr};
  Particles *sim_tracks_{nullptr};
  Particles *reco_tracks_{nullptr};
  int pid_code_=2212;
  double y_beam_=0.0;

  std::map<int, int> fields_id_;
  TH3F* gen_tracks_prim_cent_;
  std::vector<TH2F *> gen_tracks_prim_;
  std::vector<TH2F *> gen_tracks_sec_;
};
} // namespace AnalysisTree
#endif // EFFICIENCY_SRC_SIM_ACCEPTANCE_H_
