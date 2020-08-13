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

#include <centrality.h>

namespace AnalysisTree {
class RecoAcceptance : public FillTask{
public:
  RecoAcceptance() = default;
  ~RecoAcceptance() override = default;
  void Init( std::map<std::string, void*>& branch_map ) override;
  void Exec() override;
  void Finish() override;
private:
  enum FIELDS{
    HITS_TOF_RPC,
    SIM_GEANT_PID,
    RECO_GEANT_PID,
    IS_PRIMARY,
  };
  EventHeader* reco_header_{nullptr};
  Particles* sim_tracks_{nullptr};
  Particles* reco_tracks_{nullptr};
  HitDetector* meta_hits_{nullptr};
  Matching* sim_reco_matching_{nullptr};
  Matching*mdc_meta_matching_{nullptr};

  std::map<int, int> fields_id_;
  TProfile* momentum_err_{nullptr};

  std::vector<TH2F*> pdg_tracks_prim_; // Gen-PID == Reco-PID, is_primary
  std::vector<TH2F*> pid_tracks_prim_; // Gen-PID == Reco-PID, !is_primary
  std::vector<TH2F*> pdg_tracks_sec_; // is_primary
  std::vector<TH2F*> pid_tracks_sec_; // !is_primary
  std::vector<TH2F*> pid_tracks_mismatch_; // Gen-PID != Reco-PID
  std::vector<TH2F*> pid_reco_; // pid_prim+pid_sec+pid_mismatch
};
} // namespace AnalysisTree
#endif // QUALITY_ASSURANCE_SRC_TREE_READER_H_
