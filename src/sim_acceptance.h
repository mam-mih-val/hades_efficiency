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

#include <centrality.h>

namespace AnalysisTree {
class SimAcceptance : public FillTask {
  void Init(std::map<std::string, void *> &branch_map) override;
  void Exec() override;
  void Finish() override;

private:
  enum FIELDS {
    HITS_TOF_RPC,
    SIM_GEANT_PID,
    IS_PRIMARY,
    PSI_RP
  };
  EventHeader *sim_header_{nullptr};
  EventHeader *reco_header_{nullptr};
  Particles *sim_tracks_{nullptr};

  std::map<int, int> fields_id_;
  std::vector<TH2F *> gen_tracks_prim_;
  std::vector<TH2F *> gen_tracks_sec_;
  std::vector<TH2F *> gen_prim_phi_pt_midrapidity_;
  std::vector<TH2F *> gen_prim_delta_phi_pt_midrapidity_;
};
} // namespace AnalysisTree
#endif // EFFICIENCY_SRC_SIM_ACCEPTANCE_H_
