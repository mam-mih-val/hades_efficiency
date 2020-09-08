//
// Created by mikhail on 6/16/20.
//

#include <iostream>
#include <chrono>

#include <AnalysisTree/TaskManager.hpp>

#include <cuts.h>

#include "reco_acceptance.h"
#include "sim_acceptance.h"


int main(int n_args, char** args){
  if(n_args<2){
    std::cout << "Error: missing file list operand" << std::endl;
    std::cout << "Please use: ./acceptance file.list" << std::endl;
    return 1;
  }
  std::string list{args[1]};
  AnalysisTree::TaskManager manager({list}, {"hades_analysis_tree"});

  auto *reco_acceptance = new AnalysisTree::RecoAcceptance;
  auto *sim_acceptance = new AnalysisTree::SimAcceptance;

  manager.SetEventCuts(HadesUtils::Cuts::Get(HadesUtils::Cuts::BRANCH_TYPE::EVENT_HEADER,
                                             HadesUtils::DATA_TYPE::AuAu_1_23AGeV));
  manager.AddBranchCut(HadesUtils::Cuts::Get(HadesUtils::Cuts::BRANCH_TYPE::MDC_TRACKS,
                                             HadesUtils::DATA_TYPE::AuAu_1_23AGeV));
  manager.AddBranchCut(HadesUtils::Cuts::Get(HadesUtils::Cuts::BRANCH_TYPE::META_HITS,
                                             HadesUtils::DATA_TYPE::AuAu_1_23AGeV));

  manager.AddTask(reco_acceptance);
  manager.AddTask(sim_acceptance);
  manager.SetOutFileName("out.root");
  manager.Init();
  manager.Run(100000);
  manager.Finish();
  return 0;
}