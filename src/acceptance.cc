//
// Created by mikhail on 6/16/20.
//

#include <iostream>
#include <chrono>

#include <AnalysisTree/TaskManager.hpp>

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
  int pid_code = 2212;
  std::string centrality_parameter_file_name;
  TFile* centrality_parameter_file{nullptr};
  TH1* centrality_parameter_histogram{nullptr};
  if( n_args > 2 ) {
    pid_code=atoi(args[2]);
    std::cout << "Found PID-code " << pid_code << std::endl;
  }
  if( n_args > 3 ) {
    centrality_parameter_file_name = args[3];
    std::cout << "Found centrality parameter file " << centrality_parameter_file_name << std::endl;
    centrality_parameter_file = TFile::Open( centrality_parameter_file_name.c_str() );
    if( centrality_parameter_file )
      centrality_parameter_file->GetObject( "Centrality/TOFRPC_5pc_fixedCuts", centrality_parameter_histogram );
    else
      std::cerr << "Couldn't open centrality parameter file " << centrality_parameter_file_name << std::endl;
  }
  auto *reco_acceptance = new AnalysisTree::RecoAcceptance;
  reco_acceptance->SetPidCode(pid_code);
  auto *sim_acceptance = new AnalysisTree::SimAcceptance;
  sim_acceptance->SetPidCode(pid_code);
  if( centrality_parameter_histogram ){
    reco_acceptance->SetCentralityParameters(centrality_parameter_histogram);
    sim_acceptance->SetCentralityParameters(centrality_parameter_histogram);
  }
  manager.AddTask(reco_acceptance);
  manager.AddTask(sim_acceptance);
//  manager.SetEventCuts(HadesUtils::Cuts::Get(HadesUtils::Cuts::BRANCH_TYPE::EVENT_HEADER,
//                                             HadesUtils::DATA_TYPE::AuAu_1_23AGeV));
//  manager.AddBranchCut(HadesUtils::Cuts::Get(HadesUtils::Cuts::BRANCH_TYPE::MDC_TRACKS,
//                                             HadesUtils::DATA_TYPE::AuAu_1_23AGeV));
//  manager.AddBranchCut(HadesUtils::Cuts::Get(HadesUtils::Cuts::BRANCH_TYPE::META_HITS,
//                                             HadesUtils::DATA_TYPE::AuAu_1_23AGeV));
  manager.SetOutFileName("out.root");
  manager.Init();
  manager.Run(-1);
  manager.Finish();
  return 0;
}