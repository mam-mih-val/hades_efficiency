//
// Created by mikhail on 6/16/20.
//

#include <iostream>
#include <chrono>

#include "reco_acceptance.h"
#include "sim_acceptance.h"
#include <AnalysisTree/TaskManager.h>
int main(int n_args, char** args){
  if(n_args<2){
    std::cout << "Error: missing file list operand" << std::endl;
    std::cout << "Please use: ./run file.list" << std::endl;
    return 1;
  }
  std::string list{args[1]};
  AnalysisTree::TaskManager manager(list, "hades_analysis_tree");
  auto * mom_rec_eff = new AnalysisTree::RecoAcceptance;
  auto * mom_sim_eff = new AnalysisTree::SimAcceptance;
  manager.AddTask(mom_rec_eff);
  manager.AddTask(mom_sim_eff);
  manager.SetOutFileName("out.root");
  manager.Init();
  auto start = std::chrono::system_clock::now();
  manager.Run(-1);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> duration = end-start;
  manager.Finish();
  std::cout << "Elapsed time: " << duration.count() << " s" << std::endl;
  return 0;
}