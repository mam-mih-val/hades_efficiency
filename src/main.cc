//
// Created by mikhail on 6/16/20.
//

#include "momentum_reconstruction_eff.h"
#include <AnalysisTree/TaskManager.h>
#include <iostream>
int main(int n_args, char** args){
  if(n_args<2){
    std::cout << "Error: missing file list operand" << std::endl;
    std::cout << "Please use: ./run file.list" << std::endl;
  }
  std::string list{args[1]};
  AnalysisTree::TaskManager manager(list, "hades_analysis_tree");
  auto * mom_rec_eff = new AnalysisTree::MomentumReconctructionEff;
  manager.AddTask(reinterpret_cast<AnalysisTree::FillTask*>(mom_rec_eff));
  manager.SetOutFileName("out.root");
  manager.Init();
  manager.Run(-1);
  manager.Finish();
  return 0;
}