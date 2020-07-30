//
// Created by mikhail on 6/29/20.
//

#include <TFile.h>
#include <TH2F.h>
#include <iostream>
int main(int n, char** args){
  if( n < 3 ){
    std::cout << "Aim not specified" << std::endl;
    return 1;
  }
  std::string  file_name{ args[1] };
  auto *file_in = TFile::Open( file_name.data() );
  if( !file_in){
    std::cout << "No such file" << std::endl;
    return 2;
  }
  // detector efficiency
  std::vector<TH2F*> eff_prim;
  std::vector<TH2F*> eff_sec;
  // particle identification effects
  std::vector<TH2F*> pid_mismatch_prim;
  std::vector<TH2F*> pid_mismatch_sec;

  TH2F* gen_acc_prim;
  TH2F* gen_acc_sec;
  TH2F* pdg_acc_prim;
  TH2F* pdg_acc_sec;
  TH2F* pid_acc_prim;
  TH2F* pid_acc_sec;

  int percentile = 2;
  while( percentile < 40 ){
    std::string name = "gen_acceptance_prim_" + std::to_string(percentile);
    file_in->GetObject( name.data(), gen_acc_prim );
    name = "gen_acceptance_sec_" + std::to_string(percentile);
    file_in->GetObject( name.data(), gen_acc_sec );

    name = "pdg_acceptance_prim_" + std::to_string(percentile);
    file_in->GetObject( name.data(), pdg_acc_prim );
    name = "pdg_acceptance_sec_" + std::to_string(percentile);
    file_in->GetObject( name.data(), pdg_acc_sec );

    name = "pid_acceptance_prim_" + std::to_string(percentile);
    file_in->GetObject( name.data(), pid_acc_prim );
    name = "pid_acceptance_sec_" + std::to_string(percentile);
    file_in->GetObject( name.data(), pid_acc_sec );
    pid_acc_prim->Divide( pdg_acc_prim );
    pid_mismatch_prim.emplace_back( pid_acc_prim );
    pid_acc_sec->Divide( pdg_acc_sec );
    pid_mismatch_sec.emplace_back( pid_acc_sec );

    pdg_acc_prim->Divide( gen_acc_prim );
    pdg_acc_sec->Divide( gen_acc_sec );
    eff_prim.emplace_back(pdg_acc_prim);
    eff_sec.emplace_back(pdg_acc_sec);

    percentile+=5;
  }

  auto* file_out = TFile::Open( "efficiency_out.root", "recreate" );
  file_out->cd();
  percentile=2;
  int i=0;
  while (percentile<40){
    std::string name = "efficiency_prim_" + std::to_string(percentile);
    eff_prim.at(i)->Write(name.data());

    name = "efficiency_sec_" + std::to_string(percentile);
    eff_sec.at(i)->Write(name.data());

    name = "pid_mismatch_prim_" + std::to_string(percentile);
    pid_mismatch_prim.at(i)->Write(name.data());

    name = "pid_mismatch_sec_" + std::to_string(percentile);
    pid_mismatch_sec.at(i)->Write(name.data());

    i++;
    percentile+=5;
  }
  file_in->Close();
  file_out->Close();
  return 0;
}