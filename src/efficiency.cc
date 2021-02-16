//
// Created by mikhail on 6/29/20.
//

#include <TFile.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <iostream>

int main(int n, char** args){
  if( n < 1 ){
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
  std::vector<TH2F*> efficiency;
  std::vector<TH2F*> midrapidity_pt_phi_efficiency;
  std::vector<TH2F*> midrapidity_pt_delta_phi_efficiency;
  std::vector<TH2F*> contamination;
  std::vector<TH2F*> mismatch;

  TH2F* gen_prim;
  TH2F* gen_prim_pt_phi;
  TH2F* gen_prim_pt_delta_phi;
  TH3F* gen_prim_pt_y_n_tracks;
  TH2F* gen_sec;
  TH2F* pdg_prim;
  TH2F* pdg_prim_pt_phi;
  TH2F* pdg_prim_pt_delta_phi;
  TH3F* pdg_prim_pt_y_n_tracks;
  TH2F* pdg_sec;
  TH2F* pid_prim;
  TH2F* pid_sec;
  TH2F* pid_reco;
  TH2F* pid_mismatch;
  TProfile* n_tracks_sector;
  file_in->GetObject( "det_tracks_sector_vs_centrality", n_tracks_sector );
  int percentile = 2;
  while( percentile < 40 ){
    std::string name = "gen_tracks_prim_" + std::to_string(percentile);
    file_in->GetObject( name.data(), gen_prim);
    gen_prim->Sumw2();
    name = "gen_tracks_sec_" + std::to_string(percentile);
    file_in->GetObject( name.data(), gen_sec);

    name = "pdg_tracks_prim_" + std::to_string(percentile);
    file_in->GetObject( name.data(), pdg_prim);
    pdg_prim->Sumw2();
    name = "pdg_tracks_sec_" + std::to_string(percentile);
    file_in->GetObject( name.data(), pdg_sec);

    name = "pid_tracks_prim_" + std::to_string(percentile);
    file_in->GetObject( name.data(), pid_prim);
    name = "pid_tracks_sec_" + std::to_string(percentile);
    file_in->GetObject( name.data(), pid_sec);

    name = "pid_reco_" + std::to_string(percentile);
    file_in->GetObject( name.data(), pid_reco );
    name = "pid_tracks_mismatch_" + std::to_string(percentile);
    file_in->GetObject( name.data(), pid_mismatch );

    name = "pdg_prim_phi_pt_midrapidity_" + std::to_string(percentile);
    file_in->GetObject( name.data(), pdg_prim_pt_phi );

    name = "pid_prim_delta_phi_pt_midrapidity_" + std::to_string(percentile);
    file_in->GetObject( name.data(), pdg_prim_pt_delta_phi );

    name = "gen_prim_phi_pt_midrapidity_" + std::to_string(percentile);
    file_in->GetObject( name.data(), gen_prim_pt_phi );

    name = "gen_prim_delta_phi_pt_midrapidity_" + std::to_string(percentile);
    file_in->GetObject( name.data(), gen_prim_pt_delta_phi );

    pdg_sec->Add(pid_mismatch);
    pdg_sec->Divide(pid_reco);

    pdg_prim->Divide(gen_prim);
    pid_mismatch->Divide(pid_reco);

    pdg_prim_pt_phi->Divide(gen_prim_pt_phi);
//    pdg_prim_pt_delta_phi->Rebin2D(5, 4);
//    gen_prim_pt_delta_phi->Rebin2D(5, 4);
//    pdg_prim_pt_delta_phi->Divide(gen_prim_pt_delta_phi);

    efficiency.emplace_back(pdg_prim);
    contamination.emplace_back(pdg_sec);
    mismatch.emplace_back(pid_mismatch);
    midrapidity_pt_phi_efficiency.emplace_back(pdg_prim_pt_phi);
    midrapidity_pt_delta_phi_efficiency.emplace_back(pdg_prim_pt_delta_phi);
    percentile+=5;
  }
  file_in->GetObject( "gen_prim_pT_y_n_tracks_sector", gen_prim_pt_y_n_tracks );
  file_in->GetObject( "pdg_prim_pT_y_n_tracks_sector", pdg_prim_pt_y_n_tracks );

  pdg_prim_pt_y_n_tracks->Sumw2();
  gen_prim_pt_y_n_tracks->Sumw2();
  pdg_prim_pt_y_n_tracks->Divide( gen_prim_pt_y_n_tracks );
  auto* file_out = TFile::Open( "efficiency_protons.root", "recreate" );
  file_out->cd();
  pdg_prim_pt_y_n_tracks->Write("efficiency_pT_y_n_tacks_sector");
//  n_tracks_sector->Write();
  percentile=2;
  int i=0;
  while (percentile<40){
    std::string name = "efficiency_" + std::to_string(percentile);
    efficiency.at(i)->Write(name.data());

    name = "contamination_" + std::to_string(percentile);
    contamination.at(i)->Write(name.data());

    name = "pid_mismatch_" + std::to_string(percentile);
    mismatch.at(i)->Write(name.data());

    name = "pt_phi_eff_midrapidity_"+std::to_string(percentile);
//    midrapidity_pt_phi_efficiency.at(i)->Write(name.data());

    name = "pt_delta_phi_eff_midrapidity_"+std::to_string(percentile);
//    midrapidity_pt_delta_phi_efficiency.at(i)->Write(name.data());

    i++;
    percentile+=5;
  }
  file_in->Close();
  file_out->Close();
  return 0;
}