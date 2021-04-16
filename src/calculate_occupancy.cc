//
// Created by mikhail on 1/19/21.
//
#include <iostream>

#include <TFile.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TF1.h>

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
  TProfile* mult_sector{nullptr};
  file_in->GetObject("det_tracks_sector_vs_centrality", mult_sector);
  std::vector<TH2F*> efficiencies;
  int percentile = 2;
  TH2F* histo_ptr{nullptr};
  while( percentile < 40 ){
    std::string histo_name{ "efficiency_"+std::to_string(percentile) };
    file_in->GetObject(histo_name.c_str(), histo_ptr);
    efficiencies.emplace_back(histo_ptr);
    percentile+=5;
  }
  std::vector<TGraphErrors*> centrality_dependencies;
  float y_axis[16];
  for(int j=0; j<16; ++j){ y_axis[j]=-0.75f+0.1f* (float) j; }
  float pt_axis[]={0, 0.29375, 0.35625, 0.41875, 0.48125, 0.54375, 0.61875, 0.70625, 0.81875, 1.01875, 2.0};
  auto histo_par0 = new TH2F("par0", ";y_{cm};p_{T} [GeV/c];slope",
                          15, y_axis, 10, pt_axis );
  auto histo_par1 = new TH2F("par1", ";y_{cm};p_{T} [GeV/c];offset",
                          15, y_axis, 10, pt_axis );
  auto histo_par2 = new TH2F("par2", ";y_{cm};p_{T} [GeV/c];offset",
                          15, y_axis, 10, pt_axis );
  auto histo_par3 = new TH2F("par3", ";y_{cm};p_{T} [GeV/c];offset",
                          15, y_axis, 10, pt_axis );
  TF1* fit_function;
  for( int y_bin=1; y_bin <= 15; y_bin++ ){
    for (int pT_bin = 1; pT_bin <= 10; ++pT_bin) {
      auto y = efficiencies.front()->GetXaxis()->GetBinCenter(y_bin);
      auto pT = efficiencies.front()->GetYaxis()->GetBinCenter(pT_bin);
      centrality_dependencies.push_back(new TGraphErrors(8));
      std::string graph_name{
          "y="+std::to_string(y)+"_pT="+std::to_string(pT) };
      centrality_dependencies.back()->SetName(graph_name.c_str());
      centrality_dependencies.back()->SetTitle(graph_name.c_str());
      graph_name="fit_"+graph_name;
      fit_function = new TF1( graph_name.c_str(), "pol3", 0.0, 30.0 );
      for (int cc = 0; cc < 8; ++cc) {
        auto efficiency = efficiencies.at(cc);
        auto Ntr = mult_sector->GetBinContent(cc+1);
        auto e = efficiency->GetBinContent(y_bin, pT_bin);
        auto e_err = efficiency->GetBinError(y_bin, pT_bin);
        centrality_dependencies.back()->SetPoint(cc, Ntr, e);
        centrality_dependencies.back()->SetPointError(cc, 0.0, e_err);
      }
      centrality_dependencies.back()->Fit(fit_function);
      auto par0 = fit_function->GetParameter(0);
      auto par0_err = fit_function->GetParError(0);
      histo_par0->SetBinContent(y_bin, pT_bin, par0);
      histo_par0->SetBinError(y_bin, pT_bin, par0_err);

      auto par1 = fit_function->GetParameter(1);
      auto par1_err = fit_function->GetParError(1);
      histo_par1->SetBinContent(y_bin, pT_bin, par1);
      histo_par1->SetBinError(y_bin, pT_bin, par1_err);

      auto par2 = fit_function->GetParameter(2);
      auto par2_err = fit_function->GetParError(2);
      histo_par2->SetBinContent(y_bin, pT_bin, par2);
      histo_par2->SetBinError(y_bin, pT_bin, par2_err);

      auto par3 = fit_function->GetParameter(3);
      auto par3_err = fit_function->GetParError(3);
      histo_par3->SetBinContent(y_bin, pT_bin, par3);
      histo_par3->SetBinError(y_bin, pT_bin, par3_err);
    }
  }
  auto file_out = TFile::Open("occupancy.root", "recreate");
  file_out->mkdir("fits");
  file_out->cd("fits");
  for( auto graph : centrality_dependencies ){
    graph->Write();
  }
  file_out->cd("/");
  histo_par0->Write();
  histo_par1->Write();
  histo_par2->Write();
  histo_par3->Write();
  file_out->Close();
  return 0;
}