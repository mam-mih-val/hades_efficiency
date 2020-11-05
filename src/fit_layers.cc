//
// Created by mikhail on 9/27/20.
//

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TH3F.h>
#include <string>

int main(int n, char** args) {
  std::string input_file=args[1];
  auto file = TFile::Open(input_file.c_str());
  TH3F* histogram_layers_pt_phi;
  std::vector<TH2F*> matrix_first_harm{
      new TH2F( "v1_low_pT", ";in;out", 9, 0.0, 9.0, 9, 0.0, 9.0 ),
      new TH2F( "v1_mid_pT", ";in;out", 9, 0.0, 9.0, 9, 0.0, 9.0 ),
      new TH2F( "v1_high_pT", ";in;out", 9, 0.0, 9.0, 9, 0.0, 9.0 ),
  };
  std::vector<TH2F*> matrix_second_harm{
      new TH2F( "v2_low_pT", ";in;out", 9, 0.0, 9.0, 9, 0.0, 9.0 ),
      new TH2F( "v2_mid_pT", ";in;out", 9, 0.0, 9.0, 9, 0.0, 9.0 ),
      new TH2F( "v2_high_pT", ";in;out", 9, 0.0, 9.0, 9, 0.0, 9.0 ),
  };
  std::map<int, int> map_combination_position{
      {44, 1},
      {45, 2},
      {46, 3},
      {54, 4},
      {55, 5},
      {56, 6},
      {64, 7},
      {65, 8},
      {66, 9},
  };
  std::vector<std::string> bin_labels{"44", "45", "46", "54", "55", "56", "64", "65", "66"};
  for( auto histo: matrix_first_harm )
    for( int i=1; i<10; ++i ){
      histo->GetXaxis()->SetBinLabel(i, bin_labels.at(i-1).c_str());
      histo->GetYaxis()->SetBinLabel(i, bin_labels.at(i-1).c_str());
    }
  for( auto histo: matrix_second_harm )
    for( int i=1; i<10; ++i ){
      histo->GetXaxis()->SetBinLabel(i, bin_labels.at(i-1).c_str());
      histo->GetYaxis()->SetBinLabel(i, bin_labels.at(i-1).c_str());
    }
  std::string th3f_name{"pgd_prim_delta_phi_pt_layers_all_22"};
  file->GetObject(th3f_name.c_str(), histogram_layers_pt_phi);
  assert(histogram_layers_pt_phi);
  for( int n1 = 4; n1 <7; ++n1){
    for (int n2 = 4; n2 < 7; ++n2) {
      for (int n3 = 4; n3 < 7; ++n3) {
        for (int n4 = 4; n4 < 7; ++n4) {
          int combination = n1*1000 + n2*100 + n3*10 + n4*1;
          int comb_in = n1*10+n2;
          int comb_out = n3*10+n4;
          std::string projection_name =
              "delta_phi_low_pT_" + std::to_string(combination);
          auto *histogram_phi = histogram_layers_pt_phi->ProjectionY(
              projection_name.c_str(), 1, 5, combination+1, combination+1);
          histogram_phi->Scale(1.0 / histogram_phi->GetEntries());
          projection_name = "fit_delta_phi_low_pT_" + std::to_string(combination);
          auto *fit_function =
              new TF1(projection_name.c_str(), "[0]*(1+2*[1]*cos(x)+2*[2]*cos(2*x))",
                      -3.0, 3.0);
          histogram_phi->Fit(fit_function, "Q", "", -3.0, 3.0);
          auto v1 = fit_function->GetParameter(1);
          auto v1_err = fit_function->GetParError(1);
          matrix_first_harm.at(0)->SetBinContent(
              map_combination_position.at(comb_in),
              map_combination_position.at(comb_out),
              v1);
          matrix_first_harm.at(0)->SetBinError(
              map_combination_position.at(comb_in),
              map_combination_position.at(comb_out),
              v1_err);
          auto v2 = fit_function->GetParameter(2);
          auto v2_err = fit_function->GetParError(2);
          matrix_second_harm.at(0)->SetBinContent(
              map_combination_position.at(comb_in),
              map_combination_position.at(comb_out),
              v2);
          matrix_second_harm.at(0)->SetBinError(
              map_combination_position.at(comb_in),
              map_combination_position.at(comb_out),
              v2_err);
//          std::cout << comb_in << " : " << comb_out << " " << v1 <<  std::endl;
        }
      }
    }
  }
  for( int n1 =4; n1 <7; ++n1){
    for (int n2 = 4; n2 < 7; ++n2) {
      for (int n3 = 4; n3 < 7; ++n3) {
        for (int n4 = 4; n4 < 7; ++n4) {
          int combination = n1*1000 + n2*100 + n3*10 + n4*1;
          int comb_in = n1*10+n2;
          int comb_out = n3*10+n4;
          std::string projection_name =
              "delta_phi_mid_pT_" + std::to_string(combination);
          auto *histogram_phi = histogram_layers_pt_phi->ProjectionY(
              projection_name.c_str(), 6, 10, combination+1, combination+1);
          histogram_phi->Scale(1.0 / histogram_phi->GetEntries());
          projection_name = "fit_delta_mid_low_pT_" + std::to_string(combination);
          auto *fit_function =
              new TF1(projection_name.c_str(), "[0]*(1+2*[1]*cos(x)+2*[2]*cos(2*x))",
                      -3.0, 3.0);
          histogram_phi->Fit(fit_function, "Q", "", -3.0, 3.0);
          auto v1 = fit_function->GetParameter(1);
          auto v1_err = fit_function->GetParError(1);
          matrix_first_harm.at(1)->SetBinContent(
              map_combination_position.at(comb_in),
              map_combination_position.at(comb_out),
              v1);
          matrix_first_harm.at(1)->SetBinError(
              map_combination_position.at(comb_in),
              map_combination_position.at(comb_out),
              v1_err);
          auto v2 = fit_function->GetParameter(2);
          auto v2_err = fit_function->GetParError(2);
          matrix_second_harm.at(1)->SetBinContent(
              map_combination_position.at(comb_in),
              map_combination_position.at(comb_out),
              v2);
          matrix_second_harm.at(1)->SetBinError(
              map_combination_position.at(comb_in),
              map_combination_position.at(comb_out),
              v2_err);
//          std::cout << comb_in << " : " << comb_out << " " << v1 <<  std::endl;
        }
      }
    }
  }
  for( int n1 =4; n1 <7; ++n1){
    for (int n2 = 4; n2 < 7; ++n2) {
      for (int n3 = 4; n3 < 7; ++n3) {
        for (int n4 = 4; n4 < 7; ++n4) {
          int combination = n1*1000 + n2*100 + n3*10 + n4*1;
          int comb_in = n1*10+n2;
          int comb_out = n3*10+n4;
          std::string projection_name =
              "delta_phi_high_pT_" + std::to_string(combination);
          auto *histogram_phi = histogram_layers_pt_phi->ProjectionY(
              projection_name.c_str(), 11, 16, combination+1, combination+1);
          histogram_phi->Scale(1.0 / histogram_phi->GetEntries());
          projection_name = "fit_delta_mid_high_pT_" + std::to_string(combination);

          auto *fit_function =
              new TF1(projection_name.c_str(), "[0]*(1+2*[1]*cos(x)+2*[2]*cos(2*x))",
                      -3.0, 3.0);
          histogram_phi->Fit(fit_function, "Q", "", -3.0, 3.0);
          auto v1 = fit_function->GetParameter(1);
          auto v1_err = fit_function->GetParError(1);
          matrix_first_harm.at(2)->SetBinContent(
              map_combination_position.at(comb_in),
              map_combination_position.at(comb_out),
              v1);
          matrix_first_harm.at(2)->SetBinError(
              map_combination_position.at(comb_in),
              map_combination_position.at(comb_out),
              v1_err);
          auto v2 = fit_function->GetParameter(2);
          auto v2_err = fit_function->GetParError(2);
          matrix_second_harm.at(2)->SetBinContent(
              map_combination_position.at(comb_in),
              map_combination_position.at(comb_out),
              v2);
          matrix_second_harm.at(2)->SetBinError(
              map_combination_position.at(comb_in),
              map_combination_position.at(comb_out),
              v2_err);
//          std::cout << comb_in << " : " << comb_out << " " << v1 <<  std::endl;
        }
      }
    }
  }
  auto file_out = TFile::Open("fit_parameters.root", "RECREATE");
  file_out->cd();
  for(auto graph : matrix_first_harm)
    graph->Write();

  for(auto graph : matrix_second_harm)
    graph->Write();

  file_out->Close();
}