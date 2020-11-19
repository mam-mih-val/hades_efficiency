//
// Created by mikhail on 11/5/20.
//

#include <TF1.h>
#include <TFile.h>
#include <TGraph2DErrors.h>
#include <TH2F.h>
#include <TH3F.h>
#include <string>
int main(int n, char** args) {
  std::string input_file=args[1];
  auto file = TFile::Open(input_file.c_str());
  std::vector<TH3F*> pdg_delta_phi_vs_pt_rapidity;
  std::vector<TH3F*> gen_delta_phi_vs_pt_rapidity;
  std::vector<TH3F*> eff_vs_delta_phi_pt_rapidity;
  std::vector<TH2F*> v0_vs_pt_rapidity;
  std::vector<TH2F*> v1_vs_pt_rapidity;
  std::vector<TH2F*> v2_vs_pt_rapidity;

  int p=2;
  while( p<40 ){
    std::string name = "pdg_prim_delta_phi_pt_midrapidity_" + std::to_string(p);
    pdg_delta_phi_vs_pt_rapidity.emplace_back();
    file->GetObject( name.c_str(), pdg_delta_phi_vs_pt_rapidity.back() );
    pdg_delta_phi_vs_pt_rapidity.back()->Rebin3D(1, 2, 2);
    name = "gen_prim_delta_phi_pt_midrapidity_" + std::to_string(p);
    gen_delta_phi_vs_pt_rapidity.emplace_back();
    file->GetObject(name.c_str(), gen_delta_phi_vs_pt_rapidity.back());
    gen_delta_phi_vs_pt_rapidity.back()->Rebin3D(1, 2, 2);
    eff_vs_delta_phi_pt_rapidity.emplace_back();
//    pdg_delta_phi_vs_pt_rapidity.back()->Divide(gen_delta_phi_vs_pt_rapidity.back());
    p+=5;
  }
  p=2;
  auto file_out = TFile::Open("fit_parameters.root", "RECREATE");
  file_out->mkdir("fits");
  file_out->cd("/fits");
  for( size_t i=0; i<pdg_delta_phi_vs_pt_rapidity.size(); i++ ){
    size_t n_bins_y = pdg_delta_phi_vs_pt_rapidity.at(i)->GetNbinsX();
    size_t n_bins_pt = pdg_delta_phi_vs_pt_rapidity.at(i)->GetNbinsY();
    std::string name_v1 = "v1_y_pT_"+std::to_string(p);
    v1_vs_pt_rapidity.push_back( new TH2F(name_v1.c_str(), ";y;pT;v_{1}",
                                         15, -0.75, 0.75,
                                         10, 0.0, 2.0 ) );
    v1_vs_pt_rapidity.back()->GetSumw2();
    std::string name_v2 = "v2_y_pT_"+std::to_string(p);
    v2_vs_pt_rapidity.push_back( new TH2F(name_v2.c_str(), ";y;pT;v_{2}",
                                         15, -0.75, 0.75,
                                         10, 0.0, 2.0 ) );
    for( size_t y =1; y <=n_bins_y; ++y){
      for(size_t pT =1; pT <n_bins_pt; ++pT){
        auto bin_y = pdg_delta_phi_vs_pt_rapidity.at(i)->GetXaxis()->GetBinCenter(y);
        auto bin_pT = pdg_delta_phi_vs_pt_rapidity.at(i)->GetYaxis()->GetBinCenter(pT);
        name_v1 = "v1_y_pT_"+std::to_string(p) + "_"+std::to_string(bin_y)+"_"+std::to_string(bin_pT);
        auto proj_num = pdg_delta_phi_vs_pt_rapidity.at(i)->ProjectionZ(name_v1.c_str(),y, y, pT, pT );
        std::string name_gen = name_v1+"_gen";
        auto proj_den = gen_delta_phi_vs_pt_rapidity.at(i)->ProjectionZ(name_gen.c_str(),y, y, pT, pT );
        proj_num->Divide(proj_den);
        name_v1 +="_fit";
        auto fit_func = new TF1(name_v1.c_str(), "[0]*( 1+2*[1]*cos(x)+2*[2]*cos(2*x) )", -3.15, 3.15 );
        proj_num->Fit(fit_func);
        v1_vs_pt_rapidity.back()->SetBinContent( y, pT, fit_func->GetParameter(1) );
        v1_vs_pt_rapidity.back()->SetBinError( y, pT, fit_func->GetParError(1) );
        v2_vs_pt_rapidity.back()->SetBinContent( y, pT, fit_func->GetParameter(2) );
        v2_vs_pt_rapidity.back()->SetBinError( y, pT, fit_func->GetParError(2) );
        proj_num->Write();
      }
    }
    p+=5;
  }
  file_out->cd("/");
  for( auto histo : v1_vs_pt_rapidity )
    histo->Write();
  for( auto histo : v2_vs_pt_rapidity )
    histo->Write();
  file_out->Close();
}