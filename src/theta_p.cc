//
// Created by mikhail on 11/29/20.
//

#include "TCanvas.h"
#include "TFile.h"
#include <TF1.h>
#include "TMath.h"
#include <vector>

int main(){
  std::vector<TF1*> y_edges;
  std::vector<TF1*> pt_edges;
  double y=-0.01;
  while( y<=1.49 ){
    std::string y2 = std::to_string( 2*y );
//    std::string formula = "y = TMath::ACos( TMath::Sqrt( x*x+0.94*0.94 ) * (exp("+y2+") - 1)/( x*("+y2+"+1) ) )";
    std::string formula = "TMath::ACos( TMath::Sqrt( x*x+0.94*0.94 ) * (exp("+y2+") - 1)/( x*( exp("+y2+")+1) ) )";
    std::string name = "edge_y_"+std::to_string(y);
    y_edges.push_back( new TF1( name.c_str(), formula.c_str(), 0.0, 3.0 ) );
    y+=0.1;
  }
  double pT=0.0;
  while (pT < 2.0){
    std::string pt = std::to_string(pT);
    std::string formula = "TMath::ASin("+pt+"/x)";
    std::string name = "edge_pT_"+pt;
    pt_edges.push_back( new TF1( name.c_str(), formula.c_str(), 0.0, 3.0 ) );
    pT+=0.2;
  }
  auto canv = new TCanvas( "canv", "", 1000, 1100 );
  canv->cd();
  y_edges.front()->Draw();
  for( int i=1; i<y_edges.size(); i++ ){
    y_edges.at(i)->Draw("same");
  }
  canv->SaveAs("canv.png");
  auto file_out = TFile::Open( "theta_p.root", "recreate" );
  file_out->cd();
  for( auto edge : y_edges )
    edge->Write();
  for( auto edge : pt_edges )
    edge->Write();
  file_out->Close();
  return 0;
}