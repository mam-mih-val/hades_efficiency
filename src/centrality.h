
#ifndef CENTRALITY_H
#define CENTRALITY_H

#include <iostream>
#include <string>

#include "TH1F.h"

class Centrality
{
public:
  static Centrality* GetInstance(){
    if(!instance_)
      instance_ = new Centrality();
    return instance_;
  }
  int GetCentralityClass5pc(int multiplicity){
    auto bin = centrality_5pc_->FindBin(multiplicity);
    return static_cast<int>(centrality_5pc_->GetBinContent(bin) - 1);
  };
  int GetCentralityClass10pc(int multiplicity){
    auto bin = centrality_10pc_->FindBin(multiplicity);
    return static_cast<int>(centrality_10pc_->GetBinContent(bin) - 1);
  };
  float GetCentrality5pc(int multiplicity){
    return 5.0*GetCentralityClass5pc(multiplicity)+2.5;
  }
  float GetCentrality10pc(int multiplicity){
    return 10.0*GetCentralityClass10pc(multiplicity)+5.0;
  }
  int GetNumClasses5pc(){
    return centrality_5pc_->GetNbinsX()-1;
  };
  int GetNumClasses10pc(){
    return centrality_10pc_->GetNbinsX()-1;
  };
  void LoadCentrality5pc(){
    Double_t xAxis[17] = {0, 10, 17, 23, 29, 36, 45, 55, 68, 82, 99, 117, 136, 157, 180, 312, 499};
    centrality_5pc_ = new TH1F("TOFRPC_5pc_fixedCuts__1","TOFRPC_5pc_fixedCuts",16, xAxis);
    centrality_5pc_->SetBinContent(1,15);
    centrality_5pc_->SetBinContent(2,14);
    centrality_5pc_->SetBinContent(3,13);
    centrality_5pc_->SetBinContent(4,12);
    centrality_5pc_->SetBinContent(5,11);
    centrality_5pc_->SetBinContent(6,10);
    centrality_5pc_->SetBinContent(7,9);
    centrality_5pc_->SetBinContent(8,8);
    centrality_5pc_->SetBinContent(9,7);
    centrality_5pc_->SetBinContent(10,6);
    centrality_5pc_->SetBinContent(11,5);
    centrality_5pc_->SetBinContent(12,4);
    centrality_5pc_->SetBinContent(13,3);
    centrality_5pc_->SetBinContent(14,2);
    centrality_5pc_->SetBinContent(15,1);
    centrality_5pc_->SetBinError(1,1.8301);
    centrality_5pc_->SetBinError(2,5.29551);
    centrality_5pc_->SetBinError(3,5.07632);
    centrality_5pc_->SetBinError(4,4.68075);
    centrality_5pc_->SetBinError(5,4.83154);
    centrality_5pc_->SetBinError(6,5.30302);
    centrality_5pc_->SetBinError(7,4.84388);
    centrality_5pc_->SetBinError(8,5.20623);
    centrality_5pc_->SetBinError(9,4.82995);
    centrality_5pc_->SetBinError(10,5.22526);
    centrality_5pc_->SetBinError(11,4.99358);
    centrality_5pc_->SetBinError(12,4.83847);
    centrality_5pc_->SetBinError(13,4.94002);
    centrality_5pc_->SetBinError(14,5.02256);
    centrality_5pc_->SetBinError(15,5.08179);
    centrality_5pc_->SetBinError(16,0.00100254);
    centrality_5pc_->SetEntries(16);
  };
  void LoadCentrality10pc(){
    Double_t xAxis1[10] = {0, 10, 23, 36, 55, 82, 117, 157, 312, 499};

    centrality_10pc_ = new TH1F("TOFRPC_10pc_fixedCuts__1","TOFRPC_10pc_fixedCuts",9, xAxis1);
    centrality_10pc_->SetBinContent(1,8);
    centrality_10pc_->SetBinContent(2,7);
    centrality_10pc_->SetBinContent(3,6);
    centrality_10pc_->SetBinContent(4,5);
    centrality_10pc_->SetBinContent(5,4);
    centrality_10pc_->SetBinContent(6,3);
    centrality_10pc_->SetBinContent(7,2);
    centrality_10pc_->SetBinContent(8,1);
    centrality_10pc_->SetBinError(1,1.8301);
    centrality_10pc_->SetBinError(2,10.3718);
    centrality_10pc_->SetBinError(3,9.5123);
    centrality_10pc_->SetBinError(4,10.1469);
    centrality_10pc_->SetBinError(5,10.0362);
    centrality_10pc_->SetBinError(6,10.2188);
    centrality_10pc_->SetBinError(7,9.77849);
    centrality_10pc_->SetBinError(8,10.1044);
    centrality_10pc_->SetBinError(9,0.00100254);
    centrality_10pc_->SetEntries(9);
  };
private:
  Centrality() :
  centrality_5pc_(nullptr),
  centrality_10pc_(nullptr){
    LoadCentrality5pc();
    LoadCentrality10pc();
  };
  ~Centrality() = default;
  static Centrality* instance_;
  TH1F* centrality_5pc_;
  TH1F* centrality_10pc_;
};

#endif // CENTRALITY_H