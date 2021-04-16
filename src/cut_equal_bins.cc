//
// Created by mikhail on 11/11/20.
//

#include <TFile.h>
#include <TH1D.h>
#include <iostream>

int main(int n_args, char** args){
  if( n_args < 4 )
    throw std::runtime_error("2 argumets expected. Try: ./equal_bins path/to/file.root histogram_name n_bins");
  auto file = TFile::Open( args[1], "read" );
  if( !file )
    throw std::runtime_error("Couldn't open the file "+std::string(args[1]));
  TH1F* histo;
  file->GetObject( args[2], histo );
  if( !histo )
    throw std::runtime_error("Couldn't read the histogram "+std::string(args[2])+" from file "+std::string(args[1]));
  int n_bins = atoi(args[3]);
  double part = 1.0/n_bins;
  auto entries = histo->GetEntries();
  int bin_enries=0;
  std::vector<double> new_axis;
  new_axis.push_back( histo->GetXaxis()->GetBinLowEdge(1) );
  for( auto i=0; i<=histo->GetNbinsX(); ++i ){
    bin_enries+=histo->GetBinContent(i);
    if( bin_enries < part*entries )
      continue;
    else{
      new_axis.push_back(histo->GetBinLowEdge(i-1));
      bin_enries=0;
    }
  }
  new_axis.push_back( histo->GetXaxis()->GetBinUpEdge(histo->GetNbinsX()) );
  std::cout << "{ ";
  for( auto num : new_axis ){
    std::cout << num << ", " ;
  }
  std::cout << "\b\b }" << std::endl;
  return 0;
}