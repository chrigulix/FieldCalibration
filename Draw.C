// C++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <sstream>
#include <algorithm>

// C headers
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <math.h>
#include <cstring>
#include <unistd.h>

// ROOT headers
#include "TROOT.h"
#include "TObject.h"
#include "TApplication.h"
#include "TAttLine.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TProfile.h"
#include "TPolyLine3D.h"
#include "TStyle.h"
#include "TFrame.h"
#include "TFile.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TView3D.h"
#include "TTree.h"

void Draw()
{
  const float kDetector[] = {2.5604,2.325,10.368}; // Size of the detectors in m (used from docdb 1821)
  const int kResolution[] = {25,25,100};
  
  // Remove gifs
  std::remove("RecoX.gif"); std::remove("RecoY.gif"); std::remove("RecoZ.gif");
  std::remove("ResidualX.gif"); std::remove("ResidualY.gif"); std::remove("ResidualZ.gif");
  
  
  TFile* FieldFile = new TFile("RecoDispl.root","READ");
  if (FieldFile->IsZombie()) 
  {
       std::cout << "Error opening field file" << std::endl;
       exit(-1);
  }

  std::vector<TH3F*> DistortMap;
  DistortMap.push_back((TH3F*) FieldFile->Get("Distortion_Field_X"));
  DistortMap.push_back((TH3F*) FieldFile->Get("Distortion_Field_Y"));
  DistortMap.push_back((TH3F*) FieldFile->Get("Distortion_Field_Z"));
  
//   FieldFile->Close();
  
  TFile* RecoFile = new TFile("RecoDispl.root","READ");
  if (RecoFile->IsZombie()) 
  {
       std::cout << "Error opening reco file" << std::endl;
       exit(-1);
  }
  
  std::vector<TH3F*> ReconstructionMap;
  ReconstructionMap.push_back((TH3F*) RecoFile->Get("Reco_Field_X"));
  ReconstructionMap.push_back((TH3F*) RecoFile->Get("Reco_Field_Y"));
  ReconstructionMap.push_back((TH3F*) RecoFile->Get("Reco_Field_Z"));
//   ReconstructionMap.push_back((TH3F*) RecoFile->Get("Reco_Field_X"));

  //   RecoFile->Close();
  int Nxbins = ReconstructionMap[0]->GetNbinsX();
  int Nybins = ReconstructionMap[0]->GetNbinsY();
  int Nzbins = ReconstructionMap[0]->GetNbinsZ();
  
//   unsigned int MaximumBin[3];
//   unsigned int MinimumBin[3];
  
//   double Maximum[3];
//   double Minimum[3];
  
//   for(unsigned coord = 0; coord < 3; coord++)
//   {
//     MaximumBin[coord] = ReconstructionMap[coord]->GetMaximumBin();
//     MinimumBin[coord] = ReconstructionMap[coord]->GetMinimumBin();
    
//     Maximum[coord] = ReconstructionMap[coord]->GetMaximum();
//     Minimum[coord] = ReconstructionMap[coord]->GetMinimum();
    
//     Maximum[coord] = ReconstructionMap[coord]->GetBinContent(MaximumBin);
//     Minimum[coord] = ReconstructionMap[coord]->GetBinContent(MinimumBin);
//   }
  
  
  int Smoothing = 3; // order of smoothing
  int SRange = (int)(Smoothing/2.);
  int Nelements; // # pixels used for smoothing
  
//   for(int coord = 0; coord < 3; coord++)
//   {
//     ReconstructionMap[3] = ReconstructionMap[coord];
//     for(int xbin = 1; xbin <= Nxbins; xbin++) for(int ybin = 1; ybin <= Nybins; ybin++) for(int zbin = 1; zbin <= Nzbins; zbin++)
//     {
//       Nelements = 0;
//       double SmoothEntry = 0.;
//       for(int dx = -SRange; dx <= SRange; dx++) for(int dy = -SRange; dy <= SRange; dy++) for(int dz = -SRange; dz <= SRange; dz++)
//       {
// 	if(xbin+dx > 1 && xbin+dx <= Nxbins && ybin+dy > 1 && ybin+dy <= Nybins && zbin+dz > 1 && zbin+dz <= Nzbins)
// 	{
// 	  Nelements++;
// 	  SmoothEntry += ReconstructionMap[coord]->GetBinContent(xbin+dx,ybin+dy,zbin+dz);
// 	}
//       }
//       SmoothEntry /= (double)Nelements;
//       ReconstructionMap[3]->SetBinContent(xbin,ybin,zbin,SmoothEntry);
//     }
//     ReconstructionMap[coord] = ReconstructionMap[3];
//   }
    
//   delete ReconstructionMap[3];
//   ReconstructionMap.pop_back();
  
  
  std::vector<TH2F*> ProjectionXZ;
  ProjectionXZ.push_back(new TH2F("Reco Map X","Reco Map X",kResolution[2],0,kDetector[2],kResolution[0],0,kDetector[0]));
  ProjectionXZ.push_back(new TH2F("Reco Map Y","Reco Map Y",kResolution[2],0,kDetector[2],kResolution[0],0,kDetector[0]));
  ProjectionXZ.push_back(new TH2F("Reco Map Z","Reco Map Z",kResolution[2],0,kDetector[2],kResolution[0],0,kDetector[0]));
  
  int MinBin = 1+SRange;
  for (int ybin = MinBin; ybin <= kResolution[1]-SRange; ybin++) 
  {
    for (int xbin = MinBin; xbin <= kResolution[0]-SRange; xbin++) for (int zbin = MinBin; zbin <= kResolution[2]-SRange; zbin++)
    {
//       One Laser
      for (int coord = 0; coord < 3; coord++) ProjectionXZ[coord] -> SetBinContent(zbin,xbin,ReconstructionMap[coord]->GetBinContent(xbin,ybin,zbin));
//       Two laser approach in symmetric distributions
//       for (int coord = 0; coord < 3; coord++) 
// 	ProjectionXZ[coord] -> SetBinContent(zbin,xbin,(-ReconstructionMap[coord]->GetBinContent(xbin,ybin,zbin)-ReconstructionMap[coord]->GetBinContent(xbin,ybin,kResolution[3]-zbin+1))/2.0);
    }
    for (int coord = 0; coord < 3; coord++) ProjectionXZ[coord] -> SetStats(0);
    
    ProjectionXZ[0] -> SetMaximum(40);
    ProjectionXZ[0] -> SetMinimum(-40);
    ProjectionXZ[1] -> SetMaximum(40);
    ProjectionXZ[1] -> SetMinimum(-40);
    ProjectionXZ[2] -> SetMaximum(40);
    ProjectionXZ[2] -> SetMinimum(-40);
    
    TCanvas * C1 = new TCanvas("Reco Map X","Reco Map X",1000,500);
    ProjectionXZ[0] -> Draw("colz");
    C1 -> Print("RecoX.gif+5","gif+5");
    TCanvas * C2 = new TCanvas("Reco Map Y","Reco Map Y",1000,500);
    ProjectionXZ[1] -> Draw("colz");
    C2 -> Print("RecoY.gif+5","gif+5");
    TCanvas * C3 = new TCanvas("Reco Map Z","Reco Map Z",1000,500);
    ProjectionXZ[2] -> Draw("colz");
    C3 -> Print("RecoZ.gif+5","gif+5");
    
    delete C1;
    delete C2;
    delete C3;
  }
  
  for (int ybin = MinBin; ybin <= kResolution[1]-SRange; ybin++) 
  {
    for (int xbin = MinBin; xbin <= kResolution[0]-SRange; xbin++) for (int zbin = MinBin; zbin <= kResolution[2]-SRange; zbin++)
    {
      // One laser
      for (int coord = 0; coord < 3; coord++) ProjectionXZ[coord] -> SetBinContent(zbin,xbin,(ReconstructionMap[coord]->GetBinContent(xbin,ybin,zbin)-DistortMap[coord]->GetBinContent(xbin,ybin,zbin)));
      // Two lasers
//       for (int coord = 0; coord < 3; coord++) 
// 	ProjectionXZ[coord] -> SetBinContent(zbin,xbin,(ReconstructionMap[coord]->GetBinContent(xbin,ybin,zbin)+ReconstructionMap[coord]->GetBinContent(xbin,ybin,kResolution[3]-zbin+1))/2.0+DistortMap[coord]->GetBinContent(xbin,ybin,zbin));
    }
    for (int coord = 0; coord < 3; coord++) ProjectionXZ[coord] -> SetStats(0);
    
    ProjectionXZ[0] -> SetMaximum(0.01);
    ProjectionXZ[0] -> SetMinimum(-0.01);
    ProjectionXZ[1] -> SetMaximum(0.015);
    ProjectionXZ[1] -> SetMinimum(-0.015);
    ProjectionXZ[2] -> SetMaximum(0.05);
    ProjectionXZ[2] -> SetMinimum(-0.05);
    
    TCanvas * C1 = new TCanvas("Residual Map X","Residual Map X",1000,500);
    ProjectionXZ[0] -> Draw("colz");
    C1 -> Print("ResidualX.gif+5","gif+50");
    TCanvas * C2 = new TCanvas("Residual Map Y","Residual Map Y",1000,500);
    ProjectionXZ[1] -> Draw("colz");
    C2 -> Print("ResidualY.gif+5","gif+50");
    TCanvas * C3 = new TCanvas("Residual Map Z","Residual Map Z",1000,500);
    ProjectionXZ[2] -> Draw("colz");
    C3 -> Print("ResidualZ.gif+5","gif+50");
    
    delete C1;
    delete C2;
    delete C3;
  }
//   for(int coord = 0; coord < 3; coord++)
//   {
//     delete DistortMap[coord];
//     delete ReconstructionMap[coord];
//     delete ProjectionXZ[coord];
//   }
//   FieldFile->Close();
//   RecoFile->Close();
};