// C++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <algorithm>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <chrono>
#include <cmath>
#include <cstring>
#include <thread>
#include <array>

// C headers
#include <pthread.h>
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
#include "TChain.h"
#include "TVector3.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

// Own Files
#include "include/LaserTrack.hpp"
#include "include/ThreeVector.hpp"
#include "include/TPCVolumeHandler.hpp"
#include "include/Interpolation3D.hpp"
#include "include/Matrix3x3.hpp"
#include "include/Laser.hpp"

int main(int argc, char** argv);
Laser ReadRecoTracks(int argc, char** argv);
void LaserInterpThread(Laser&, const Laser&, const Delaunay&);
void WriteRootFile(std::vector<ThreeVector<float>>&, TPCVolumeHandler&);
void WriteTextFile(std::vector<ThreeVector<float>>&);

// Main function
int main(int argc, char** argv)
{
    // Start timer, just because it's nice to know how long shit takes
    time_t timer;
    std::time(&timer);

    std::cout << argc << std::endl;
    
    // If there are to few input arguments abord
    if(argc < 2)
    {
        std::cerr << "ERROR: Too few arguments, use ./LaserCal <input file names>" << std::endl;
        return -1;
    }

    // Coose detector dimensions, coordinate system offset and resolutions
    ThreeVector<float> DetectorSize = {256.04,232.5,1036.8};
    ThreeVector<float> DetectorOffset = {0.0,-DetectorSize[1]/static_cast<float>(2.0),0.0};
    ThreeVector<unsigned long> DetectorResolution = {26,26,101};
   
    // Create the detector volume
    TPCVolumeHandler Detector(DetectorSize, DetectorOffset, DetectorResolution);
    
    // Read data and store it to a Laser object
    std::cout << "Reading data..." << std::endl;
    Laser LaserTrackSet = ReadRecoTracks(argc,argv);
    
    // Calculate track displacement
    std::cout << "Find track displacements.. " << std::endl;
    LaserTrackSet.CorrectTrackSet();
    
    // Create delaunay mesh
    std::cout << "Generate mesh..." << std::endl;
    Delaunay Mesh = TrackMesher(LaserTrackSet.GetTrackSet());
    
    // Interpolate Displacement Map (regularly spaced grid)
    std::vector<ThreeVector<float>> DisplacementMap = InterpolateMap(LaserTrackSet.GetTrackSet(),Mesh,Detector);
    
    // Fill displacement map into TH3 histograms and write them to file
    WriteRootFile(DisplacementMap,Detector);
    
//     Detector.GetMapMaximum().print();
//     Detector.GetMapMinimum().print();
    
//   double HistRange[3][2];
//   for (int coord = 0; coord < 3; coord++)
//   {
//       for(int minmax = 0; minmax < 2; minmax++)
//       {
//           HistRange[coord][minmax] = (DetectorSize[coord]+pow(2,-25))*( minmax - pow(-1,minmax)/((double)DetectorResolution[coord]-1.0)/2.0 );
//           std::cout << DetectorOffset[coord]+HistRange[coord][minmax] << " ";
//       }
//       std::cout << std::endl;
//   }
  
//   time_t timer;
//   std::time(&timer);
// 
//   std::vector<ThreeVector<float>> Displacement;
//   
//   ThreeVector<float> Location;
//   for(unsigned xbin = 0; xbin < DetectorResolution[0]; xbin++) 
//   {
//     std::cout << xbin << std::endl;
//     Location[0] = DetectorOffset[0]+(HistRange[0][1]-HistRange[0][0])/(float)DetectorResolution[0]*xbin;
//     for(unsigned ybin = 0; ybin < DetectorResolution[1]; ybin++) 
//     {
//       Location[1] = DetectorOffset[1]+(HistRange[1][1]-HistRange[1][0])/(float)DetectorResolution[1]*ybin;
//       for(unsigned zbin = 0; zbin < DetectorResolution[2]; zbin++)
//       {
// 	Location[2] = DetectorOffset[2]+(HistRange[2][1]-HistRange[2][0])/(float)DetectorResolution[2]*zbin;
// 	Displacement.push_back(InterpolateCGAL(LaserTrackSet.GetTrackSet(),Mesh,Location));
//       }
//     }
//   }
    std::cout << "End of program after "<< std::difftime(std::time(NULL),timer) << " s" << std::endl;
}

Laser ReadRecoTracks(int argc, char** argv)
{
    // Create Laser (collection of laser tracks) this will be the returned object
    Laser TrackSelection;
    
    // Initialize read variables, the pointers for more complex data structures 
    // are very important for Root. Rene Brun in hell (do you see what I did there?)
    int EventNumber;
    
    std::vector<TVector3> TrackSamples;
    std::vector<TVector3>* pTrackSamples = &TrackSamples;
    
    TVector3 EntryPoint;
    TVector3* pEntryPoint = &EntryPoint;
    TVector3 ExitPoint;
    TVector3* pExitPoint =&ExitPoint;
    
    // Open TChains to store all trees
    TChain* LaserInfoTree = new TChain("lasers");
    TChain* RecoTrackTree = new TChain("tracks");
    
    // Loop through all input files and add them to the TChain
    for(int arg = 1; arg < argc; arg++)
    {
        // Open input file and add to TChains
        LaserInfoTree->Add(argv[arg]);
        RecoTrackTree->Add(argv[arg]);
    }
    
    // Assigne branch addresses
    LaserInfoTree->SetBranchAddress("entry",&pEntryPoint);
    LaserInfoTree->SetBranchAddress("exit", &pExitPoint);
    RecoTrackTree->SetBranchAddress("track",&pTrackSamples);
    RecoTrackTree->SetBranchAddress("event", &EventNumber);
    
    // Only start read out when both trees have the same amount of entries 
    if(LaserInfoTree->GetEntries() == RecoTrackTree->GetEntries())
    {
        // Loop over all tree entries
        for(Size_t tree_index = 0; tree_index < RecoTrackTree->GetEntries(); tree_index++)
        {
            // Get tree entries of both trees
            LaserInfoTree->GetEntry(tree_index);
            RecoTrackTree->GetEntry(tree_index);
            

            // This here sorts the tracks by their distance to the EntryPoint. The algorithm uses a lambda
            // It will compare the distance to the EntryPoint of two vector entries A & B
            std::sort(TrackSamples.begin(),TrackSamples.end(), [&EntryPoint](TVector3 A, TVector3 B)
                {
                    A -= EntryPoint;
                    B -= EntryPoint;
                    // Her only the squared distance was used to avoid costly sqrt operations
                    return A.Mag2() > B.Mag2();
                }    
            );
            
            // This step will erase all double entries. First std::unique shifts every double to the end
            // of the vector and gives back the new end point of the data set. After that we erase the HistRange
            // between this new end and the real end of the vector
            TrackSamples.erase( std::unique(TrackSamples.begin(),TrackSamples.end()), TrackSamples.end());
            
            // Add new track to Laser TrackSelection
            TrackSelection.AppendTrack(LaserTrack(EntryPoint,ExitPoint,TrackSamples));
        }
    }
    else // If the trees don't have the same amount of entries, through error (I know not propper error handling)
    {
        std::cerr << "ERROR: Two TTrees don't have the same ammount of entries!" << std::endl;
    }
    
//     delete pEntryPoint;
//     delete pExitPoint;
//     delete pTrackSamples;
    
    gDirectory->GetList()->Delete();
    
    delete LaserInfoTree;
    delete RecoTrackTree;
    
    return TrackSelection;
}

void LaserInterpThread(Laser& LaserTrackSet, const Laser& InterpolationLaser, const Delaunay& InterpolationMesh)
{
  LaserTrackSet.InterpolateTrackSet(InterpolationLaser, InterpolationMesh);
}

void WriteRootFile(std::vector<ThreeVector<float>>& InterpolationData, TPCVolumeHandler& TPCVolume)
{ 
    // Store TPC properties which are important for the TH3 generation
    ThreeVector<unsigned long> Resolution = TPCVolume.GetDetectorResolution();
    ThreeVector<float> MinimumCoord = TPCVolume.GetMapMinimum();
    ThreeVector<float> MaximumCoord = TPCVolume.GetMapMaximum();
  
    // Initialize all TH3F
    std::vector<TH3F> RecoField;
    RecoField.push_back(TH3F("Reco_Field_X","Reco Field X",Resolution[0],MinimumCoord[0],MaximumCoord[0],Resolution[1],MinimumCoord[1],MaximumCoord[1],Resolution[2],MinimumCoord[2],MaximumCoord[2]));
    RecoField.push_back(TH3F("Reco_Field_Y","Reco Field Y",Resolution[0],MinimumCoord[0],MaximumCoord[0],Resolution[1],MinimumCoord[1],MaximumCoord[1],Resolution[2],MinimumCoord[2],MaximumCoord[2]));
    RecoField.push_back(TH3F("Reco_Field_Z","Reco Field Z",Resolution[0],MinimumCoord[0],MaximumCoord[0],Resolution[1],MinimumCoord[1],MaximumCoord[1],Resolution[2],MinimumCoord[2],MaximumCoord[2]));
  

    // Loop over all xbins
    for(unsigned xbin = 0; xbin < TPCVolume.GetDetectorResolution()[0]; xbin++)
    {
        // Loop over all ybins
        for(unsigned ybin = 0; ybin < TPCVolume.GetDetectorResolution()[1]; ybin++) 
        {
            // Loop over all zbins
            for(unsigned zbin = 0; zbin < TPCVolume.GetDetectorResolution()[2]; zbin++)
            {
                // Loop over all coordinates
                for(unsigned coord = 0; coord < 3; coord++)
                {
                    // Fill interpolated grid points into histograms
                    RecoField[coord].SetBinContent(xbin+1,ybin+1,zbin+1, InterpolationData[zbin+( TPCVolume.GetDetectorResolution()[2]*(ybin+TPCVolume.GetDetectorResolution()[1]*xbin) )][coord]);
                } // end coordinate loop
            } // end zbin loop
        } // end ybin loop
    } // end zbin loop
  
  TFile *OutputFile = new TFile("RecoField.root", "recreate");
  for(unsigned coord = 0; coord < RecoField.size(); coord++)
  {
    RecoField[coord].Write();
  }
  OutputFile -> Close();
  delete OutputFile;
//   gDirectory->GetList()->Delete();
}

void WriteTextFile(std::vector<ThreeVector<float>>& InterpolationData)
{
  std::ofstream OutputFile;
  OutputFile.open("Reco.txt", std::ios::out);
  for(unsigned entry = 0; entry < InterpolationData.size(); entry++)
  {
    OutputFile << InterpolationData[entry][0] << InterpolationData[entry][1] << InterpolationData[entry][2];
  }
  OutputFile.close();
}

// void DrawSpaceCharge()
// {
//   TH2F * hSpaceCharge = new TH2F("Space Charge","Space Charge",kResolution[2],0,kDetector[2],kResolution[0],0,kDetector[0]);
//   for (int ybin = 0; ybin < kResolution[1]; ybin++) 
//   {
//     for (int xbin = 0; xbin < kResolution[0]; xbin++) for (int zbin = 0; zbin < kResolution[2]; zbin++)
//     {
//       hSpaceCharge -> SetBinContent(zbin+1,xbin+1,-fChargeDistribution[xbin][ybin][zbin]);
//     }
//     hSpaceCharge -> SetStats(0);
//     hSpaceCharge -> SetMaximum(1e-8);
//     hSpaceCharge -> SetMinimum(0);
//     TCanvas * C0 = new TCanvas("Space Charge","Space Charge",1000,500);
//     hSpaceCharge -> Draw("colz");
//     C0 -> Print("SpaceCharge.gif+5","gif+5");
//   }
//   delete hSpaceCharge;
// }
// 
// void DrawField()
// {
//   TH2F * hFieldX = new TH2F("Field Map X","Field Map X",kResolution[2],0,kDetector[2],kResolution[0],0,kDetector[0]);
//   TH2F * hFieldY = new TH2F("Field Map Y","Field Map Y",kResolution[2],0,kDetector[2],kResolution[0],0,kDetector[0]);
//   TH2F * hFieldZ = new TH2F("Field Map Z","Field Map Z",kResolution[2],0,kDetector[2],kResolution[0],0,kDetector[0]);
//   for (int ybin = 0; ybin < kResolution[1]; ybin++) 
//   {
//     for (int xbin = 0; xbin < kResolution[0]; xbin++) for (int zbin = 0; zbin < kResolution[2]; zbin++)
//     {
//       hFieldX -> SetBinContent(zbin+1,xbin+1,fRecoField[0][xbin][ybin][zbin]);
//       hFieldY -> SetBinContent(zbin+1,xbin+1,fRecoField[1][xbin][ybin][zbin]);
//       hFieldZ -> SetBinContent(zbin+1,xbin+1,fRecoField[2][xbin][ybin][zbin]);
//     }
//     hFieldX -> SetStats(0);
//     hFieldY -> SetStats(0);
//     hFieldZ -> SetStats(0);
//     
//     hFieldX -> SetMaximum(0.02);
//     hFieldX -> SetMinimum(-0.02);
//     hFieldY -> SetMaximum(0.07);
//     hFieldY -> SetMinimum(-0.07);
//     hFieldZ-> SetMaximum(0.05);
//     hFieldZ -> SetMinimum(-0.05);
//     
//     TCanvas * C1 = new TCanvas("Field Map X","Field Map X",1000,500);
//     hFieldX -> Draw("colz");
//     C1 -> Print("FieldX.gif+5","gif+5");
//     TCanvas * C2 = new TCanvas("Field Map Y","Field Map Y",1000,500);
//     hFieldY -> Draw("colz");
//     C2 -> Print("FieldY.gif+5","gif+5");
//     TCanvas * C3 = new TCanvas("Field Map Z","Field Map Z",1000,500);
//     hFieldZ -> Draw("colz");
//     C3 -> Print("FieldZ.gif+5","gif+5");
//     
//     delete C1;
//     delete C2;
//     delete C3;
//   }
//   delete hFieldX;
//   delete hFieldY;
//   delete hFieldZ;
// }
// 
// void FillHisto ()
// { 
//   TH3F *DistMapX = new TH3F("Distortion_Field_X","Distortion Field X",kResolution[0],0,kDetector[0],kResolution[1],0,kDetector[1],kResolution[2],0,kDetector[2]);
//   TH3F *DistMapY = new TH3F("Distortion_Field_Y","Distortion Field Y",kResolution[0],0,kDetector[0],kResolution[1],0,kDetector[1],kResolution[2],0,kDetector[2]);
//   TH3F *DistMapZ = new TH3F("Distortion_Field_Z","Distortion Field Z",kResolution[0],0,kDetector[0],kResolution[1],0,kDetector[1],kResolution[2],0,kDetector[2]);
//   
//   for (int xbin = 0; xbin < kResolution[0]; xbin++) for (int ybin = 0; ybin < kResolution[1]; ybin++) for(int zbin = 0; zbin < kResolution[2]; zbin++)
//   {
//     DistMapX -> SetBinContent(xbin+1,ybin+1,zbin+1,fRecoField[0][xbin][ybin][zbin]);
//     DistMapY -> SetBinContent(xbin+1,ybin+1,zbin+1,fRecoField[1][xbin][ybin][zbin]);
//     DistMapZ -> SetBinContent(xbin+1,ybin+1,zbin+1,fRecoField[2][xbin][ybin][zbin]);
//   }
//   TFile *OutputFile = new TFile("Field.root", "recreate");
//   DistMapX -> Write();
//   DistMapY -> Write();
//   DistMapZ -> Write(); 
//   OutputFile -> Close();
//   
//   delete DistMapX;
//   delete DistMapY;
//   delete DistMapZ;
// }