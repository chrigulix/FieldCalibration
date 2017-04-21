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

// #include "TTreeReader.h"
// #include "TTreeReaderValue.h"

// Own Files
// #include "Geometry.hpp"
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
// void DrawSpaceCharge();
// void DrawField();

const double kPi = 3.14159265358979323846;
const float kDetector[] = {2.5604,2.325,10.368}; // Size of the detectors in m (used from docdb 1821)
const int kResolution[] = {26,26,101}; // Field resolution in every coordinate
const int kTrackRes = 100; // Number of samples in every laser track
const int kSqrtNumberOfTracks = 101;
const float kLaserOffset[3] = {1.5,1.1625,-0.21}; // Laser head offset from field cage in [m]

// Main function
int main(int argc, char** argv)
{
    std::cout << argc << std::endl;
    
    // If there are to few input arguments abord
    if(argc < 2)
    {
        std::cerr << "ERROR: Too fiew arguments, use ./LaserCal <input file names>" << std::endl;
        return -1;
    }

    // Coose detector dimensions, coordinate system offset and resolutions
    ThreeVector<float> DetectorSize = {2.5604,2.325,10.368};
    ThreeVector<float> DetectorOffset = {0.0,-DetectorSize[1]/(float)2.0,0.0};
    ThreeVector<int> DetectorResolution = {26,26,101};
   
    // Create the detector volume
    TPCVolumeHandler Detector(DetectorSize, DetectorOffset, DetectorResolution);
    
    // Read data and store it to a Laser object
    std::cout << "Reading data..." << std::endl;
    ReadRecoTracks(argc,argv);
    
  double HistRange[3][2];
  for (int coord = 0; coord < 3; coord++) for(int minmax = 0; minmax < 2; minmax++)
  {
    HistRange[coord][minmax] = (DetectorSize[coord]+pow(2,-25))*( minmax - pow(-1,minmax)/((double)DetectorResolution[coord]-1.0)/2.0 );
  }
  
//   auto start = std::chrono::high_resolution_clock::now();
//   auto elapsed = std::chrono::high_resolution_clock::now() - start;
//   std::cout << "Function Time [ns]: " << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count() << std::endl;
    
  time_t timer;
  std::time(&timer);
  // Enter Functions here:
  unsigned BeamBins = 101;
  std::array<float,2> LaserAngles_1;
  std::array<float,2> LaserAngles_2;
  ThreeVector<float> LaserPosition_1 = {1.3,0.0,-0.21};
  ThreeVector<float> LaserPosition_2 = {1.3,0.0,Detector.GetDetectorSize().at(2)+ (float)0.21};
  
//   std::vector<Laser> LaserTrackSets = ReadMooneyTracks("laserDataSCE.root",Detector);
  
  std::vector<Laser> LaserTrackSets;
  LaserTrackSets.resize(2);
  
//   std::vector<LaserTrack> TrackVector;
  for(unsigned theta_entry = 0; theta_entry < BeamBins; theta_entry++)
  {
    LaserAngles_1[0] = -45.0 + 90.0/(float)(BeamBins-1)*theta_entry;
    LaserAngles_2[0] = -45.0 + 90.0/(float)(BeamBins-1)*theta_entry;
    for(unsigned phi_entry = 0; phi_entry < BeamBins; phi_entry++)
    {
      LaserAngles_1[1] = -45 + 90.0/(float)(BeamBins-1)*phi_entry;
      LaserAngles_2[1] = -225 + 90.0/(float)(BeamBins-1)*phi_entry;
      
//       LaserTrackSets.at(0).AppendTrack(LaserTrack(50, LaserAngles_1, LaserPosition_1, Detector));
//       LaserTrackSets.at(1).AppendTrack(LaserTrack(50, LaserAngles_2, LaserPosition_2, Detector));
    }
  }
  
  for(auto& laser : LaserTrackSets)
  {
//     laser.DistortTrackSet("Field.root",Detector);
  }
  
  std::vector<Delaunay> MeshVector;
  
  for(unsigned set_no = 0; set_no < LaserTrackSets.size(); set_no++)
  {
    LaserTrackSets.at(set_no).CorrectTrackSet();
    MeshVector.push_back( TrackMesher(LaserTrackSets.at(set_no).GetTrackSet()) );
  }
  
//   LaserTrackSets.front().DrawTrack(1000);
  
  Laser TempLaser_1 = LaserTrackSets.front();
  Laser TempLaser_2 = LaserTrackSets.back();
  
//   std::cout << "Start 1st track interpolation" << std::endl;
//   LaserInterpThread(LaserTrackSets.front(),TempLaser_2,MeshVector.back());
//   std::thread Thread_0( LaserInterpThread, std::ref(LaserTrackSets.front()), std::ref(TempLaser_2), std::ref(MeshVector.back()) );
  
//   std::cout << "Start 2nd track interpolation" << std::endl;
//   LaserInterpThread(LaserTrackSets.back(),TempLaser_1,MeshVector.front());
//   std::thread Thread_1 ( LaserInterpThread, std::ref(LaserTrackSets.back()), std::ref(TempLaser_1), std::ref(MeshVector.front()) );
  
//   Thread_0.join();
//   Thread_1.join();
  
//   for(unsigned track_no = 0; track_no < LaserTrackSets.front().GetNumberOfTracks(); track_no++)
//   {
//     InterpolateTrack(LaserTrackSets.front().GetTrack(track_no),LaserTrackSets.back().GetTrackSet(),MeshVector.back());
//   }
  
  
  
//   for(unsigned track_no = 0; track_no < LaserTrackSets.back().GetNumberOfTracks(); track_no++)
//   {
//     InterpolateTrack(LaserTrackSets.back().GetTrack(track_no),TempLaser.GetTrackSet(),MeshVector.front());
//   }
  
  std::cout << "Cleanup!" << std::endl;
  
  MeshVector.clear();
  
  std::vector<LaserTrack> TrackVector = LaserTrackSets.back().GetTrackSet();
  LaserTrackSets.pop_back();
  TrackVector.insert(TrackVector.begin(),LaserTrackSets.front().begin(), LaserTrackSets.front().end());
  LaserTrackSets.clear();
  
  
  std::cout << "Start field map interpolation" << std::endl;
  
  Delaunay Mesh = TrackMesher(TrackVector);

  
//   for(unsigned iter = 0; iter < TrackVector.size(); iter++)
//   {
//     for(unsigned coord = 0; coord < 3; coord++)
//       std::cout << TrackVector[iter].LaserCorrection[20][coord] << " " << TrackVector[iter].LaserCorrection[20][coord] << " ";
//     std::cout << std::endl;
//   }

  std::vector<ThreeVector<float>> Displacement;
  
//   std::cout << "Create Mesh..." << std::endl;
//   Delaunay Mesh_1 = TrackMesher(TrackVector);
//   std::cout << "Mesh done!" << std::endl;
  
//   Point bla(1.0,1.0,3.0);
//   ThreeVector<float> blabla = {1.0,1.0,3.0};
//   std::vector<ThreeVector<float>> shit;
//   shit.push_back(blabla);
//   Delaunay::Cell_handle Cell =  Mesh.locate(bla);
  
//   for(unsigned vertex_no = 0; vertex_no < 4; vertex_no++)
//   {
//     std::cout << Cell->vertex(vertex_no)->point()[0] << " " << Cell->vertex(vertex_no)->point()[1] << " " << Cell->vertex(vertex_no)->point()[2] << std::endl;
//   }
  
//   InterpolateCGAL(TrackVector,Mesh,Location);
  
//   ThreeVector<float> Location = {1.2,1.0,5.1};
  
  ThreeVector<float> Location;
  for(unsigned xbin = 0; xbin < DetectorResolution[0]; xbin++) 
  {
    std::cout << xbin << std::endl;
    Location[0] = DetectorOffset[0]+(HistRange[0][1]-HistRange[0][0])/(float)DetectorResolution[0]*xbin;
    for(unsigned ybin = 0; ybin < DetectorResolution[1]; ybin++) 
    {
      Location[1] = DetectorOffset[1]+(HistRange[1][1]-HistRange[1][0])/(float)DetectorResolution[1]*ybin;
      for(unsigned zbin = 0; zbin < DetectorResolution[2]; zbin++)
      {
	Location[2] = DetectorOffset[2]+(HistRange[2][1]-HistRange[2][0])/(float)DetectorResolution[2]*zbin;
	Displacement.push_back(InterpolateCGAL(TrackVector,Mesh,Location));
      }
    }
  }
  
  WriteRootFile(Displacement,Detector);
//   WriteTextFile(Displacement);
  
//   Displacement = InterpolateCGAL(TrackVector,Mesh,Location);
//   ThreeVector<float> Location = {1.2,1.0,5.1};
//   Displacement = InterpolateCGAL(TrackVector,Mesh,Location);
//   GetClosestTracksInfo(TrackVector,4);
  
  std::cout << "End of program after "<< std::difftime(std::time(NULL),timer) << " s" << std::endl;
//   App->Run();
}

Laser ReadRecoTracks(int argc, char** argv)
{
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
        }
        
        EntryPoint.Print();
        for(const auto& sample : TrackSamples) sample.Print();
    }
    else // If the trees don't have the same amount of entries, through error (I know not propper error handling)
    {
        std::cerr << "ERROR: Two TTrees don't have the same ammount of entries!" << std::endl;
    }
    
//     LaserInfoTree->Print();
//     RecoTrackTree->Print();
}

void LaserInterpThread(Laser& LaserTrackSet, const Laser& InterpolationLaser, const Delaunay& InterpolationMesh)
{
  LaserTrackSet.InterpolateTrackSet(InterpolationLaser, InterpolationMesh);
}

void WriteRootFile(std::vector<ThreeVector<float>>& InterpolationData, TPCVolumeHandler& TPCVolume)
{ 
  double HistRange[3][2];
  for (int coord = 0; coord < 3; coord++) for(int minmax = 0; minmax < 2; minmax++)
  {
    HistRange[coord][minmax] = (TPCVolume.GetDetectorSize()[coord]+pow(2,-25))*( minmax - pow(-1,minmax)/((double)TPCVolume.GetDetectorResolution()[coord]-1.0)/2.0 );
  }
  
  ThreeVector<float> DetectorOffset = {0.0,-TPCVolume.GetDetectorSize()[1]/(float)2.0,0.0};
  
  std::vector<TH3F*> RecoField;
  RecoField.push_back(new TH3F("Reco_Field_X","Reco Field X",TPCVolume.GetDetectorResolution()[0],HistRange[0][0],HistRange[0][1],TPCVolume.GetDetectorResolution()[1],HistRange[1][0],HistRange[1][1],TPCVolume.GetDetectorResolution()[2],HistRange[2][0],HistRange[2][1]));
  RecoField.push_back(new TH3F("Reco_Field_Y","Reco Field Y",TPCVolume.GetDetectorResolution()[0],HistRange[0][0],HistRange[0][1],TPCVolume.GetDetectorResolution()[1],HistRange[1][0],HistRange[1][1],TPCVolume.GetDetectorResolution()[2],HistRange[2][0],HistRange[2][1]));
  RecoField.push_back(new TH3F("Reco_Field_Z","Reco Field Z",TPCVolume.GetDetectorResolution()[0],HistRange[0][0],HistRange[0][1],TPCVolume.GetDetectorResolution()[1],HistRange[1][0],HistRange[1][1],TPCVolume.GetDetectorResolution()[2],HistRange[2][0],HistRange[2][1]));
  
  
  for(unsigned xbin = 0; xbin < TPCVolume.GetDetectorResolution()[0]; xbin++) for(unsigned ybin = 0; ybin < TPCVolume.GetDetectorResolution()[1]; ybin++) for(unsigned zbin = 0; zbin < TPCVolume.GetDetectorResolution()[2]; zbin++)
  {
    for(unsigned coord = 0; coord < 3; coord++)
    {
      RecoField[coord]->SetBinContent(xbin+1,ybin+1,zbin+1, InterpolationData[zbin+( TPCVolume.GetDetectorResolution()[2]*(ybin+TPCVolume.GetDetectorResolution()[1]*xbin) )][coord]);
      std::cout << InterpolationData[zbin+( TPCVolume.GetDetectorResolution()[2]*(ybin+TPCVolume.GetDetectorResolution()[1]*xbin) )][coord] << " ";
    }
    std::cout << std::endl;
  }
  
  TFile *OutputFile = new TFile("RecoField.root", "recreate");
  for(unsigned coord = 0; coord < RecoField.size(); coord++)
  {
    RecoField[coord]->Write();
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