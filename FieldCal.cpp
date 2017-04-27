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

// Initialize functions defined below
Laser ReadRecoTracks(int argc, char** argv);
void WriteRootFile(std::vector<ThreeVector<float>>&, TPCVolumeHandler&);
void WriteTextFile(std::vector<ThreeVector<float>>&);
void LaserInterpThread(Laser&, const Laser&, const Delaunay&);

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
    std::cout << "Find track displacements... " << std::endl;
    // Choose displacement algorithm (available so far: TrackDerivative or ClosestPoint)
    LaserTrackSet.CalcDisplacement(LaserTrack::ClosestPoint);
    
    // Add displacement to reconstructed track to change to detector coordinates (only for map generation)
    LaserTrackSet.AddDisplToReco();
    
    // Create delaunay mesh
    std::cout << "Generate mesh..." << std::endl;
    Delaunay Mesh = TrackMesher(LaserTrackSet.GetTrackSet());
    
    // Interpolate Displacement Map (regularly spaced grid)
    std::cout << "Start interpolation..." << std::endl;
    std::vector<ThreeVector<float>> DisplacementMap = InterpolateMap(LaserTrackSet.GetTrackSet(),Mesh,Detector);
    
    // Fill displacement map into TH3 histograms and write them to file
    std::cout << "Write to File ..." << std::endl;
    WriteRootFile(DisplacementMap,Detector);
    
    std::cout << "End of program after "<< std::difftime(std::time(NULL),timer) << " s" << std::endl;
} // end main


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
} // end ReadRecoTracks


void WriteRootFile(std::vector<ThreeVector<float>>& InterpolationData, TPCVolumeHandler& TPCVolume)
{ 
    // Store TPC properties which are important for the TH3 generation
    ThreeVector<unsigned long> Resolution = TPCVolume.GetDetectorResolution();
    ThreeVector<float> MinimumCoord = TPCVolume.GetMapMinimum();
    ThreeVector<float> MaximumCoord = TPCVolume.GetMapMaximum();
  
    // Initialize all TH3F
    std::vector<TH3F> RecoDisplacement;
    RecoDisplacement.push_back(TH3F("Reco_Displacement_X","Reco Displacement X",Resolution[0],MinimumCoord[0],MaximumCoord[0],Resolution[1],MinimumCoord[1],MaximumCoord[1],Resolution[2],MinimumCoord[2],MaximumCoord[2]));
    RecoDisplacement.push_back(TH3F("Reco_Displacement_Y","Reco Displacement Y",Resolution[0],MinimumCoord[0],MaximumCoord[0],Resolution[1],MinimumCoord[1],MaximumCoord[1],Resolution[2],MinimumCoord[2],MaximumCoord[2]));
    RecoDisplacement.push_back(TH3F("Reco_Displacement_Z","Reco Displacement Z",Resolution[0],MinimumCoord[0],MaximumCoord[0],Resolution[1],MinimumCoord[1],MaximumCoord[1],Resolution[2],MinimumCoord[2],MaximumCoord[2]));
  

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
                    RecoDisplacement[coord].SetBinContent(xbin+1,ybin+1,zbin+1, InterpolationData[zbin+( TPCVolume.GetDetectorResolution()[2]*(ybin+TPCVolume.GetDetectorResolution()[1]*xbin) )][coord]);
                } // end coordinate loop
            } // end zbin loop
        } // end ybin loop
    } // end zbin loop
    
    // Open and recreate output file
    TFile OutputFile("RecoDispl.root", "recreate");
    
    // Loop over space coordinates
    for(unsigned coord = 0; coord < RecoDisplacement.size(); coord++)
    {
        // Write every TH3 map into file
        RecoDisplacement[coord].Write();
    }
    
    // Close output file and clean up
    OutputFile.Close();
    gDirectory->GetList()->Delete();
}


void WriteTextFile(std::vector<ThreeVector<float>>& InterpolationData)
{
    // Initialize stream to file
    std::ofstream OutputFile;
  
    // Open output file
    OutputFile.open("Reco.txt", std::ios::out);
  
    // Loop over all interpolated data points
    for(unsigned entry = 0; entry < InterpolationData.size(); entry++)
    {
        // Write every point into a seperate line
        OutputFile << InterpolationData[entry][0] << InterpolationData[entry][1] << InterpolationData[entry][2];
    }

    // Close file
    OutputFile.close();
} // WriteRootFile


// This is the multi-threading interpolation function. Just hangs out here, for legacy purposes 
void LaserInterpThread(Laser& LaserTrackSet, const Laser& InterpolationLaser, const Delaunay& InterpolationMesh)
{
  LaserTrackSet.InterpolateTrackSet(InterpolationLaser, InterpolationMesh);
} // LaserInterpThread