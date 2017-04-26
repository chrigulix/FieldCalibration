#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <tuple>
#include <array>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <string>
#include <iomanip>

#include <TFile.h>
#include <TH3.h>

#include "ThreeVector.hpp"
#include "TPCVolumeHandler.hpp"
#include <TVector3.h>

#ifndef LASERTRACK_H
#define LASERTRACK_H

class LaserTrack
{
// public:
private:
  unsigned NumberOfTracksegments;
  std::array<float,2> TrackAngles;
  ThreeVector<float> LaserPosition;
  ThreeVector<float> EntryPoint;
  ThreeVector<float> ExitPoint;
  
//   std::vector<ThreeVector<float>> LaserTrue;
  std::vector<ThreeVector<float>> LaserReco;
  std::vector<ThreeVector<float>> LaserDisplacement;
  
  void FillTrack();
  void FindBoundaries(const TPCVolumeHandler&);
  
public:
  LaserTrack();
  
  // Constructor with reco data input with entry point, exit point, and reco track 
  LaserTrack(const TVector3& InEntryPoint, const TVector3& InExitPoint, const std::vector<TVector3>& RecoTrack);
  LaserTrack(std::array<float,2>&, ThreeVector<float>&, const TPCVolumeHandler&);
  LaserTrack(const unsigned int,std::array<float,2>&, ThreeVector<float>&, const TPCVolumeHandler&);
  
  ThreeVector<float> GetPoyntingVector();
  void DistortTrack(std::string, const TPCVolumeHandler&);
  void CalcDisplacement();
  void AddDisplToReco();
  void AddToDisplacement(ThreeVector<float>&, unsigned long);
  
  ThreeVector<float> GetLaserPosition();
  std::array<float,2> GetAngles();
  unsigned long GetNumberOfSamples() const;
  ThreeVector<float> GetSamplePosition(const unsigned int&) const;
  ThreeVector<float> GetDisplacement(const unsigned int&) const;
  ThreeVector<float> GetEntryPoint();
  ThreeVector<float> GetExitPoint();
  void AppendSample(ThreeVector<float>& SamplePosition, ThreeVector<float>& SampleDisplacement);
  void AppendSample(ThreeVector<float>&);
  void AppendSample(float SamplePos_x, float SamplePos_y, float SamplePos_z, float SampleCorr_x, float SampleCorr_y, float SampleCorr_z);
  void AppendSample(float,float,float);
  
  static void DistortTracks(std::vector<LaserTrack>&, const std::string&, const TPCVolumeHandler&);  
};

#endif