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
  std::vector<ThreeVector<float>> LaserCorrection;
  
  void FillTrack();
  void FindBoundaries(const TPCVolumeHandler&);
  
public:
  LaserTrack();
  LaserTrack(std::array<float,2>&, ThreeVector<float>&, const TPCVolumeHandler&);
  LaserTrack(const unsigned int,std::array<float,2>&, ThreeVector<float>&, const TPCVolumeHandler&);
  
  ThreeVector<float> GetPoyntingVector();
  void DistortTrack(std::string, const TPCVolumeHandler&);
  void CorrectTrack();
  void AddToCorrection(ThreeVector<float>&, unsigned long);
  
  ThreeVector<float> GetLaserPosition();
  std::array<float,2> GetAngles();
  unsigned long GetNumberOfSamples() const;
  ThreeVector<float> GetSamplePosition(const unsigned int&) const;
  ThreeVector<float> GetCorrection(const unsigned int&) const;
  ThreeVector<float> GetEntryPoint();
  ThreeVector<float> GetExitPoint();
  void AppendSample(ThreeVector<float>& SamplePosition, ThreeVector<float>& SampleCorrection);
  void AppendSample(ThreeVector<float>&);
  void AppendSample(float SamplePos_x, float SamplePos_y, float SamplePos_z, float SampleCorr_x, float SampleCorr_y, float SampleCorr_z);
  void AppendSample(float,float,float);
  
  static void DistortTracks(std::vector<LaserTrack>&, const std::string&, const TPCVolumeHandler&);  
};

#endif