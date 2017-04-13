#include <array>
#include <vector>
#include <cmath>
#include "ThreeVector.hpp"

#ifndef TPCVOLUMEHANDLER_H
#define TPCVOLUMEHANDLER_H

class TPCVolumeHandler
{
public:
  TPCVolumeHandler();
  TPCVolumeHandler(const ThreeVector<float>&, const ThreeVector<float>&, const ThreeVector<int>&);
  TPCVolumeHandler(const std::array<float,3>&, const std::array<float,3>&, const std::array<int,3>&);
  TPCVolumeHandler(float[3],float[3],int[3]);
  
  std::vector<ThreeVector<float>> GetNormalVectors();
  std::vector<ThreeVector<float>> GetNormVectorOffset();
  
  ThreeVector<float> GetDetectorSize();
  ThreeVector<float> GetDetectorOffset();
  ThreeVector<int> GetDetectorResolution();
  
private:
  ThreeVector<float> DetectorSize;
  ThreeVector<float> DetectorOffset;
  ThreeVector<int> DetectorResolution;
  
  void CalcNormal(ThreeVector<float>&, ThreeVector<float>&);
  
  std::vector<ThreeVector<float>> NormalVector;
  std::vector<ThreeVector<float>> NormalVectorOffSet;
  
};

#endif