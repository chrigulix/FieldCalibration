#include "../include/LaserTrack.hpp"

// Default constructor
LaserTrack::LaserTrack()
{

}

LaserTrack::LaserTrack(std::array<float,2>& Angles, ThreeVector<float>& Position, TPCVolumeHandler& TPCVolume) : TrackAngles(Angles) , LaserPosition(Position)
{  
  FindBoundaries(TPCVolume);
}


// This Constructor initializes the laser track vectors and fills them according to the angles and segments
LaserTrack::LaserTrack(const unsigned int NSegments, std::array<float,2>& Angles, ThreeVector<float>& Position, TPCVolumeHandler& TPCVolume) : NumberOfTracksegments(NSegments) , TrackAngles(Angles) , LaserPosition(Position)
{
  // Convert degrees to radian
  for(unsigned angle = 0; angle<TrackAngles.size(); angle++)
    TrackAngles[angle] *= M_PI/180.0;

//   LaserTrue.resize(NumberOfTracksegments+1);
  LaserReco.resize(NumberOfTracksegments+1);
  LaserCorrection.resize(NumberOfTracksegments+1);
  
  FindBoundaries(TPCVolume);
  FillTrack();
}

ThreeVector<float> LaserTrack::GetPoyntingVector()
{
  // Create unit vector of Poynting vector from spherical coordinates
  ThreeVector<float> vec_res;
  vec_res[0] = (float)(std::cos(TrackAngles[0])*std::sin(TrackAngles[1])); // cos(theta)*cos(phi)
  vec_res[1] = (float)std::sin(TrackAngles[0]); // sin(theta)
  vec_res[2] = (float)(std::cos(TrackAngles[0])*std::cos(TrackAngles[1])); // cos(theta)*cos(phi)
  
  return vec_res;
}

void LaserTrack::FillTrack()
{
  ThreeVector<float> Poynting = GetPoyntingVector();
  
  LaserReco.front() = EntryPoint;
  LaserReco.back() = ExitPoint;
  
  // Parameter length of a single track segment
  float SegmentParameter = ThreeVector<float>::VectorNorm(ExitPoint - EntryPoint) / (float)NumberOfTracksegments;
  
  // Fill segments
  for(unsigned segment = 1; segment < NumberOfTracksegments; segment++)
  {
    LaserReco[segment] = EntryPoint + (SegmentParameter*segment)*Poynting;
  }
}

void LaserTrack::FindBoundaries(TPCVolumeHandler& TPCVolume)
{
  ThreeVector<float> Poynting = GetPoyntingVector();
  
  // This is the absolute error allowed for the boundary conditions of the laser tracks
  float epsilon_abs = 1e-6;
  
//   TPCVolumeHandler TPCVolume(DetectorSize, DetectorOffset);
  auto SurfaceNormVec = TPCVolume.GetNormalVectors();
  auto NormVecOffset = TPCVolume.GetNormVectorOffset();
  
  std::vector<float> BoundaryParameter;
  std::vector<ThreeVector<float>> BoundaryPoint;
  
  for(unsigned NSurface = 0; NSurface < SurfaceNormVec.size(); NSurface++)
  {
    float IntersecParameter;
    if(ThreeVector<float>::DotProduct(SurfaceNormVec[NSurface],Poynting) != 0.0)
      IntersecParameter = -ThreeVector<float>::DotProduct(SurfaceNormVec[NSurface],LaserPosition-NormVecOffset[NSurface]) / ThreeVector<float>::DotProduct(SurfaceNormVec[NSurface],Poynting);
    else
      IntersecParameter = -0xDEADBEEF;
    
    // Calculate the intersection point (UnitVector is added and suptracted to prevent rounding errors)
    ThreeVector<float> IntersecPoint = LaserPosition + Poynting * IntersecParameter;
    bool IntersectionFlag = true;
    
    for(unsigned coord = 0; coord < IntersecPoint.size(); coord++)
    {
      // Check if Intersection Point is within the surface area with an error of epsilon_abs
      if( IntersecPoint[coord] > TPCVolume.GetDetectorSize()[coord]+TPCVolume.GetDetectorOffset()[coord]+epsilon_abs || IntersecPoint[coord] < TPCVolume.GetDetectorOffset()[coord]-epsilon_abs || IntersecParameter < 0.0 )
      {
	IntersectionFlag = false;
	break;
      }
    }
    if(IntersectionFlag)
    {
      BoundaryParameter.push_back(IntersecParameter);
      BoundaryPoint.push_back(IntersecPoint);
    }
  }
  
  // If there are no hits on the TPC the laser track will be initialized empty with an error
  if(!BoundaryParameter.size())
  {
//     LaserTrue.clear();
    LaserReco.clear();
    std::cout << "ERROR: Laser beam with start position = (" << LaserPosition[0] << ", " << LaserPosition[1] << ", " << LaserPosition[2] << ") "
	      << "and track angles (theta, phi) = (" << TrackAngles[0] << ", " << TrackAngles[1] << ") "
	      << "is missing the chamber!" << std::endl;
    return;
  }
  else if(BoundaryParameter.size() < 2) // If the laser beam is produced in the volume add the laser position and the parameter as boundary conditions
  {
    BoundaryParameter.push_back(0.0);
    BoundaryPoint.push_back(LaserPosition);
    std::cout << "Warning: Laser beam with start position = (" << LaserPosition[0] << ", " << LaserPosition[1] << ", " << LaserPosition[2] << ") "
	      << "and track angles (theta, phi) = (" << TrackAngles[0]/M_PI*180.0 << ", " << TrackAngles[1]/M_PI*180.0 << ") "
	      << "is placed in the chamber!" << std::endl;
  }
  
//   std::cout << BoundaryPoint.size() << std::endl;
//   for(auto iter : BoundaryPoint)
//     std::cout << iter[0] << " " << iter[1] << " " << iter[2] << std::endl;
  
  // Check order of boundary points and swap if they are in wrong order
  if(BoundaryParameter[0] > BoundaryParameter[1])
  {
    std::swap(BoundaryPoint.front(),BoundaryPoint.back());
    std::swap(BoundaryParameter.front(),BoundaryParameter.back());
  }
  
  // Fill the entry and exit point
  EntryPoint = BoundaryPoint.front();
  ExitPoint = BoundaryPoint.back();
}


void LaserTrack::DistortTrack(std::string MapFileName, TPCVolumeHandler& TPCVolume)
{
  float epsilon_abs = 1e-6;
  
  TFile* FieldFile = new TFile(MapFileName.c_str());
  if (FieldFile->IsZombie()) 
  {
       std::cout << "ERROR: " << MapFileName << " can not be opened!" << std::endl;
       return;
  }
  
  std::vector<TH3F*> DistortMap;
  DistortMap.push_back((TH3F*) FieldFile->Get("Distortion_Field_X"));
  DistortMap.push_back((TH3F*) FieldFile->Get("Distortion_Field_Y"));
  DistortMap.push_back((TH3F*) FieldFile->Get("Distortion_Field_Z"));
  
  float InterpolDistortion = 0.0;
  
  // Add field distortion to every track segment point
  for(unsigned segment = 0; segment < NumberOfTracksegments+1; segment++)
  {
    for(unsigned coord = 0; coord < LaserReco[segment].size(); coord++)
    {
      if(InterpolDistortion = DistortMap[coord]->Interpolate(LaserReco[segment][0],LaserReco[segment][1]+TPCVolume.GetDetectorSize()[1]/(float)2.0,LaserReco[segment][2]))
      {
	LaserReco[segment][coord] += InterpolDistortion;
      }
    }
  }
  for (int coord = 0; coord < LaserReco[0].size(); coord++) delete DistortMap[coord];
  FieldFile->Close();
  delete FieldFile;
  gDirectory->GetList()->Delete();
  
//   for(unsigned segment = 0; segment < NumberOfTracksegments+1; segment++)
//   {
//     std::cout << LaserTrue[segment][0] << " " << LaserTrue[segment][1] << " " << LaserTrue[segment][2] << std::endl;
//     std::cout << LaserReco[segment][0] << " " << LaserReco[segment][1] << " " << LaserReco[segment][2] << std::endl;
//   }
}

void LaserTrack::CorrectTrack()
{
  float TrueTrackPara;
  
//   std::cout << LaserReco.size() << std::endl;
  
  if(LaserReco.size() > 2)
  {
    // Correct all intermediate track samples with the perpendicularity method
    for(unsigned sample_no = 1; sample_no < LaserReco.size(); sample_no++)
    {
      TrueTrackPara = ThreeVector<float>::DotProduct(LaserReco[sample_no-1]-LaserReco[sample_no+1],LaserReco[sample_no]-EntryPoint)
		    / ThreeVector<float>::DotProduct(LaserReco[sample_no-1]-LaserReco[sample_no+1],ExitPoint-EntryPoint);

      LaserCorrection[sample_no] = EntryPoint - LaserReco[sample_no] + TrueTrackPara*(ExitPoint-EntryPoint);
    }
  }
  
  if(LaserReco.size() > 1)
  {
    // Correct first and last track samples
    LaserCorrection.front() =  EntryPoint - LaserReco.front();
    LaserCorrection.back() =  ExitPoint - LaserReco.back();
  }
}

void LaserTrack::AddToCorrection(ThreeVector<float>& AdditionVector, unsigned long sample)
{
  LaserCorrection[sample] += AdditionVector;
}

ThreeVector<float> LaserTrack::GetLaserPosition()
{
  return LaserPosition;
}

std::array<float,2> LaserTrack::GetAngles()
{
  return TrackAngles;
}

unsigned long LaserTrack::GetNumberOfSamples() const
{
  return LaserReco.size();  
}

ThreeVector<float> LaserTrack::GetSamplePosition(const unsigned int& SampleNumber) const
{
  return LaserReco[SampleNumber];
}

ThreeVector< float > LaserTrack::GetCorrection(const unsigned int& SampleNumber) const
{
  return LaserCorrection[SampleNumber];
}

ThreeVector<float> LaserTrack::GetEntryPoint()
{
  return EntryPoint;
}

ThreeVector<float> LaserTrack::GetExitPoint()
{
  return ExitPoint;
}



void LaserTrack::AppendSample(ThreeVector<float>& SamplePosition)
{
  LaserReco.push_back(SamplePosition);
  LaserCorrection.push_back(ThreeVector<float>(0.0,0.0,0.0));
}

void LaserTrack::AppendSample(float SamplePos_x, float SamplePos_y, float SamplePos_z)
{
  LaserReco.push_back( ThreeVector<float>(SamplePos_x,SamplePos_y,SamplePos_z) );
  LaserCorrection.push_back(ThreeVector<float>(0.0,0.0,0.0));
}

void LaserTrack::DistortTracks(std::vector<LaserTrack>& LaserTracks, const std::string& MapFileName, TPCVolumeHandler& TPCVolume)
{
  float epsilon_abs = 1e-6;
  
  TFile* FieldFile = new TFile(MapFileName.c_str());
  if (FieldFile->IsZombie()) 
  {
       std::cout << "ERROR: " << MapFileName << " can not be opened!" << std::endl;
       return;
  }
  
  std::vector<TH3F*> DistortMap;
  DistortMap.push_back((TH3F*) FieldFile->Get("Distortion_Field_X"));
  DistortMap.push_back((TH3F*) FieldFile->Get("Distortion_Field_Z"));
  DistortMap.push_back((TH3F*) FieldFile->Get("Distortion_Field_Y"));
  
  float InterpolDistortion = 0.0;
  
  for(unsigned track_no = 0; track_no < LaserTracks.size(); track_no++)
  {
//     // Add field distortion to every track segment point
    for(unsigned segment = 0; segment < LaserTracks[track_no].NumberOfTracksegments + 1; segment++)
    {
      for(unsigned coord = 0; coord < LaserTracks[track_no].LaserReco[segment].size(); coord++)
      {
	LaserTracks[track_no].LaserReco[segment][coord] += DistortMap[coord]->Interpolate(
	LaserTracks[track_no].LaserReco[segment][0]-TPCVolume.GetDetectorOffset().at(0),
	LaserTracks[track_no].LaserReco[segment][1]-TPCVolume.GetDetectorOffset().at(1),
	LaserTracks[track_no].LaserReco[segment][2]-TPCVolume.GetDetectorOffset().at(2)
	);
// 	  std::cout << LaserTracks[track_no].LaserTrue[segment][coord] << " " << LaserTracks[track_no].LaserReco[segment][coord] << " ";
// 	  std::cout << DistortMap[coord]->Interpolate(LaserTracks[track_no].LaserTrue[segment][0],LaserTracks[track_no].LaserTrue[segment][1]+DetectorSize[1]/(float)2.0,LaserTracks[track_no].LaserTrue[segment][2]) << " ";
      }
    } 
//       std::cout << std::endl;
  }
//   
  for(unsigned coord = 0; coord < DistortMap.size(); coord++)
  {
    delete DistortMap[coord];
  }
  FieldFile->Close();
  delete FieldFile;
  gDirectory->GetList()->Delete();
}
