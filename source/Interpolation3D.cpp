#include "../include/Interpolation3D.hpp"
#include "../include/Laser.hpp"

std::vector<std::pair<unsigned,float>> GetClosestTracksInfo(std::vector<LaserTrack>& LaserTracks, const unsigned NumberOfClosestTracks)
{
//   std::vector<LaserTrack> ClosestTracks;
//   ClosestTracks.resize(NumberOfClosestTracks);
  
  std::vector<std::pair<unsigned,float>> ClosestTracksInfo;
  for(unsigned NumberOfTracks = 0; NumberOfTracks < NumberOfClosestTracks; NumberOfTracks++)
  {
    ClosestTracksInfo.push_back(std::make_pair(0,0xDEADBEEF));
  }
  
  
  ThreeVector<float> PoyntingVector = {1.0,1.0,3.0};
  PoyntingVector -= LaserTracks.front().GetLaserPosition();
  PoyntingVector /= PoyntingVector.GetNorm();
  std::array<float,2> AnglesOfPoint = AnglesFromPoynting(PoyntingVector);
  
  std::cout << AnglesOfPoint[0]*180/M_PI << " " << AnglesOfPoint[1]*180/M_PI << std::endl;
  
  for(unsigned track = 0; track < LaserTracks.size(); track++)
  {
    float Radius = 0;
    for(unsigned angle = 0; angle < AnglesOfPoint.size(); angle++)
      Radius += std::pow(AnglesOfPoint[angle]-LaserTracks[track].GetAngles()[angle],2);
    
    if(ClosestTracksInfo.back().second > Radius)
    {
      ClosestTracksInfo.back().first = track;
      ClosestTracksInfo.back().second = Radius;
      
      std::sort(ClosestTracksInfo.begin(),ClosestTracksInfo.end(),PairSortFunction);
//       std::cout << ClosestTracksInfo[0].first << " " << ClosestTracksInfo[0].second << " " << ClosestTracksInfo[1].first  << " " << ClosestTracksInfo[1].second << " " 
// 		<< ClosestTracksInfo[2].first << " " << ClosestTracksInfo[2].second /*<< " " << ClosestTracksInfo[3].first  << " " << ClosestTracksInfo[3].second */<< std::endl; 
    }
  }
  
  return ClosestTracksInfo;
}

std::vector<std::pair<unsigned int, unsigned int>> GetClosestLaserSample(std::vector<LaserTrack>& LaserTracks, const unsigned int NumberOfClosestSamples)
{
  ThreeVector<float> InterpolPoint = {1.0,1.0,3.0};
  
  std::vector< std::pair<unsigned int, unsigned int>> ClosestLaserSample;
  ClosestLaserSample.resize(NumberOfClosestSamples);
  
  std::vector<std::pair<unsigned int, float>> ClosestSampleInfo;
  for(unsigned info = 0; info < NumberOfClosestSamples; info++)
    ClosestSampleInfo.push_back(std::make_pair(0,0xDEADBEEF));
  
  std::vector<std::pair<unsigned int,float>> ClosestTracksInfo = GetClosestTracksInfo(LaserTracks,3);
  
  float Radius = 0;
  for(unsigned info = 0; info < ClosestTracksInfo.size(); info++)
  {
    for(unsigned tracksample = 0; tracksample < LaserTracks[ClosestTracksInfo[info].first].GetNumberOfSamples(); tracksample++)
    {
      Radius = (InterpolPoint - LaserTracks[ClosestTracksInfo[info].first].GetSamplePosition(tracksample)).GetNorm();
      if(!info && ClosestSampleInfo[1].second > Radius)
      {
	ClosestSampleInfo[1].first = tracksample;
	ClosestSampleInfo[1].second = Radius;
	
	std::sort(ClosestSampleInfo.begin(),ClosestSampleInfo.begin()+1,PairSortFunction);
// 	std::cout << ClosestSampleInfo.front().first << " " << ClosestSampleInfo.front().second << " " << ClosestSampleInfo[1].first <<" " << ClosestSampleInfo[1].second << std::endl;
      }
      else if(info && ClosestSampleInfo.back().second > Radius)
      {
	ClosestSampleInfo.back().first = tracksample;
	ClosestSampleInfo.back().second = Radius;
	
	std::sort(ClosestSampleInfo.end()-1,ClosestSampleInfo.end(),PairSortFunction);
      }
    }   
  }
  return ClosestLaserSample;
}

std::array<float,2> AnglesFromPoynting(ThreeVector<float>& Poynting)
{
  Poynting /= Poynting.GetNorm();
  
  std::array<float,2> ang_res;
  ang_res[0] = std::asin(Poynting[1]);
  ang_res[1] = std::asin(Poynting[0]/std::cos(ang_res[0]));
  return ang_res;
}

ThreeVector<float> PoyntingFromAngles(const std::array<float,2>& Angles)
{
  // Create unit vector of Poynting vector from spherical coordinates
  ThreeVector<float> vec_res;
  vec_res[0] = (float)(std::cos(Angles[0])*std::sin(Angles[1])); // cos(theta)*cos(phi)
  vec_res[1] = (float)std::sin(Angles[0]); // sin(theta)
  vec_res[2] = (float)(std::cos(Angles[0])*std::cos(Angles[1])); // cos(theta)*cos(phi)
  
  return vec_res;
}

bool PairSortFunction(std::pair<unsigned,float> left_pair, std::pair<unsigned,float> right_pair)
{
  return (left_pair.second < right_pair.second);
}

Delaunay TrackMesher(const std::vector<LaserTrack>& LaserTracks)
{ 
  std::vector< std::pair<Point,std::pair<unsigned long, unsigned long>> > Points;
  
  for(unsigned long track = 0; track < LaserTracks.size(); track++)
  {
    for(unsigned long sample = 0; sample < LaserTracks[track].GetNumberOfSamples(); sample++)
    {
        Point SamplePoint = Point(LaserTracks[track].GetSamplePosition(sample)[0],LaserTracks[track].GetSamplePosition(sample)[1],LaserTracks[track].GetSamplePosition(sample)[2]);
        
        std::pair<unsigned long, unsigned long> SamplePointIndex = std::make_pair(track,sample);
        
        Points.push_back( std::make_pair(SamplePoint,SamplePointIndex) );
    }
  }
  
  Delaunay DelaunayMesh(Points.begin(), Points.end());
  
  return DelaunayMesh;
}

ThreeVector<float> PointToVector(Point& InputPoint)
{
  ThreeVector<float> vec_res = {(float)InputPoint[0], (float)InputPoint[1], (float)InputPoint[2]};
  return vec_res;
}

Point VectorToPoint(ThreeVector<float>& InputVector)
{
  Point point_res(InputVector[0],InputVector[1],InputVector[2]);
  return point_res;
}

ThreeVector<float> InterpolateCGAL(const std::vector<LaserTrack>& LaserTracks, const Delaunay& Mesh, ThreeVector<float> Location)
{
//   std::cout << "Start CGAL" << std::endl;
  
//   ThreeVector<float> Location = {0.0,0.0,0.0};
  std::array<std::pair<unsigned long, unsigned long>,4> PointIndex;
  ThreeVector<float> InterpolatedDispl = {0.0,0.0,0.0};
  
//   if(Test == Location) std::cout << "cool" << std::endl;
//   else std::cout << "shit" << std::endl;
//   Delaunay Mesh = TrackMesher(LaserTracks);
  
//   std::cout << "Mesh done!" << std::endl;
  
//   auto start = std::chrono::highresolution_clock::now();

  Delaunay::Cell_handle Cell =  Mesh.locate(VectorToPoint(Location));
  
  ThreeVector<float> VertexPoint;
  for(unsigned vertex_no = 0; vertex_no < PointIndex.size(); vertex_no++)
  {
    PointIndex[vertex_no] = Cell->vertex(vertex_no)->info();
  }
  
  Matrix3x3 TransMatrix = {{0,0,0},{0,0,0},{0,0,0}};
  
  for(unsigned row = 0; row < 3; row++)
  {
    for(unsigned column = 0; column < 3; column++)
    {
      TransMatrix[row][column] = LaserTracks[PointIndex[column].first].GetSamplePosition(PointIndex[column].second)[row] - LaserTracks[PointIndex.back().first].GetSamplePosition(PointIndex.back().second)[row];
    }
  }

  Location -= LaserTracks[PointIndex.back().first].GetSamplePosition(PointIndex.back().second);
//   std::cout << TransMatrix.Determinant() << std::endl;
  
  std::vector<float> BaryCoord;
  if(TransMatrix.Invert())
  {
    BaryCoord = (TransMatrix * Location).GetStdVector();
    BaryCoord.push_back(1-BaryCoord[0]-BaryCoord[1]-BaryCoord[2]);
//     std::cout << BaryCoord[0] << " " << BaryCoord[1] << " " << BaryCoord[2] << " " << BaryCoord[3] << std::endl;
  }
  else
  {
//     BaryCoord = {0.0,0.0,0.0,0.0};
    InterpolatedDispl = {0.0,0.0,0.0};
    return InterpolatedDispl;
  }
  
  if(BaryCoord[0]< 0.0 || BaryCoord[1]< 0.0 || BaryCoord[2]< 0.0 || BaryCoord[3]< 0.0)
  {
//     std::cout << BaryCoord[0] << " " << BaryCoord[1] << " " << BaryCoord[2] << " " << BaryCoord[3] << std::endl;
    InterpolatedDispl = {0.0,0.0,0.0};
    return InterpolatedDispl;
  }
  
  for(unsigned vertex_no = 0; vertex_no < 4; vertex_no++)
  {
    InterpolatedDispl += (LaserTracks[PointIndex[vertex_no].first].GetCorrection(PointIndex[vertex_no].second) * BaryCoord[vertex_no]);
//     std::cout << LaserTracks[PointIndex[vertex_no].first].GetCorrection(PointIndex[vertex_no].second)[0] << " " << LaserTracks[PointIndex[vertex_no].first].GetCorrection(PointIndex[vertex_no].second)[1] << " " << LaserTracks[PointIndex[vertex_no].first].GetCorrection(PointIndex[vertex_no].second)[2] << std::endl;
//     for(unsigned coord = 0; coord < 3; coord++)  std::cout << (LaserTracks[PointIndex[vertex_no].first].GetCorrection(PointIndex[vertex_no].second)[coord] * BaryCoord[vertex_no]) << " ";
//     std::cout << std::endl;
  }
//   std::cout << InterpolatedDispl[0] << " " << InterpolatedDispl[0] << " " << InterpolatedDispl[0] << std::endl;
//   std::cout << InterpolatedDispl[0] << " " << InterpolatedDispl[1] << " " << InterpolatedDispl[2] << std::endl;
//   for(unsigned row = 0; row < 3; row++)
//   {
//     for(unsigned column = 0; column < 3; column++)
//     {
//       std::cout << TransMatrix[row][column] << " ";
//     }
//     std::cout << std::endl;
//   }
  
//   auto elapsed = std::chrono::high_resolution_clock::now() - start;
  
//   std::cout << "Function Time [ms]: " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() << std::endl;
  
  return InterpolatedDispl;
//   Triangulation::Cell_handle c = D.;
//   CGAL::Delaunay_triangulation_3<InterpKernel, Tds, CGAL::Fast_location> a(Points.begin(), Points.end());
}

void InterpolateTrack(LaserTrack& Track, const std::vector<LaserTrack>& LaserTracks, const Delaunay& Mesh)
{
  std::array<std::pair<unsigned long, unsigned long>,4> PointIndex;
  ThreeVector<float> InterpolatedDispl = {0.0,0.0,0.0};
  ThreeVector<float> Location;
  
  for(unsigned sample_no = 0; sample_no < Track.GetNumberOfSamples(); sample_no++)
  {
    Location = Track.GetSamplePosition(sample_no);
    Delaunay::Cell_handle Cell =  Mesh.locate(VectorToPoint(Location));
    
    ThreeVector<float> VertexPoint;
  
    for(unsigned vertex_no = 0; vertex_no < 4; vertex_no++)  
    {
      PointIndex[vertex_no] = Cell->vertex(vertex_no)->info();
    }
    
    Matrix3x3 TransMatrix = {{0,0,0},{0,0,0},{0,0,0}};
  
    for(unsigned row = 0; row < 3; row++)
    {
      for(unsigned column = 0; column < 3; column++)
      {
	TransMatrix[row][column] = LaserTracks[PointIndex[column].first].GetSamplePosition(PointIndex[column].second)[row] - LaserTracks[PointIndex.back().first].GetSamplePosition(PointIndex.back().second)[row];
      }
    }
    
    Location -= LaserTracks[PointIndex.back().first].GetSamplePosition(PointIndex.back().second);
//     std::cout << TransMatrix.Determinant() << std::endl;
  
    std::vector<float> BaryCoord;
    if(TransMatrix.Invert())
    {
    BaryCoord = (TransMatrix * Location).GetStdVector();
    BaryCoord.push_back(1-BaryCoord[0]-BaryCoord[1]-BaryCoord[2]);
//     std::cout << BaryCoord[0] << " " << BaryCoord[1] << " " << BaryCoord[2] << " " << BaryCoord[3] << std::endl;
    }
    else
    {
//       BaryCoord = {0.0,0.0,0.0,0.0};
      InterpolatedDispl = {0.0,0.0,0.0};
      return;
    }
    
    if(BaryCoord[0]< 0.0 || BaryCoord[1]< 0.0 || BaryCoord[2]< 0.0 || BaryCoord[3]< 0.0)
    {
//       std::cout << BaryCoord[0] << " " << BaryCoord[1] << " " << BaryCoord[2] << " " << BaryCoord[3] << std::endl;
      InterpolatedDispl = {0.0,0.0,0.0};
      return;
    }
    
    for(unsigned vertex_no = 0; vertex_no < 4; vertex_no++)
    {
      InterpolatedDispl += (LaserTracks[PointIndex[vertex_no].first].GetCorrection(PointIndex[vertex_no].second) * BaryCoord[vertex_no]);
//       std::cout << LaserTracks[PointIndex[vertex_no].first].GetCorrection(PointIndex[vertex_no].second)[0] << " " << LaserTracks[PointIndex[vertex_no].first].GetCorrection(PointIndex[vertex_no].second)[1] << " " << LaserTracks[PointIndex[vertex_no].first].GetCorrection(PointIndex[vertex_no].second)[2] << std::endl;
//       for(unsigned coord = 0; coord < 3; coord++)  std::cout << (LaserTracks[PointIndex[vertex_no].first].GetCorrection(PointIndex[vertex_no].second)[coord] * BaryCoord[vertex_no]) << " ";
//       std::cout << std::endl;
    }
    
    InterpolatedDispl = ThreeVector<float>::DotProduct(InterpolatedDispl,Track.GetPoyntingVector()) * Track.GetPoyntingVector();
    Track.AddToCorrection(InterpolatedDispl,sample_no);
  }
}
