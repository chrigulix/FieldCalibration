#include "../include/Laser.hpp"

Laser::Laser() {}

Laser::Laser(std::vector<LaserTrack> InputSet)
{
  LaserTrackSet = InputSet;
}

std::vector<LaserTrack>::iterator Laser::begin()
{
  return LaserTrackSet.begin();
}

std::vector<LaserTrack>::iterator Laser::end()
{
  return LaserTrackSet.end();
}

void Laser::AppendTrack(const LaserTrack& InputTrack)
{
  LaserTrackSet.push_back(InputTrack);
}

LaserTrack Laser::GetTrack(const long unsigned int& TrackNumber)
{
  LaserTrackSet.at(TrackNumber);
}

LaserTrack Laser::GetLastTrack()
{
  return LaserTrackSet.back();
}

LaserTrack Laser::GetFirstTrack()
{
  return LaserTrackSet.front();
}

long unsigned int Laser::GetNumberOfTracks() const
{
  return LaserTrackSet.size();
}

std::vector<LaserTrack> Laser::GetTrackSet() const
{
  return LaserTrackSet;
}

void Laser::DistortTrackSet(std::string MapFileName, TPCVolumeHandler& TPCVolume)
{
  LaserTrack::DistortTracks(LaserTrackSet,MapFileName,TPCVolume);
}

// Applies correction algorithm on all tracks of a laser track set
void Laser::CorrectTrackSet()
{
    // Loop over all tracks, correct every single one
    for(auto& Track : LaserTrackSet)
    {
        Track.CorrectTrack();
    }
}

void Laser::InterpolateTrackSet(const std::vector<LaserTrack>& LaserTracks, const Delaunay& Mesh)
{
  for(auto& Track : LaserTrackSet)
  {
    InterpolateTrack(Track,LaserTracks,Mesh);
  }
}

void Laser::InterpolateTrackSet(const Laser& LaserTracks, const Delaunay& Mesh)
{
  InterpolateTrackSet(LaserTracks.GetTrackSet(), Mesh);
}

void Laser::DrawTrack(const long unsigned int& TrackNumber)
{
  
  TApplication theApp("App",0,0);
  
//   std::cout << LaserTrackSet.at(2).GetNumberOfSamples() << std::endl;
  TGraph2D* TrueTrack = new TGraph2D(2);
  TPolyLine3D* DistortedTrack = new TPolyLine3D(LaserTrackSet.at(TrackNumber).GetNumberOfSamples());
  TPolyLine3D* CorrectedTrack = new TPolyLine3D(LaserTrackSet.at(TrackNumber).GetNumberOfSamples());
  
  
  for(unsigned sample_no = 0; sample_no < LaserTrackSet.at(TrackNumber).GetNumberOfSamples(); sample_no++)
  {
    DistortedTrack->SetPoint(sample_no,
			     LaserTrackSet.at(TrackNumber).GetSamplePosition(sample_no).at(0),
			     LaserTrackSet.at(TrackNumber).GetSamplePosition(sample_no).at(1),
			     LaserTrackSet.at(TrackNumber).GetSamplePosition(sample_no).at(2));
    CorrectedTrack->SetPoint(sample_no,
			     LaserTrackSet.at(TrackNumber).GetSamplePosition(sample_no).at(0)+LaserTrackSet.at(TrackNumber).GetCorrection(sample_no).at(0),
			     LaserTrackSet.at(TrackNumber).GetSamplePosition(sample_no).at(1)+LaserTrackSet.at(TrackNumber).GetCorrection(sample_no).at(1),
			     LaserTrackSet.at(TrackNumber).GetSamplePosition(sample_no).at(2)+LaserTrackSet.at(TrackNumber).GetCorrection(sample_no).at(2));
  }
  
  TrueTrack->SetPoint(0,LaserTrackSet.at(TrackNumber).GetEntryPoint().at(0),LaserTrackSet.at(TrackNumber).GetEntryPoint().at(1),LaserTrackSet.at(TrackNumber).GetEntryPoint().at(2));
  TrueTrack->SetPoint(1,LaserTrackSet.at(TrackNumber).GetExitPoint().at(0),LaserTrackSet.at(TrackNumber).GetExitPoint().at(1),LaserTrackSet.at(TrackNumber).GetExitPoint().at(2));
  
//   std::cout << LaserTrackSet.at(TrackNumber).GetEntryPoint().at(0) << " " << LaserTrackSet.at(TrackNumber).GetEntryPoint().at(1) << " " << LaserTrackSet.at(TrackNumber).GetEntryPoint().at(2) << std::endl;
//   std::cout << LaserTrackSet.at(TrackNumber).GetSamplePosition(0).at(0)+LaserTrackSet.at(TrackNumber).GetCorrection(0).at(0) << " " << 
// 	       LaserTrackSet.at(TrackNumber).GetSamplePosition(0).at(1)+LaserTrackSet.at(TrackNumber).GetCorrection(0).at(1) << " " << 
// 	       LaserTrackSet.at(TrackNumber).GetSamplePosition(0).at(2)+LaserTrackSet.at(TrackNumber).GetCorrection(0).at(2) << std::endl;
  
  std::string PictureName = "Track_" + std::to_string(TrackNumber);
  
  TrueTrack->SetLineColor(kBlue);
  DistortedTrack->SetLineColor(kRed);
  CorrectedTrack->SetLineColor(kGreen);
  
  TCanvas * C0 = new TCanvas(PictureName.c_str(),PictureName.c_str(),1000,700);
  TrueTrack->Draw("LINE");
  DistortedTrack->Draw("same");
  CorrectedTrack->Draw("same");
  C0 -> Print((PictureName+".png").c_str(),"png");
  
  theApp.Run();
  
  delete DistortedTrack;
  delete CorrectedTrack;
  delete TrueTrack;
  delete C0;
  gDirectory->GetList()->Delete();
}



