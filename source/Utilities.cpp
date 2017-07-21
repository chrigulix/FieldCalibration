//
// Created by matthias on 28.04.17.
//

#include "../include/Utilities.hpp"

std::vector<Laser> ReachedExitPoint(const Laser& LaserSet, float ExitBoundary) {

    std::vector<Laser> Selection;
    Selection.resize(2);

    for(auto& Track : LaserSet.GetTrackSet())
    {
        auto ExitToLast = Track.GetExitPoint() - Track.GetBack();
        auto d = ExitToLast.GetNorm();

        if (d < ExitBoundary) {
            Selection.front().AppendTrack(Track);
        }
        else {
            Selection.back().AppendTrack(Track);
        }
    }
    return Selection;
}

std::vector<Laser> SplitTrackSet(const Laser& LaserSet, unsigned int Downsample) {
    std::vector<Laser> Sets;
    Sets.resize(Downsample);

    for(auto& Track : LaserSet.GetTrackSet())
    {

        auto SourceTrack = Track.GetReco();
        unsigned long TrackSize = Track.GetNumberOfSamples();

        for (unsigned long offset=0; offset < Sets.size(); offset++){
            std::vector<ThreeVector<float>> SampledRecoTrack;

            for (unsigned long idx=offset; idx < SourceTrack.size(); idx += Downsample) {
                SampledRecoTrack.push_back(SourceTrack[idx]);
            }
            LaserTrack SampledTrack(Track.GetEntryPoint(), Track.GetExitPoint(), SampledRecoTrack);
            Sets[offset].AppendTrack(SampledTrack);
            SampledRecoTrack.clear();

        }
    }
    return Sets;
}