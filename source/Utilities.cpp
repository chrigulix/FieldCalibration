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