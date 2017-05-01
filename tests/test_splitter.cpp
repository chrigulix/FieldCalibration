//
// Created by matthias on 28.04.17.
//

#include <TVector3.h>

#include <gtest/gtest.h>
#include "../include/Utilities.hpp"

// Tests factorial of 0.
TEST(TestSplitter, CheckSelection) {

    TVector3 entry(1.,1.,1.);
    TVector3 exit(0.,0.,0.);

    std::vector<TVector3> RecobLaserTrack;

    // fill the dummy vector in the right order. Tracks are assumed to be filled ordered,
    // starting at the entry point and ending closest to the exit point.
    for (int i=9; i > 0; i--){
        const float pt = i / 10.0;
        RecobLaserTrack.push_back(TVector3(pt,pt,pt));
    }

    // Make a laser track out of entry, exit and track points. Then add it to the Laser collection
    LaserTrack Track1 = LaserTrack(entry, exit, RecobLaserTrack);
    std::vector<LaserTrack> PewPew{ Track1 };
    Laser Las(PewPew);

    // Now the actual testing of the function follows
    // Check if track gets into the first selection:
    auto res = ReachedExitPoint(PewPew, 0.5);
    auto InSelection = res.front();
    auto OutSelection = res.back();
    ASSERT_EQ(InSelection.GetNumberOfTracks(), 1);
    ASSERT_EQ(OutSelection.GetNumberOfTracks(), 0);
    res.clear();

    // Check if track gets into the second selection under new condition:
    res = ReachedExitPoint(PewPew, 0.001);
    InSelection = res.front();
    OutSelection = res.back();
    ASSERT_EQ(InSelection.GetNumberOfTracks(), 0);
    ASSERT_EQ(OutSelection.GetNumberOfTracks(), 1);

    //float a = 1;
    ASSERT_TRUE(true);
    //EXPECT_EQ(1, ReachedExitPoint(Test, 1));
}

int main(int ac, char* av[])
{
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}