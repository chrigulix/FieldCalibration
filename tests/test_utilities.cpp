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

TEST(TestDownsampler, SplitSingleEvenTrack) {

    TVector3 entry(1.,1.,1.);
    TVector3 exit(0.,0.,0.);

    std::vector<TVector3> RecobLaserTrack;

    // fill the dummy vector in the right order. Tracks are assumed to be filled ordered,
    // starting at the entry point and ending closest to the exit point.
    for (int i=0; i < 10; i++){
        const float pt = i;
        RecobLaserTrack.push_back(TVector3(pt,pt,pt));
    }

    // Make a laser track out of entry, exit and track points. Then add it to the Laser collection
    LaserTrack Track1 = LaserTrack(entry, exit, RecobLaserTrack);
    std::vector<LaserTrack> PewPew{ Track1 };
    Laser Las(PewPew);

    std::vector<Laser> LaserSets = SplitTrackSet(Las, 2);

    auto first_set = LaserSets.front();
    auto first_track = first_set.GetFirstTrack();

    ASSERT_TRUE(first_track.GetReco().front() == ThreeVector<float>(0., 0., 0.));
    ASSERT_TRUE(first_track.GetReco().back() == ThreeVector<float>(8., 8., 8.));
    ASSERT_EQ(first_track.GetNumberOfSamples(), 5);

    auto second_set = LaserSets.back();
    auto second_track = second_set.GetFirstTrack();

    ASSERT_TRUE(second_track.GetReco().front() == ThreeVector<float>(1., 1., 1.));
    ASSERT_TRUE(second_track.GetReco().back() == ThreeVector<float>(9., 9., 9.));
    ASSERT_EQ(second_track.GetNumberOfSamples(), 5);
}

TEST(TestDownsampler, SplitSingleOddTrack) {

    TVector3 entry(1.,1.,1.);
    TVector3 exit(0.,0.,0.);

    std::vector<TVector3> RecobLaserTrack;

    // fill the dummy vector in the right order. Tracks are assumed to be filled ordered,
    // starting at the entry point and ending closest to the exit point.
    for (int i=0; i < 11; i++){
        const float pt = i;
        RecobLaserTrack.push_back(TVector3(pt,pt,pt));
    }

    // Make a laser track out of entry, exit and track points. Then add it to the Laser collection
    LaserTrack Track1 = LaserTrack(entry, exit, RecobLaserTrack);
    std::vector<LaserTrack> PewPew{ Track1 };
    Laser Las(PewPew);

    std::vector<Laser> LaserSets = SplitTrackSet(Las, 2);

    auto first_set = LaserSets.front();
    auto first_track = first_set.GetFirstTrack();

    ASSERT_TRUE(first_track.GetReco().front() == ThreeVector<float>(0., 0., 0.));
    ASSERT_TRUE(first_track.GetReco().back() == ThreeVector<float>(10., 10., 10.));
    ASSERT_EQ(first_track.GetNumberOfSamples(), 6);

    auto second_set = LaserSets.back();
    auto second_track = second_set.GetFirstTrack();

    ASSERT_TRUE(second_track.GetReco().front() == ThreeVector<float>(1., 1., 1.));
    ASSERT_TRUE(second_track.GetReco().back() == ThreeVector<float>(9., 9., 9.));
    ASSERT_EQ(second_track.GetNumberOfSamples(), 5);
}

TEST(TestDownsampler, HigherOrder) {

    TVector3 entry(1.,1.,1.);
    TVector3 exit(0.,0.,0.);

    std::vector<TVector3> RecobLaserTrack;

    // fill the dummy vector in the right order. Tracks are assumed to be filled ordered,
    // starting at the entry point and ending closest to the exit point.
    for (int i=1; i < 11; i++){
        const float pt = i;
        RecobLaserTrack.push_back(TVector3(pt,pt,pt));
    }

    // Make a laser track out of entry, exit and track points. Then add it to the Laser collection
    LaserTrack Track1 = LaserTrack(entry, exit, RecobLaserTrack);
    std::vector<LaserTrack> PewPew{ Track1 };
    Laser Las(PewPew);

    std::vector<Laser> LaserSets = SplitTrackSet(Las, 4);

    ASSERT_EQ(LaserSets.size(), 4);

    auto first_set = LaserSets[0];
    auto first_track = first_set.GetFirstTrack();

    ASSERT_TRUE(first_track.GetReco().front() == ThreeVector<float>(1., 1., 1.));
    ASSERT_TRUE(first_track.GetReco()[1] == ThreeVector<float>(5., 5., 5.));
    ASSERT_TRUE(first_track.GetReco().back() == ThreeVector<float>(9., 9., 9.));
    ASSERT_EQ(first_track.GetNumberOfSamples(), 3);

    auto second_set = LaserSets.back();
    auto second_track = second_set.GetFirstTrack();

    ASSERT_TRUE(second_track.GetReco().front() == ThreeVector<float>(4., 4., 4.));
    ASSERT_TRUE(second_track.GetReco().back() == ThreeVector<float>(8., 8., 8.));
    ASSERT_EQ(second_track.GetNumberOfSamples(), 2);

    ASSERT_EQ(LaserSets[1].GetFirstTrack().GetNumberOfSamples(), 3);
    ASSERT_EQ(LaserSets[2].GetFirstTrack().GetNumberOfSamples(), 2);
}

TEST(TestDownsampler, TwoTracks) {

    TVector3 entry(1.,1.,1.);
    TVector3 exit(0.,0.,0.);

    std::vector<TVector3> RecobLaserTrack1;
    std::vector<TVector3> RecobLaserTrack2;

    // fill the dummy vector in the right order. Tracks are assumed to be filled ordered,
    // starting at the entry point and ending closest to the exit point.
    for (int i=1; i < 11; i++){
        const float pt = i;
        RecobLaserTrack1.push_back(TVector3(pt,pt,pt));
        RecobLaserTrack1.push_back(TVector3(-pt,-pt,-pt));
    }

    // Make a laser track out of entry, exit and track points. Then add it to the Laser collection
    LaserTrack Track1 = LaserTrack(entry, exit, RecobLaserTrack1);
    LaserTrack Track2 = LaserTrack(entry, exit, RecobLaserTrack2);
    std::vector<LaserTrack> PewPew{ Track1, Track2 };
    Laser Las(PewPew);

    std::vector<Laser> LaserSets = SplitTrackSet(Las, 2);

    ASSERT_EQ(LaserSets.size(), 2);

    auto first_set = LaserSets[0];
    auto second_set = LaserSets[1];

    ASSERT_EQ( first_set.GetNumberOfTracks(), 2);
    ASSERT_EQ(second_set.GetNumberOfTracks(), 2);

    auto reco_set0_track0 = first_set.GetTrack(0);
    auto reco_set0_track1 = first_set.GetTrack(1);

    auto reco_set1_track0 = second_set.GetTrack(0);
    auto reco_set1_track1 = second_set.GetTrack(1);


    ASSERT_TRUE(reco_set0_track0.GetReco().front() == ThreeVector<float>(1., 1., 1.));
    ASSERT_TRUE(reco_set1_track0.GetReco().front() == ThreeVector<float>(2., 2., 2.));

    ASSERT_TRUE(reco_set0_track1.GetReco().front() == ThreeVector<float>(-1., -1., -1.));
    ASSERT_TRUE(reco_set1_track1.GetReco().front() == ThreeVector<float>(-2., -2., -2.));

}


int main(int ac, char* av[])
{
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}